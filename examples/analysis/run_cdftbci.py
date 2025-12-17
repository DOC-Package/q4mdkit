#!/usr/bin/env python3
"""
CDFTB-CI analysis script.

Compute the CDFTB-CI Hamiltonian and transfer integrals from CDFTB output.

Usage:
    python run_cdftbci.py [config.yaml]

If no config file is specified, defaults to 'cdftb_settings.yaml' in the current directory.
"""

import argparse
import numpy as np
from pathlib import Path
import yaml

from qm4d4crystal.analysis import (
    UnrestrictedOrbitalData,
    build_cdftbci_hamiltonian_unrestricted,
    solve_cdftbci_unrestricted,
    compute_transfer_integral_unrestricted,
)
from qm4d4crystal.analysis.cdftb_result_reader import (
    load_spin_polarized_calculation,
    get_orbital_info_from_eigenvec,
    get_atom_orbital_map,
    build_fragment_weight_matrix,
)


def run_cdftbci_analysis(config_path: Path) -> None:
    """
    Run CDFTB-CI analysis on existing CDFTB output.
    
    Parameters
    ----------
    config_path : Path
        Path to YAML configuration file.
    """
    config_path = Path(config_path).resolve()
    base_dir = config_path.parent
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Get output directory
    output_cfg = config.get('output', {})
    output_dir = base_dir / output_cfg.get('directory', 'output')
    
    # Get fragment info
    fragments = config.get('fragments', [])
    if len(fragments) != 2:
        raise ValueError("CDFTB-CI requires exactly 2 fragments")
    
    frag1_name = fragments[0]['name']
    frag2_name = fragments[1]['name']
    
    print("=" * 70)
    print("CDFTB-CI Analysis")
    print("=" * 70)
    print(f"Config file: {config_path}")
    print(f"Output directory: {output_dir}")
    print(f"Fragment A: {frag1_name}")
    print(f"Fragment B: {frag2_name}")
    print("=" * 70)
    
    # Read energies from energies.dat
    energies_file = output_dir / "energies.dat"
    if not energies_file.exists():
        raise FileNotFoundError(
            f"Energies file not found: {energies_file}\n"
            "Please run run_cdftb.py first to generate CDFTB output."
        )
    
    energies_data = np.loadtxt(energies_file, comments='#')
    if energies_data.ndim == 1:
        energies_data = energies_data.reshape(1, -1)
    
    n_frames = energies_data.shape[0]
    print(f"\nFound {n_frames} frames in energies.dat")
    
    # Prepare output file
    ci_output_file = output_dir / "cdftbci_results.dat"
    
    with open(ci_output_file, 'w') as f:
        f.write("# CDFTB-CI Results\n")
        f.write("# Frame  Time(fs)  E_A(Ha)  E_B(Ha)  V_A(Ha)  V_B(Ha)  S_AB  "
                "H_AB(Ha)  J_direct(meV)  J_lowdin(meV)  E1(Ha)  E2(Ha)  dE(eV)\n")
    
    # Process each frame
    for i in range(n_frames):
        frame_id = int(energies_data[i, 0])
        time_fs = energies_data[i, 1]
        E_A = energies_data[i, 2]
        E_B = energies_data[i, 3]
        
        print(f"\nFrame {frame_id:05d} (t = {time_fs:.1f} fs)")
        
        # Fragment directories
        frag1_dir = output_dir / f"frame_{frame_id:05d}" / frag1_name
        frag2_dir = output_dir / f"frame_{frame_id:05d}" / frag2_name
        
        if not frag1_dir.exists() or not frag2_dir.exists():
            print(f"  WARNING: Fragment directories not found, skipping frame")
            continue
        
        try:
            # Load spin-polarized data
            data_A = load_spin_polarized_calculation(frag1_dir)
            data_B = load_spin_polarized_calculation(frag2_dir)
            
            # Create UnrestrictedOrbitalData objects
            orb_A = UnrestrictedOrbitalData(
                C_alpha=data_A['C_alpha'],
                C_beta=data_A['C_beta'],
                occ_alpha=data_A['occ_alpha'],
                occ_beta=data_A['occ_beta'],
                n_alpha=data_A['n_alpha'],
                n_beta=data_A['n_beta']
            )
            orb_B = UnrestrictedOrbitalData(
                C_alpha=data_B['C_alpha'],
                C_beta=data_B['C_beta'],
                occ_alpha=data_B['occ_alpha'],
                occ_beta=data_B['occ_beta'],
                n_alpha=data_B['n_alpha'],
                n_beta=data_B['n_beta']
            )
            
            # Read AO overlap
            S_AO = data_A['S']
            n_orbitals = data_A['n_orbitals']
            
            # Build weight matrices
            eigenvec_file = frag1_dir / "eigenvec.out"
            orbital_info = get_orbital_info_from_eigenvec(eigenvec_file)
            atom_to_orbitals = get_atom_orbital_map(orbital_info)
            n_atoms = len(atom_to_orbitals)
            n_atoms_half = n_atoms // 2
            
            frag_A_atoms = list(range(1, n_atoms_half + 1))
            frag_B_atoms = list(range(n_atoms_half + 1, n_atoms + 1))
            
            w_A = build_fragment_weight_matrix(S_AO, frag_A_atoms, n_atoms,
                                               atom_to_orbitals=atom_to_orbitals)
            w_B = build_fragment_weight_matrix(S_AO, frag_B_atoms, n_atoms,
                                               atom_to_orbitals=atom_to_orbitals)
            
            # Get constraint potentials
            V_A = data_A['Vc']
            V_B = data_B['Vc']
            
            # Target populations (Mulliken population for constrained fragment)
            # Read from config or use default
            ci_cfg = config.get('cdftb_ci', {})
            N_A = ci_cfg.get('N_A', 101)  # Default for pentacene cation
            N_B = ci_cfg.get('N_B', 101)
            
            # Debug output
            print(f"  n_alpha_A={orb_A.n_alpha}, n_beta_A={orb_A.n_beta}")
            print(f"  V_A={V_A:.6f} Ha, V_B={V_B:.6f} Ha")
            print(f"  N_A={N_A}, N_B={N_B}")
            
            # Build Hamiltonian
            ham = build_cdftbci_hamiltonian_unrestricted(
                orb_A, orb_B, S_AO, w_A, w_B,
                E_A, E_B, V_A, V_B, N_A, N_B
            )
            
            # Solve eigenvalue problem
            eigenvalues, eigenvectors = solve_cdftbci_unrestricted(ham.H, ham.S)
            
            # Compute transfer integrals
            J_direct = compute_transfer_integral_unrestricted(ham.H, ham.S, method="direct")
            J_lowdin = compute_transfer_integral_unrestricted(ham.H, ham.S, method="lowdin")
            
            # Convert to meV
            J_direct_meV = J_direct * 27211.386
            J_lowdin_meV = J_lowdin * 27211.386
            dE_eV = (eigenvalues[1] - eigenvalues[0]) * 27.211386
            
            print(f"  S_AB = {ham.S_AB:.6f} (α: {ham.S_AB_alpha:.6f}, β: {ham.S_AB_beta:.6f})")
            print(f"  W_BA = {ham.W_BA:.6f}, W_AB = {ham.W_AB:.6f}")
            print(f"  H_AB = {ham.H_AB:.6f} Ha")
            print(f"  J (direct) = {J_direct_meV:.2f} meV")
            print(f"  J (Löwdin) = {J_lowdin_meV:.2f} meV")
            print(f"  E1 = {eigenvalues[0]:.6f} Ha, E2 = {eigenvalues[1]:.6f} Ha")
            print(f"  ΔE = {dE_eV:.4f} eV")
            
            # Write results
            with open(ci_output_file, 'a') as f:
                f.write(f"{frame_id:6d}  {time_fs:8.2f}  {E_A:14.10f}  {E_B:14.10f}  "
                        f"{V_A:12.8f}  {V_B:12.8f}  {ham.S_AB:12.8f}  {ham.H_AB:14.10f}  "
                        f"{J_direct_meV:10.4f}  {J_lowdin_meV:10.4f}  "
                        f"{eigenvalues[0]:14.10f}  {eigenvalues[1]:14.10f}  {dE_eV:10.6f}\n")
            
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print("\n" + "=" * 70)
    print(f"Results saved to {ci_output_file}")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description="Run CDFTB-CI analysis on MD trajectory"
    )
    parser.add_argument(
        "config",
        nargs="?",
        default="cdftb_settings.yaml",
        help="Path to YAML configuration file (default: cdftb_settings.yaml)"
    )
    args = parser.parse_args()
    
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"Error: Configuration file not found: {config_path}")
        exit(1)
    
    print(f"Running CDFTB-CI analysis with config: {config_path}")
    run_cdftbci_analysis(config_path)


if __name__ == "__main__":
    main()
