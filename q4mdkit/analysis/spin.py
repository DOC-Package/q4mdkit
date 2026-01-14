"""
Compute spin populations for fragments from CDFTB calculations.

This module reads the α and β Mulliken populations from DFTB+ calculations
and computes spin expectation values for each fragment.

Spin expectation value: ⟨Sz⟩_frag = (N_α_frag - N_β_frag) / 2
"""

import numpy as np
from pathlib import Path
from typing import Dict, Tuple
import yaml

from .cdftb_result_reader import (
    load_spin_polarized_calculation,
    get_orbital_info_from_eigenvec,
    get_atom_orbital_map,
)


def compute_fragment_spin_population(
    C_alpha: np.ndarray,
    C_beta: np.ndarray,
    S_AO: np.ndarray,
    n_alpha: int,
    n_beta: int,
    fragment_atoms: list,
    atom_to_orbitals: dict
) -> Dict[str, float]:
    """
    Compute spin population for a fragment.
    
    Parameters
    ----------
    C_alpha : np.ndarray
        α MO coefficients (n_AO x n_MO).
    C_beta : np.ndarray
        β MO coefficients (n_AO x n_MO).
    S_AO : np.ndarray
        AO overlap matrix (n_AO x n_AO).
    n_alpha : int
        Number of α electrons.
    n_beta : int
        Number of β electrons.
    fragment_atoms : list
        List of atom indices (1-indexed) in the fragment.
    atom_to_orbitals : dict
        Mapping from atom index to list of orbital indices.
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'n_alpha_frag': Number of α electrons on fragment (Mulliken)
        - 'n_beta_frag': Number of β electrons on fragment (Mulliken)
        - 'spin_density': Spin density (n_α - n_β)
    """
    # Extract occupied orbitals
    C_alpha_occ = C_alpha[:, :n_alpha]
    C_beta_occ = C_beta[:, :n_beta]
    
    # Density matrices: P = C @ C^T
    P_alpha = C_alpha_occ @ C_alpha_occ.T
    P_beta = C_beta_occ @ C_beta_occ.T
    
    # Mulliken population on fragment: N_frag = Tr(P @ S @ W)
    # where W_μν = δ_{μ∈frag} for fragment projection
    # For Mulliken: N_frag = Σ_{μ∈frag} Σ_ν P_μν S_νμ = Σ_{μ∈frag} (P @ S)_μμ
    
    # Get orbital indices for fragment
    fragment_orbitals = []
    for atom_idx in fragment_atoms:
        if atom_idx in atom_to_orbitals:
            fragment_orbitals.extend(atom_to_orbitals[atom_idx])
    
    # Compute PS matrices
    PS_alpha = P_alpha @ S_AO
    PS_beta = P_beta @ S_AO
    
    # Sum over fragment orbitals
    n_alpha_frag = sum(PS_alpha[orb_idx, orb_idx] for orb_idx in fragment_orbitals)
    n_beta_frag = sum(PS_beta[orb_idx, orb_idx] for orb_idx in fragment_orbitals)
    
    # Spin density
    spin_density = n_alpha_frag - n_beta_frag
    
    return {
        'n_alpha_frag': n_alpha_frag,
        'n_beta_frag': n_beta_frag,
        'spin_density': spin_density,
    }


def compute_spin_populations_from_cdftb(config_path: Path) -> None:
    """
    Compute spin populations for all fragments from CDFTB output.
    
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
        raise ValueError("This function requires exactly 2 fragments")
    
    frag1_name = fragments[0]['name']
    frag2_name = fragments[1]['name']
    
    print("=" * 70)
    print("Spin Population Analysis")
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
            "Please run CDFTB calculations first."
        )
    
    energies_data = np.loadtxt(energies_file, comments='#')
    if energies_data.ndim == 1:
        energies_data = energies_data.reshape(1, -1)
    
    n_frames = energies_data.shape[0]
    print(f"\nFound {n_frames} frames in energies.dat\n")
    
    # Prepare output files
    spin_output_file_A = output_dir / "spin1.dat"
    spin_output_file_B = output_dir / "spin2.dat"
    
    with open(spin_output_file_A, 'w') as f:
        f.write("# Spin Populations from CDFTB Calculations (State A)\n")
        f.write("# Spin density = N_α - N_β on each fragment\n")
        f.write(f"# Frame  Time(fs)  ")
        f.write(f"{frag1_name}_n_alpha  {frag1_name}_n_beta  {frag1_name}_spin  ")
        f.write(f"{frag2_name}_n_alpha  {frag2_name}_n_beta  {frag2_name}_spin\n")
    
    with open(spin_output_file_B, 'w') as f:
        f.write("# Spin Populations from CDFTB Calculations (State B)\n")
        f.write("# Spin density = N_α - N_β on each fragment\n")
        f.write(f"# Frame  Time(fs)  ")
        f.write(f"{frag1_name}_n_alpha  {frag1_name}_n_beta  {frag1_name}_spin  ")
        f.write(f"{frag2_name}_n_alpha  {frag2_name}_n_beta  {frag2_name}_spin\n")
    
    # Process each frame
    for i in range(n_frames):
        frame_id = int(energies_data[i, 0])
        time_fs = energies_data[i, 1]
        
        print(f"Frame {frame_id:05d} (t = {time_fs:.1f} fs)")
        
        # Fragment directories
        frag1_dir = output_dir / f"frame_{frame_id:05d}" / frag1_name
        frag2_dir = output_dir / f"frame_{frame_id:05d}" / frag2_name
        
        if not frag1_dir.exists() or not frag2_dir.exists():
            print(f"  WARNING: Fragment directories not found, skipping frame")
            continue
        
        try:
            # Load spin-polarized data for both fragments
            data_A = load_spin_polarized_calculation(frag1_dir)
            data_B = load_spin_polarized_calculation(frag2_dir)
            
            # Get orbital-to-atom mapping (same for both calculations)
            eigenvec_file = frag1_dir / "eigenvec.out"
            orbital_info = get_orbital_info_from_eigenvec(eigenvec_file)
            atom_to_orbitals = get_atom_orbital_map(orbital_info)
            n_atoms = len(atom_to_orbitals)
            n_atoms_half = n_atoms // 2
            
            # Define fragment atoms (1-indexed)
            frag_A_atoms = list(range(1, n_atoms_half + 1))
            frag_B_atoms = list(range(n_atoms_half + 1, n_atoms + 1))
            
            # Compute spin populations for both fragments in state A
            spin_A_frag1 = compute_fragment_spin_population(
                data_A['C_alpha'], data_A['C_beta'], data_A['S'],
                data_A['n_alpha'], data_A['n_beta'],
                frag_A_atoms, atom_to_orbitals
            )
            spin_A_frag2 = compute_fragment_spin_population(
                data_A['C_alpha'], data_A['C_beta'], data_A['S'],
                data_A['n_alpha'], data_A['n_beta'],
                frag_B_atoms, atom_to_orbitals
            )
            
            # Compute spin populations for both fragments in state B
            spin_B_frag1 = compute_fragment_spin_population(
                data_B['C_alpha'], data_B['C_beta'], data_B['S'],
                data_B['n_alpha'], data_B['n_beta'],
                frag_A_atoms, atom_to_orbitals
            )
            spin_B_frag2 = compute_fragment_spin_population(
                data_B['C_alpha'], data_B['C_beta'], data_B['S'],
                data_B['n_alpha'], data_B['n_beta'],
                frag_B_atoms, atom_to_orbitals
            )
            
            print(f"  State A: {frag1_name} spin={spin_A_frag1['spin_density']:+.4f}, "
                  f"{frag2_name} spin={spin_A_frag2['spin_density']:+.4f}")
            print(f"  State B: {frag1_name} spin={spin_B_frag1['spin_density']:+.4f}, "
                  f"{frag2_name} spin={spin_B_frag2['spin_density']:+.4f}")
            
            # Write to files
            with open(spin_output_file_A, 'a') as f:
                f.write(f"{frame_id:5d}  {time_fs:8.2f}  ")
                f.write(f"{spin_A_frag1['n_alpha_frag']:16.8f}  {spin_A_frag1['n_beta_frag']:16.8f}  "
                        f"{spin_A_frag1['spin_density']:16.8f}  ")
                f.write(f"{spin_A_frag2['n_alpha_frag']:16.8f}  {spin_A_frag2['n_beta_frag']:16.8f}  "
                        f"{spin_A_frag2['spin_density']:16.8f}\n")
            
            with open(spin_output_file_B, 'a') as f:
                f.write(f"{frame_id:5d}  {time_fs:8.2f}  ")
                f.write(f"{spin_B_frag1['n_alpha_frag']:16.8f}  {spin_B_frag1['n_beta_frag']:16.8f}  "
                        f"{spin_B_frag1['spin_density']:16.8f}  ")
                f.write(f"{spin_B_frag2['n_alpha_frag']:16.8f}  {spin_B_frag2['n_beta_frag']:16.8f}  "
                        f"{spin_B_frag2['spin_density']:16.8f}\n")
            
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print("\n" + "=" * 70)
    print(f"Results saved to:")
    print(f"  State A: {spin_output_file_A}")
    print(f"  State B: {spin_output_file_B}")
    print("=" * 70)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Compute spin populations for fragments from CDFTB calculations"
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
    
    compute_spin_populations_from_cdftb(config_path)
