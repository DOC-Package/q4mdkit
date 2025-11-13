from ase.io import read
import sys
from pathlib import Path

def build_topology_file(gro_file, itp_file, system_name, atoms_per_molecule,
                        forcefield="oplsaa.ff/forcefield.itp",
                        molecule_name="MOL",
                        output_file="topol.top"):
    """
    Build a GROMACS topology file from a GRO structure file for a molecular crystal.
    """
    # Read GRO file
    atoms = read(gro_file)
    
    # Calculate number of molecules
    n_atoms = len(atoms)
    n_molecules = n_atoms // atoms_per_molecule
    
    if n_atoms % atoms_per_molecule != 0:
        print(f"⚠️ WARNING: Total atoms ({n_atoms}) is not evenly divisible by atoms per molecule ({atoms_per_molecule})")
        print(f"   Remainder: {n_atoms % atoms_per_molecule} atoms")
    
    print(f"GRO file: {gro_file}")
    print(f"Atoms: {n_atoms}")
    print(f"Molecules: {n_molecules}")
    
    # Create topology file content
    top_content = f"""; Topology file for {system_name}
    ; Generated from {gro_file}

    ; Include forcefield
    #include "{forcefield}"

    ; Include molecule topology
    #include "{itp_file}"

    [ system ]
    {system_name}

    [ molecules ]
    ; Compound        #mols
    {molecule_name:<16s}  {n_molecules}
    """
    
    # Save topology file
    with open(output_file, "w") as f:
        f.write(top_content)
    
    print(f"✅ Created: {output_file}")

    return {'n_atoms': n_atoms, 'n_molecules': n_molecules, 'output_file': output_file, 'gro_file': gro_file, 'itp_file': itp_file}

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Create GROMACS topology file from GRO structure")
    parser.add_argument("--gro", help="Input GRO file")
    parser.add_argument("--itp", help="Molecule topology ITP file")
    parser.add_argument("--system-name", help="System name")
    parser.add_argument("--atoms-per-mol", type=int, help="Number of atoms per molecule")
    parser.add_argument("--ff", "--forcefield", default="oplsaa.ff/forcefield.itp", help="Forcefield ITP file")
    parser.add_argument("--mol-name", default="MOL", help="Molecule name in ITP file")
    parser.add_argument("-o", "--output", default="topol.top", help="Output topology file")
    args = parser.parse_args()
    build_topology_file(gro_file=args.gro, itp_file=args.itp, system_name=args.system_name, atoms_per_molecule=args.atoms_per_mol, forcefield=args.ff, molecule_name=args.mol_name, output_file=args.output)

if __name__ == "__main__":
    main()