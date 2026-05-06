#!/usr/bin/env python3
"""Generate qmatoms list from PDB file for QM/MM calculations.

Usage:
    python make_qmatoms.py input.pdb [options]

Options:
    --resid N       Select residue by number (default: 1)
    --resname NAME  Select residue by name (e.g., QUC)
    --output FILE   Output file (default: qmatoms.txt)
"""

import argparse
from pathlib import Path


def get_qmatoms(pdb_file, resid=None, resname=None):
    """Extract atom indices for QM region from PDB file.
    
    Args:
        pdb_file: Path to PDB file
        resid: Residue number to select (1-indexed)
        resname: Residue name to select
    
    Returns:
        List of atom indices (0-indexed for ASH)
    """
    atoms = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                # PDB serial is 1-indexed; convert to 0-indexed for ASH
                atom_idx = int(line[6:11]) - 1
                res_name = line[17:20].strip()
                res_num = int(line[22:26])
                
                if resid is not None and res_num == resid:
                    atoms.append(atom_idx)
                elif resname is not None and res_name == resname:
                    atoms.append(atom_idx)
    
    return atoms


def main():
    parser = argparse.ArgumentParser(description='Generate qmatoms for QM/MM')
    parser.add_argument('pdb', help='Input PDB file')
    parser.add_argument('--resid', type=int, default=1, help='Residue number (default: 1)')
    parser.add_argument('--resname', type=str, help='Residue name (overrides resid)')
    parser.add_argument('--output', '-o', default='qmatoms', help='Output file')
    args = parser.parse_args()
    
    resid = None if args.resname else args.resid
    resname = args.resname
    
    atoms = get_qmatoms(args.pdb, resid=resid, resname=resname)
    
    if not atoms:
        print(f"No atoms found for resid={resid} or resname={resname}")
        return
    
    # Write qmatoms file (indices only, space-separated)
    with open(args.output, 'w') as f:
        f.write(' '.join(map(str, atoms)))
    
    print(f"Written {len(atoms)} atoms to {args.output}")


if __name__ == '__main__':
    main()
