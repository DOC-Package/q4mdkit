#!/usr/bin/env python3
"""Convert CHARMM CRD to PDB format using pyCHARMM."""
import argparse
import os
import pycharmm
from pycharmm import lingo


def crd_to_pdb(psf_path: str, crd_path: str, pdb_path: str):
    """
    Convert CRD to PDB using pyCHARMM.
    
    Parameters
    ----------
    psf_path : str
        Path to PSF file
    crd_path : str
        Path to input CRD file
    pdb_path : str
        Path to output PDB file
    """
    # Convert to absolute paths
    psf_path = os.path.abspath(psf_path)
    crd_path = os.path.abspath(crd_path)
    pdb_path = os.path.abspath(pdb_path)
    
    print(f"Converting CRD to PDB...")
    print(f"  PSF: {psf_path}")
    print(f"  CRD: {crd_path}")
    print(f"  PDB: {pdb_path}")
    
    # Load base CHARMM force field (if present) before reading PSF so atom types exist
    base_rtf = "/home/takahashi/charmm/charmm/toppar/top_all36_cgenff.rtf"
    base_prm = "/home/takahashi/charmm/charmm/toppar/par_all36_cgenff.prm"
    if os.path.exists(base_rtf) and os.path.exists(base_prm):
        print(f"Loading base CHARMM36 CGenFF force field from {base_rtf} and {base_prm} ...")
        lingo.charmm_script(f"BOMLEV -2\nopen read unit 1 card name {base_rtf}\nread rtf card unit 1\nclose unit 1\nopen read unit 1 card name {base_prm}\nread param card unit 1\nclose unit 1\nBOMLEV 0\n")
    else:
        print("Warning: base CHARMM toppar not found at expected location; PSF read may fail if it requires these types.")

    # Read PSF
    lingo.charmm_script(f"""
BOMLEV -2
open read unit 10 card name {psf_path}
read psf card unit 10
close unit 10
BOMLEV 0
""")
    
    # Read CRD
    lingo.charmm_script(f"""
open read unit 20 card name {crd_path}
read coor card unit 20
close unit 20
""")
    
    # Write PDB
    lingo.charmm_script(f"""
open write unit 30 card name {pdb_path}
write coor pdb unit 30
close unit 30
""")
    
    print(f"Successfully converted: {crd_path} -> {pdb_path}")


def main():
    parser = argparse.ArgumentParser(description="Convert CHARMM CRD to PDB using pyCHARMM")
    parser.add_argument("--psf", required=True, help="Input PSF file")
    parser.add_argument("--crd", required=True, help="Input CRD file")
    parser.add_argument("--pdb", required=True, help="Output PDB file")
    
    args = parser.parse_args()
    crd_to_pdb(args.psf, args.crd, args.pdb)


if __name__ == "__main__":
    main()
