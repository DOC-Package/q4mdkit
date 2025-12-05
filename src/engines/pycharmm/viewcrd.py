#!/usr/bin/env python3
"""Visualize CHARMM CRD file using various methods."""
import argparse
import os


def crd_to_pdb_mda(psf_path: str, crd_path: str, pdb_path: str):
    """
    Convert CRD to PDB using MDAnalysis.
    
    Parameters
    ----------
    psf_path : str
        Path to PSF file
    crd_path : str
        Path to input CRD file
    pdb_path : str
        Path to output PDB file
    """
    try:
        import MDAnalysis as mda
    except ImportError:
        print("Error: MDAnalysis not installed")
        print("Install with: conda install -c conda-forge MDAnalysis")
        print("         or: pip install MDAnalysis")
        return False
    
    print(f"Converting CRD to PDB with MDAnalysis...")
    print(f"  PSF: {psf_path}")
    print(f"  CRD: {crd_path}")
    print(f"  PDB: {pdb_path}")
    
    u = mda.Universe(psf_path, crd_path, format='CRD')
    u.atoms.write(pdb_path)
    print(f"Successfully converted: {crd_path} -> {pdb_path}")
    return True


def view_crd_nglview(psf_path: str, crd_path: str):
    """
    Open interactive viewer with nglview (Jupyter only).
    
    Parameters
    ----------
    psf_path : str
        Path to PSF file
    crd_path : str
        Path to CRD file
    
    Returns
    -------
    nglview widget (Jupyter only)
    """
    try:
        import MDAnalysis as mda
        import nglview as nv
    except ImportError:
        print("Error: MDAnalysis and nglview required")
        print("Install with: conda install -c conda-forge MDAnalysis nglview")
        print("         or: pip install MDAnalysis nglview")
        return None
    
    print(f"Loading structure...")
    print(f"  PSF: {psf_path}")
    print(f"  CRD: {crd_path}")
    
    u = mda.Universe(psf_path, crd_path, format='CRD')
    view = nv.show_mdanalysis(u)
    print("Interactive viewer created (use in Jupyter)")
    return view


def open_with_vmd(psf_path: str, crd_path: str):
    """
    Open CRD file with VMD.
    
    Parameters
    ----------
    psf_path : str
        Path to PSF file
    crd_path : str
        Path to CRD file
    """
    import subprocess
    
    print(f"Opening with VMD...")
    print(f"  PSF: {psf_path}")
    print(f"  CRD: {crd_path}")
    
    try:
        subprocess.run(['vmd', psf_path, crd_path])
    except FileNotFoundError:
        print("Error: VMD not found in PATH")
        print("Make sure VMD is installed and accessible")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize or convert CHARMM CRD files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert to PDB using MDAnalysis
  viewcrd --psf system.psf --crd coords.crd --pdb output.pdb
  
  # Open with VMD
  viewcrd --psf system.psf --crd coords.crd --vmd
  
  # Interactive viewer in Jupyter
  viewcrd --psf system.psf --crd coords.crd --view
"""
    )
    parser.add_argument("--psf", required=True, help="PSF file")
    parser.add_argument("--crd", required=True, help="CRD file")
    parser.add_argument("--pdb", help="Output PDB file (convert mode)")
    parser.add_argument("--view", action="store_true", help="Open nglview (Jupyter only)")
    parser.add_argument("--vmd", action="store_true", help="Open with VMD")
    
    args = parser.parse_args()
    
    # Check files exist
    if not os.path.exists(args.psf):
        print(f"Error: PSF file not found: {args.psf}")
        return
    if not os.path.exists(args.crd):
        print(f"Error: CRD file not found: {args.crd}")
        return
    
    # Execute requested action
    if args.pdb:
        crd_to_pdb_mda(args.psf, args.crd, args.pdb)
    elif args.vmd:
        open_with_vmd(args.psf, args.crd)
    elif args.view:
        view_crd_nglview(args.psf, args.crd)
    else:
        print("Error: Specify --pdb, --vmd, or --view")
        parser.print_help()


if __name__ == "__main__":
    main()
