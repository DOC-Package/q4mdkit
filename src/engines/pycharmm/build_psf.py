import argparse
import os
import tempfile
import re
import pycharmm
import pycharmm.lingo as lingo


def build_psf(
    pdb_path: str,
    crd_path: str,
    toppar_stream: str,
    resname: str = "MOL",
    nmol: int = 54,
    segid: str = "SYS",
    write_psf_path: str = "system.psf",
    write_pdb_path: str | None = None
):
    """
    Build CHARMM PSF topology using pyCHARMM.
    
    See docs/build_psf.md for detailed documentation.
    
    Parameters
    ----------
    pdb_path : str
        Template PDB file (single molecule)
    crd_path : str
        CRD file with all molecules (extended format)
    toppar_stream : str
        CGenFF stream file with topology/parameters
    resname : str, optional
        Residue name (default: "MOL")
    nmol : int, optional
        Number of molecules (default: 54)
    segid : str, optional
        Segment ID (default: "SYS")
    write_psf_path : str, optional
        Output PSF path (default: "system.psf")
    write_pdb_path : str, optional
        Output PDB path (default: None)
    """
    
    print(f"Building PSF with pyCHARMM...")
    print(f"  PDB (sequence): {pdb_path}")
    print(f"  CRD (coords):   {crd_path}")
    print(f"  Topology:       {toppar_stream}")
    print(f"  PSF output:     {write_psf_path}")
    
    # Save current directory and change to directory containing files
    # This avoids path case-sensitivity issues in pyCHARMM
    original_dir = os.getcwd()
    work_dir = os.path.dirname(os.path.abspath(pdb_path))
    os.chdir(work_dir)
    
    # Use relative filenames
    pdb_filename = os.path.basename(pdb_path)
    crd_filename = os.path.basename(crd_path)
    toppar_filename = os.path.basename(toppar_stream)
    psf_filename = os.path.basename(write_psf_path)
    if write_pdb_path:
        pdb_out_filename = os.path.basename(write_pdb_path)
    
    try:
        # Read molecule-specific topology/parameters
        with open(toppar_filename, 'r') as f:
            content = f.read()
        
        # Check if base CHARMM36 CGenFF files exist
        base_rtf = "/home/takahashi/charmm/charmm/toppar/top_all36_cgenff.rtf"
        base_prm = "/home/takahashi/charmm/charmm/toppar/par_all36_cgenff.prm"
        
        if os.path.exists(base_rtf) and os.path.exists(base_prm):
            print(f"Loading base CHARMM36 CGenFF force field...")
            # Load base CHARMM36 CGenFF force field
            # Set BOMLEV to -2 to ignore NBFIX warnings (expected when reading CGenFF without protein/nucleic acid FF)
            lingo.charmm_script(f"""
BOMLEV -2

open unit 10 read card name {base_rtf}
read rtf card unit 10
close unit 10

open unit 10 read card name {base_prm}
read param card unit 10
close unit 10

BOMLEV 0
""")
        else:
            print("Warning: Base CHARMM36 CGenFF files not found")
            print("  RTF:", base_rtf)
            print("  PRM:", base_prm)
            print("Attempting to use minimal parameters (may be incomplete)")
            
            # Fallback: Extract atom types and create minimal base RTF/PRM
            atom_types = set()
            for match in re.finditer(r'ATOM\s+\S+\s+(\S+)', content):
                atom_types.add(match.group(1))
            
            minimal_masses = """* Minimal mass definitions
*

read rtf card
* Minimal RTF
*
36 1

AUTO ANGLES DIHE

"""
            for atype in sorted(atom_types):
                if atype.startswith('C'):
                    minimal_masses += f"MASS -1 {atype:8s} 12.01100  ! carbon\n"
                elif atype.startswith('H'):
                    minimal_masses += f"MASS -1 {atype:8s}  1.00800  ! hydrogen\n"
                elif atype.startswith('N'):
                    minimal_masses += f"MASS -1 {atype:8s} 14.00700  ! nitrogen\n"
                elif atype.startswith('O'):
                    minimal_masses += f"MASS -1 {atype:8s} 15.99900  ! oxygen\n"
                else:
                    minimal_masses += f"MASS -1 {atype:8s} 12.01100  ! other\n"
            
            minimal_masses += "\nEND\n\nread param card\n* Minimal parameters\n*\n\nEND\n"
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.str', delete=False) as tmp:
                tmp_rtf_path = tmp.name
                tmp.write(minimal_masses)
            
            lingo.charmm_script(f"stream {tmp_rtf_path}")
            os.unlink(tmp_rtf_path)
        
        # Modify CGenFF stream file (remove 'flex' option)
        content_modified = content.replace('read param card flex append', 'read param card append')
        content_modified = content_modified.replace('READ PARAM CARD FLEX APPEND', 'READ PARAM CARD APPEND')
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.str', delete=False) as tmp:
            tmp_str_path = tmp.name
            tmp.write(content_modified)
        
        # Execute CHARMM workflow
        # 1. Load molecule-specific topology/parameters (append to base)
        lingo.charmm_script(f"stream {tmp_str_path}")
        
        # 2. Set extended I/O format
        lingo.charmm_script("IOFORMAT EXTENDED")
        
        # 3. Define sequence explicitly (nmol copies of resname)
        # Split into multiple lines (10 residues per line) for readability
        print(f"Defining sequence: {nmol} x {resname}...")
        seq_lines = []
        for i in range(0, nmol, 10):
            chunk = [resname] * min(10, nmol - i)
            seq_lines.append(" ".join(chunk))
        seq_text = "\n".join(seq_lines)
        
        lingo.charmm_script(f"""
read sequence card
* Molecular crystal sequence
*
{nmol}
{seq_text}

""")
        
        # 4. Generate PSF structure
        print(f"Generating PSF structure...")
        lingo.charmm_script(f"generate {segid} setup warn")
        
        # 5. Autogenerate angles, dihedrals, and impropers
        lingo.charmm_script("autogenerate angles dihedrals impropers")
        
        # 6. Read coordinates from CRD (high precision)
        print(f"Reading coordinates from {crd_filename}...")
        lingo.charmm_script(f"""
open unit 11 read card name {crd_filename}
read coor card unit 11
close unit 11
""")
        
        # 7. Write PSF
        lingo.charmm_script(f"""
open unit 20 write card name {psf_filename}
write psf card unit 20
* 
""")
        print(f"Successfully created {psf_filename}")
        
        # 8. Optionally write PDB
        if write_pdb_path:
            lingo.charmm_script(f"""
open unit 21 write card name {pdb_out_filename}
write coor pdb unit 21
close unit 21
""")
            print(f"Successfully wrote {pdb_out_filename}")
        
        # Clean up temporary stream file
        if os.path.exists(tmp_str_path):
            os.unlink(tmp_str_path)
    
    finally:
        # Return to original directory
        os.chdir(original_dir)


def main():
    """Command-line interface for building PSF files.
    
    See docs/build_psf.md for detailed documentation and examples.
    """
    parser = argparse.ArgumentParser(
        description="Build CHARMM PSF from PDB sequence + CRD coordinates (pyCHARMM)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="See docs/build_psf.md for detailed documentation and examples."
    )
    
    parser.add_argument("--pdb", required=True, help="Input PDB file (for reference)")
    parser.add_argument("--crd", required=True, help="Input CRD file (for coordinates)")
    parser.add_argument("--toppar", required=True, help="CHARMM topology/parameter stream file")
    parser.add_argument("--resname", default="MOL", help="Residue name (default: MOL)")
    parser.add_argument("--nmol", type=int, required=True, help="Number of molecules in crystal")
    parser.add_argument("--segid", default="SYS", help="Segment ID (default: SYS)")
    parser.add_argument("--psf", default="system.psf", help="Output PSF file (default: system.psf)")
    parser.add_argument("--out-pdb", default=None, help="Output PDB file (optional)")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.pdb):
        raise FileNotFoundError(f"PDB file not found: {args.pdb}")
    if not os.path.exists(args.crd):
        raise FileNotFoundError(f"CRD file not found: {args.crd}")
    if not os.path.exists(args.toppar):
        raise FileNotFoundError(f"Toppar file not found: {args.toppar}")
    
    # Build PSF
    build_psf(pdb_path=args.pdb, crd_path=args.crd, toppar_stream=args.toppar, resname=args.resname,
        nmol=args.nmol, segid=args.segid, write_psf_path=args.psf, write_pdb_path=args.out_pdb)


if __name__ == "__main__":
    main()
