"""Extract box vectors from PDB CRYST1 record."""

import numpy as np
from pathlib import Path


def cryst1_to_vectors(a, b, c, alpha, beta, gamma):
    """Convert cell parameters to 3x3 box vectors (row vectors, Angstrom).
    
    Standard convention: a along x, b in xy plane, c general.
    
    Args:
        a, b, c: Cell lengths in Angstrom
        alpha, beta, gamma: Cell angles in degrees
    
    Returns:
        numpy.ndarray: 3x3 matrix of box vectors (row vectors)
    """
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)
    
    cos_alpha = np.cos(alpha_rad)
    cos_beta = np.cos(beta_rad)
    cos_gamma = np.cos(gamma_rad)
    sin_gamma = np.sin(gamma_rad)
    
    # Vector a along x
    ax, ay, az = a, 0.0, 0.0
    
    # Vector b in xy plane
    bx = b * cos_gamma
    by = b * sin_gamma
    bz = 0.0
    
    # Vector c general
    cx = c * cos_beta
    cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cz = np.sqrt(c**2 - cx**2 - cy**2)
    
    return np.array([
        [ax, ay, az],
        [bx, by, bz],
        [cx, cy, cz]
    ])


def parse_cryst1(pdb_file):
    """Parse CRYST1 record from PDB file.
    
    Args:
        pdb_file: Path to PDB file
    
    Returns:
        tuple: (a, b, c, alpha, beta, gamma)
    
    Raises:
        ValueError: If no CRYST1 record found
    """
    with open(pdb_file) as f:
        for line in f:
            if line.startswith('CRYST1'):
                a = float(line[6:15])
                b = float(line[15:24])
                c = float(line[24:33])
                alpha = float(line[33:40])
                beta = float(line[40:47])
                gamma = float(line[47:54])
                return a, b, c, alpha, beta, gamma
    raise ValueError(f"No CRYST1 record found in {pdb_file}")


def write_box(box_file, vectors, a, b, c, alpha, beta, gamma):
    """Write box vectors to file.
    
    Args:
        box_file: Output file path
        vectors: 3x3 numpy array of box vectors
        a, b, c: Cell lengths in Angstrom
        alpha, beta, gamma: Cell angles in degrees
    """
    with open(box_file, 'w') as f:
        f.write("# PBC Box vectors (3x3 matrix, each row is a vector in Angstrom)\n")
        f.write(f"# Cell parameters: a={a:.6f} b={b:.6f} c={c:.6f} "
                f"alpha={alpha:.4f} beta={beta:.4f} gamma={gamma:.4f}\n")
        for row in vectors:
            f.write(f"{row[0]:.10f} {row[1]:.10f} {row[2]:.10f}\n")


def pdb_to_box(pdb_file, output=None):
    """Extract box vectors from PDB and write to .box file.
    
    Args:
        pdb_file: Input PDB file path
        output: Output .box file path (default: same name with .box extension)
    
    Returns:
        dict: Box information with keys 'a', 'b', 'c', 'alpha', 'beta', 'gamma',
              'vectors', 'output'
    """
    pdb_path = Path(pdb_file)
    output = output or str(pdb_path.with_suffix('.box'))
    
    a, b, c, alpha, beta, gamma = parse_cryst1(pdb_file)
    vectors = cryst1_to_vectors(a, b, c, alpha, beta, gamma)
    
    write_box(output, vectors, a, b, c, alpha, beta, gamma)
    
    return {
        'a': a, 'b': b, 'c': c,
        'alpha': alpha, 'beta': beta, 'gamma': gamma,
        'vectors': vectors,
        'output': output
    }


def read_box(box_file):
    """Read box vectors from .box file.
    
    Args:
        box_file: Path to .box file
    
    Returns:
        numpy.ndarray: 3x3 matrix of box vectors (row vectors)
    """
    vectors = []
    with open(box_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) == 3:
                vectors.append([float(x) for x in parts])
    return np.array(vectors)
