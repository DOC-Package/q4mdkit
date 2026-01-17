"""
MM (Force Field) Configuration Module

Load settings from YAML configuration file and provide factory functions
for creating ASH OpenMM theory objects for classical MD simulations.

This module is a simplified version of qmmm_config that uses only
force field (OpenMM) for MD, without any QM region.
"""

import os
import yaml
from pathlib import Path
import numpy as np


class MMConfig:
    """Configuration manager for MM (Force Field) calculations."""
    
    def __init__(self, config_file=None):
        """
        Initialize configuration from YAML file.
        
        Args:
            config_file: Path to YAML config file. 
                        If None, uses 'mm_settings.yaml' in the same directory.
        """
        if config_file is None:
            config_file = Path(__file__).parent / "mm_settings.yaml"
        
        self.config_file = Path(config_file)
        self._load_config()
    
    def _load_config(self):
        """Load configuration from YAML file."""
        with open(self.config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # File paths
        paths = config.get('paths', {})
        # AMBER topology files
        self.amber_prmtop = paths.get('amber_prmtop', '')
        self.amber_inpcrd = paths.get('amber_inpcrd', '')
        # PDB file for coordinates (optional, can override inpcrd coordinates)
        self.pdbfile = paths.get('pdbfile', '')
        # Box file for PBC dimensions
        self.boxfile = paths.get('boxfile', None)
        
        # Parallel settings
        parallel = config.get('parallel', {})
        self.numcores = parallel.get('numcores', 1)
        
        # OpenMM settings
        openmm = config.get('openmm', {})
        self.periodic = openmm.get('periodic', True)
        self.periodic_nonbonded_cutoff = openmm.get('periodic_nonbonded_cutoff', 9.0)
        self.periodic_cell_dimensions = openmm.get('periodic_cell_dimensions', None)
        self.use_boxfile = openmm.get('use_boxfile', False)
        self.autoconstraints = openmm.get('autoconstraints', None)
        self.rigidwater = openmm.get('rigidwater', False)
        self.hydrogenmass = openmm.get('hydrogenmass', 1.5)
        self.platform = openmm.get('platform', 'CPU')
        
        # Validate: both periodic_cell_dimensions and use_boxfile cannot be set
        if self.periodic_cell_dimensions is not None and self.use_boxfile:
            raise ValueError(
                "Cannot specify both 'periodic_cell_dimensions' and 'use_boxfile: true'. "
                "Please use only one method to define the PBC box."
            )
        
        # Load box vectors from file if use_boxfile is True
        if self.use_boxfile:
            if not self.boxfile:
                raise ValueError(
                    "'use_boxfile' is true but 'boxfile' path is not set in paths section."
                )
            self.pbc_vectors = self._load_boxfile(self.boxfile)
        else:
            self.pbc_vectors = None
    
    def _load_boxfile(self, boxfile):
        """
        Load PBC box vectors from file.
        
        Args:
            boxfile: Path to box file (format: 3x3 matrix, each row is a vector)
        
        Returns:
            list: [[ax, ay, az], [bx, by, bz], [cx, cy, cz]]
        """
        vectors = []
        with open(boxfile, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    values = [float(x) for x in line.split()]
                    if len(values) == 3:
                        vectors.append(values)
                    elif len(values) == 6:
                        # Old format: a b c alpha beta gamma - convert to vectors
                        a, b, c, alpha, beta, gamma = values
                        return self._cell_to_box_vectors(a, b, c, alpha, beta, gamma)
        if len(vectors) == 3:
            return vectors
        return None
    
    @staticmethod
    def _cell_to_box_vectors(a, b, c, alpha, beta, gamma):
        """
        Convert cell parameters to box vectors for triclinic cells.
        
        Args:
            a, b, c: Cell lengths in Angstrom
            alpha, beta, gamma: Cell angles in degrees
        
        Returns:
            list: 3x3 box vectors [[ax, ay, az], [bx, by, bz], [cx, cy, cz]]
        """
        # Convert angles to radians
        alpha_rad = np.radians(alpha)
        beta_rad = np.radians(beta)
        gamma_rad = np.radians(gamma)
        
        # Box vector a along x-axis
        ax, ay, az = a, 0.0, 0.0
        
        # Box vector b in xy-plane
        bx = b * np.cos(gamma_rad)
        by = b * np.sin(gamma_rad)
        bz = 0.0
        
        # Box vector c
        cx = c * np.cos(beta_rad)
        cy = c * (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad)
        cz = np.sqrt(c**2 - cx**2 - cy**2)
        
        return [[ax, ay, az], [bx, by, bz], [cx, cy, cz]]
    
    def create_openmm_theory(self, prmtopfile=None, inpcrdfile=None, pbc_vectors=None):
        """
        Create OpenMMTheory object using AMBER topology files.
        
        Args:
            prmtopfile: Path to AMBER prmtop file. If None, uses the path from config.
            inpcrdfile: Path to AMBER inpcrd file. If None, uses the path from config.
            pbc_vectors: [[ax,ay,az], [bx,by,bz], [cx,cy,cz]] in Angstrom.
                         If None, uses the value from config.
        
        Returns:
            OpenMMTheory: Configured OpenMM theory object.
        """
        from ash import OpenMMTheory
        
        if prmtopfile is None:
            prmtopfile = self.amber_prmtop
        if inpcrdfile is None:
            inpcrdfile = self.amber_inpcrd
        if pbc_vectors is None:
            pbc_vectors = self.pbc_vectors
        
        # Print PBC vectors info
        if pbc_vectors is not None:
            print("Using PBC box vectors:")
            for i, vec in enumerate(pbc_vectors):
                print(f"  v{i+1}: [{vec[0]:.4f}, {vec[1]:.4f}, {vec[2]:.4f}]")
        
        return OpenMMTheory(
            Amberfiles=True,
            amberprmtopfile=prmtopfile,
            periodic=self.periodic,
            periodic_nonbonded_cutoff=self.periodic_nonbonded_cutoff,
            PBCvectors=pbc_vectors,
            autoconstraints=self.autoconstraints,
            rigidwater=self.rigidwater,
            hydrogenmass=self.hydrogenmass,
            platform=self.platform,
            numcores=self.numcores
        )
    
    def print_system_info(self, frag):
        """
        Print system information summary.
        
        Args:
            frag: ASH Fragment object.
        """
        print(f"\nSystem Information:")
        print(f"  Total atoms:     {frag.numatoms}")
        print(f"  MM cutoff:       {self.periodic_nonbonded_cutoff} Å")
        print(f"  H mass:          {self.hydrogenmass} Da")
        print(f"  Platform:        {self.platform}")
    
    def print_config(self):
        """Print current configuration summary."""
        print(f"\nMM Configuration ({self.config_file}):")
        print(f"  AMBER prmtop:    {self.amber_prmtop}")
        print(f"  AMBER inpcrd:    {self.amber_inpcrd}")
        if self.pdbfile:
            print(f"  PDB file:        {self.pdbfile}")
        print(f"  Cores:           {self.numcores}")
        print(f"  Platform:        {self.platform}")
        print(f"  Periodic:        {self.periodic}")
        if self.pbc_vectors:
            print(f"  PBC Box vectors:")
            for i, vec in enumerate(self.pbc_vectors):
                print(f"    v{i+1}: [{vec[0]:.4f}, {vec[1]:.4f}, {vec[2]:.4f}]")


# Default configuration instance
_default_config = None

def get_config(config_file=None):
    """
    Get configuration instance.
    
    Args:
        config_file: Path to YAML config file.
                    If None and first call, uses default path.
    
    Returns:
        MMConfig: Configuration instance.
    """
    global _default_config
    
    if config_file is not None:
        return MMConfig(config_file)
    
    if _default_config is None:
        _default_config = MMConfig()
    
    return _default_config


# Convenience functions using default config
def create_openmm_theory(prmtopfile=None, inpcrdfile=None, pbc_vectors=None):
    """Create OpenMMTheory object using AMBER topology."""
    return get_config().create_openmm_theory(prmtopfile, inpcrdfile, pbc_vectors)

def print_system_info(frag):
    """Print system information."""
    return get_config().print_system_info(frag)


if __name__ == "__main__":
    # Test configuration loading
    config = get_config()
    config.print_config()
