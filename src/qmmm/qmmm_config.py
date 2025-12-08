"""
QM/MM Configuration Module

Load settings from YAML configuration file and provide factory functions
for creating ASH theory objects.
"""

import os
import yaml
from pathlib import Path


class QMMMConfig:
    """Configuration manager for QM/MM calculations."""
    
    def __init__(self, config_file=None):
        """
        Initialize configuration from YAML file.
        
        Args:
            config_file: Path to YAML config file. 
                        If None, uses 'qmmm_settings.yaml' in the same directory.
        """
        if config_file is None:
            config_file = Path(__file__).parent / "qmmm_settings.yaml"
        
        self.config_file = Path(config_file)
        self._load_config()
    
    def _load_config(self):
        """Load configuration from YAML file."""
        with open(self.config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # File paths
        paths = config.get('paths', {})
        self.sk_dir = paths.get('sk_dir', '')
        # AMBER topology files
        self.amber_prmtop = paths.get('amber_prmtop', '')
        self.amber_inpcrd = paths.get('amber_inpcrd', '')
        # PDB file for coordinates (optional, can override inpcrd coordinates)
        self.pdbfile = paths.get('pdbfile', '')
        # Box file for PBC dimensions
        self.boxfile = paths.get('boxfile', None)
        # Atom selection files
        self.qatoms_file = paths.get('qatoms', '')
        self.actatoms_file = paths.get('actatoms', '')
        
        # QM region settings
        qm = config.get('qm', {})
        self.qm_charge = qm.get('charge', 0)
        self.qm_mult = qm.get('mult', 1)
        
        # Parallel settings
        parallel = config.get('parallel', {})
        self.numcores_qm = parallel.get('numcores_qm', 1)
        self.numcores_mm = parallel.get('numcores_mm', 1)
        
        # DFTB settings
        dftb = config.get('dftb', {})
        sk_files = dftb.get('slater_koster_files', {})
        self.slater_koster_files = {
            k: v.format(sk_dir=self.sk_dir) for k, v in sk_files.items()
        }
        self.hubbard_derivs = dftb.get('hubbard_derivs', {})
        self.hcorrection_zeta = dftb.get('hcorrection_zeta', 4.0)
        
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
        import numpy as np
        
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
    
    def load_qmatoms(self, filename=None):
        """
        Load QM atom indices from file.
        
        Args:
            filename: Path to file containing atom indices.
                     If None, uses the path from config.
        
        Returns:
            list: List of atom indices (0-indexed).
        """
        if filename is None:
            filename = self.qatoms_file
        
        with open(filename, 'r') as f:
            return [int(x) for x in f.read().strip().split()]
    
    def load_active_atoms(self, filename=None):
        """
        Load active atom indices from file.
        
        Args:
            filename: Path to file containing atom indices.
                     If None, uses the path from config.
        
        Returns:
            list: List of atom indices (0-indexed).
        """
        if filename is None:
            filename = self.actatoms_file
        
        with open(filename, 'r') as f:
            return [int(x) for x in f.read().strip().split()]
    
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
            numcores=self.numcores_mm
        )
    
    def create_dftb_theory(self):
        """
        Create DFTBTheory object with DFTB3/3ob settings.
        
        Returns:
            DFTBTheory: Configured DFTB theory object.
        """
        from ash import DFTBTheory
        
        return DFTBTheory(
            hamiltonian="DFTB",
            SCC=True,
            ThirdOrderFull=True,
            slaterkoster_dict=self.slater_koster_files,
            hubbard_derivs_dict=self.hubbard_derivs,
            hcorrection_zeta=self.hcorrection_zeta,
            numcores=self.numcores_qm,
            printlevel=2
        )
    
    def create_qmmm_theory(self, frag, qmatoms, omm=None, qm_dftb=None):
        """
        Create QMMMTheory object.
        
        Args:
            frag: ASH Fragment object.
            qmatoms: List of QM atom indices.
            omm: OpenMMTheory object. If None, creates a new one.
            qm_dftb: DFTBTheory object. If None, creates a new one.
        
        Returns:
            QMMMTheory: Configured QM/MM theory object.
        """
        from ash import QMMMTheory
        
        if omm is None:
            omm = self.create_openmm_theory()
        if qm_dftb is None:
            qm_dftb = self.create_dftb_theory()
        
        return QMMMTheory(
            qm_theory=qm_dftb,
            mm_theory=omm,
            fragment=frag,
            qmatoms=qmatoms,
            embedding="elstat",
            qm_charge=self.qm_charge,
            qm_mult=self.qm_mult,
            printlevel=2
        )
    
    def print_system_info(self, frag, qmatoms):
        """
        Print system information summary.
        
        Args:
            frag: ASH Fragment object.
            qmatoms: List of QM atom indices.
        """
        mm_atoms = len([i for i in range(frag.numatoms) if i not in qmatoms])
        print(f"\nSystem Information:")
        print(f"  Total atoms:     {frag.numatoms}")
        print(f"  QM atoms:        {len(qmatoms)}")
        print(f"  MM atoms:        {mm_atoms}")
        print(f"  QM charge/mult:  {self.qm_charge}/{self.qm_mult}")
        print(f"  DFTB method:     DFTB3/3ob-3-1")
        print(f"  MM cutoff:       {self.periodic_nonbonded_cutoff} Å")
        print(f"  H mass:          {self.hydrogenmass} Da")
    
    def print_config(self):
        """Print current configuration summary."""
        print(f"\nQM/MM Configuration ({self.config_file}):")
        print(f"  SK directory:    {self.sk_dir}")
        print(f"  AMBER prmtop:    {self.amber_prmtop}")
        print(f"  AMBER inpcrd:    {self.amber_inpcrd}")
        if self.pdbfile:
            print(f"  PDB file:        {self.pdbfile}")
        print(f"  QM atoms file:   {self.qatoms_file}")
        print(f"  QM charge/mult:  {self.qm_charge}/{self.qm_mult}")
        print(f"  QM cores:        {self.numcores_qm}")
        print(f"  MM cores:        {self.numcores_mm}")
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
        QMMMConfig: Configuration instance.
    """
    global _default_config
    
    if config_file is not None:
        return QMMMConfig(config_file)
    
    if _default_config is None:
        _default_config = QMMMConfig()
    
    return _default_config


# Convenience functions using default config
def load_qmatoms(filename=None):
    """Load QM atom indices from file."""
    return get_config().load_qmatoms(filename)

def load_active_atoms(filename=None):
    """Load active atom indices from file."""
    return get_config().load_active_atoms(filename)

def create_openmm_theory(prmtopfile=None, inpcrdfile=None, pbc_vectors=None):
    """Create OpenMMTheory object using AMBER topology."""
    return get_config().create_openmm_theory(prmtopfile, inpcrdfile, pbc_vectors)

def create_dftb_theory():
    """Create DFTBTheory object."""
    return get_config().create_dftb_theory()

def create_qmmm_theory(frag, qmatoms, omm=None, qm_dftb=None):
    """Create QMMMTheory object."""
    return get_config().create_qmmm_theory(frag, qmatoms, omm, qm_dftb)

def print_system_info(frag, qmatoms):
    """Print system information."""
    return get_config().print_system_info(frag, qmatoms)


if __name__ == "__main__":
    # Test configuration loading
    config = get_config()
    config.print_config()
