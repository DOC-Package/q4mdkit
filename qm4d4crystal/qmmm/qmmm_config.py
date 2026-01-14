"""
QM/MM Configuration Module

Load settings from YAML configuration file and provide factory functions
for creating ASH theory objects.
"""

import os
import re
import shutil
import yaml
from pathlib import Path


class DFTBTheory_LogSCC:
    """
    DFTBTheory wrapper that logs SCC convergence information.
    
    Extracts SCC iteration count, error, and optionally total energy from 
    detailed.out after each gradient calculation and appends to log files.
    """
    
    def __init__(self, *args, scc_logfile="scc_error.dat",
                 keep_detailed=False, output_dir="output",
                 energy_log=False, energy_logfile="qm_energy.dat", **kwargs):
        """
        Initialize DFTBTheory with SCC logging.
        
        Args:
            scc_logfile: Filename for SCC error log (saved in output_dir).
            keep_detailed: If True, save detailed.out for each step.
            output_dir: Directory to save log files.
            energy_log: If True, also log total energy from detailed.out.
            energy_logfile: Filename for energy log (saved in output_dir).
            *args, **kwargs: Passed to DFTBTheory.__init__
        """
        from ash import DFTBTheory
        self._dftb = DFTBTheory(*args, **kwargs)
        self._callidx = 0
        # Use absolute path to avoid issues when working directory changes
        self._output_dir = Path(output_dir).resolve()
        self._output_dir.mkdir(parents=True, exist_ok=True)
        self._log = self._output_dir / scc_logfile
        self._keep_detailed = keep_detailed
        self._energy_log = energy_log
        self._energy_logfile = self._output_dir / energy_logfile
        # Overwrite log files at start of each run
        self._log.write_text("# call_index  iSCC  SCC_error(a.u.)\n")
        if self._energy_log:
            self._energy_logfile.write_text("# call_index  Electronic(H)         Repulsive(H)          Total(H)\n")
    
    def __getattr__(self, name):
        """Delegate attribute access to wrapped DFTBTheory."""
        return getattr(self._dftb, name)
    
    def _read_scc_from_detailed(self, detailed="detailed.out"):
        """
        Read SCC iteration info from DFTB+ detailed.out file.
        
        Args:
            detailed: Path to detailed.out file.
        
        Returns:
            tuple: (iSCC, error) or None if not found.
        """
        p = Path(detailed)
        if not p.exists():
            return None
        txt = p.read_text(errors="ignore")
        
        # Parse SCC convergence lines from detailed.out
        # Format: " iSCC Total electronic   Diff electronic      SCC error"
        # Example: "   12   -0.87075037E+02   -0.61469052E-09    0.74172173E-05"
        pat = re.compile(
            r"^\s*(\d+)\s+([-+]?\d+\.\d+E[+-]\d+)\s+([-+]?\d+\.\d+E[+-]\d+)\s+([-+]?\d+\.\d+E[+-]\d+)\s*$",
            re.M
        )
        matches = pat.findall(txt)
        if not matches:
            return None
        iSCC, _, _, err = matches[-1]
        err = float(err)
        return int(iSCC), err
    
    def _read_energy_from_detailed(self, detailed="detailed.out"):
        """
        Read energies from DFTB+ detailed.out file.
        
        Args:
            detailed: Path to detailed.out file.
        
        Returns:
            tuple: (electronic_energy, repulsive_energy, total_energy) in Hartree,
                   or None if not found.
        """
        p = Path(detailed)
        if not p.exists():
            return None
        txt = p.read_text(errors="ignore")
        
        # Parse energy lines from detailed.out (Hartree only)
        # Format: "Total Electronic energy:           -87.4133454695 H        -2378.6382 eV"
        # Format: "Repulsive energy:                    2.0023495898 H           54.4867 eV"
        # Format: "Total energy:                      -85.4109958798 H        -2324.1514 eV"
        
        pat_elec = re.compile(r"Total Electronic energy:\s+([-+]?\d+\.\d+)\s+H")
        pat_rep = re.compile(r"Repulsive energy:\s+([-+]?\d+\.\d+)\s+H")
        pat_total = re.compile(r"Total energy:\s+([-+]?\d+\.\d+)\s+H")
        
        m_elec = pat_elec.search(txt)
        m_rep = pat_rep.search(txt)
        m_total = pat_total.search(txt)
        
        if not (m_elec and m_rep and m_total):
            return None
        
        elec = float(m_elec.group(1))
        rep = float(m_rep.group(1))
        total = float(m_total.group(1))
        return elec, rep, total
    
    def run(self, *args, **kwargs):
        """
        Run DFTB calculation and log SCC convergence if Grad=True.
        
        Args:
            *args, **kwargs: Passed to DFTBTheory.run()
        
        Returns:
            Result from DFTBTheory.run()
        """
        Grad = kwargs.get("Grad", False)
        res = self._dftb.run(*args, **kwargs)
        
        if Grad:
            # Log SCC error
            info = self._read_scc_from_detailed("detailed.out")
            if info is not None:
                iSCC, err = info
                with self._log.open("a") as f:
                    f.write(f"{self._callidx:8d}  {iSCC:4d}  {err:.12e}\n")
            
            # Log energy if enabled
            if self._energy_log:
                energy_info = self._read_energy_from_detailed("detailed.out")
                if energy_info is not None:
                    elec, rep, total = energy_info
                    with self._energy_logfile.open("a") as f:
                        f.write(f"{self._callidx:8d}  {elec:20.10f}  {rep:20.10f}  {total:20.10f}\n")
            
            # Keep detailed.out if requested
            if self._keep_detailed:
                dest = self._output_dir / f"detailed_{self._callidx:06d}.out"
                shutil.copy("detailed.out", dest)
            
            self._callidx += 1
        
        return res


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
        self.max_scc_iterations = dftb.get('max_scc_iterations', 300)
        self.third_order_full = dftb.get('third_order_full', True)
        # SCC logging settings
        self.scc_log = dftb.get('scc_log', False)
        self.scc_logfile = dftb.get('scc_logfile', 'scc_error.dat')
        self.keep_detailed = dftb.get('keep_detailed', False)
        self.energy_log = dftb.get('energy_log', False)
        self.energy_logfile = dftb.get('energy_logfile', 'qm_energy.dat')
        
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
        
        If scc_log is True, returns DFTBTheory_LogSCC which logs
        SCC convergence information from detailed.out at each step.
        
        Returns:
            DFTBTheory or DFTBTheory_LogSCC: Configured DFTB theory object.
        """
        if self.scc_log:
            return DFTBTheory_LogSCC(
                hamiltonian="DFTB",
                SCC=True,
                ThirdOrderFull=self.third_order_full,
                slaterkoster_dict=self.slater_koster_files,
                hubbard_derivs_dict=self.hubbard_derivs,
                hcorrection_zeta=self.hcorrection_zeta,
                MaxSCCIterations=self.max_scc_iterations,
                numcores=self.numcores_qm,
                printlevel=2,
                scc_logfile=self.scc_logfile,
                keep_detailed=self.keep_detailed,
                output_dir="output",
                energy_log=self.energy_log,
                energy_logfile=self.energy_logfile
            )
        else:
            from ash import DFTBTheory
            return DFTBTheory(
                hamiltonian="DFTB",
                SCC=True,
                ThirdOrderFull=self.third_order_full,
                slaterkoster_dict=self.slater_koster_files,
                hubbard_derivs_dict=self.hubbard_derivs,
                hcorrection_zeta=self.hcorrection_zeta,
                MaxSCCIterations=self.max_scc_iterations,
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
