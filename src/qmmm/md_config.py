"""
MD Configuration Module

Load MD settings from YAML configuration file and provide methods
for running molecular dynamics simulations.
"""

import os
import yaml
from pathlib import Path


class MDConfig:
    """Configuration manager for MD simulations."""
    
    def __init__(self, config_file=None):
        """
        Initialize configuration from YAML file.
        
        Args:
            config_file: Path to YAML config file. 
                        If None, uses 'md_settings.yaml' in the same directory.
        """
        if config_file is None:
            config_file = Path(__file__).parent / "md_settings.yaml"
        
        self.config_file = Path(config_file)
        self._load_config()
    
    def _load_config(self):
        """Load configuration from YAML file."""
        with open(self.config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # MD settings
        md_settings = config.get('md', {})
        self.timestep = md_settings.get('timestep', 0.001)  # ps
        self.traj_frequency = md_settings.get('traj_frequency', 50)
        
        # Integrators for each ensemble (with defaults)
        integrators = md_settings.get('integrator', {})
        if isinstance(integrators, str):
            # Backward compatibility: single integrator for all
            self.nvt_integrator = integrators
            self.npt_integrator = integrators
            self.nve_integrator = 'VelocityVerletIntegrator'
        else:
            self.nvt_integrator = integrators.get('nvt', 'LangevinMiddleIntegrator')
            self.npt_integrator = integrators.get('npt', 'LangevinMiddleIntegrator')
            self.nve_integrator = integrators.get('nve', 'VelocityVerletIntegrator')
        
        # NVT settings (temperature, coupling_frequency)
        nvt_settings = md_settings.get('nvt', {})
        if 'temperature' not in nvt_settings:
            raise ValueError("Temperature must be specified in md.nvt.temperature")
        self.temperature = nvt_settings['temperature']  # K
        self.coupling_frequency = nvt_settings.get('coupling_frequency', 1)  # 1/ps
        
        # NPT settings (barostat, pressure, barostat_frequency)
        npt_settings = md_settings.get('npt', {})
        self.barostat = npt_settings.get('barostat', 'MonteCarloBarostat')
        self.pressure = npt_settings.get('pressure', 1)  # bar
        self.barostat_frequency = npt_settings.get('barostat_frequency', 25)  # steps
        
        # MD simulation times for each phase
        md_times = md_settings.get('simulation_time', {})
        self.nvt_time = md_times.get('nvt', 100.0)  # ps
        self.npt_time = md_times.get('npt', 100.0)  # ps
        self.nve_time = md_times.get('nve', 100.0)  # ps
        
        # Output settings
        output = config.get('output', {})
        self.output_dir = output.get('directory', 'output')
        self.save_gro = output.get('save_gro', True)
    
    def print_config(self):
        """Print MD configuration summary."""
        print(f"\nMD Configuration ({self.config_file}):")
        print(f"  Timestep:           {self.timestep} ps ({self.timestep * 1000} fs)")
        print(f"  Traj frequency:     {self.traj_frequency}")
        print(f"  NVT:")
        print(f"    Temperature:      {self.temperature} K")
        print(f"    Coupling freq:    {self.coupling_frequency} /ps")
        print(f"    Integrator:       {self.nvt_integrator}")
        print(f"  NPT:")
        print(f"    Pressure:         {self.pressure} bar")
        print(f"    Barostat:         {self.barostat}")
        print(f"    Barostat freq:    {self.barostat_frequency}")
        print(f"    Integrator:       {self.npt_integrator}")
        print(f"  NVE:")
        print(f"    Integrator:       {self.nve_integrator}")
        print(f"  NVT time:        {self.nvt_time} ps")
        print(f"  NPT time:        {self.npt_time} ps")
        print(f"  NVE time:        {self.nve_time} ps")
    
    def run_nvt(self, frag, theory, output_dir=None, simulation_time=None):
        """
        Run NVT equilibration.
        
        Args:
            frag: ASH Fragment object.
            theory: Theory object (QMMMTheory or OpenMMTheory).
            output_dir: Output directory path. If None, uses config value.
            simulation_time: Simulation time in ps. If None, uses config value.
        
        Returns:
            None
        """
        from ash import OpenMM_MD
        
        if output_dir is None:
            output_dir = self.output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        if simulation_time is None:
            simulation_time = self.nvt_time
        
        OpenMM_MD(
            fragment=frag,
            theory=theory,
            timestep=self.timestep,
            simulation_time=simulation_time,
            temperature=self.temperature,
            coupling_frequency=self.coupling_frequency,
            integrator=self.nvt_integrator,
            traj_frequency=self.traj_frequency,
            trajfilename=f"{output_dir}/nvt",
            datafilename=f"{output_dir}/nvt.txt"
        )
    
    def run_npt(self, frag, theory, output_dir=None, simulation_time=None):
        """
        Run NPT equilibration.
        
        Args:
            frag: ASH Fragment object.
            theory: Theory object (QMMMTheory or OpenMMTheory).
            output_dir: Output directory path. If None, uses config value.
            simulation_time: Simulation time in ps. If None, uses config value.
        
        Returns:
            None
        """
        from ash import OpenMM_MD
        
        if output_dir is None:
            output_dir = self.output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        if simulation_time is None:
            simulation_time = self.npt_time
        
        OpenMM_MD(
            fragment=frag,
            theory=theory,
            timestep=self.timestep,
            simulation_time=simulation_time,
            temperature=self.temperature,
            barostat=self.barostat,
            pressure=self.pressure,
            barostat_frequency=self.barostat_frequency,
            integrator=self.npt_integrator,
            traj_frequency=self.traj_frequency,
            trajfilename=f"{output_dir}/npt",
            datafilename=f"{output_dir}/npt.txt"
        )
    
    def run_nve(self, frag, theory, output_dir=None, simulation_time=None):
        """
        Run NVE production.
        
        Args:
            frag: ASH Fragment object.
            theory: Theory object (QMMMTheory or OpenMMTheory).
            output_dir: Output directory path. If None, uses config value.
            simulation_time: Simulation time in ps. If None, uses config value.
        
        Returns:
            None
        """
        from ash import OpenMM_MD
        
        if output_dir is None:
            output_dir = self.output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        if simulation_time is None:
            simulation_time = self.nve_time
        
        OpenMM_MD(
            fragment=frag,
            theory=theory,
            timestep=self.timestep,
            simulation_time=simulation_time,
            integrator=self.nve_integrator,
            traj_frequency=self.traj_frequency,
            trajfilename=f"{output_dir}/nve",
            datafilename=f"{output_dir}/nve.txt"
        )
    
    def save_final_structure(self, output_dir=None, prefix="md"):
        """
        Save final structure in GRO format.
        
        Args:
            output_dir: Output directory path. If None, uses config value.
            prefix: Prefix for output files (nvt, npt, nve, etc.)
        """
        if not self.save_gro:
            return
        
        if output_dir is None:
            output_dir = self.output_dir
        
        try:
            import mdtraj as md
            traj = md.load(f"{output_dir}/{prefix}_lastframe.pdb")
            traj.save_gro(f"{output_dir}/{prefix}_lastframe.gro")
            print(f"Saved: {output_dir}/{prefix}_lastframe.gro")
        except Exception as e:
            print(f"Warning: Could not save GRO file: {e}")


# Default configuration instance
_default_md_config = None

def get_md_config(config_file=None):
    """
    Get MD configuration instance.
    
    Args:
        config_file: Path to YAML config file.
                    If None and first call, uses default path.
    
    Returns:
        MDConfig: Configuration instance.
    """
    global _default_md_config
    
    if config_file is not None:
        return MDConfig(config_file)
    
    if _default_md_config is None:
        _default_md_config = MDConfig()
    
    return _default_md_config
