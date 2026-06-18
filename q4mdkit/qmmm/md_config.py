"""
MD Configuration Module

Load MD settings from YAML configuration file and provide methods
for running molecular dynamics simulations.
"""

import os
import shutil
import xml.etree.ElementTree as ET
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
        self.special_wrapping = md_settings.get('special_wrapping', False)
        self.special_wrapping_updatepos = md_settings.get('special_wrapping_updatepos', False)
        wrapping_atoms_value = md_settings.get('wrapping_atoms', None)
        wrapping_atoms_file = md_settings.get('wrapping_atoms_file', None)
        self.wrapping_atoms = self._resolve_atom_indices(
            wrapping_atoms_value,
            wrapping_atoms_file,
        )
        
        # Integrators for each ensemble (with defaults)
        integrators = md_settings.get('integrator', {})
        if isinstance(integrators, str):
            # Backward compatibility: single integrator for all
            self.nvt_integrator = integrators
            self.npt_integrator = integrators
            self.nve_integrator = 'VerletIntegrator'
        else:
            self.nvt_integrator = integrators.get('nvt', 'LangevinMiddleIntegrator')
            self.npt_integrator = integrators.get('npt', 'LangevinMiddleIntegrator')
            self.nve_integrator = integrators.get('nve', 'VerletIntegrator')
        
        # NVT settings (temperature, coupling_frequency)
        nvt_settings = md_settings.get('nvt', {})
        self.temperature = nvt_settings.get('temperature', None)  # K (required for NVT/NPT)
        self.coupling_frequency = nvt_settings.get('coupling_frequency', 1)  # 1/ps
        
        # NPT settings (barostat, pressure, barostat_frequency)
        npt_settings = md_settings.get('npt', {})
        self.barostat = npt_settings.get('barostat', 'MonteCarloBarostat')
        self.pressure = npt_settings.get('pressure', 1)  # bar
        self.barostat_frequency = npt_settings.get('barostat_frequency', 25)  # steps
        
        # MD simulation times for each phase
        md_times = md_settings.get('simulation_time', {})
        if isinstance(md_times, (int, float)):
            # Single value: use for all ensembles
            self.nvt_time = float(md_times)
            self.npt_time = float(md_times)
            self.nve_time = float(md_times)
        else:
            self.nvt_time = md_times.get('nvt', 100.0)  # ps
            self.npt_time = md_times.get('npt', 100.0)  # ps
            self.nve_time = md_times.get('nve', 100.0)  # ps
        
        # Output settings
        output = config.get('output', {})
        self.output_dir = output.get('directory', 'output')
        self.save_gro = output.get('save_gro', False)
        postprocess_wrapped = md_settings.get('postprocess_wrapped_trajectory', None)
        self.postprocess_wrapped = (
            self.special_wrapping
            if postprocess_wrapped is None
            else postprocess_wrapped
        )
        wrap_final_pdb = md_settings.get('wrap_final_pdb', None)
        self.wrap_final_pdb = (
            self.special_wrapping
            if wrap_final_pdb is None
            else wrap_final_pdb
        )
        wrap_final_state = md_settings.get('wrap_final_state', None)
        self.wrap_final_state = (
            self.special_wrapping
            if wrap_final_state is None
            else wrap_final_state
        )

    def _resolve_atom_indices(self, atom_indices, atom_indices_file=None):
        """Resolve atom indices from a list or from a file path."""
        if atom_indices_file is not None:
            return self._load_atom_indices_file(atom_indices_file)
        if isinstance(atom_indices, str):
            return self._load_atom_indices_file(atom_indices)
        return atom_indices

    def _load_atom_indices_file(self, filename):
        """Load whitespace-separated 0-index atom indices from file."""
        filepath = Path(filename)
        if not filepath.is_absolute():
            filepath = self.config_file.parent / filepath
        with open(filepath, 'r') as f:
            return [int(x) for x in f.read().strip().split()]
    
    def print_config(self):
        """Print MD configuration summary."""
        print(f"\nMD Configuration ({self.config_file}):")
        print(f"  Timestep:           {self.timestep} ps ({self.timestep * 1000} fs)")
        print(f"  Traj frequency:     {self.traj_frequency}")
        print(f"  Special wrapping:   {self.special_wrapping}")
        print(f"  Wrap update pos:    {self.special_wrapping_updatepos}")
        print(f"  Wrapping atoms:     {self.wrapping_atoms}")
        print(f"  Postprocess wrap:   {self.postprocess_wrapped}")
        print(f"  Wrap final PDB:     {self.wrap_final_pdb}")
        print(f"  Wrap final state:   {self.wrap_final_state}")
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
    
    def run_nvt(self, frag, theory, output_dir=None, simulation_time=None, statefile=None):
        """
        Run NVT equilibration.
        
        Args:
            frag: ASH Fragment object.
            theory: Theory object (QMMMTheory or OpenMMTheory).
            output_dir: Output directory path. If None, uses config value.
            simulation_time: Simulation time in ps. If None, uses config value.
            statefile: Path to OpenMM state XML file to load initial state from.
                      If provided, positions and velocities are loaded from this file.
        
        Returns:
            None
        """
        from ash import OpenMM_MD
        
        if self.temperature is None:
            raise ValueError("Temperature must be specified in md.nvt.temperature for NVT simulation")
        
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
            special_wrapping=self.special_wrapping,
            special_wrapping_updatepos=self.special_wrapping_updatepos,
            wrapping_atoms=self.wrapping_atoms,
            trajfilename=f"{output_dir}/nvt",
            datafilename=f"{output_dir}/nvt.csv",
            statefile=statefile
        )
        self.postprocess_wrapped_trajectory(output_dir=output_dir, prefix="nvt")
        self.postprocess_wrapped_pdb(output_dir=output_dir, prefix="nvt")
        self.postprocess_wrapped_state(output_dir=output_dir, prefix="nvt")
    
    def run_npt(self, frag, theory, output_dir=None, simulation_time=None, statefile=None):
        """
        Run NPT equilibration.
        
        Args:
            frag: ASH Fragment object.
            theory: Theory object (QMMMTheory or OpenMMTheory).
            output_dir: Output directory path. If None, uses config value.
            simulation_time: Simulation time in ps. If None, uses config value.
            statefile: Path to OpenMM state XML file to load initial state from.
                      If provided, positions and velocities are loaded from this file.
        
        Returns:
            None
        """
        from ash import OpenMM_MD
        
        if self.temperature is None:
            raise ValueError("Temperature must be specified in md.nvt.temperature for NPT simulation")
        
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
            special_wrapping=self.special_wrapping,
            special_wrapping_updatepos=self.special_wrapping_updatepos,
            wrapping_atoms=self.wrapping_atoms,
            trajfilename=f"{output_dir}/npt",
            datafilename=f"{output_dir}/npt.csv",
            statefile=statefile
        )
        self.postprocess_wrapped_trajectory(output_dir=output_dir, prefix="npt")
        self.postprocess_wrapped_pdb(output_dir=output_dir, prefix="npt")
        self.postprocess_wrapped_state(output_dir=output_dir, prefix="npt")
    
    def run_nve(self, frag, theory, output_dir=None, simulation_time=None, statefile=None):
        """
        Run NVE production.
        
        Args:
            frag: ASH Fragment object.
            theory: Theory object (QMMMTheory or OpenMMTheory).
            output_dir: Output directory path. If None, uses config value.
            simulation_time: Simulation time in ps. If None, uses config value.
            statefile: Path to OpenMM state XML file to load initial state from.
                      If provided, positions and velocities are loaded from this file.
                      Note: For NVE from NPT, use remove_montecarlo_params() first.
        
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
            special_wrapping=self.special_wrapping,
            special_wrapping_updatepos=self.special_wrapping_updatepos,
            wrapping_atoms=self.wrapping_atoms,
            trajfilename=f"{output_dir}/nve",
            datafilename=f"{output_dir}/nve.csv",
            statefile=statefile
        )
        self.postprocess_wrapped_trajectory(output_dir=output_dir, prefix="nve")
        self.postprocess_wrapped_pdb(output_dir=output_dir, prefix="nve")
        self.postprocess_wrapped_state(output_dir=output_dir, prefix="nve")

    def _anchor_molecules_for_mdtraj(self, topology):
        """Build mdtraj anchor_molecules from configured atom indices."""
        if self.wrapping_atoms is None:
            return None
        return [
            set(topology.atom(int(i)) for i in self.wrapping_atoms)
        ]

    def postprocess_wrapped_trajectory(self, output_dir=None, prefix="md",
                                       trajectory_file=None, topology_file=None,
                                       output_file=None, save_pdb_snapshots=False):
        """
        Save a visualization trajectory with the wrapping anchor centered.

        ASH's special_wrapping_updatepos is applied in the step-by-step
        QM/MM-style loops, but pure OpenMM MM runs write DCD frames via
        OpenMM reporters. Re-imaging after the run keeps the raw trajectory
        intact and provides a centered trajectory for visualization.
        """
        if not self.special_wrapping or not self.postprocess_wrapped:
            return None

        if output_dir is None:
            output_dir = self.output_dir
        output_path = Path(output_dir)

        if trajectory_file is None:
            trajectory_file = output_path / f"{prefix}.dcd"
        else:
            trajectory_file = Path(trajectory_file)

        if topology_file is None:
            topology_file = output_path / f"{prefix}_firstframe.pdb"
        else:
            topology_file = Path(topology_file)

        if output_file is None:
            output_file = output_path / f"{prefix}_wrapped.dcd"
        else:
            output_file = Path(output_file)

        if not trajectory_file.exists():
            print(f"Warning: trajectory not found, skipping wrapped postprocess: {trajectory_file}")
            return None
        if not topology_file.exists():
            print(f"Warning: topology not found, skipping wrapped postprocess: {topology_file}")
            return None

        try:
            import mdtraj as md
        except ImportError:
            print("Warning: mdtraj is not available; could not write wrapped trajectory")
            return None

        print(f"Postprocessing wrapped trajectory: {trajectory_file}")
        try:
            traj = md.load(str(trajectory_file), top=str(topology_file))
            anchor_molecules = self._anchor_molecules_for_mdtraj(traj.topology)
            imaged = traj.image_molecules(anchor_molecules=anchor_molecules)
            imaged.save_dcd(str(output_file))

            if save_pdb_snapshots:
                firstframe_file = output_file.with_name(f"{output_file.stem}_firstframe.pdb")
                lastframe_file = output_file.with_name(f"{output_file.stem}_lastframe.pdb")
                imaged[0].save_pdb(str(firstframe_file))
                imaged[-1].save_pdb(str(lastframe_file))
        except Exception as exc:
            print(f"Warning: could not write wrapped trajectory: {exc}")
            return None

        print(f"Saved wrapped trajectory: {output_file}")
        if save_pdb_snapshots:
            print(f"Saved wrapped first/last frames: {firstframe_file}, {lastframe_file}")
        return output_file

    def postprocess_wrapped_pdb(self, output_dir=None, prefix="md",
                                pdb_file=None, backup_file=None,
                                wrapped_copy_file=None):
        """Rewrite the final PDB with the wrapping anchor centered."""
        if not self.special_wrapping or not self.wrap_final_pdb:
            return None

        if output_dir is None:
            output_dir = self.output_dir
        output_path = Path(output_dir)

        if pdb_file is None:
            pdb_file = output_path / f"{prefix}_lastframe.pdb"
        else:
            pdb_file = Path(pdb_file)

        if backup_file is None:
            backup_file = pdb_file.with_name(f"{pdb_file.stem}_raw{pdb_file.suffix}")
        else:
            backup_file = Path(backup_file)

        if wrapped_copy_file is not None:
            wrapped_copy_file = Path(wrapped_copy_file)

        if not pdb_file.exists():
            print(f"Warning: final PDB not found, skipping wrapped PDB postprocess: {pdb_file}")
            return None

        try:
            import mdtraj as md
        except ImportError:
            print("Warning: mdtraj is not available; could not wrap final PDB")
            return None

        print(f"Wrapping final PDB for next calculation: {pdb_file}")
        try:
            traj = md.load(str(pdb_file))
            anchor_molecules = self._anchor_molecules_for_mdtraj(traj.topology)
            imaged = traj.image_molecules(anchor_molecules=anchor_molecules)

            shutil.copy2(pdb_file, backup_file)
            imaged.save_pdb(str(pdb_file))
            if wrapped_copy_file is not None and wrapped_copy_file != pdb_file:
                shutil.copy2(pdb_file, wrapped_copy_file)
        except Exception as exc:
            print(f"Warning: could not wrap final PDB: {exc}")
            return None

        print(f"Saved raw final PDB backup: {backup_file}")
        print(f"Saved wrapped final PDB: {pdb_file}")
        if wrapped_copy_file is not None:
            print(f"Saved wrapped final PDB copy: {wrapped_copy_file}")
        return pdb_file

    def postprocess_wrapped_state(self, output_dir=None, prefix="md",
                                  state_file=None, topology_file=None,
                                  backup_file=None):
        """Rewrite OpenMM final State XML positions with wrapped coordinates."""
        if not self.special_wrapping or not self.wrap_final_state:
            return None

        if output_dir is None:
            output_dir = self.output_dir
        output_path = Path(output_dir)

        if state_file is None:
            state_file = Path("OpenMM_MD_final_state.xml")
        else:
            state_file = Path(state_file)

        if topology_file is None:
            topology_file = output_path / f"{prefix}_firstframe.pdb"
        else:
            topology_file = Path(topology_file)

        if backup_file is None:
            backup_file = state_file.with_name(f"{state_file.stem}_raw{state_file.suffix}")
        else:
            backup_file = Path(backup_file)

        if not state_file.exists():
            print(f"Warning: final state XML not found, skipping wrapped state postprocess: {state_file}")
            return None
        if not topology_file.exists():
            print(f"Warning: topology not found, skipping wrapped state postprocess: {topology_file}")
            return None

        try:
            import mdtraj as md
        except ImportError:
            print("Warning: mdtraj is not available; could not wrap final state XML")
            return None

        print(f"Wrapping final OpenMM state for next calculation: {state_file}")
        try:
            tree = ET.parse(state_file)
            root = tree.getroot()
            positions_node = root.find("Positions")
            box_node = root.find("PeriodicBoxVectors")
            if positions_node is None:
                raise ValueError("State XML has no Positions node")
            if box_node is None:
                raise ValueError("State XML has no PeriodicBoxVectors node")

            position_nodes = list(positions_node.findall("Position"))
            coords_nm = [
                [
                    float(pos.attrib["x"]),
                    float(pos.attrib["y"]),
                    float(pos.attrib["z"]),
                ]
                for pos in position_nodes
            ]

            box_vectors_nm = []
            for vector_name in ("A", "B", "C"):
                vector_node = box_node.find(vector_name)
                if vector_node is None:
                    raise ValueError(f"State XML has no {vector_name} box vector")
                box_vectors_nm.append([
                    float(vector_node.attrib.get("x", 0.0)),
                    float(vector_node.attrib.get("y", 0.0)),
                    float(vector_node.attrib.get("z", 0.0)),
                ])

            topology = md.load(str(topology_file)).topology
            if len(coords_nm) != topology.n_atoms:
                raise ValueError(
                    f"State/topology atom count mismatch: {len(coords_nm)} vs {topology.n_atoms}"
                )

            import numpy as np
            traj = md.Trajectory(np.array(coords_nm, dtype=float).reshape(1, -1, 3), topology)
            traj.unitcell_vectors = np.array(box_vectors_nm, dtype=float).reshape(1, 3, 3)
            anchor_molecules = self._anchor_molecules_for_mdtraj(traj.topology)
            imaged = traj.image_molecules(anchor_molecules=anchor_molecules)

            for pos_node, xyz in zip(position_nodes, imaged.xyz[0]):
                pos_node.set("x", f"{float(xyz[0]):.16g}")
                pos_node.set("y", f"{float(xyz[1]):.16g}")
                pos_node.set("z", f"{float(xyz[2]):.16g}")

            shutil.copy2(state_file, backup_file)
            tree.write(state_file, encoding="unicode", xml_declaration=True)
        except Exception as exc:
            print(f"Warning: could not wrap final state XML: {exc}")
            return None

        print(f"Saved raw final state backup: {backup_file}")
        print(f"Saved wrapped final state: {state_file}")
        return state_file
    
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
