import sys
from types import SimpleNamespace


class _FakeDftbPlus:
    pass


sys.modules.setdefault("dftbplus", SimpleNamespace(DftbPlus=_FakeDftbPlus))
sys.modules.setdefault("hsd", SimpleNamespace(load=None, dump=None))
sys.modules.setdefault("mdtraj", SimpleNamespace(Topology=object, Trajectory=object))

from q4mdkit.analysis.cdftbci import (  # noqa: E402
    _coarse_frame_selected,
    _fine_frame_indices_between,
    _phase_adaptive_refinement_indices,
    _write_interval_ci_output,
    load_cdftbci_config,
)


def test_load_cdftbci_config_reads_phase_adaptive_refinement(tmp_path):
    config_path = tmp_path / "cdftbci.yaml"
    config_path.write_text(
        """
input:
  trajectory: traj/nve.dcd
  topology: traj/top.pdb
  qm_atoms: traj/qmatoms
  pccharges_template: input/PCcharges.dat
  hsd_template: input/dftb_in.hsd
output:
  directory: output
frames:
  start: 10
  n_frames: 3
  dt_fs: 4.0
  stride: 4
fragments:
  - name: fragment1
    atom_range: "1:36"
  - name: fragment2
    atom_range: "37:72"
phase_tracking:
  enabled: true
  adaptive_refinement:
    enabled: true
    trajectory_1fs: traj/nve-1fs.dcd
    sigma_min_threshold: 0.7
    fine_dt_fs: 1.0
    log_file: adaptive_refinement.dat
""",
        encoding="utf-8",
    )

    config = load_cdftbci_config(config_path)

    assert config.phase_adaptive_refinement_enabled is True
    assert config.phase_adaptive_refinement_traj_path == tmp_path / "traj/nve-1fs.dcd"
    assert config.phase_adaptive_refinement_sigma_min_threshold == 0.7
    assert config.phase_adaptive_refinement_fine_dt_fs == 1.0
    assert config.phase_adaptive_refinement_log_file == tmp_path / "output/adaptive_refinement.dat"
    assert config.start_frame == 10
    assert config.n_frames == 3
    assert config.dt_fs == 4.0
    assert config.frame_stride == 4


def test_coarse_frame_selected_uses_stride_independent_of_dt_fs():
    assert _coarse_frame_selected(frame_id=10, start_frame=10, stride=4)
    assert not _coarse_frame_selected(frame_id=11, start_frame=10, stride=4)
    assert not _coarse_frame_selected(frame_id=13, start_frame=10, stride=4)
    assert _coarse_frame_selected(frame_id=14, start_frame=10, stride=4)


def test_write_interval_ci_output_skips_intermediate_refinement_frames():
    assert _write_interval_ci_output(frame_is_fine=False, remaining_refinement_frames=0)
    assert not _write_interval_ci_output(frame_is_fine=True, remaining_refinement_frames=3)
    assert not _write_interval_ci_output(frame_is_fine=True, remaining_refinement_frames=1)
    assert _write_interval_ci_output(frame_is_fine=True, remaining_refinement_frames=0)


def test_phase_adaptive_refinement_defaults_to_input_trajectory_and_dt(tmp_path):
    config_path = tmp_path / "cdftbci.yaml"
    config_path.write_text(
        """
input:
  trajectory: traj/nve-1fs.dcd
  topology: traj/top.pdb
  qm_atoms: traj/qmatoms
  pccharges_template: input/PCcharges.dat
  hsd_template: input/dftb_in.hsd
output:
  directory: output
frames:
  dt_fs: 1.0
  stride: 4
fragments:
  - name: fragment1
    atom_range: "1:36"
  - name: fragment2
    atom_range: "37:72"
phase_tracking:
  enabled: true
  adaptive_refinement:
    enabled: true
    sigma_min_threshold: 0.7
""",
        encoding="utf-8",
    )

    config = load_cdftbci_config(config_path)

    assert config.phase_adaptive_refinement_enabled is True
    assert config.phase_adaptive_refinement_traj_path == config.traj_path
    assert config.phase_adaptive_refinement_fine_dt_fs == config.dt_fs


def test_fine_frame_indices_between_coarse_phase_frames():
    assert _fine_frame_indices_between(
        prev_time_fs=0.0,
        curr_time_fs=4.0,
        fine_dt_fs=1.0,
    ) == [1, 2, 3]


def test_phase_adaptive_refinement_indices_trigger_only_below_threshold():
    assert _phase_adaptive_refinement_indices(
        enabled=True,
        previous_time_fs=0.0,
        current_time_fs=4.0,
        sigma_min_a=0.8,
        sigma_min_b=0.6,
        threshold=0.7,
        fine_dt_fs=1.0,
    ) == [1, 2, 3, 4]
    assert _phase_adaptive_refinement_indices(
        enabled=True,
        previous_time_fs=0.0,
        current_time_fs=4.0,
        sigma_min_a=0.8,
        sigma_min_b=0.75,
        threshold=0.7,
        fine_dt_fs=1.0,
    ) == []
    assert _phase_adaptive_refinement_indices(
        enabled=False,
        previous_time_fs=0.0,
        current_time_fs=4.0,
        sigma_min_a=0.1,
        sigma_min_b=0.1,
        threshold=0.7,
        fine_dt_fs=1.0,
    ) == []
