import numpy as np
import pytest

from q4mdkit.analysis.nwchem_cdftci import (
    NWChemCDFTCIConfig,
    NWChemCDFTCIState,
    build_nwchem_orbital_data,
    compute_nwchem_cdftci_from_outputs,
    load_cdftci_config,
    parse_nwchem_mo_coefficients,
    reconstruct_ao_overlap,
)
from q4mdkit.analysis.cdft_ci_core import compute_transfer_integral_unrestricted


def test_parse_nwchem_mo_coefficients():
    C_alpha, C_beta, occ_alpha, occ_beta = parse_nwchem_mo_coefficients(_MO_OUTPUT)

    expected = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]])
    assert C_alpha.shape == (3, 2)
    assert C_beta.shape == (3, 2)
    assert np.allclose(C_alpha, expected)
    assert np.allclose(C_beta, expected)
    assert np.allclose(occ_alpha, [1.0, 0.0])
    assert np.allclose(occ_beta, [1.0, 0.0])


def test_build_nwchem_orbital_data():
    orb = build_nwchem_orbital_data(_MO_OUTPUT)
    assert orb.n_alpha == 1
    assert orb.n_beta == 1
    assert orb.C_alpha.shape == (3, 2)


def test_parse_mo_coefficients_missing_section_raises():
    with pytest.raises(ValueError, match="Final MO vectors"):
        parse_nwchem_mo_coefficients("no molecular orbital analysis here")


def test_reconstruct_ao_overlap_identity():
    C = np.eye(4)
    assert np.allclose(reconstruct_ao_overlap(C), np.eye(4))


def _render_mo_section(C: np.ndarray, occ: list[float]) -> str:
    lines = [
        " DFT Final Alpha Molecular Orbital Analysis",
        " ------------------------------------------",
        "",
    ]
    n_bf = C.shape[0]
    n_mo = C.shape[1]
    for j in range(n_mo):
        lines.append(
            f" Vector    {j + 1}  Occ={occ[j]:.6f}D+00  E=-1.0D+00  Symmetry=a"
        )
        lines.append("              MO Center= 0.0D+00, 0.0D+00, 0.0D+00, r^2= 0.0D+00")
        lines.append(
            "   Bfn.  Coefficient  Atom+Function         Bfn.  Coefficient  Atom+Function"
        )
        lines.append(
            "  ----- ------------  ---------------      ----- ------------  ---------------"
        )
        for i in range(0, n_bf, 2):
            left = f"     {i + 1}      {C[i, j]:.6f}  {i + 1} C  s"
            if i + 1 < n_bf:
                right = f"     {i + 2}      {C[i + 1, j]:.6f}  {i + 2} C  s"
            else:
                right = ""
            lines.append(left + "          " + right)
        lines.append("")
    return "\n".join(lines)


def _render_evecs(C: np.ndarray, spin: str) -> str:
    n_bf, n_mo = C.shape
    lines = [f" global array: {spin} evecs[1:{n_bf},1:{n_mo}],  handle: 0", ""]
    for c0 in range(0, n_mo, 6):
        cols = list(range(c0, min(c0 + 6, n_mo)))
        lines.append("  " + "".join(f"{c + 1:>12d}" for c in cols))
        lines.append("  " + "".join("------------" for _ in cols))
        for r in range(n_bf):
            lines.append("  " + f"{r + 1:>3d}" + "".join(f"{C[r, c]:12.6f}" for c in cols))
        lines.append("")
    return "\n".join(lines)


def _render_output(C: np.ndarray, occ: list[float], energy: float, mult: list[float]) -> str:
    alpha = _render_mo_section(C, occ)
    beta = alpha.replace("Alpha", "Beta")
    alpha_e = _render_evecs(C, "alpha")
    beta_e = _render_evecs(C, "beta")
    mult_lines = "\n".join(f"      {i + 1}        {m:.10f}" for i, m in enumerate(mult))
    return "\n".join(
        [
            alpha,
            "",
            beta,
            "",
            " Final MO vectors",
            " ----------------",
            "",
            alpha_e,
            beta_e,
            " alpha - beta orbital overlaps",
            "  CDFT final multipliers",
            mult_lines,
            "",
            f"         Total DFT energy =    {energy:.10f}",
        ]
    )


_MO_OUTPUT = _render_output(
    np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]), [1.0, 0.0], 0.0, [0.0]
)


def test_compute_nwchem_cdftci_from_outputs():
    C_A = np.array([[1.0, 0.0], [0.0, 1.0]])
    C_B = np.array([[0.6, 0.8], [0.8, -0.6]])
    occ = [1.0, 0.0]

    out_A = _render_output(C_A, occ, -1.0, [0.5])
    out_B = _render_output(C_B, occ, -2.0, [0.7])

    ham = compute_nwchem_cdftci_from_outputs(
        out_A,
        out_B,
        fragment_A_atoms=[1],
        fragment_B_atoms=[2],
        N_A=1.0,
        N_B=1.0,
    )

    # S_AB = 0.6 * 0.6 = 0.36, W_BA = 0.72, W_AB = 0
    assert np.isclose(ham.S_AB, 0.36, atol=1e-10)
    # H_AB = 0.5*(-1-2+0.5+0.7)*0.36 - 0.5*(0.5*0.72) = -0.504
    assert np.isclose(ham.H_AB, -0.504, atol=1e-10)

    J = compute_transfer_integral_unrestricted(ham.H, ham.S, method="direct")
    assert np.isclose(J, 0.036, atol=1e-10)


def test_compute_nwchem_cdftci_from_outputs_applies_net_charge_correction():
    C_A = np.array([[1.0, 0.0], [0.0, 1.0]])
    C_B = np.array([[0.6, 0.8], [0.8, -0.6]])
    occ = [1.0, 0.0]

    out_A = _render_output(C_A, occ, -1.0, [0.5])
    out_B = _render_output(C_B, occ, -2.0, [0.7])

    ham = compute_nwchem_cdftci_from_outputs(
        out_A,
        out_B,
        fragment_A_atoms=[1],
        fragment_B_atoms=[2],
        N_A=1.0,
        N_B=1.0,
        Z_A=1.0,
        Z_B=1.0,
    )

    assert np.isclose(ham.S_AB, 0.36, atol=1e-10)
    assert np.isclose(ham.W_BA, -0.36, atol=1e-10)
    assert np.isclose(ham.W_AB, 0.36, atol=1e-10)
    assert np.isclose(ham.H_AB, -0.36, atol=1e-10)

    J = compute_transfer_integral_unrestricted(ham.H, ham.S, method="direct")
    assert np.isclose(J, 0.18, atol=1e-10)


def test_state_build_constraints_charge_only_and_with_spin():
    charge_only = NWChemCDFTCIState("a", "1:36", 1.0)
    assert [c.kind for c in charge_only.build_constraints()] == ["charge"]

    with_spin = NWChemCDFTCIState("a", "1:36", 1.0, spin=1.0)
    assert [c.kind for c in with_spin.build_constraints()] == ["charge", "spin"]


def test_load_cdftci_config(tmp_path):
    (tmp_path / "traj.dcd").write_text("", encoding="utf-8")
    (tmp_path / "top.pdb").write_text("", encoding="utf-8")
    (tmp_path / "qmatoms").write_text("1 2 3 4\n", encoding="utf-8")
    config_path = tmp_path / "nwchem_cdftci.yaml"
    config_path.write_text(
        """
input:
  trajectory: traj.dcd
  topology: top.pdb
  qm_atoms: qmatoms
output:
  directory: out
frames:
  start: 0
  n_frames: 1
cdftci:
  state_a:
    name: donor
    atoms: "1:2"
    charge: 1.0
    population: mulliken
  state_b:
    name: acceptor
    atoms: "3:4"
    charge: -1.0
    spin: 1.0
""",
        encoding="utf-8",
    )

    config = load_cdftci_config(config_path)

    assert isinstance(config, NWChemCDFTCIConfig)
    assert config.state_A.name == "donor"
    assert config.state_A.atoms == "1:2"
    assert config.state_A.charge == 1.0
    assert config.state_B.name == "acceptor"
    assert config.state_B.atoms == "3:4"
    assert config.state_B.charge == -1.0
    assert config.state_B.spin == 1.0
    assert config.coupling_file == tmp_path / "out" / "cdftci.dat"
