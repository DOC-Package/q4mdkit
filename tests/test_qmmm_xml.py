"""Read and validate OpenMM restart cells without ASH or OpenMM installed."""

import xml.etree.ElementTree as ET

import pytest

from q4mdkit.qmmm import xml


@pytest.fixture
def restart(tmp_path):
    path = tmp_path / "state.xml"
    state = ET.Element("State")
    box = ET.SubElement(state, "PeriodicBoxVectors")
    ET.SubElement(box, "A", x="4.4255050604395985", y="0", z="0")
    ET.SubElement(box, "B", x="0.2", y="4.5", z="0")
    ET.SubElement(box, "C", x="0.1", y="0.3", z="4.6")
    for group, particle in (("Positions", "Position"), ("Velocities", "Velocity")):
        node = ET.SubElement(state, group)
        for _ in range(2):
            ET.SubElement(node, particle, x="0.1", y="0.2", z="0.3")
    ET.ElementTree(state).write(path)
    return path


def test_restart_box_is_reusable_and_preserves_precision_and_source(restart):
    original = restart.read_bytes()
    vectors = xml.read_restart_box(restart)
    assert vectors[0] == pytest.approx([44.255050604395985, 0, 0], abs=1e-12)
    assert vectors[1] == pytest.approx([2, 45, 0])
    assert vectors[2] == pytest.approx([1, 3, 46])
    assert restart.read_bytes() == original


def test_restart_box_checks_explicit_system_atom_count(restart):
    with pytest.raises(ValueError, match="Positions must contain 3 atoms"):
        xml.read_restart_box(restart, expected_atoms=3)


def test_restart_box_checks_velocity_count_when_atom_count_is_inferred(restart):
    tree = ET.parse(restart)
    velocities = tree.getroot().find("Velocities")
    velocities.remove(velocities[0])
    tree.write(restart)
    with pytest.raises(ValueError, match="Velocities must contain 2 atoms"):
        xml.read_restart_box(restart)


@pytest.mark.parametrize("group,particle", [("Positions", "Position"), ("Velocities", "Velocity")])
@pytest.mark.parametrize("value", [None, "nan", "invalid"])
def test_restart_box_rejects_incomplete_or_nonfinite_particles(restart, group, particle, value):
    tree = ET.parse(restart)
    atom = tree.getroot().find(f"{group}/{particle}")
    if value is None:
        del atom.attrib["x"]
    else:
        atom.set("x", value)
    tree.write(restart)
    with pytest.raises(ValueError, match=group):
        xml.read_restart_box(restart, expected_atoms=2)


@pytest.mark.parametrize("problem", ["missing_vector", "missing_component", "nan", "negative"])
def test_restart_box_rejects_invalid_cell(restart, problem):
    tree = ET.parse(restart)
    box = tree.getroot().find("PeriodicBoxVectors")
    if problem == "missing_vector":
        box.remove(box.find("B"))
    elif problem == "missing_component":
        del box.find("B").attrib["x"]
    elif problem == "nan":
        box.find("B").set("y", "nan")
    else:
        box.find("B").set("y", "-1")
    tree.write(restart)
    with pytest.raises(ValueError, match="box"):
        xml.read_restart_box(restart)
