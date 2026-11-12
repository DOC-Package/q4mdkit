import math
import xml.etree.ElementTree as ET


def read_restart_box(statefile, expected_atoms=None):
    """Validate an OpenMM State XML and return its box vectors in angstrom.

    Positions and velocities must contain finite x/y/z values for the same
    number of atoms. Pass ``expected_atoms`` to check a particular system;
    otherwise the atom count is inferred from the State's positions.
    OpenMM stores box vectors in nm, so all vector components are multiplied
    by 10 without rounding. The source file is left unchanged.

    Args:
        statefile: Path to an OpenMM State XML file.
        expected_atoms: Optional expected number of atoms in the restart.

    Returns:
        The three box vectors as a 3-by-3 list of floats in angstrom.

    Raises:
        ValueError: Particle arrays or periodic box vectors are invalid.
    """
    state = ET.parse(statefile).getroot()
    if expected_atoms is None:
        expected_atoms = len(state.findall("Positions/Position"))
    if expected_atoms <= 0:
        raise ValueError(f"{statefile}: Positions must contain at least one atom")
    for group, particle in (("Positions", "Position"), ("Velocities", "Velocity")):
        particles = state.findall(f"{group}/{particle}")
        if len(particles) != expected_atoms:
            raise ValueError(f"{statefile}: {group} must contain {expected_atoms} atoms")
        try:
            valid = all(math.isfinite(float(atom.attrib[axis]))
                        for atom in particles for axis in ("x", "y", "z"))
        except (KeyError, ValueError):
            valid = False
        if not valid:
            raise ValueError(f"{statefile}: {group} require finite numeric x/y/z values")
    vectors = []
    for name in ("A", "B", "C"):
        node = state.find(f"PeriodicBoxVectors/{name}")
        if node is None:
            raise ValueError(f"{statefile}: missing periodic box vector {name}")
        try:
            vector = [10.0 * float(node.attrib[axis]) for axis in ("x", "y", "z")]
        except (KeyError, ValueError) as error:
            raise ValueError(f"{statefile}: box vector {name} requires numeric x/y/z values") from error
        if not all(math.isfinite(value) for value in vector):
            raise ValueError(f"{statefile}: non-finite periodic box vector {name}")
        vectors.append(vector)
    if any(vectors[i][i] <= 0 for i in range(3)):
        raise ValueError(f"{statefile}: periodic box lengths must be positive")
    return vectors


def remove_montecarlo_xml(input_xml, output_xml):
    """
    Remove MonteCarlo barostat parameters from OpenMM state XML file.
    """
    # Parse XML
    tree = ET.parse(input_xml)
    root = tree.getroot()
    
    # Find State element
    state = root.find('State') if root.tag != 'State' else root
    if state is None:
        state = root  # Root might be State itself
    
    # Find and remove MonteCarlo-related attributes from Parameters element
    params = state.find('Parameters')
    if params is not None:
        # Remove MonteCarlo-related attributes
        attribs_to_remove = []
        for attrib in params.attrib:
            if 'MonteCarlo' in attrib:
                attribs_to_remove.append(attrib)
        
        for attrib in attribs_to_remove:
            del params.attrib[attrib]
            print(f"  Removed parameter: {attrib}")
        
        # If Parameters element is now empty, remove it entirely
        if len(params.attrib) == 0:
            state.remove(params)
            print("  Removed empty Parameters element")
    
    # Write to new file
    tree.write(output_xml, xml_declaration=True, encoding='unicode')
    print(f"  Saved to: {output_xml}")
