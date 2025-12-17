import xml.etree.ElementTree as ET

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