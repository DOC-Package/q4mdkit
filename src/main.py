"""
Main entry point
"""
import logging
import sys
from pathlib import Path

# Add parent directory to path (for module imports)
sys.path.insert(0, str(Path(__file__).parent))

from config import (
    SIMULATION_CONFIG,
    SIMULATION_ENGINE,
    LOG_LEVEL,
    LOG_FORMAT,
    ensure_directories,
    get_output_path,
)
from common import (
    setup_logging,
    get_timestamp,
    print_header,
    print_dict,
    CrystalStructure,
    ResultAnalyzer,
)
from engines.openmm import OpenMMSimulator
from engines.pycharmm import PycharmmSimulator

logger = logging.getLogger(__name__)


def get_simulator(engine_name: str, config: dict):
    """Return a simulator instance for the requested engine.

    Args:
        engine_name: Engine name ('openmm' or 'pycharmm')
        config: Simulation configuration

    Returns:
        A simulator instance
    """
    if engine_name.lower() == "openmm":
        return OpenMMSimulator(config)
    elif engine_name.lower() == "pycharmm":
        return PycharmmSimulator(config)
    else:
        raise ValueError(f"Unsupported engine: {engine_name}")


def main():
    """Main execution function."""
    # Configure logging
    setup_logging(LOG_LEVEL, LOG_FORMAT)
    
    # Ensure directories exist
    ensure_directories()
    
    print_header(f"Crystal MD Simulation ({SIMULATION_ENGINE.upper()})")
    
    # 1. Generate crystal structure
    logger.info("Generating crystal structure...")
    crystal = CrystalStructure(
        lattice_constant=0.4,  # nm
        num_cells=(3, 3, 3)
    )
    positions = crystal.generate_fcc_lattice()
    box_vectors = crystal.get_box_vectors()
    logger.info(f"Number of atoms: {crystal.get_num_atoms()}")
    logger.info(f"Box size: {box_vectors[0, 0]:.2f} x {box_vectors[1, 1]:.2f} x {box_vectors[2, 2]:.2f} nm")
    
    # 2. Select simulation engine
    logger.info("\n" + "="*60)
    logger.info(f"Simulation engine: {SIMULATION_ENGINE}")
    simulator = get_simulator(SIMULATION_ENGINE, SIMULATION_CONFIG)
    
    # 3. Set up and run simulation
    simulator.setup(positions, box_vectors)
    simulator.run()
    
    # 4. Analyze results
    logger.info("\n" + "="*60)
    logger.info("Analyzing results...")
    results = simulator.get_results()
    analyzer = ResultAnalyzer(results)
    
    # Compute statistics
    stats = analyzer.calculate_statistics()
    print("\nStatistics:")
    print_dict(stats, indent=2)
    
    # 5. Save results
    logger.info("\n" + "="*60)
    timestamp = get_timestamp()
    output_base = get_output_path(f"simulation_{SIMULATION_ENGINE}_{timestamp}")
    
    analyzer.save_data(output_base)
    analyzer.plot_results(output_base)
    
    print_header("Done")
    logger.info(f"Results saved to {output_base.parent}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("\nInterrupted by user")
        sys.exit(0)
    except Exception as e:
        logger.error(f"An error occurred: {e}", exc_info=True)
        sys.exit(1)
