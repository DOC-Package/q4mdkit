"""
Tests for analysis module
"""
import sys
from pathlib import Path
import tempfile

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from common.analysis import ResultAnalyzer


def test_statistics_calculation():
    """Test statistics calculation."""
    # Dummy data
    results = {
        "steps": [0, 100, 200, 300],
        "energy": [-1000.0, -1005.0, -1002.0, -1003.0],
        "temperature": [300.0, 298.0, 302.0, 301.0],
    }
    
    analyzer = ResultAnalyzer(results)
    stats = analyzer.calculate_statistics()
    
    # Check presence of statistics
    assert "energy_mean" in stats
    assert "energy_std" in stats
    assert "temperature_mean" in stats
    assert "temperature_std" in stats
    
    # Check value ranges
    assert -1010 < stats["energy_mean"] < -1000
    assert 295 < stats["temperature_mean"] < 305
    
    print("✓ Statistics calculation test passed")


def test_data_saving():
    """Test data saving."""
    results = {
        "steps": [0, 100, 200],
        "energy": [-1000.0, -1005.0, -1002.0],
        "temperature": [300.0, 298.0, 302.0],
    }
    
    analyzer = ResultAnalyzer(results)
    
    # Save to temporary file
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "test_output"
        analyzer.save_data(output_path)
        
        csv_path = output_path.with_suffix(".csv")
        assert csv_path.exists(), "CSV file not created"
        
        # Check file contents
        with open(csv_path, "r") as f:
            lines = f.readlines()
            assert len(lines) == 4, "Invalid number of data rows"  # header + 3 rows
            assert "step,energy,temperature" in lines[0], "Invalid header"
    
    print("✓ Data saving test passed")


if __name__ == "__main__":
    test_statistics_calculation()
    test_data_saving()
    print("\nAll tests passed!")
