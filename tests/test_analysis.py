"""
解析モジュールのテスト
"""
import sys
from pathlib import Path
import tempfile

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from modules.analysis import ResultAnalyzer


def test_statistics_calculation():
    """統計量計算のテスト"""
    # ダミーデータ
    results = {
        "steps": [0, 100, 200, 300],
        "energy": [-1000.0, -1005.0, -1002.0, -1003.0],
        "temperature": [300.0, 298.0, 302.0, 301.0],
    }
    
    analyzer = ResultAnalyzer(results)
    stats = analyzer.calculate_statistics()
    
    # 統計量の存在確認
    assert "energy_mean" in stats
    assert "energy_std" in stats
    assert "temperature_mean" in stats
    assert "temperature_std" in stats
    
    # 値の範囲確認
    assert -1010 < stats["energy_mean"] < -1000
    assert 295 < stats["temperature_mean"] < 305
    
    print("✓ 統計量計算テスト合格")


def test_data_saving():
    """データ保存のテスト"""
    results = {
        "steps": [0, 100, 200],
        "energy": [-1000.0, -1005.0, -1002.0],
        "temperature": [300.0, 298.0, 302.0],
    }
    
    analyzer = ResultAnalyzer(results)
    
    # 一時ファイルに保存
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "test_output"
        analyzer.save_data(output_path)
        
        csv_path = output_path.with_suffix(".csv")
        assert csv_path.exists(), "CSVファイルが作成されていません"
        
        # ファイル内容の確認
        with open(csv_path, "r") as f:
            lines = f.readlines()
            assert len(lines) == 4, "データ行数が不正"  # ヘッダー + 3行
            assert "step,energy,temperature" in lines[0], "ヘッダーが不正"
    
    print("✓ データ保存テスト合格")


if __name__ == "__main__":
    test_statistics_calculation()
    test_data_saving()
    print("\n全テスト合格!")
