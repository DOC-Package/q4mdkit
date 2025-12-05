#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pyCHARMM でエネルギー最小化 → NVT → NPT を行う最小実装。

✔ 前提:
  - 既に PSF + 座標 (CRD/PDB) がある想定。
  - 周期境界 (CRYSTAL/IMAGE) + PME を使用。
  - 温度制御は Nosé–Hoover、圧力制御は Langevin piston。


if __name__ == "__main__":
 以降の設定を自分の系に合わせて修正
  2) python run_pycharmm_md.py などの名前で実行

注: pyCHARMM には Python API もありますが、ここでは最も移植性が高い
    CHARMM コマンド文字列 (lingo.charmm_script) を中心に使います。
"""
from __future__ import annotations
import os
import sys
import subprocess

# --- pyCHARMM のロードと簡易ラッパ ---
try:
    import pycharmm
    from pycharmm import lingo
    cs = lingo.charmm_script  # CHARMM コマンド1行を実行
except Exception as e:  # 旧/別名フォールバック
    from pycharmm import charmm_script as cs  # type: ignore


def open_unit(unit: int, mode: str, path: str, file_type: str = "CARD"):
    """Open a file unit in CHARMM.
    
    Parameters
    ----------
    unit : int
        Unit number
    mode : str
        'READ' or 'WRITE'
    path : str
        File path
    file_type : str
        'CARD' for text files, 'FILE' for binary files (e.g., DCD)
    """
    path = os.path.abspath(path)
    cs(f"OPEN {mode.upper()} {file_type.upper()} UNIT {unit:d} NAME {path}")


def find_valid_fft_dimension(n: int) -> int:
    """Find the smallest integer >= n that has only 2, 3, 5 as prime factors (for PME FFT)."""
    if n <= 1:
        return 2
    candidate = n
    while True:
        temp = candidate
        # Remove all factors of 2, 3, 5
        while temp % 2 == 0:
            temp //= 2
        while temp % 3 == 0:
            temp //= 3
        while temp % 5 == 0:
            temp //= 5
        # If only 2, 3, 5 factors, temp should be 1
        if temp == 1:
            return candidate
        candidate += 1


# --- PDBのCRYST1からセルを読む & 結晶型を推定 ---
def _read_pdb_cryst1(path: str) -> tuple[tuple[float,float,float], tuple[float,float,float]] | None:
    """PDBのCRYST1行から (a,b,c),(alpha,beta,gamma) を返す。見つからなければNone。"""
    try:
        with open(path, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("CRYST1"):
                    a = float(line[6:15]); b = float(line[15:24]); c = float(line[24:33])
                    alpha = float(line[33:40]); beta = float(line[40:47]); gamma = float(line[47:54])
                    return (a,b,c), (alpha,beta,gamma)
    except Exception:
        pass
    return None


def _read_crd_box(path: str) -> tuple[tuple[float,float,float], tuple[float,float,float]] | None:
    """CRDファイルの最終行からボックス情報 (a,b,c),(alpha,beta,gamma) を返す。"""
    try:
        with open(path, "rt", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
        # 最終行をチェック（ボックス情報は6つの浮動小数点数）
        last_line = lines[-1].strip()
        parts = last_line.split()
        if len(parts) == 6:
            a, b, c, alpha, beta, gamma = map(float, parts)
            # 妥当性チェック（ボックスサイズが正、角度が0-180度）
            if all(x > 0 for x in [a, b, c]) and all(0 < ang < 180 for ang in [alpha, beta, gamma]):
                return (a, b, c), (alpha, beta, gamma)
    except Exception:
        pass
    return None


def _guess_crystal_type(abc: tuple[float,float,float], ang: tuple[float,float,float]) -> str:
    """角度と格子長から CHARMM の型名を推定 (CUBIC/TETRagonal/HEXagonal/RHDO/ORTHorhombic/MONOclinic/TRIClinic)。"""
    a,b,c = abc; alpha,beta,gamma = ang
    def eq(x,y,tol=1e-3):
        return abs(x-y) <= tol*max(1.0,abs(x),abs(y))
    a90 = all(abs(x-90.0) < 1e-2 for x in (alpha,beta,gamma))
    if a90 and eq(a,b) and eq(b,c):
        return "CUBIC"
    if a90 and eq(a,b) and not eq(b,c):
        return "TETRagonal"
    if abs(alpha-90.0)<1e-2 and abs(beta-90.0)<1e-2 and abs(gamma-120.0)<1e-2 and eq(a,b):
        return "HEXAgonal"
    if eq(a,b) and eq(b,c) and (abs(alpha-90.0)>1e-2 or abs(beta-90.0)>1e-2 or abs(gamma-90.0)>1e-2) and abs(alpha-beta)<1e-2 and abs(beta-gamma)<1e-2:
        return "RHDO"
    if a90:
        return "ORTHorhombic"
    if abs(alpha-90.0)<1e-2 and abs(gamma-90.0)<1e-2 and abs(beta-90.0)>1e-2:
        return "MONOclinic"
    return "TRIClinic"


def setup_system(
    psf_path: str,
    coor_path: str,
    toppar_stream: str | None,
    box: tuple[float, float, float] | None = None,  # (A, B, C) in Angstrom; None→PDBから自動
    pbc_type: str = "AUTO",   # AUTO→PDBのCRYST1から推定 / それ以外は CHARMM 型名
    box_from_pdb: str | None = None,  # CRD使用時でも別PDBのCRYST1から取得したい場合
    use_pme: bool = True,
    cutnb: float = 14.0,
    ctofnb: float = 12.0,
    ctonnb: float = 10.0,
    pme_kappa: float = 0.34,
    pme_order: int = 6,
    pme_grid: tuple[int, int, int] | None = None,
):
    """PSF/座標を読み、PBC+PME をセットアップし、初期エネルギーを評価。
    - pbc_type="AUTO" かつ coor_path が PDBなら CRYST1 を読んで自動設定
    - 角度が 90 度でない斜方晶/単斜/三斜格子にも対応
    """
    cs("BOMLEV -2")  # Suppress warnings during topology reading
    cs("PRNLEV 5")

    # 1) トポロジ/パラメータ
    # Load base CGenFF force field if available
    base_rtf = "/home/takahashi/charmm/charmm/toppar/top_all36_cgenff.rtf"
    base_prm = "/home/takahashi/charmm/charmm/toppar/par_all36_cgenff.prm"
    
    if os.path.exists(base_rtf) and os.path.exists(base_prm):
        print(f"Loading base CHARMM36 CGenFF force field...")
        cs(f"open read unit 1 card name {base_rtf}")
        cs("read rtf card unit 1")
        cs("close unit 1")
        cs(f"open read unit 1 card name {base_prm}")
        cs("read param card unit 1")
        cs("close unit 1")
    
    # Load molecule-specific parameters
    if toppar_stream:
        print(f"Loading molecule-specific parameters: {toppar_stream}")
        # Read and modify stream file to remove 'flex' keyword (causes issues with append)
        import tempfile
        with open(toppar_stream, 'r') as f:
            content = f.read()
        content_fixed = content.replace('read param card flex append', 'read param card append')
        content_fixed = content_fixed.replace('READ PARAM CARD FLEX APPEND', 'READ PARAM CARD APPEND')
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.str', delete=False) as tmp:
            tmp_path = tmp.name
            tmp.write(content_fixed)
        
        cs(f"stream {tmp_path}")
        os.unlink(tmp_path)
    
    cs("BOMLEV 0")  # Reset to normal warning level

    # 2) PSF & 座標の読込
    open_unit(10, "READ", psf_path)
    cs("READ PSF CARD UNIT 10")
    cs("CLOSE UNIT 10")

    open_unit(20, "READ", coor_path)
    is_pdb = coor_path.lower().endswith((".pdb", ".ent"))
    if is_pdb:
        cs("READ COOR PDB UNIT 20")
    else:
        cs("READ COOR CARD UNIT 20")
    cs("CLOSE UNIT 20")
    # Don't use COOR ORIENT with PBC - it breaks the box alignment
    # cs("COOR ORIENT")
    cs("COOR STAT")

    # 3) ボックス情報の取得
    abc = None
    ang = None
    
    if pbc_type.upper() == "AUTO":
        # 自動検出: PDBならCRYST1、CRDなら最終行、またはbox_from_pdbから
        if is_pdb:
            cryst = _read_pdb_cryst1(coor_path)
        else:
            # CRDからボックス情報を読む
            cryst = _read_crd_box(coor_path)
            # CRDにない場合、box_from_pdbがあればそこから
            if cryst is None and box_from_pdb:
                cryst = _read_pdb_cryst1(box_from_pdb)
        
        if cryst is not None:
            abc, ang = cryst
            pbc_type_eff = _guess_crystal_type(abc, ang)
            print(f"Auto-detected box: {abc[0]:.3f} {abc[1]:.3f} {abc[2]:.3f} Å, {ang[0]:.2f}° {ang[1]:.2f}° {ang[2]:.2f}°")
            print(f"Crystal type: {pbc_type_eff}")
        else:
            pbc_type_eff = "ORTHorhombic"
            print("Warning: Could not auto-detect box, will use manual 'box' parameter")
    else:
        pbc_type_eff = pbc_type

    # ボックスサイズの決定
    if abc is None:
        if box is None:
            raise ValueError("box を指定するか、PDBのCRYST1/CRDのボックス情報から自動取得してください (pbc_type='AUTO')")
        ax, by, cz = box
        ang = (90.0, 90.0, 90.0)
    else:
        ax, by, cz = abc

    alpha,beta,gamma = ang if ang is not None else (90.0,90.0,90.0)
    cs(f"CRYSTAL DEFINE {pbc_type_eff} {ax:.3f} {by:.3f} {cz:.3f} {alpha:.2f} {beta:.2f} {gamma:.2f}")
    cs("CRYSTAL BUILD NOPER 0")
    cs("IMAGE BYRES XCEN 0.0 YCEN 0.0 ZCEN 0.0 SELECT ALL END")

    # 4) 非結合 + PME 設定
    if use_pme:
        if pme_grid is None:
            # PME grid must have only 2, 3, 5 as prime factors
            gx = find_valid_fft_dimension(max(32, int(round(ax))))
            gy = find_valid_fft_dimension(max(32, int(round(by))))
            gz = find_valid_fft_dimension(max(32, int(round(cz))))
            print(f"PME grid: {gx} x {gy} x {gz}")
        else:
            gx, gy, gz = pme_grid
        cs(
            "NBOND ATOM CDIEL SWITCH VATOM VDIST VSWITCH "
            f"CUTNB {cutnb:.2f} CTOFNB {ctofnb:.2f} CTONNB {ctonnb:.2f} "
            f"EWALD PMEWALD KAPPA {pme_kappa:.3f} ORDER {pme_order:d} "
            f"FFTX {gx:d} FFTY {gy:d} FFTZ {gz:d}"
        )
    else:
        cs(
            "NBOND ATOM CDIEL SWITch VATOM VDIST VSWITCH "
            f"CUTNB {cutnb:.2f} CTOFNB {ctofnb:.2f} CTONNB {ctonnb:.2f}"
        )

    cs("ENERGY")


def setup_shake(shake_hyd: bool = True, shake_fast: bool = False, tol: float = 1e-8):
    """SHAKE制約を設定する。
    
    Parameters
    ----------
    shake_hyd : bool
        水素を含む結合にSHAKEを適用 (推奨: True)
    shake_fast : bool
        高速原子 (e.g. CH3) にもSHAKEを適用 (デフォルト: False)
    tol : float
        SHAKE収束許容誤差 (デフォルト: 1e-8)
    """
    shake_cmd = "SHAKE"
    
    if shake_hyd:
        shake_cmd += " BOND"  # 水素を含む結合を固定
    
    if shake_fast:
        shake_cmd += " FAST"  # 高速原子も固定
    
    shake_cmd += f" TOL {tol:.2e}"
    
    cs(shake_cmd)
    print(f"SHAKE constraints enabled: {shake_cmd}")


def minimize_energy(nsd: int = 1000, nabnr: int = 2000):
    """Steepest Descent → ABNR の 2段最小化"""
    cs(f"MINIMIZE SD NSTEP {nsd:d} NPRINt 100")
    cs(f"MINIMIZE ABNR NSTEP {nabnr:d} NPRINt 100")
    cs("ENERGY")


def run_nvt(
    temp: float,
    nsteps: int,
    dt_ps: float,
    traj_path: str = "nvt.dcd",
    rst_path: str = "nvt.restart",
    nsavc: int = 500,
    nsavr: int = 5000,
    iasvel: int = 1,
    seed: int = 314159,
    tmass: float = 1000.0,
    use_shake: bool = True,
):
    """Nosé–Hoover で NVT。
    DYNAmics CPT の temperature-spec のみを与える。
    """
    shake_flag = "SHAKE" if use_shake else ""

    # Open trajectory/restart units only if paths are provided
    dyn_io = ""
    if traj_path:
        # DCD files are binary, use UNFORM (unformatted) type
        cs(f"OPEN WRITE UNFORM UNIT 51 NAME {os.path.abspath(traj_path)}")
        dyn_io += f"IUNCRD 51 NSAVC {nsavc:d} "
    
    # Restart file handling - let CHARMM handle it automatically
    if rst_path:
        dyn_io += f"IUNWRI 52 NSAVR {nsavr:d} IUNREA -1 "

    cs(
        f"DYNAMICS CPT START {shake_flag} "
        f"NSTEPS {nsteps:d} TIMESTEP {dt_ps:.6f} "
        f"{dyn_io}NTRFRQ 100 "
        f"IASORS 1 IASVEL {iasvel:d} ISEED {seed:d} "
        f"HOOVER TMASS {tmass:.1f} REFT {temp:.1f} "
    )


def run_npt(
    temp: float,
    press_atm: float,
    nsteps: int,
    dt_ps: float,
    traj_path: str = "npt.dcd",
    rst_path: str = "npt.restart",
    nsavc: int = 500,
    nsavr: int = 5000,
    tmass: float = 1000.0,
    pmass: float = 400.0,
    pgamma: float = 20.0,
    fbeta: float = 5.0,
    iasvel: int = 0,
    seed: int = 12345,
    use_shake: bool = True,
):
    """Nosé–Hoover(温度) + Langevin piston(圧力) で NPT。
    DYNAmics CPT の pressure-spec + temperature-spec を併用。
    """
    # Disable energy check for NPT with large pressure changes during equilibration
    cs("BOMLEV -2")
    #cs("WRNLEV -")
    # Set friction coefficient for Langevin piston
    cs(f"SCALAR FBETA SET {fbeta:.2f} SELECT ALL END")
    
    shake_flag = "SHAKE" if use_shake else ""

    dyn_io = ""
    if traj_path:
        # DCD files are binary, use UNFORM (unformatted) type
        cs(f"OPEN WRITE UNFORM UNIT 61 NAME {os.path.abspath(traj_path)}")
        dyn_io += f"IUNCRD 61 NSAVC {nsavc:d} "
    
    # Restart file handling - let CHARMM handle it automatically
    if rst_path:
        dyn_io += f"IUNWRI 62 NSAVR {nsavr:d} IUNREA -1 "

    cs(
        f"DYNAMICS CPT START {shake_flag} "
        f"NSTEPS {nsteps:d} TIMESTEP {dt_ps:.6f} "
        f"{dyn_io}NTRFRQ 100 "
        f"IASORS 1 IASVEL {iasvel:d} ISEED {seed:d} "
        # 圧力 (Langevin piston)
        f"PCONST LANGEVIN PMASS {pmass:.1f} PGAMMA {pgamma:.1f} TBATH {temp:.1f} "
        f"PREF {press_atm:.3f} "
        # 温度 (Hoover)
        f"HOOVER TMASS {tmass:.1f} REFT {temp:.1f} "
    )


def write_snapshots(prefix: str = "out", psf_path: str = "pentacene.psf"):
    """Write snapshot as CRD and PDB."""
    os.makedirs("out", exist_ok=True)
    
    # Write CRD using multi-line script to avoid EOF issues
    crd_file = os.path.join("out", f"{prefix}.coor.crd")
    crd_path_abs = os.path.abspath(crd_file)
    
    # Use multi-line string to keep all commands together in one script block
    script = f"""
    OPEN WRITE UNIT 70 CARD NAME {crd_path_abs}
    WRITE COOR CARD UNIT 70
    CLOSE UNIT 70
    """
    cs(script)
    print(f"  Wrote CRD: {crd_file}")
    
    # Convert CRD to PDB using crd2pdb tool (workaround for Fortran EOF issue)
    pdb_file = os.path.join("out", f"{prefix}.coor.pdb")
    crd2pdb_script = os.path.join(
        os.path.dirname(__file__), "..", "..", "..", "src", "engines", "pycharmm", "crd2pdb.py"
    )
    if os.path.exists(crd2pdb_script):
        try:
            result = subprocess.run(
                [sys.executable, crd2pdb_script, "--psf", psf_path, "--crd", crd_file, "--pdb", pdb_file],
                check=True,
                capture_output=True,
                text=True
            )
            print(f"  Wrote PDB: {pdb_file}")
        except subprocess.CalledProcessError as e:
            print(f"  Warning: PDB conversion failed: {e.stderr}")


if __name__ == "__main__":
    # ======= ここを自分の系に合わせて設定 =======
    PSF   = "pentacene.psf"        # 既存 PSF
    COOR  = "pentacene.crd"        # Wrapped coordinates (recommended for PBC)
    #PDB   = "pentacene.pdb"         # CRD使用時にセル情報を取るPDB (必要なら指定)

    TOPPAR_STREAM = "pentacene.str"  # CGenFF topology/parameter stream file

    # ボックス (Å) - AUTOならCRD/PDBから自動読み込み
    BOX = None  # None = 自動検出
    PBC = "AUTO"  # AUTO = CRD/PDBから自動判定

    # 最小化
    MIN_SD   = 200   # Quick test
    MIN_ABNR = 300

    # NVT 条件
    TEMP_K   = 300.0
    NVT_STEPS = 20000  # Quick test: 0.4 ps @ 2 fs
    DT_PS     = 0.001

    # NPT 条件
    PRESS_ATM = 1.01325  # ~1 atm
    NPT_STEPS = 10  # Minimal test to check if NPT runs at all   # Quick test: 0.4 ps @ 2 fs

    # ======= 実行シーケンス =======
    setup_system(
        psf_path=PSF,
        coor_path=COOR,
        toppar_stream=TOPPAR_STREAM,
        box=BOX,
        pbc_type=PBC,
        use_pme=True,
        cutnb=14.0, ctofnb=12.0, ctonnb=10.0,
        pme_kappa=0.34, pme_order=6,
        pme_grid=None,  # 自動見積もり (用途により明示指定を推奨)
    )

    print("==> Setting up SHAKE constraints for hydrogen bonds")
    setup_shake(shake_hyd=True, shake_fast=False, tol=1e-8)

    print("==> Minimization: SD -> ABNR")
    minimize_energy(nsd=MIN_SD, nabnr=MIN_ABNR)
    
    print("==> Writing minimized structure")
    try:
        write_snapshots(prefix="min", psf_path=PSF)
        print("==> Successfully wrote min.coor.crd and min.coor.pdb")
    except Exception as e:
        print(f"Warning: Error writing snapshots: {e}")
        print("Continuing with dynamics...")

    print("==> NVT (Nosé–Hoover)")
    run_nvt(
        temp=TEMP_K,
        nsteps=NVT_STEPS,
        dt_ps=DT_PS,
        traj_path="out/nvt.dcd",
        rst_path="out/nvt.restart",
        nsavc=50,  # Save every 50 steps
        nsavr=5000,
        iasvel=1,
        seed=20251111,
        tmass=1000.0,
        use_shake=True,
    )
    
    print("==> Writing NVT structure")
    try:
        write_snapshots(prefix="after_nvt", psf_path=PSF)
        print("==> Successfully wrote after_nvt.coor.crd and after_nvt.coor.pdb")
    except Exception as e:
        print(f"Warning: Error writing snapshots: {e}")
        print("Continuing with NPT...")

    print("==> NPT (Hoover + Langevin piston)")
    run_npt(
        temp=TEMP_K,
        press_atm=PRESS_ATM,
        nsteps=NPT_STEPS,
        dt_ps=DT_PS,
        traj_path="out/npt.dcd",
        rst_path="out/npt.restart",
        nsavc=1,  # Save every step for debugging
        nsavr=5000,
        tmass=1000.0,
        pmass=400.0,
        pgamma=20.0,
        iasvel=0,  # Don't assign initial velocities (continue from NVT)
        seed=20251114,
        use_shake=True,
    )
    
    print("==> Writing NPT structure")
    try:
        write_snapshots(prefix="after_npt", psf_path=PSF)
        print("==> Successfully wrote after_npt.coor.crd and after_npt.coor.pdb")
    except Exception as e:
        print(f"Warning: Error writing snapshots: {e}")

    print("All done.")
