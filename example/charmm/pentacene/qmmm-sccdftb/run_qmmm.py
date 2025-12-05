#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pyCHARMM で QM/MM シミュレーションを実行する実装。

✔ QM領域: 中央の2分子 (residue 27, 28)
✔ QMプログラム: SCC-DFTB (Self-Consistent Charge DFTB)
✔ MM領域: 残りの52分子 (CHARMM36 CGenFF)
✔ 手順:
  1) システムセットアップ (PBC + PME)
  2) QM領域の定義 (SCC-DFTB)
  3) エネルギー最小化
  4) QM/MM NVT シミュレーション

run_md.py をベースに QM/MM 機能を追加。

Note: SCC-DFTB はCHARMM内蔵で、PBC + PME と完全互換。
      BOMLEV -2 で大きな分子のGROUP extent警告を抑制。
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
except Exception as e:
    from pycharmm import charmm_script as cs  # type: ignore


def open_unit(unit: int, mode: str, path: str, file_type: str = "CARD"):
    """Open a file unit in CHARMM."""
    path = os.path.abspath(path)
    cs(f"OPEN {mode.upper()} {file_type.upper()} UNIT {unit:d} NAME {path}")


def find_valid_fft_dimension(n: int) -> int:
    """Find the smallest integer >= n that has only 2, 3, 5 as prime factors (for PME FFT)."""
    if n <= 1:
        return 2
    candidate = n
    while True:
        temp = candidate
        for factor in [2, 3, 5]:
            while temp % factor == 0:
                temp //= factor
        if temp == 1:
            return candidate
        candidate += 1


def _read_crd_box(path: str) -> tuple[tuple[float,float,float], tuple[float,float,float]] | None:
    """CRDファイルの最終行からボックス情報を返す。"""
    try:
        with open(path, "rt", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
        last_line = lines[-1].strip()
        parts = last_line.split()
        if len(parts) == 6:
            a, b, c, alpha, beta, gamma = map(float, parts)
            if all(x > 0 for x in [a, b, c]) and all(0 < ang < 180 for ang in [alpha, beta, gamma]):
                return (a, b, c), (alpha, beta, gamma)
    except Exception:
        pass
    return None


def _guess_crystal_type(abc: tuple[float,float,float], ang: tuple[float,float,float]) -> str:
    """角度と格子長から CHARMM の型名を推定。"""
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
    box: tuple[float, float, float] | None = None,
    pbc_type: str = "AUTO",
    use_pme: bool = True,
    cutnb: float = 14.0,
    ctofnb: float = 12.0,
    ctonnb: float = 10.0,
    pme_kappa: float = 0.34,
    pme_order: int = 6,
    pme_grid: tuple[int, int, int] | None = None,
):
    """PSF/座標を読み、PBC+PME をセットアップ。"""
    cs("BOMLEV -2")
    cs("PRNLEV 5")

    # Load base CGenFF force field
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
    
    cs("BOMLEV 0")

    # Read PSF & coordinates
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
    cs("COOR STAT")

    # Box information
    cryst = _read_crd_box(coor_path) if not is_pdb else None
    
    if cryst is not None:
        abc, ang = cryst
        pbc_type_eff = _guess_crystal_type(abc, ang)
        print(f"Auto-detected box: {abc[0]:.3f} {abc[1]:.3f} {abc[2]:.3f} Å, {ang[0]:.2f}° {ang[1]:.2f}° {ang[2]:.2f}°")
        print(f"Crystal type: {pbc_type_eff}")
    else:
        if box is None:
            raise ValueError("box parameter required or use AUTO detection")
        abc = box
        ang = (90.0, 90.0, 90.0)
        pbc_type_eff = pbc_type if pbc_type != "AUTO" else "ORTHorhombic"

    ax, by, cz = abc
    alpha, beta, gamma = ang
    cs(f"CRYSTAL DEFINE {pbc_type_eff} {ax:.3f} {by:.3f} {cz:.3f} {alpha:.2f} {beta:.2f} {gamma:.2f}")
    cs("CRYSTAL BUILD NOPER 0")
    cs("IMAGE BYRES XCEN 0.0 YCEN 0.0 ZCEN 0.0 SELECT ALL END")

    # PME setup
    if use_pme:
        if pme_grid is None:
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


def setup_nbond_for_qmmm(cutnb: float = 14.0, ctofnb: float = 12.0, ctonnb: float = 10.0):
    """Setup nonbond parameters for QM/MM with PBC.
    
    Use GROUP-based cutoff for QM/MM with IMAGE.
    BOMLEV -2 suppresses group extent warnings for large molecules.
    """
    print("==> Setting up GROUP-based cutoff for QM/MM with PBC")
    cs("BOMLEV -2")  # Suppress group extent warnings
    cs(
        "NBOND GROUP CDIEL SWITCH VGROUP VDIST VSWITCH "
        f"CUTNB {cutnb:.2f} CTOFNB {ctofnb:.2f} CTONNB {ctonnb:.2f}"
    )
    cs("ENERGY")
    cs("BOMLEV 0")  # Restore normal warning level


def setup_qmmm(
    qm_residues: list[int],
    qm_charge: int = 0,
    qm_mult: int = 1,
):
    print("==> Setting up SCC-DFTB QM/MM")

    # QM領域 selection（resid 27 .or. resid 28）
    qm_selection = " .or. ".join([f"resid {resid}" for resid in qm_residues])

    # QMREG 定義（あとで使う用）
    cs(f"DEFINE QMREG SELE {qm_selection} END")

    # SCCDFTB 初期化
    # chrg: QM電荷, spin: 2S (= 2*mult-2; 一重項なら 0)
    spin = 2*qm_mult - 2
    cmd = (
        f"SCCDFTB CHRG {qm_charge} SPIN {spin} REMO "
        f"SELE {qm_selection} END "
        f"CUTF PME UPDT 1"
    )
    # CUTF PME UPDT 1: QM/MM静電をPMEで扱い、セル更新あり（必要に応じて）

    print("  CHARMM command:")
    print("  ", cmd)
    cs(cmd)

    # QM領域の確認
    cs("COOR PRINT SELE QMREG END")

    print("==> SCC-DFTB QM/MM setup complete")
    cs("ENERGY")  # 初期 QM/MM エネルギーを計算



def minimize_energy_qmmm(nsd: int = 100, nabnr: int = 200):
    """QM/MM系でのエネルギー最小化。
    
    Note: QM/MM計算は重いため、ステップ数は通常のMMより少なめに設定。
    """
    print(f"==> QM/MM minimization: SD ({nsd} steps) -> ABNR ({nabnr} steps)")
    cs(f"MINIMIZE SD NSTEP {nsd:d} NPRINt 50")
    cs(f"MINIMIZE ABNR NSTEP {nabnr:d} NPRINt 50")
    cs("ENERGY")


def run_nvt_qmmm(
    temp: float,
    nsteps: int,
    dt_ps: float,
    traj_path: str = "qmmm_nvt.dcd",
    nsavc: int = 100,
    iasvel: int = 1,
    seed: int = 314159,
    tmass: float = 1000.0,
    use_shake: bool = False,  # SHAKEはQM/MMでは通常使わない
):
    """QM/MM NVT dynamics。
    
    Note: QM/MM計算は重いため、保存頻度は低め (nsavc大きめ) を推奨。
          SHAKEはQM領域に適用できないため、デフォルトOFF。
    """
    print(f"==> QM/MM NVT dynamics: {nsteps} steps, dt={dt_ps} ps")
    
    shake_flag = "SHAKE" if use_shake else ""

    dyn_io = ""
    if traj_path:
        cs(f"OPEN WRITE UNFORM UNIT 51 NAME {os.path.abspath(traj_path)}")
        dyn_io += f"IUNCRD 51 NSAVC {nsavc:d} "

    cs(
        f"DYNAMICS CPT START {shake_flag} "
        f"NSTEPS {nsteps:d} TIMESTEP {dt_ps:.6f} "
        f"{dyn_io}NTRFRQ 100 "
        f"IASORS 1 IASVEL {iasvel:d} ISEED {seed:d} "
        f"HOOVER TMASS {tmass:.1f} REFT {temp:.1f} "
    )
    
    print(f"==> QM/MM NVT complete")


def write_snapshots(prefix: str = "out", psf_path: str = "pentacene.psf"):
    """Write snapshot as CRD and PDB."""
    os.makedirs("out", exist_ok=True)
    
    crd_file = os.path.join("out", f"{prefix}.coor.crd")
    crd_path_abs = os.path.abspath(crd_file)
    
    script = f"""
    OPEN WRITE UNIT 70 CARD NAME {crd_path_abs}
    WRITE COOR CARD UNIT 70
    CLOSE UNIT 70
    """
    cs(script)
    print(f"  Wrote CRD: {crd_file}")
    
    # Convert to PDB
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
    # ======= QM/MM設定 =======
    PSF   = "pentacene.psf"
    COOR  = "pentacene.crd"
    TOPPAR_STREAM = "pentacene.str"

    # QM領域: 中央の2分子 (54分子中の27, 28番目)
    QM_RESIDUES = [27, 28]
    QM_METHOD = "SCCDFTB"       # Self-Consistent Charge DFTB
    QM_PROGRAM = "SCCDFTB"      # QM program (CHARMM internal SCC-DFTB)
    QM_CHARGE = 0            # 中性
    QM_MULT = 1              # 一重項

    # ボックス設定
    BOX = None  # Auto-detect from CRD
    PBC = "AUTO"

    # 最小化 (QM/MMは重いので少なめ)
    MIN_SD   = 50
    MIN_ABNR = 100

    # NVT条件 (QM/MMは重いので短め)
    TEMP_K    = 300.0
    NVT_STEPS = 500   # 1 ps @ 2 fs
    DT_PS     = 0.002

    print("="*70)
    print("QM/MM Simulation with pyCHARMM + SCC-DFTB")
    print("="*70)
    print(f"QM region: {len(QM_RESIDUES)} molecules (residues {QM_RESIDUES})")
    print(f"QM method: {QM_METHOD} (Self-Consistent Charge DFTB)")
    print(f"MM region: {54 - len(QM_RESIDUES)} molecules (CHARMM36 CGenFF)")
    print("="*70)
    print("Note: SCC-DFTB is built into CHARMM - no external QM program needed")
    print("      SCC-DFTB provides high accuracy with PBC+PME support")
    print("="*70)

    # ======= 実行シーケンス =======
    
    # 1) システムセットアップ
    print("\n### Step 1: System setup (MM force field + PBC + PME)")
    print("Note: SCC-DFTB supports PBC and PME")
    setup_system(
        psf_path=PSF,
        coor_path=COOR,
        toppar_stream=TOPPAR_STREAM,
        box=BOX,
        pbc_type=PBC,
        use_pme=True,  # PME enabled for SCC-DFTB
        cutnb=14.0, ctofnb=12.0, ctonnb=10.0,
        pme_kappa=0.34, pme_order=6,
        pme_grid=None,
    )

    # 2) QM/MM用のNBOND再設定 (GROUP-based cutoff with PBC/IMAGE)
    print("\n### Step 2a: Reconfigure nonbonds for QM/MM")
    setup_nbond_for_qmmm(cutnb=14.0, ctofnb=12.0, ctonnb=10.0)
    
    # 3) QM/MMセットアップ
    print("\n### Step 2b: QM/MM setup")
    """
    setup_qmmm(
        qm_residues=QM_RESIDUES,
        qm_method=QM_METHOD,
        qm_program=QM_PROGRAM,
        qm_charge=QM_CHARGE,
        qm_mult=QM_MULT,
        boundary="DIV",  # Dividing scheme for QM/MM boundary
    )
    """
    setup_qmmm(
        qm_residues=QM_RESIDUES,
        qm_charge=QM_CHARGE,
        qm_mult=QM_MULT,
    )

    # 4) エネルギー最小化
    print("\n### Step 3: QM/MM energy minimization")
    minimize_energy_qmmm(nsd=MIN_SD, nabnr=MIN_ABNR)
    
    print("\n### Step 3b: Writing minimized structure")
    try:
        write_snapshots(prefix="qmmm_min", psf_path=PSF)
        print("==> Successfully wrote qmmm_min files")
    except Exception as e:
        print(f"Warning: Error writing snapshots: {e}")

    # 5) NVT dynamics
    print("\n### Step 4: QM/MM NVT dynamics")
    run_nvt_qmmm(
        temp=TEMP_K,
        nsteps=NVT_STEPS,
        dt_ps=DT_PS,
        traj_path="out/qmmm_nvt.dcd",
        nsavc=100,  # Save every 100 steps (0.2 ps intervals)
        iasvel=1,
        seed=20251114,
        tmass=1000.0,
        use_shake=False,  # Don't use SHAKE with QM/MM
    )
    
    print("\n### Step 4b: Writing final structure")
    try:
        write_snapshots(prefix="qmmm_after_nvt", psf_path=PSF)
        print("==> Successfully wrote qmmm_after_nvt files")
    except Exception as e:
        print(f"Warning: Error writing snapshots: {e}")

    print("\n" + "="*70)
    print("QM/MM simulation complete!")
    print("="*70)
    print("Output files:")
    print("  - out/qmmm_min.coor.{crd,pdb}        : Minimized structure")
    print("  - out/qmmm_nvt.dcd                   : NVT trajectory")
    print("  - out/qmmm_after_nvt.coor.{crd,pdb}  : Final structure")
    print("="*70)
