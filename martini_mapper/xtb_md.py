import subprocess
from pathlib import Path
import shutil
from typing import Dict, Any
import os

def run_xtb_md(
    xyz_path: str,
    workdir: str = "xtb_md_run",
    temp_k: float = 298.15,
    time_ps: float = 20.0,
    *,
    step_fs: float = 0.1,
    shake: int = 0,
    hmass: int = 4,
    nvt: bool = True,
    dump: int = 50,
) -> Dict[str, Any]:

    xtb_exe = shutil.which("xtb")
    if xtb_exe is None:
        raise RuntimeError("xtb executable not found in PATH.")

    xyz_path = Path(xyz_path).resolve()
    if not xyz_path.exists():
        raise FileNotFoundError(f"XYZ not found: {xyz_path}")

    wd = Path(workdir).resolve()
    wd.mkdir(parents=True, exist_ok=True)

    # copy input xyz into workdir
    local_xyz = wd / xyz_path.name
    local_xyz.write_text(xyz_path.read_text())

    # write md.inp (always)
    md_inp = wd / "md.inp"
    md_inp.write_text(
        "\n".join([
            "$md",
            f"  temp={temp_k}",
            f"  time={time_ps}",
            f"  dump={dump}",
            f"  step={step_fs}",
            f"  shake={shake}",
            f"  hmass={hmass}",
            f"  nvt={'true' if nvt else 'false'}",
            "$end",
            "",
        ])
    )

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_STACKSIZE"] = "2G" 
    
    # --------------------------
    # 1) OPTIMIZATION: --opt normal
    # --------------------------
    opt_cmd = [
        xtb_exe,
        local_xyz.name,
        "--gfnff",
        "--chrg", "0",
        "--uhf", "0",
        "--opt", "normal",
    ]

    opt_out = wd / "xtb_opt_stdout.log"
    opt_err = wd / "xtb_opt_stderr.log"
    with open(opt_out, "w") as logf, open(opt_err, "w") as errf:
        opt_res = subprocess.run(opt_cmd, cwd=wd, stdout=logf, stderr=errf, text=True, env=env)

    if opt_res.returncode != 0:
        raise RuntimeError(
            f"xtb optimization failed (return code {opt_res.returncode}).\n"
            f"See:\n  {opt_out}\n  {opt_err}"
        )

    # xTB writes optimized geometry to xtbopt.xyz for xyz input
    opt_xyz = wd / "xtbopt.xyz"
    if not opt_xyz.exists():
        raise RuntimeError(
            "Optimization finished but xtbopt.xyz not found.\n"
            f"Expected: {opt_xyz}\n"
            f"See:\n  {opt_out}\n  {opt_err}"
        )

    # --------------------------
    # 2) MD: --md using xtbopt.xyz
    # --------------------------
    md_cmd = [
        xtb_exe,
        opt_xyz.name,              # start MD from optimized geometry
        "--gfnff",
        "--norestart",
        "--chrg", "0",
        "--uhf", "0",
        "--input", md_inp.name,
        "--md",
    ]

    md_out = wd / "xtb_md_stdout.log"
    md_err = wd / "xtb_md_stderr.log"
    with open(md_out, "w") as logf, open(md_err, "w") as errf:
        md_res = subprocess.run(md_cmd, cwd=wd, stdout=logf, stderr=errf, text=True, env=env)

    if md_res.returncode != 0:
        raise RuntimeError(
            f"xtb MD failed (return code {md_res.returncode}).\n"
            f"See:\n  {md_out}\n  {md_err}"
        )

    trj = wd / "xtb.trj"
    if not trj.exists():
        raise RuntimeError(
            f"xtb finished but trajectory not found: {trj}\n"
            f"See:\n  {md_out}\n  {md_err}"
        )

    return {
        "workdir": str(wd),
        "opt_xyz": str(opt_xyz),
        "trajectory_xyz": str(trj),
        "stdout_log_opt": str(opt_out),
        "stderr_log_opt": str(opt_err),
        "stdout_log_md": str(md_out),
        "stderr_log_md": str(md_err),
        "input_xyz": str(local_xyz),
        "md_input": str(md_inp),
    }
