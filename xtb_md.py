import subprocess
from pathlib import Path
import shutil
from typing import Dict, Any

def _write_md_input(
    workdir: Path,
    *,
    temp_k: float,
    time_ps: float,
    include_knobs: bool = True,
    step_fs: float = 1.0,
    shake: int = 2,
    hmass: int = 4,
    sccacc: float = 2.0,
    nvt: bool = True,
    dump: int = 50,
) -> Path:
    """Create an xTB $md input file (md.inp).

    Notes
    -----
    - Many xTB versions read MD controls like temperature/time from the $md block.
      Passing --time/--temp on the command line may be ignored or treated as unknown options.
    - When include_knobs=False, we still write temp/time (and dump) but omit extra stability knobs.
    """
    md_inp = workdir / "md.inp"

    lines = [
        "$md",
        f"  temp={temp_k}",
        f"  time={time_ps}",
        f"  dump={dump}",
    ]

    if include_knobs:
        lines.extend([
            f"  step={step_fs}",
            f"  shake={shake}",
            f"  hmass={hmass}",
            f"  sccacc={sccacc}",
            f"  nvt={'true' if nvt else 'false'}",
        ])

    lines.extend(["$end", ""])
    md_inp.write_text("\n".join(lines))
    return md_inp

def run_xtb_md(
    xyz_path: str,
    workdir: str = "xtb_md_run",
    temp_k: float = 300.0,
    time_ps: float = 10.0,
    charge: int = 0,
    uhf: int = 0,
    *,
    use_md_input: bool = True,
    step_fs: float = 0.5,
    shake: int = 2,
    hmass: int = 4,
    sccacc: float = 1.0,
    nvt: bool = True,
    dump: int = 200,
    preopt: bool = False,
) -> Dict[str, Any]:
   """Run xTB MD starting from an XYZ file.

   Outputs are written into *workdir* (e.g., xtb.trj, xtb.log, xtb_stdout.log).

   Important
   ---------
   - This runner controls MD settings via an $md input file (md.inp). In many xTB builds,
     temperature and simulation time are *not* accepted as --temp/--time flags and must be
     provided in $md instead.

   Parameters
   ----------
   xyz_path:
       Input XYZ file.
   workdir:
       Output directory for xTB.
   temp_k, time_ps:
       Thermostat temperature (K) and simulated time (ps). Written into md.inp.
   use_md_input:
       If True, include extra stability knobs in md.inp (step/shake/hmass/sccacc/nvt).
       If False, still writes md.inp for temp/time but omits extra knobs.
   preopt:
       If True, run `xtb <xyz> --opt` first and start MD from xtbopt.xyz.

   Returns
   -------
   dict with paths to outputs (trajectory, logs, inputs).
   """
   xtb_exe = shutil.which("xtb")
   if xtb_exe is None:
       raise RuntimeError(
           "xtb executable not found in PATH. "
           "Install via conda: conda install -c conda-forge xtb"
       )

   xyz_path = Path(xyz_path).resolve()
   if not xyz_path.exists():
       raise FileNotFoundError(f"XYZ not found: {xyz_path}")

   wd = Path(workdir).resolve()
   wd.mkdir(parents=True, exist_ok=True)

   # Copy the input xyz into workdir to keep xtb outputs local
   local_xyz = wd / xyz_path.name
   local_xyz.write_text(xyz_path.read_text())

   # Always write md.inp so temp/time are honored; optionally include extra knobs.
   md_inp = _write_md_input(
       wd,
       temp_k=temp_k,
       time_ps=time_ps,
       include_knobs=use_md_input,
       step_fs=step_fs,
       shake=shake,
       hmass=hmass,
       sccacc=sccacc,
       nvt=nvt,
       dump=dump,
   )

   # Optional pre-optimization step
   md_start_xyz = local_xyz
   if preopt:
       opt_cmd = [
           xtb_exe,
           local_xyz.name,
           "--chrg", str(charge),
           "--uhf", str(uhf),
           "--opt",
       ]
       opt_log_path = wd / "xtb_opt_stdout.log"
       opt_err_path = wd / "xtb_opt_stderr.log"
       with open(opt_log_path, "w") as logf, open(opt_err_path, "w") as errf:
           opt_res = subprocess.run(opt_cmd, cwd=wd, stdout=logf, stderr=errf, text=True)
       if opt_res.returncode != 0:
           raise RuntimeError(
               f"xtb optimization failed (return code {opt_res.returncode}). "
               f"See logs:\n  {opt_log_path}\n  {opt_err_path}"
           )
       xtbopt = wd / "xtbopt.xyz"
       if not xtbopt.exists():
           raise RuntimeError(
               f"xtb optimization finished but xtbopt.xyz not found: {xtbopt}\n"
               f"Check logs:\n  {opt_log_path}\n  {opt_err_path}"
           )
       md_start_xyz = xtbopt

   # MD command (time/temp are controlled via md.inp)
   cmd = [
       xtb_exe,
       md_start_xyz.name,
       "--chrg", str(charge),
       "--uhf", str(uhf),
       "--input", md_inp.name,
       "--md",
   ]

   # Run and capture stdout/stderr into files for debugging
   log_path = wd / "xtb_stdout.log"
   err_path = wd / "xtb_stderr.log"
   with open(log_path, "w") as logf, open(err_path, "w") as errf:
       result = subprocess.run(
           cmd,
           cwd=wd,
           stdout=logf,
           stderr=errf,
           text=True
       )

   if result.returncode != 0:
       raise RuntimeError(
           f"xtb MD failed (return code {result.returncode}). "
           f"See logs:\n  {log_path}\n  {err_path}"
       )

   # Typical xtb outputs
   trj = wd / "xtb.trj"
   main_log = wd / "xtb.log"

   if not trj.exists():
       raise RuntimeError(
           f"xtb finished but trajectory not found: {trj}\n"
           f"Check logs:\n  {log_path}\n  {err_path}"
       )

   return {
       "workdir": str(wd),
       "trajectory_xyz": str(trj),
       "xtb_log": str(main_log) if main_log.exists() else None,
       "stdout_log": str(log_path),
       "stderr_log": str(err_path),
       "input_xyz": str(local_xyz),
       "md_start_xyz": str(md_start_xyz),
       "md_input": str(md_inp),
   }
