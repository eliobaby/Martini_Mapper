import os
import subprocess

base_dir = os.getcwd()
mapper_dir = os.path.abspath("..")

failed_report = os.path.join(base_dir, "not_mapped.txt")

with open(failed_report, "w") as report:

    for mol_name in sorted(os.listdir(base_dir)):
        mol_path = os.path.join(base_dir, mol_name)

        if not os.path.isdir(mol_path):
            continue

        smiles_path = os.path.join(mol_path, f"{mol_name}.smiles")
        if not os.path.exists(smiles_path):
            continue

        with open(smiles_path, "r") as f:
            smiles = f.readline().strip()

        if not smiles:
            report.write(f"{mol_name} : empty SMILES\n")
            continue

        print(f"Running martini_mapper for {mol_name}: {smiles}")

        log_file = os.path.join(mol_path, "martini_mapper.log")

        try:
            with open(log_file, "w") as log:
                subprocess.run(
                    [
                        "python",
                        "-m",
                        "martini_mapper",
                        mol_name,
                        smiles,
                    ],
                    cwd=mol_path,
                    env={**os.environ, "PYTHONPATH": mapper_dir},
                    stdout=log,
                    stderr=log,
                    check=True
                )

        except subprocess.CalledProcessError:
            report.write(
                f"{mol_name} : not mapped, for error look at {mol_name}/martini_mapper.log\n"
            )
            print(f"⚠️  {mol_name} failed — logged")
            continue
