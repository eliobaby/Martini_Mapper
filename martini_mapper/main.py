from __future__ import annotations

import argparse
import copy
from pathlib import Path
from typing import Optional, Tuple

from .errors import ExternalToolError, InvalidSmilesError, MappingError, OutputError
from .utils import (
    clean_smiles,
    ensure_nonempty,
    organize_outputs,
    safe_mol_from_smiles,
)

# Core mapping imports (keep these eager; they are the library's main feature)
from .setup_mapping import get_atom_properties, connectivity_matrix
from .martini_3_dictionary import get_m3_dict
from .mapping_scheme import map_molecule, parse_smiles
from .algorithm import map_martini_beads


def run_mapping(
    compound_name: str,
    smiles: str,
    *,
    run_xtb: bool = True,
    write_files: bool = True,
    out_dir: Optional[Path] = None,
    keep_intermediates: bool = False,
) -> Tuple[list[str], list]:
    """Run the full Martini mapping pipeline.

    Returns (final_beads, mapping_structure).

    Parameters
    ----------
    run_xtb:
        If True, attempt to run the xtb-based coordinate pipeline.
        If xtb is not installed, raises ExternalToolError.
    write_files:
        If True, write .txt/.gro/.itp files via text_file helpers.
    out_dir:
        If provided, move the kept output files into this directory (or into
        ./<compound_name> if None).
    keep_intermediates:
        If False (default), remove the intermediate files/directories that
        xtb/MD steps tend to generate.
    """
    compound_name = ensure_nonempty(compound_name, field="name")
    smiles = ensure_nonempty(smiles, field="smiles")

    smiles_clean = clean_smiles(smiles)

    mol = safe_mol_from_smiles(smiles_clean)

    atom_properties = get_atom_properties(mol)
    if not atom_properties:
        raise MappingError("No atoms detected after RDKit parsing; cannot map.")

    final = ["" for _ in atom_properties]
    conn_mat = connectivity_matrix(mol, len(atom_properties))

    bead_dict = get_m3_dict()
    smiles_token = parse_smiles(smiles_clean)
    mapping = map_molecule(smiles_token, conn_mat, atom_properties)
    mapping_copy = copy.deepcopy(mapping)

    try:
        final = map_martini_beads(mapping, final, bead_dict)
    except Exception as e:
        # normalize to MappingError for callers
        if isinstance(e, MappingError):
            raise
        raise MappingError(str(e)) from e

    # Write mapping output (optional)
    if write_files:
        try:
            from .text_file import export_bead_mapping, gro_maker, itp_maker
        except Exception as e:  # pragma: no cover
            raise OutputError(f"Failed to import output writers: {e}") from e

        beads_data, connections_data = export_bead_mapping(final, mapping_copy, smiles_clean, compound_name)
        gro_maker(beads_data, compound_name)
        itp_maker(beads_data, connections_data, compound_name)

    # xtb pipeline (optional)
    if run_xtb:
        try:
            from .xtb import smiles_to_ref
            from .cg_from_xtb_pipeline import build_cg_from_xtb
        except Exception as e:
            raise ExternalToolError(
                "xtb pipeline requested but required modules could not be imported. "
                "Install xtb-python (conda-forge) or run with --no-xtb."
            ) from e

        try:
            smiles_to_ref(compound_name, smiles_clean)
            build_cg_from_xtb(
                smiles=smiles_clean,
                final=final,
                mapping=mapping,
                ref_gro=f"ref_{compound_name}.gro",
                ref_xtc=f"ref_{compound_name}.xtc",
                compound=compound_name,
                out_prefix=compound_name,
                T=300.0,
            )
        except Exception as e:
            raise ExternalToolError(f"xtb pipeline failed: {e}") from e

    # organize outputs (optional)
    if write_files:
        try:
            organize_outputs(
                compound_name,
                out_dir=out_dir,
                keep_intermediates=keep_intermediates,
            )
        except Exception as e:
            raise OutputError(f"Failed to organize outputs: {e}") from e

    return final, mapping


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="martini_mapper")
    p.add_argument("name", nargs="?", help="Compound name (e.g., Benzene)")
    p.add_argument("smiles", nargs="?", help="SMILES string (e.g. c1ccccc1)")
    p.add_argument("--no-xtb", action="store_true", help="Skip xtb coordinate/CG generation.")
    p.add_argument("--no-files", action="store_true", help="Do not write .txt/.gro/.itp outputs (library mode).")

    p.add_argument(
        "--out-dir",
        default=None,
        help="Directory to place final outputs (default: ./<name>/)",
    )
    p.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Do not delete intermediate xtb/MD files.",
    )
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    compound_name = args.name if args.name is not None else input("Name: ")
    smiles = args.smiles if args.smiles is not None else input("Smiles: ")

    out_dir = Path(args.out_dir) if args.out_dir else None

    try:
        run_mapping(
            compound_name,
            smiles,
            run_xtb=not args.no_xtb,
            write_files=not args.no_files,
            out_dir=out_dir,
            keep_intermediates=args.keep_intermediates,
        )
    except (InvalidSmilesError, MappingError, ExternalToolError, OutputError) as e:
        # Human-friendly error for CLI usage
        raise SystemExit(str(e)) from e

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
