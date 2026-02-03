from __future__ import annotations

from pathlib import Path
import shutil
from typing import Iterable, Sequence, Optional

from rdkit import Chem

from .errors import InvalidSmilesError, MappingError, OutputError


def ensure_nonempty(value: str, *, field: str = "value") -> str:
    if value is None:
        raise InvalidSmilesError(f"{field} cannot be None")
    if not isinstance(value, str):
        raise InvalidSmilesError(f"{field} must be a string, got {type(value)!r}")
    v = value.strip()
    if not v:
        raise InvalidSmilesError(f"{field} cannot be empty")
    return v


def clean_smiles(smiles: str) -> str:
    """Apply project-specific SMILES cleanup.

    - removes '*' wildcards used in some datasets
    - fixes duplicated ring numbers (project-specific helper)
    """
    s = ensure_nonempty(smiles, field="smiles")
    try:
        from .setup_mapping import duplicate_ring_number_handler, remove_stars
    except Exception as e:  # pragma: no cover
        raise OutputError(f"Could not import SMILES cleanup utilities: {e}") from e

    s = remove_stars(s)
    s = duplicate_ring_number_handler(s)
    return s


def safe_mol_from_smiles(smiles: str) -> Chem.Mol:
    """Parse SMILES safely and raise a clear exception on failure."""
    s = ensure_nonempty(smiles, field="smiles")
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        raise InvalidSmilesError(
            "Invalid SMILES string: RDKit could not parse it. "
            "Check for typos and ensure you are passing a canonical SMILES."
        )
    return mol


def validate_final_mapping(final: Sequence[str]) -> None:
    """Validate that every atom has a bead assigned."""
    if final is None:
        raise MappingError("final mapping is None")
    if not isinstance(final, (list, tuple)):
        raise MappingError(f"final mapping must be a list/tuple, got {type(final)!r}")
    missing = [i for i, v in enumerate(final) if not isinstance(v, str) or v.strip() == ""]
    if missing:
        preview = ", ".join(map(str, missing[:20]))
        extra = "" if len(missing) <= 20 else f" (and {len(missing)-20} more)"
        raise MappingError(f"Unmapped atom indices: {preview}{extra}")


def organize_outputs(
    compound_name: str,
    *,
    out_dir: Optional[Path] = None,
    keep_intermediates: bool = False,
) -> Path:
    """Move final outputs into a folder and optionally delete intermediates."""
    compound_name = ensure_nonempty(compound_name, field="name")
    out_dir = out_dir or Path(compound_name)
    out_dir.mkdir(exist_ok=True)

    keep = {
        f"{compound_name}.gro",
        f"{compound_name}.itp",
        f"{compound_name}.txt",
    }

    for fname in keep:
        src = Path(fname)
        if src.exists():
            shutil.move(str(src), str(out_dir / src.name))

    if keep_intermediates:
        return out_dir

    trash_files = [
        f"{compound_name}.xyz",
        f"{compound_name}.pdb",
        f"ref_{compound_name}.gro",
        f"ref_{compound_name}.xtc",
        f"{compound_name}_cg_ref.gro",
        f"{compound_name}_cg_ref.xtc",
        f"{compound_name}_bonded_report.txt",
    ]
    for f in trash_files:
        p = Path(f)
        if p.exists():
            try:
                p.unlink()
            except Exception as e:  # pragma: no cover
                raise OutputError(f"Failed to delete intermediate file {p}: {e}") from e

    trash_dirs = [f"{compound_name}_md"]
    for d in trash_dirs:
        p = Path(d)
        if p.exists() and p.is_dir():
            try:
                shutil.rmtree(p)
            except Exception as e:  # pragma: no cover
                raise OutputError(f"Failed to remove intermediate dir {p}: {e}") from e

    return out_dir
