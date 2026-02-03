import pytest

from martini_mapper.errors import InvalidSmilesError
from martini_mapper.utils import ensure_nonempty, safe_mol_from_smiles


def test_ensure_nonempty_rejects_empty():
    with pytest.raises(InvalidSmilesError):
        ensure_nonempty("", field="smiles")
    with pytest.raises(InvalidSmilesError):
        ensure_nonempty("   ", field="smiles")


def test_safe_mol_from_smiles_rejects_invalid():
    with pytest.raises(InvalidSmilesError):
        safe_mol_from_smiles("this_is_not_smiles")


def test_safe_mol_from_smiles_accepts_valid():
    mol = safe_mol_from_smiles("c1ccccc1")
    assert mol.GetNumAtoms() == 6
