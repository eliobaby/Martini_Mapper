import pytest

from martini_mapper.errors import InvalidSmilesError
from martini_mapper.main import run_mapping


@pytest.mark.parametrize(
    "name,smiles",
    [
        ("Benzene", "c1ccccc1"),
        ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
    ],
)
def test_run_mapping_library_mode(name, smiles):
    final, mapping = run_mapping(
        name,
        smiles,
        run_xtb=False,
        write_files=False,
    )
    assert isinstance(final, list)
    assert len(final) > 0
    assert all(isinstance(x, str) and x for x in final)  # every atom mapped to a bead string
    assert isinstance(mapping, list)


def test_run_mapping_rejects_empty_smiles():
    with pytest.raises(InvalidSmilesError):
        run_mapping("X", "", run_xtb=False, write_files=False)


def test_run_mapping_rejects_invalid_smiles():
    with pytest.raises(InvalidSmilesError):
        run_mapping("X", "not_a_smiles", run_xtb=False, write_files=False)
