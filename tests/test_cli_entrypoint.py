from martini_mapper.main import main


def test_cli_no_xtb_no_files():
    rc = main(["Benzene", "c1ccccc1", "--no-xtb", "--no-files"])
    assert rc == 0
