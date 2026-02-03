from pathlib import Path
import shutil
from rdkit import Chem
from setup_mapping import get_atom_properties, connectivity_matrix, duplicate_ring_number_handler, remove_stars
from martini_3_dictionary import get_m3_dict
from mapping_scheme import map_molecule, parse_smiles
from algorithm import map_martini_beads
from cg_from_xtb_pipeline import build_cg_from_xtb
from xtb import smiles_to_ref
import copy
import argparse

def main(argv=None):
    parser = argparse.ArgumentParser(prog="martini_mapper")
    parser.add_argument("name", nargs="?", help="Compound name (e.g., Benzene)")
    parser.add_argument("smiles", nargs="?", help="SMILES string (e.g., c1ccccc1)")
    args = parser.parse_args(argv)
    if args.name is None:
        compound_name = input("Name: ")
    else:
        compound_name = args.name
    if args.smiles is None:
        smiles = input("Smiles: ")
    else:
        smiles = args.smiles
    smiles = remove_stars(smiles)
    smiles_fix = duplicate_ring_number_handler(smiles)
    mol = Chem.MolFromSmiles(smiles)
    atom_properties = get_atom_properties(mol)
    final = ["" for _ in atom_properties]
    conn_mat = connectivity_matrix(mol, len(atom_properties))
    bead_dict = get_m3_dict()
    smiles_token = parse_smiles(smiles_fix)
    mapping = map_molecule(smiles_token, conn_mat, atom_properties)
    mapping_copy = copy.deepcopy(mapping)
    final = map_martini_beads(mapping, final, bead_dict)
    '''
    from text_file import export_bead_mapping, gro_maker, itp_maker
    beads_data, connections_data = export_bead_mapping(final, mapping_copy, smiles, compound_name)
    gro_maker(beads_data, compound_name)
    itp_maker(beads_data, connections_data, compound_name)
    '''
    from text_file import export_bead_mapping, gro_maker
    beads_data, connections_data = export_bead_mapping(final, mapping_copy, smiles, compound_name)
    gro_maker(beads_data, compound_name)
    
    smiles_to_ref(compound_name, smiles)
    build_cg_from_xtb(
        smiles=smiles,
        final=final,
        mapping=mapping,
        ref_gro=f"ref_{compound_name}.gro",
        ref_xtc=f"ref_{compound_name}.xtc",
        compound = compound_name,
        out_prefix = compound_name,
        T=300.0,
    )
    
    # -----------------------------
    # Keep only: name/name.gro, name/name.itp, name/name.txt
    # -----------------------------
    out_dir = Path(compound_name)
    out_dir.mkdir(exist_ok=True)
    
    keep = {
        f"{compound_name}.gro",
        f"{compound_name}.itp",
        f"{compound_name}.txt",
    }
    
    # Move the kept files into ./name/
    for fname in keep:
        src = Path(fname)
        if src.exists():
            shutil.move(str(src), str(out_dir / src.name))
    
    # Delete common intermediates produced by smiles_to_ref + build_cg_from_xtb
    trash_files = [
        f"{compound_name}.xyz",
        f"{compound_name}.pdb",
        f"ref_{compound_name}.gro",
        f"ref_{compound_name}.xtc",
        f"{compound_name}_cg_ref.gro",
        f"{compound_name}_cg_ref.xtc",
        f"{compound_name}_bonded_report.txt",  # only if it ever exists
    ]
    
    for f in trash_files:
        p = Path(f)
        if p.exists():
            p.unlink()
    
    trash_dirs = [
        f"{compound_name}_md",  # xtb MD working directory
    ]
    
    for d in trash_dirs:
        p = Path(d)
        if p.exists() and p.is_dir():
            shutil.rmtree(p)
            
if __name__ == "__main__":
    main()