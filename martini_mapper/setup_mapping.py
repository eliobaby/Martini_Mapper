import numpy as np
import re

def get_atom_properties(mol):
    """
    Returns per-atom properties:
      [symbol, is_in_ring, is_edge, num_H]
    where num_H = total number of H attached (implicit + explicit).
    """
    atom_properties = []

    for atom in mol.GetAtoms():
        neighbors = atom.GetNeighbors()

        is_edge = (any(n.GetMass() > 1.1 for n in neighbors) and atom.GetDegree() == 1)

        # Total Hs attached to this atom (implicit + explicit)
        num_H = atom.GetTotalNumHs(includeNeighbors=True)

        atom_data = [
            atom.GetSymbol(),   # Atom type (e.g., C, O, N)
            atom.IsInRing(),    # Part of a ring?
            is_edge,            # Edge node?
            num_H               # How many H attached to this atom
        ]
        atom_properties.append(atom_data)

    return atom_properties

def connectivity_matrix(mol, length):
    # Create the NxN connectivity matrix
    matrix = np.zeros((length, length), dtype=float)
    for bond in mol.GetBonds():
        # Check if the bond is aromatic and assign a bond order accordingly
        if bond.GetIsAromatic():
            # We can choose how to represent aromatic bonds.
            # For instance, 1.5 is a common representation.
            bond_order = 1.5
        else:
            bond_order = bond.GetBondTypeAsDouble()
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        matrix[idx1][idx2] = bond_order
        matrix[idx2][idx1] = bond_order
    return matrix

def duplicate_ring_number_handler(smiles: str) -> str:
    """
    Ensure all ring-closure numbers in a SMILES string are unique by
    renumbering extra pairs beyond the first for any repeated digit.

    - Detects tokens: '%dd' for 10+, and single digits '1'–'9'.
    - Records each token's positions and groups by number.
    - For each number with >2 occurrences, reassigns extra pairs to new sequential numbers.
    """
    # Regex to find %dd or single-digit closures
    token_pattern = re.compile(r'%\d{2}|[1-9]')
    tokens = []  # list of (num:int, start:int, end:int)

    for m in token_pattern.finditer(smiles):
        txt = m.group()
        if txt.startswith('%'):
            num = int(txt[1:])
        else:
            num = int(txt)
        tokens.append({'num': num, 'start': m.start(), 'end': m.end(), 'txt': txt})

    # Group tokens by ring number
    groups = {}
    for t in tokens:
        groups.setdefault(t['num'], []).append(t)

    # Find highest existing number
    max_num = max(groups.keys(), default=0)
    next_num = max_num + 1

    # Prepare remapping for extra pairs
    remap = {}  # index in tokens list -> new number
    for num, lst in groups.items():
        count = len(lst)
        # must be even
        if count <= 2:
            continue
        # sort by appearance
        lst_sorted = sorted(lst, key=lambda x: x['start'])
        # process extra pairs
        for i in range(2, count, 2):
            t1, t2 = lst_sorted[i], lst_sorted[i+1]
            remap[id(t1)] = next_num
            remap[id(t2)] = next_num
            next_num += 1

    # If nothing to remap, return original
    if not remap:
        return smiles

    # Build new SMILES by replacing tokens
    parts = []
    last = 0
    for t in sorted(tokens, key=lambda x: x['start']):
        parts.append(smiles[last:t['start']])
        if id(t) in remap:
            newn = remap[id(t)]
            # format: single digit or %dd
            rep = str(newn) if newn < 10 else f"%{newn:02d}"
            parts.append(rep)
        else:
            parts.append(t['txt'])
        last = t['end']
    parts.append(smiles[last:])
    return ''.join(parts)

def remove_stars(s: str) -> str:
    """
    Return a copy of s with all '*' characters removed,
    and then all literal '()' substrings removed.
    """
    # Step 1: remove asterisks
    s = s.replace('*', '')
    s = s.replace('/', '')
    s = s.replace('\\', '')
    # Step 2: remove empty parentheses
    cleaned = s.replace('()', '')
    return cleaned