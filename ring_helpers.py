"""
Helper utilities for TASK 3.e (non-ring section 1–bead mapping).

These were extracted from algorithm.py to keep map_ring_...()
small and readable.
"""

from __future__ import annotations

import string
import random
from typing import Any, Dict, List, Tuple, Optional

def qualifies_foreign_non_CO(foreign_sec) -> bool:
    """Return True if foreign section has len>=3 OR len==2 but not (C,O)."""
    if len(foreign_sec) >= 3:
        return True
    if len(foreign_sec) == 2:
        elems = {a[1].upper() for a in foreign_sec}
        return elems != {"C", "O"}
    return False

def generate_random_string(length=6):
    """
    Generates a random string of a given length using all printable characters.
    
    Parameters:
        length (int): The length of the generated string (default is 6).
    
    Returns:
        str: A random string of the specified length.
    """
    allowed_chars = string.printable.strip()  # Excludes whitespace characters
    return ''.join(random.choices(allowed_chars, k=length))

def assign_bead(final: List[str], bead_key: str, *atom_ids: int) -> str:
    """
    Universal bead assignment helper.

    Usage:
        assign_bead(final, bead_key, a, b)
        assign_bead(final, bead_key, a, b, c)
        assign_bead(final, bead_key, *some_tuple)
        assign_bead(final, bead_key, *some_list)

    - If rstr not provided, generates one.
    - Assigns the same bead to every atom id provided.
    - Returns the bead string.
    """
    if bead_key is None:
        raise ValueError("assign_bead: bead_key is None")
    rstr = generate_random_string()
    bead = bead_key + rstr
    for gid in atom_ids:
        final[gid] = bead
    return bead

def pick_bead_key(martini_dict: Dict[str, List[Any]],
                  element: str,
                  bond_order: Optional[int] = None,
                  kind: str = 'S') -> Optional[str]:
    """
    Pick a bead‐type key from martini_dict whose
      • key starts with `kind` ('S' or 'T'),
      • val[0] == 2,
      • and whose val[2] (the "(...)" part) matches our element + optional '='.

    element: e.g. 'O', 'N', 'S', 'CL', 'I', 'C', 'BR', F
    bond_order: if element == 'O' and bond_order == 2, we look for '(=O)' else '(O)'
    kind: 'S' for SN6*, 'SC6'…, 'SN6a', etc.; 'T' for TN6*, 'TC4'…, 'TX2', etc.
    """
    # build the pattern to look for inside the parentheses
    el = element.capitalize()  # so 'CL' -> 'Cl'
    if element == 'O' and bond_order == 2:
        pat = '(=O)'
    else:
        pat = f'({el})'

    for key, val in martini_dict.items():
        if not key.startswith(kind):
            if key.startswith('X') and kind == "S":
                key
            else:
                continue
        if val[0] != 2:
            continue
        # val[2] is something like 'CC(O)' or 'C(Cl)' or 'CC(=O)'
        if val[2].endswith(pat):
            return key
    return None

def outer_info(full_mapping, outer_tup):
    """
    outer_tup is (foreign_section_index, foreign_atom_local_index, bond_order)
    Returns: (foreign_section, foreign_atom, bond_order)
    """
    sec_idx, local_idx, bond_order = outer_tup
    foreign_sec = full_mapping[sec_idx]
    foreign_atom = foreign_sec[local_idx]
    return foreign_sec, foreign_atom, bond_order

def benzene_collect_array3_candidates(
    *,
    section: List[List[Any]],
    atoms_by_gid: Dict[int, List[Any]],
    final: List[str],
    full_mapping: List[List[List[Any]]],
    array4: List[int],
    atom_global_idx: int
) -> List[Tuple[int, float]]:
    """
    Same logic as Step 6 nested collect_array3_candidates().
    Returns [(candidate_gid, bond_order), ...]
    """
    if len(array4) != 0:
        return []

    atom = atoms_by_gid[atom_global_idx]
    candidates: List[Tuple[int, float]] = []

    for nbr_local, bond_order in atom[4]:
        if not (0 <= nbr_local < len(section)):
            continue

        neighbor = section[nbr_local]
        g_nbr = neighbor[0]

        if final[g_nbr] != "":
            continue

        if neighbor[1].upper() != "C":
            continue

        qualifies = False
        if len(neighbor[3]) == 0:
            qualifies = True
        else:
            for outer in neighbor[3]:
                foreign_sec = full_mapping[outer[0]]
                if qualifies_foreign_non_CO(foreign_sec):
                    qualifies = True
                    break

        if qualifies:
            candidates.append((g_nbr, bond_order))

    return candidates


def benzene_assign_S_triplet(
    *,
    center_gid: int,
    chosen_candidate_gid: int,
    atoms_by_gid: Dict[int, List[Any]],
    final: List[str],
    martini_dict: Dict[str, List[Any]],
    full_mapping: List[List[List[Any]]],
) -> None:
    atom = atoms_by_gid[center_gid]
    foreign_tup = atom[3][0]
    foreign_sec, foreign_atom, bond_order = outer_info(full_mapping, foreign_tup)

    key = pick_bead_key(martini_dict, foreign_atom[1].upper(), bond_order, kind='S')
    if key is None:
        raise ValueError(f"No S-key for foreign {foreign_atom[1]} bond={bond_order}")

    assign_bead(final, key, atom[0], chosen_candidate_gid, foreign_atom[0])


def benzene_assign_T_pair(
    *,
    center_gid: int,
    atoms_by_gid: Dict[int, List[Any]],
    final: List[str],
    martini_dict: Dict[str, List[Any]],
    full_mapping: List[List[List[Any]]],
) -> None:
    atom = atoms_by_gid[center_gid]
    foreign_tup = atom[3][0]
    foreign_sec, foreign_atom, bond_order = outer_info(full_mapping, foreign_tup)

    if atom[1].upper() == "C":
        key = pick_bead_key(martini_dict, foreign_atom[1].upper(), bond_order, kind='T')
    elif foreign_atom[1].upper() == "C" and atom[1].upper() == "N":
        key = find_bead(martini_dict, 2, "N(C)")
    else:
        raise ValueError("atom and foreign are both not C")

    if key is None:
        raise ValueError(f"No T-key found for bond={bond_order}")

    assign_bead(final, key, atom[0], foreign_atom[0])


def benzene_candidate_adjacent_to_mapped(
    *,
    section: List[List[Any]],
    atoms_by_gid: Dict[int, List[Any]],
    final: List[str],
    candidate_gid: int,
) -> bool:
    cand_atom = atoms_by_gid[candidate_gid]
    for nbr_local, _ in cand_atom[4]:
        if 0 <= nbr_local < len(section):
            g = section[nbr_local][0]
            if final[g] != "":
                return True
    return False


def benzene_process_array3_two_passes(
    *,
    section: List[List[Any]],
    atoms_by_gid: Dict[int, List[Any]],
    final: List[str],
    martini_dict: Dict[str, List[Any]],
    full_mapping: List[List[List[Any]]],
    array3: List[int],
    array4: List[int],
) -> None:
    """
    Entire Step 6 logic moved out. Mutates final + array3.
    """

    seen_candidates_pass1 = set()

    # PASS 1
    for a_idx in list(array3):
        if final[a_idx] != "":
            array3.remove(a_idx)
            continue

        nbr_candidates = benzene_collect_array3_candidates(
            section=section,
            atoms_by_gid=atoms_by_gid,
            final=final,
            full_mapping=full_mapping,
            array4=array4,
            atom_global_idx=a_idx
        )

        for cand_g, _ in nbr_candidates:
            seen_candidates_pass1.add(cand_g)

        if len(nbr_candidates) == 1:
            chosen = nbr_candidates[0][0]
            benzene_assign_S_triplet(
                center_gid=a_idx,
                chosen_candidate_gid=chosen,
                atoms_by_gid=atoms_by_gid,
                final=final,
                martini_dict=martini_dict,
                full_mapping=full_mapping,
            )
            array3.remove(a_idx)

        elif len(nbr_candidates) == 2:
            # pass 1 skip
            continue

        else:
            benzene_assign_T_pair(
                center_gid=a_idx,
                atoms_by_gid=atoms_by_gid,
                final=final,
                martini_dict=martini_dict,
                full_mapping=full_mapping,
            )
            array3.remove(a_idx)

    # PASS 2
    for a_idx in list(array3):
        if final[a_idx] != "":
            array3.remove(a_idx)
            continue

        if len(seen_candidates_pass1) == 4:
            benzene_assign_T_pair(
                center_gid=a_idx,
                atoms_by_gid=atoms_by_gid,
                final=final,
                martini_dict=martini_dict,
                full_mapping=full_mapping,
            )
            array3.remove(a_idx)
            continue

        nbr_candidates = benzene_collect_array3_candidates(
            section=section,
            atoms_by_gid=atoms_by_gid,
            final=final,
            full_mapping=full_mapping,
            array4=array4,
            atom_global_idx=a_idx
        )

        if len(nbr_candidates) == 1:
            chosen = nbr_candidates[0][0]

        elif len(nbr_candidates) == 2:
            chosen = None
            for cand_g, _ in nbr_candidates:
                if benzene_candidate_adjacent_to_mapped(
                    section=section,
                    atoms_by_gid=atoms_by_gid,
                    final=final,
                    candidate_gid=cand_g
                ):
                    chosen = cand_g
                    break
            if chosen is None:
                chosen = nbr_candidates[0][0]

        else:
            benzene_assign_T_pair(
                center_gid=a_idx,
                atoms_by_gid=atoms_by_gid,
                final=final,
                martini_dict=martini_dict,
                full_mapping=full_mapping,
            )
            array3.remove(a_idx)
            continue

        benzene_assign_S_triplet(
            center_gid=a_idx,
            chosen_candidate_gid=chosen,
            atoms_by_gid=atoms_by_gid,
            final=final,
            martini_dict=martini_dict,
            full_mapping=full_mapping,
        )
        array3.remove(a_idx)

def find_bead(martini_dict: Dict[str, List[Any]], sect: int, pattern: str, tag: str | None = None) -> Optional[str]:
    pat = pattern.upper()
    rev = pat[::-1]
    #print(pattern, sect)
    for k, v in martini_dict.items():
        if tag is not None and (len(v) <= 4 or v[4] != tag):
            continue
        if v[0] == sect and v[2].upper() in (pat, rev):
            return k
    return None

def merge_phenol_to_diol(section, final, martini_dict, full_mapping, tn_prefix="TN6+"):
    """
    Post-pass for ring sections:
    - Find ring atoms whose final[gid] starts with tn_prefix (default "TN6") => phenol
    - Determine actual ring cycle order from inner connectivity (local indices)
    - Pair adjacent phenols into diols (including wrap-around adjacency)
    - Odd phenol left unchanged
    - For each paired phenol, find its outer O (size-1 non-ring section) and map 4 atoms into diol bead
    """

    n = len(section)
    if n < 3:
        return

    # ---------- Build ring cycle order using LOCAL indices ----------
    # inner connections are (neighbor_local, bond_order)
    nbrs = []
    for i in range(n):
        nbrs_i = [nbr for (nbr, _bo) in section[i][4] if 0 <= nbr < n]
        nbrs.append(nbrs_i)

    # walk the ring as a cycle
    start = 0
    prev = None
    curr = start
    order = []

    while True:
        order.append(curr)

        choices = nbrs[curr]
        if not choices:
            # not a ring
            return

        if prev is None:
            nxt = choices[0]
        else:
            if len(choices) == 1:
                nxt = choices[0]
            else:
                nxt = choices[0] if choices[0] != prev else choices[1]

        prev, curr = curr, nxt

        if curr == start:
            break
        if len(order) > n + 1:
            # failed to close -> not a simple cycle
            return

    ring_len = len(order)
    if ring_len != n:
        # section isn't a simple 1-cycle ring (rare for your ring tasks) -> don't touch
        return

    # map cycle position -> local index, and local -> cycle position
    pos_of_local = {loc: p for p, loc in enumerate(order)}

    # ---------- Collect phenol positions on the ring ----------
    phenol_pos = []
    for loc in range(n):
        gid = section[loc][0]
        if isinstance(final[gid], str) and final[gid].startswith(tn_prefix):
            phenol_pos.append(pos_of_local[loc])

    if len(phenol_pos) < 2:
        return

    phenol_pos.sort()

    # ---------- Break phenols into consecutive runs (circular-aware) ----------
    runs = []
    run = [phenol_pos[0]]
    for p in phenol_pos[1:]:
        if p == run[-1] + 1:
            run.append(p)
        else:
            runs.append(run)
            run = [p]
    runs.append(run)

    # merge wraparound run if needed: e.g., [0,1] and [4,5] in a 6-ring should merge if 5 adjacent to 0
    if len(runs) > 1 and runs[0][0] == 0 and runs[-1][-1] == ring_len - 1:
        merged = runs[-1] + runs[0]
        runs = [merged] + runs[1:-1]

    # ---------- Find DIOL bead key (val[0]=2, val[2]=C(O)(CO)) ----------
    bead_key = find_bead(martini_dict, 2, "C(O)(CO)")
    
    # ---------- Helper: find the phenolic outer oxygen gid for a ring local index ----------
    def phenol_outer_O_gid(ring_local: int):
        atom = section[ring_local]
        # outer tuples: (foreign_section_index, foreign_atom_local_index, bond_order)
        for outer_tup in atom[3]:
            foreign_sec, foreign_atom, _bo = outer_info(full_mapping, outer_tup)
            # phenol O should be a size-1 non-ring section (type 0) and element O
            if len(foreign_sec) == 1 and foreign_atom[2] == 0 and foreign_atom[1].upper() == "O":
                return foreign_atom[0]
        return None

    # ---------- Pair within each run: (0,1), (2,3), ...; leave odd leftover ----------
    for run in runs:
        # run is a list of cycle positions
        for i in range(0, len(run) - 1, 2):
            p1, p2 = run[i], run[i + 1]
            loc1 = order[p1]
            loc2 = order[p2]
            gid1 = section[loc1][0]
            gid2 = section[loc2][0]

            # sanity: still phenols
            if not (final[gid1].startswith(tn_prefix) and final[gid2].startswith(tn_prefix)):
                continue

            o1 = phenol_outer_O_gid(loc1)
            o2 = phenol_outer_O_gid(loc2)
            if o1 is None or o2 is None:
                continue

            # "rip everything off" for these 4 atoms (optional but matches your intent)
            for g in (gid1, gid2, o1, o2):
                final[g] = ""

            # assign diol bead to 4 atoms
            assign_bead(final, bead_key, gid1, gid2, o1, o2)
