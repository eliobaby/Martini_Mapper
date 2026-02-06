"""
Helper utilities for TASK 3.a
"""

from __future__ import annotations

import string
import random
from typing import Any, Dict, List, Tuple, Optional, Set

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
    atom_global_idx: int
) -> List[Tuple[int, float]]:
    """
    Same logic as Step 6 nested collect_array3_candidates().
    Returns [(candidate_gid, bond_order), ...]
    """
    atom = atoms_by_gid[atom_global_idx]
    candidates: List[Tuple[int, float]] = []
    seen_neighbors: Set[int] = set()

    for nbr_local, bond_order in atom[4]:
        if not (0 <= nbr_local < len(section)):
            continue

        neighbor = section[nbr_local]
        g_nbr = neighbor[0]

        if final[g_nbr] != "":
            continue
        
        has_len1_foreign = any(len(full_mapping[outer[0]]) == 1 for outer in neighbor[3])
        if has_len1_foreign:
            continue
        
        # seen_neighbors counts ANY qualifying neighbor next to the atom
        seen_neighbors.add(g_nbr)

        # candidates remain carbon-only
        if neighbor[1].upper() == "C":
            candidates.append((g_nbr, bond_order))

    return candidates, seen_neighbors


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
) -> None:
    """
    Entire Step 6 logic moved out. Mutates final + array3.
    """

    seen_candidates_pass1 = set()
    nbr_candidates = []

    # Edge case (must run before PASS 1)
    benzene_handle_array3_size3_unmapped_edgecase(
        section=section,
        atoms_by_gid=atoms_by_gid,
        final=final,
        martini_dict=martini_dict,
        array3=array3,
    )

    # PASS 1
    for a_idx in list(array3):
        if final[a_idx] != "":
            array3.remove(a_idx)
            continue

        _, seen_neighbors = benzene_collect_array3_candidates(
            section=section,
            atoms_by_gid=atoms_by_gid,
            final=final,
            full_mapping=full_mapping,
            atom_global_idx=a_idx
        )

        seen_candidates_pass1.update(seen_neighbors)
        
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
        
        nbr_candidates, _ = benzene_collect_array3_candidates(
            section=section,
            atoms_by_gid=atoms_by_gid,
            final=final,
            full_mapping=full_mapping,
            atom_global_idx=a_idx
        )

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

    # PASS 3
    for a_idx in list(array3):
        if final[a_idx] != "":
            array3.remove(a_idx)
            continue

        nbr_candidates, _ = benzene_collect_array3_candidates(
            section=section,
            atoms_by_gid=atoms_by_gid,
            final=final,
            full_mapping=full_mapping,
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

def benzene_step8_pair_remaining_array1(section, atoms_by_gid, gid_to_local, array1, final, martini_dict):
    """
    Step 8 extracted into a helper.
    Builds pairs from array1 and assigns beads to those pairs.

    CHANGE: if len(array1)==2, only pair if they are directly connected in the ring.
    """
    pairs1 = []

    def are_connected(g1, g2) -> bool:
        a = atoms_by_gid[g1]
        # a[4] = [(nbr_local, bo), ...]
        for nbr_local, _bo in a[4]:
            if 0 <= nbr_local < len(section):
                if section[nbr_local][0] == g2:
                    return True
        return False

    if len(array1) == 2:
        if are_connected(array1[0], array1[1]):
            pairs1.append((array1[0], array1[1]))
        # else: do nothing (leave them for later logic/fallback)

    elif len(array1) == 4:
        candidate_pair = None
        for a_global in array1:
            a_atom = atoms_by_gid[a_global]
            connected_neighbors = []
            for nbr_local, _ in a_atom[4]:
                if 0 <= nbr_local < len(section):
                    neighbor = section[nbr_local]
                    if neighbor[0] in array1:
                        connected_neighbors.append(neighbor[0])
            if len(connected_neighbors) == 1:
                candidate_pair = (a_global, connected_neighbors[0])
                break

        if candidate_pair is not None:
            pairs1.append(candidate_pair)
            remaining = [x for x in array1 if x not in candidate_pair]
            if len(remaining) == 2:
                pairs1.append((remaining[0], remaining[1]))
            else:
                for i in range(0, len(remaining), 2):
                    pairs1.append((remaining[i], remaining[i+1]))
        else:
            sorted_array1 = sorted(array1, key=gid_to_local.get)
            for i in range(0, len(sorted_array1), 2):
                pairs1.append((sorted_array1[i], sorted_array1[i+1]))

    elif len(array1) == 6:
        candidate_pair = None

        # Prefer N–N pairing first (your existing behavior)
        n_nodes = [g for g in array1 if atoms_by_gid[g][1].upper() == "N"]
        if len(n_nodes) >= 2:
            for a_global in n_nodes:
                a_atom = atoms_by_gid[a_global]
                for nbr_local, _ in a_atom[4]:
                    if 0 <= nbr_local < len(section):
                        neighbor = section[nbr_local]
                        if neighbor[0] in n_nodes:
                            candidate_pair = (a_global, neighbor[0])
                            break
                if candidate_pair is not None:
                    break

        if candidate_pair is None:
            for a_global in array1:
                a_atom = atoms_by_gid[a_global]
                connected_neighbors = []
                for nbr_local, _ in a_atom[4]:
                    if 0 <= nbr_local < len(section):
                        neighbor = section[nbr_local]
                        if neighbor[0] in array1:
                            connected_neighbors.append(neighbor[0])
                if connected_neighbors:
                    candidate_pair = (a_global, connected_neighbors[0])
                    break

        if candidate_pair is not None:
            pairs1.append(candidate_pair)
            remaining = [x for x in array1 if x not in candidate_pair]

            candidate_pair_4 = None
            for a_global in remaining:
                a_atom = next(atom for atom in section if atom[0] == a_global)
                connected_neighbors = []
                for nbr_local, _ in a_atom[4]:
                    if 0 <= nbr_local < len(section):
                        neighbor = section[nbr_local]
                        if neighbor[0] in remaining:
                            connected_neighbors.append(neighbor[0])
                if len(connected_neighbors) == 1:
                    candidate_pair_4 = (a_global, connected_neighbors[0])
                    break

            if candidate_pair_4 is not None:
                pairs1.append(candidate_pair_4)
                remaining = [x for x in remaining if x not in candidate_pair_4]
                if len(remaining) == 2:
                    pairs1.append((remaining[0], remaining[1]))
                else:
                    for i in range(0, len(remaining), 2):
                        pairs1.append((remaining[i], remaining[i+1]))
            else:
                sorted_remaining = sorted(remaining, key=gid_to_local.get)
                for i in range(0, len(sorted_remaining), 2):
                    pairs1.append((sorted_remaining[i], sorted_remaining[i+1]))
        else:
            sorted_array1 = sorted(array1, key=gid_to_local.get)
            for i in range(0, len(sorted_array1), 2):
                pairs1.append((sorted_array1[i], sorted_array1[i+1]))

    # Assign beads for the resulting pairs (your existing logic)
    for pair in pairs1:
        atom1 = atoms_by_gid[pair[0]]
        atom2 = atoms_by_gid[pair[1]]
        e1 = atom1[1].upper()
        e2 = atom2[1].upper()

        path = e1 + e2
        bead_key = find_bead(martini_dict, 1, path)
        assign_bead(final, bead_key, pair[0], pair[1])
        
def benzene_handle_array3_size3_unmapped_edgecase(
    *,
    section: List[List[Any]],
    atoms_by_gid: Dict[int, List[Any]],
    final: List[str],
    martini_dict: Dict[str, List[Any]],
    array3: List[int],
) -> bool:
    """
    Edge case handler (runs BEFORE PASS 1):

    Trigger if:
      - len(array3) == 3
      - NONE of the atoms in this ring are mapped yet (final[gid] == "" for all gids in section)
      - within array3, one atom connects to the other two (a 'hub')

    Action:
      - consider the other 3 atoms not in array3 (the complement)
      - find the 'hub' in the complement (connects to the other two)
      - group that hub with ONE of the other two, prioritizing a non-carbon partner
      - map the chosen pair using SECTION 1:
            key = find_bead(martini_dict, 1, "{type1}{type2}")
            assign_bead(final, key, hub_gid, partner_gid)

    Returns True if it applied a mapping, else False.
    """

    if len(array3) != 3:
        return False

    ring_gids = [a[0] for a in section]

    # none of the atoms in this ring mapped
    if any(final[g] != "" for g in ring_gids):
        return False

    array3_set = set(array3)

    def neighbors_in_set(gid: int, allowed: set[int]) -> List[int]:
        out = []
        atom = atoms_by_gid[gid]
        for nbr_local, _bo in atom[4]:
            if 0 <= nbr_local < len(section):
                nbr_gid = section[nbr_local][0]
                if nbr_gid in allowed:
                    out.append(nbr_gid)
        return out

    def find_hub(gids: List[int]) -> int | None:
        s = set(gids)
        for g in gids:
            if len(neighbors_in_set(g, s)) == 2:
                return g
        return None

    # array3 must have a hub connected to the other two
    hub_in_array3 = find_hub(array3)
    if hub_in_array3 is None:
        return False

    # complement (the other 3 atoms)
    other_gids = [g for g in ring_gids if g not in array3_set]
    if len(other_gids) != 3:
        return False

    hub_other = find_hub(other_gids)
    if hub_other is None:
        return False

    partners = [g for g in other_gids if g != hub_other]  # two atoms
    if len(partners) != 2:
        return False

    # pick partner: prioritize non-carbon, else any
    p1, p2 = partners
    t1 = atoms_by_gid[p1][1].upper()
    t2 = atoms_by_gid[p2][1].upper()

    if t1 != "C" and t2 == "C":
        partner = p1
    elif t2 != "C" and t1 == "C":
        partner = p2
    else:
        partner = p1  # both carbon or both non-carbon -> pick any

    hub_type = atoms_by_gid[hub_other][1].upper()
    partner_type = atoms_by_gid[partner][1].upper()

    key = find_bead(martini_dict, 1, f"{hub_type}{partner_type}")
    if key is None:
        # optional: try reversed order if your dict stores the other orientation
        key = find_bead(martini_dict, 1, f"{partner_type}{hub_type}")
    if key is None:
        raise ValueError(f"No section-1 bead for pair {hub_type}{partner_type} (or reversed)")

    assign_bead(final, key, hub_other, partner)
    return True
