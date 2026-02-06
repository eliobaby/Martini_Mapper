"""
Helper utilities for TASK 3.b/c
"""

from typing import Dict, List, Optional, Set, Tuple, Any

def build_cn_candidate_graph(
    section: List[List],
    final: List[str],
    *,
    elements: Tuple[str, ...] = ("C", "N"),
    bond_orders: Tuple[float, ...] = (2, 1.5),
) -> Tuple[List[int], Dict[int, str], Dict[int, Set[int]]]:
    """
    Build a graph on unmapped C/N atoms that participate in (2 or 1.5) inner bonds.

    Returns:
        cand: list of global atom ids in this section
        elem_of: {gid: "C"|"N"}
        neighbors: {gid: set(gid)} adjacency among cand via inner bonds (any bond order)
    """
    elems = {e.upper() for e in elements}
    cand: List[int] = []
    elem_of: Dict[int, str] = {}

    # global->local for this section
    g2local = {a[0]: i for i, a in enumerate(section)}

    for atom in section:
        g = atom[0]
        if final[g] != "":
            continue
        e = str(atom[1]).upper()
        if e not in elems:
            continue
        if any(bo in bond_orders for _, bo in atom[4]):
            cand.append(g)
            elem_of[g] = e

    cand_set = set(cand)
    neighbors: Dict[int, Set[int]] = {g: set() for g in cand}

    for g in cand:
        atom = section[g2local[g]]
        for nbr_local, _bo in atom[4]:
            if 0 <= nbr_local < len(section):
                g_nbr = section[nbr_local][0]
                if g_nbr in cand_set:
                    neighbors[g].add(g_nbr)

    return cand, elem_of, neighbors


def adjacent_in_allowed_set(
    section: List[List],
    gid_to_local: Dict[int, int],
    *,
    global_idx: int,
    allowed: Set[int],
) -> List[int]:
    """
    Return global ids in `allowed` adjacent (inner-bonded) to `global_idx` inside `section`.
    """
    out: List[int] = []
    a_atom = section[gid_to_local[global_idx]]
    for nbr_local, _ in a_atom[4]:
        if 0 <= nbr_local < len(section):
            g2 = section[nbr_local][0]
            if g2 in allowed:
                out.append(g2)
    return out


def choose_pairs_for_nonbenzene_5_ring(
    *,
    cand: List[int],
    elem_of: Dict[int, str],
    neighbors: Dict[int, Set[int]],
    final: List[str],
) -> List[Tuple[int, int]]:
    """
    Choose C/N pairings for Step A of map_nonbenzene_5_ring_section.

    Rules:
      - If exactly 5 candidates and matches the special N/C compositions described,
        enforce the required CC/CN pattern with a single N leftover.
      - Otherwise, choose a maximum-size matching (target = floor(n/2), or 2 for n==5),
        then among matchings with the same leftover count, prefer more CC pairs,
        and only then more CN pairs.

    Returns:
        list of (a_gid, b_gid) pairs.
    """

    def pair_type(a: int, b: int) -> str:
        ea, eb = elem_of[a], elem_of[b]
        if ea == "C" and eb == "C":
            return "CC"
        if ea == "N" and eb == "N":
            return "NN"
        return "CN"

    def required_pattern_for_5ring(candidates: List[int]) -> Optional[Dict[str, int]]:
        nN = sum(1 for g in candidates if elem_of[g] == "N")
        nC = len(candidates) - nN
        if len(candidates) != 5:
            return None
        if (nN, nC) == (1, 4):
            return {"CC": 2, "CN": 0, "NN": 0}
        if (nN, nC) == (2, 3):
            return {"CC": 1, "CN": 1, "NN": 0}
        if (nN, nC) == (3, 2):
            return {"CC": 0, "CN": 2, "NN": 0}
        return None

    def gen_matchings_any_size(
        nodes: List[int],
        target_pairs: int,
    ) -> List[Tuple[List[Tuple[int, int]], List[int]]]:
        """
        Generate all matchings with exactly target_pairs pairs.
        Returns (pairs, leftover_nodes).
        """
        nodes = list(nodes)

        def rec(remaining: List[int], pairs: List[Tuple[int, int]]):
            if len(pairs) == target_pairs:
                return [(pairs[:], remaining[:])]
            if len(remaining) < 2:
                return []
            out = []
            first = remaining[0]

            # Option 1: leave 'first' unmatched (only if we can still reach target_pairs)
            rest = remaining[1:]
            needed = target_pairs - len(pairs)
            if len(rest) >= 2 * needed:
                out.extend(rec(rest, pairs))

            # Option 2: match 'first' with any neighbor in remaining
            for i in range(1, len(remaining)):
                second = remaining[i]
                if second not in neighbors[first]:
                    continue
                new_remaining = remaining[1:i] + remaining[i + 1 :]
                pairs.append((first, second))
                out.extend(rec(new_remaining, pairs))
                pairs.pop()

            return out

        return rec(nodes, [])

    def gen_perfect_matchings(nodes: List[int]) -> List[List[Tuple[int, int]]]:
        """All perfect matchings on nodes (nodes length must be even)."""
        if not nodes:
            return [[]]
        out: List[List[Tuple[int, int]]] = []
        first = nodes[0]
        for i in range(1, len(nodes)):
            second = nodes[i]
            if second not in neighbors[first]:
                continue
            rest = nodes[1:i] + nodes[i + 1 :]
            for sub in gen_perfect_matchings(rest):
                out.append([(first, second)] + sub)
        return out

    def counts_of(pairs_list: List[Tuple[int, int]]) -> Dict[str, int]:
        c = {"CC": 0, "CN": 0, "NN": 0}
        for a, b in pairs_list:
            c[pair_type(a, b)] += 1
        return c

    def is_feasible(pairs_list: List[Tuple[int, int]]) -> bool:
        for a, b in pairs_list:
            if final[a] != "" or final[b] != "":
                return False
        return True

    if not cand:
        return []

    req = required_pattern_for_5ring(cand)

    # Special forced patterns (5-ring)
    if req is not None:
        n_globals = [g for g in cand if elem_of[g] == "N"]
        forced_solution = None

        # leftover must be a single N; try each N leftover, perfect match remaining 4
        for n_left in n_globals:
            remaining = [g for g in cand if g != n_left]
            for pm in gen_perfect_matchings(sorted(remaining)):
                if not is_feasible(pm):
                    continue
                cts = counts_of(pm)
                if all(cts[k] == req[k] for k in req):
                    forced_solution = (pm, [n_left])
                    break
            if forced_solution is not None:
                break

        if forced_solution is None:
            raise ValueError(
                "5-ring double bond pattern not mappable (cannot satisfy CC/CN pattern with leftover N)"
            )

        return forced_solution[0]

    # Generic matching selection
    target_pairs = 2 if len(cand) == 5 else (len(cand) // 2)
    all_matchings = gen_matchings_any_size(sorted(cand), target_pairs)

    feasible: List[Tuple[List[Tuple[int, int]], List[int]]] = []
    for ps, leftover in all_matchings:
        if is_feasible(ps):
            feasible.append((ps, leftover))

    if not feasible:
        return []

    # Prefer: (1) fewer leftovers (i.e. more pairs), then (2) more CC, then (3) more CN.
    def score(ps: List[Tuple[int, int]], leftover: List[int]) -> Tuple[int, int, int]:
        cts = counts_of(ps)
        return (-len(leftover), cts["CC"], cts["CN"])

    best = max(feasible, key=lambda x: score(x[0], x[1]))
    return best[0]

# =============================================================================
# Non-benzene 6-ring fallback helpers (count==1/2/3)
# =============================================================================

def nb6_assign_lone_carbon(final: List[str], section: List[List[Any]], lone: List[Any]) -> None:
    """
    Implements the non-benzene 6-ring fallback for a single remaining atom (count == 1).

    - lone must be a carbon.
    - If any neighbor already has a T* bead, convert that bead to an S* bead (special-case TN6 -> SN4),
      propagate to all atoms already carrying that bead tag, and assign to the lone atom.
    - Else, if any neighbor has an S* bead, copy it to lone.
    - Else, raise.
    """
    if final[lone[0]] != "":
        return

    if lone[1].upper() != "C":
        raise ValueError("non benzene 6-ring has a lone non-Carbon left")

    nbrs = [section[t[0]] for t in lone[4]]
    t_beads = [final[n[0]] for n in nbrs if final[n[0]].startswith("T")]

    if t_beads:
        old = t_beads[0]
        new = old
        # Keep the current special-case: TN6 -> SN4
        if old.startswith("T"):
            new = "S" + new[1:]
        elif not new.startswith("S"):
            new = "S" + new[1:]

        for i, v in enumerate(final):
            if v == old:
                final[i] = new
        final[lone[0]] = new
        return

    s_beads = [final[n[0]] for n in nbrs if final[n[0]].startswith("S")]
    if s_beads:
        final[lone[0]] = s_beads[0]
        return

def nb6_atoms_are_connected(section: List[List[Any]], a: List[Any], b: List[Any]) -> bool:
    """Return True if a and b are mutually connected via inner connections."""
    a_to_b = any(section[t[0]][0] == b[0] for t in a[4])
    b_to_a = any(section[t[0]][0] == a[0] for t in b[4])
    return a_to_b and b_to_a


def nb6_try_assign_connected_pair_T(section: List[List[Any]],
                                   final: List[str],
                                   martini_dict: Dict[str, List[Any]],
                                   a1: List[Any],
                                   a2: List[Any],
                                   assign_bead) -> bool:
    """
    Try the non-benzene 6-ring count==2 connected-pair logic.

    If a1/a2 are connected and both unmapped, attempt to assign a T* bead where val[0] in (3,5)
    and val[2] matches the 2-atom element string in either direction.

    Returns True if it assigned a bead, else False (including the not-connected case).
    """
    if final[a1[0]] != "" or final[a2[0]] != "":
        return False

    if not nb6_atoms_are_connected(section, a1, a2):
        return False

    types = a1[1].upper() + a2[1].upper()
    for key, val in martini_dict.items():
        if key.startswith("T") and val[0] in (3, 5) and (val[2] == types or val[2] == types[::-1]):
            assign_bead(final, key, a1[0], a2[0])
            return True

    # Connected but no key is a real failure (old logic would silently leave them unmapped).
    raise ValueError(f"non benzene 6-ring: no T-bead candidate for connected pair {types}")


def nb6_handle_count2(final: List[str],
                     section: List[List[Any]],
                     martini_dict: Dict[str, List[Any]],
                     remaining: List[List[Any]],
                     assign_bead) -> None:
    """
    Implements count==2 behavior:
    - If the two atoms are connected -> assign a T* bead (val[0] in (3,5)).
    - If not connected -> treat as two separate count==1 cases.
    """
    a1, a2 = remaining

    if nb6_try_assign_connected_pair_T(section, final, martini_dict, a1, a2, assign_bead):
        return

    # Not connected -> two independent lone assignments
    nb6_assign_lone_carbon(final, section, a1)
    nb6_assign_lone_carbon(final, section, a2)


def nb6_handle_count3(final: List[str],
                     section: List[List[Any]],
                     martini_dict: Dict[str, List[Any]],
                     full_mapping: List[List[List[Any]]],
                     remaining: List[List[Any]],
                     assign_bead) -> None:
    """
    Implements count==3 behavior (the 'actual connecting' version):

    1) Prefer the center-with-two-unmapped-neighbors case, assign an S* bead for (nbr-center-nbr).
       If still unmapped, run the 'merge' logic (T->S promotion or borrow S).

    2) If no valid center exists, try to find a connected pair among the 3 and apply count==2 to that pair,
       leaving the last atom for count==1.

    3) If none are connected, treat them as 3 independent count==1 cases.
    """
    remaining_gids = {a[0] for a in remaining}

    def unmapped_neighbors_in_remaining(a: List[Any]) -> List[List[Any]]:
        return [section[t[0]] for t in a[4]
                if final[section[t[0]][0]] == "" and section[t[0]][0] in remaining_gids]

    candidates = [a for a in remaining if len(unmapped_neighbors_in_remaining(a)) == 2]

    if not candidates:
        # Fallback A: find a connected pair among the remaining 3
        pair = None
        for i in range(len(remaining)):
            for j in range(i + 1, len(remaining)):
                a1, a2 = remaining[i], remaining[j]
                if nb6_atoms_are_connected(section, a1, a2):
                    pair = (a1, a2)
                    break
            if pair:
                break

        if pair:
            nb6_handle_count2(final, section, martini_dict, list(pair), assign_bead)
            return

        # Fallback B: none connected -> do 3 lone assignments
        for lone in remaining:
            nb6_assign_lone_carbon(final, section, lone)
        return

    center = candidates[0]
    nbrs = unmapped_neighbors_in_remaining(center)  # exactly 2
    s = nbrs[0][1].upper() + center[1].upper() + nbrs[1][1].upper()

    for key, val in martini_dict.items():
        if key.startswith("S") and val[0] in (3, 5) and (val[2] == s or val[2] == s[::-1]):
            assign_bead(final, key, center[0], nbrs[0][0], nbrs[1][0])
            break

    if final[center[0]] != "":
        return

    # --- merge logic (as in algorithm.py) ---
    merged = False

    # 1) look for a C-neighbor with a T-bead neighbor
    for nbr in nbrs:
        if nbr[1].upper() == "C":
            inner_idxs = [section[t[0]][0] for t in nbr[4]]
            outer_idxs = [full_mapping[t[0]][t[1]][0] for t in nbr[3]]
            all_idxs = inner_idxs + outer_idxs

            t_beads = [final[i] for i in all_idxs if final[i].startswith("T")]
            if t_beads:
                old = t_beads[0]
                new = "S" + old[1:]
                for j, v in enumerate(final):
                    if v == old:
                        final[j] = new
                final[nbr[0]] = new
                merged = True
                break

    # 2) if still not merged, fall back to any S-bead neighbor
    if not merged:
        for nbr in nbrs:
            inner_idxs = [section[t[0]][0] for t in nbr[4]]
            outer_idxs = [full_mapping[t[0]][t[1]][0] for t in nbr[3]]
            all_idxs = inner_idxs + outer_idxs

            s_beads = [final[i] for i in all_idxs if final[i].startswith("S")]
            if s_beads:
                final[nbr[0]] = s_beads[0]
                merged = True
                break

    if not merged:
        raise ValueError("Merging is not possible; mapping failed for count=3")
