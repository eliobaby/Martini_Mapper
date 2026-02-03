"""
Helper utilities for TASK 3.e (non-ring section 1–bead mapping).

These were extracted from algorithm.py to keep map_non_ring_section_1bead()
small and readable, and to centralize repeated candidate-selection logic
(section 7 / section 11 lookups).
"""

from __future__ import annotations

from collections import deque
from typing import Any, Dict, List, Optional, Sequence


# -------------------------------
# Path / branch tracing helpers
# -------------------------------

def trace_linear_path(section: List[List[Any]], start_local: int) -> str:
    """
    Trace a 2-edge non-ring section (a simple path) starting from an edge atom,
    returning a SMILES-like string (e.g., 'C=CNO').
    """
    visited = set()
    curr_local = start_local

    # Start with the atom’s type.
    path = section[curr_local][1]
    visited.add(curr_local)

    while True:
        # Find unvisited neighbors (inner connections only; section graph is local).
        next_candidates = [nbr for (nbr, _bo) in section[curr_local][4] if nbr not in visited]
        if not next_candidates:
            break
        # Linear path assumption: pick the only unvisited neighbor
        nxt = next_candidates[0]
        # Determine bond symbol
        bond_order = next(bo for (nbr, bo) in section[curr_local][4] if nbr == nxt)
        bond_sym = "=" if bond_order == 2 else "#" if bond_order == 3 else ""
        path += bond_sym + section[nxt][1]
        visited.add(nxt)
        curr_local = nxt

    return path


def trace_branch(section: List[List[Any]], center_local: int, next_local: int) -> str:
    """
    Trace one branch out of a center atom until an edge is reached (or dead end),
    returning a SMILES-like branch string (without the center atom).
    """
    visited = {center_local}
    curr = next_local
    branch = ""
    prev = center_local

    while True:
        # Bond symbol between prev -> curr
        bond_order = next(bo for (nbr, bo) in section[prev][4] if nbr == curr)
        bond_sym = "=" if bond_order == 2 else "#" if bond_order == 3 else ""
        branch += bond_sym + section[curr][1]
        visited.add(curr)

        # Continue forward (exclude the node we came from)
        neighbors = [nbr for (nbr, _bo) in section[curr][4] if nbr not in visited]
        if not neighbors:
            break
        prev, curr = curr, neighbors[0]

    # remove a leading bond symbol if it exists (e.g. "=O" is okay, but we keep it)
    # Keep as-is; matching code compares exact strings.
    return branch


def get_final_edge(section: List[List[Any]], center_local: int, next_local: int) -> Optional[int]:
    """
    Follow a branch from center_local via next_local and return the last local index.
    """
    visited = {center_local}
    prev = center_local
    curr = next_local
    last = curr
    while True:
        visited.add(curr)
        last = curr
        nxts = [nbr for (nbr, _bo) in section[curr][4] if nbr != prev and nbr not in visited]
        if not nxts:
            return last
        prev, curr = curr, nxts[0]


# -------------------------------
# Martini dict candidate helpers
# -------------------------------

def find_linear_candidates(martini_dict: Dict[str, List[Any]], section_id: int, path_str: str) -> List[str]:
    """
    Find candidate bead keys for a linear path in section 7/11 (val[1] == 0).
    Matching is direction-independent.
    """
    p = path_str.upper()
    rp = p[::-1]
    keys = [
        k for k, v in martini_dict.items()
        if v[0] == section_id and v[1] == 0 and (v[2].upper() == p or v[2].upper() == rp)
    ]
    # unique while preserving order-ish
    return list(dict.fromkeys(keys))


def find_branch_candidates(martini_dict: Dict[str, List[Any]], section_id: int, center_symbol: str,
                           branch_paths: Sequence[str]) -> List[str]:
    """
    Find candidate beads for a 3-edge branched structure in section 7/11.
    This matches entries with parentheses patterns like 'C(=O)(N)(C)'.
    Matching is order-independent across branches.
    """
    import re

    candidates: List[str] = []
    for key, val in martini_dict.items():
        if not (val[0] == section_id and val[1] == 1 and "(" in str(val[2])):
            continue

        prefix, _remainder = str(val[2]).split("(", 1)
        if prefix.strip() != center_symbol:
            continue

        branch_segments = re.findall(r"\(([^)]+)\)", str(val[2]))
        if len(branch_segments) != len(branch_paths):
            continue

        # order-independent multiset match
        unmatched = list(branch_segments)
        matched_all = True
        for bp in branch_paths:
            for bs in list(unmatched):
                if bp == bs:
                    unmatched.remove(bs)
                    break
            else:
                matched_all = False
                break

        if matched_all:
            candidates.append(key)

    return list(dict.fromkeys(candidates))


# -------------------------------
# Ambiguity resolution helper
# -------------------------------

def _hydrogen_filter_candidates(candidates: Sequence[str],
                               martini_dict: Dict[str, List[Any]],
                               section: List[List[Any]]) -> List[str]:
    """
    Apply the val[5] atom-type + H-count rule to prune candidates down.
    """
    filtered: List[str] = []
    for key in candidates:
        val = martini_dict.get(key)
        if not val or len(val) < 6:
            # No rule -> keep candidate
            filtered.append(key)
            continue

        rule = val[5]
        if not rule or not isinstance(rule, (list, tuple)) or len(rule) < 2:
            filtered.append(key)
            continue

        atom_type = rule[0]
        target_h = list(rule[1:])  # length 1 or 2

        matching_atoms = [atom for atom in section if len(atom) > 6 and atom[1] == atom_type]

        if len(matching_atoms) == 1:
            if matching_atoms[0][6] == target_h[0]:
                filtered.append(key)
        elif len(matching_atoms) == 2:
            if len(target_h) >= 2:
                hs = sorted([matching_atoms[0][6], matching_atoms[1][6]])
                if sorted(target_h[:2]) == hs:
                    filtered.append(key)
        else:
            # if we don't find 1-2 atoms of that type, we cannot validate -> drop
            continue

    return list(dict.fromkeys(filtered))


def resolve_unique_candidate(candidates: Sequence[str],
                             martini_dict: Dict[str, List[Any]],
                             section: List[List[Any]],
                             *,
                             path_str: Optional[str] = None,
                             foreign_bead: Optional[str] = None) -> str:
    """
    Return exactly 1 bead candidate. If multiple exist, apply the task3e val[5]
    hydrogen-rule pruning. Special handling:
      - For CO/OC: if foreign bead starts with SX4e -> TP1d; if starts with TC5 -> TN2a;
        otherwise, apply hydrogen-rule pruning.
    """
    cand = list(dict.fromkeys(candidates))

    if len(cand) == 1:
        return cand[0]
    if not cand:
        raise ValueError("No candidates to resolve.")

    # Special case for CO / OC ambiguity
    if path_str and path_str.upper() in ("CO", "OC") and foreign_bead:
        if foreign_bead.startswith("SX4e"):
            return "TP1d"
        if foreign_bead.startswith("TC5"):
            return "TN2a"
        # else fall through to hydrogen-rule pruning

    pruned = _hydrogen_filter_candidates(cand, martini_dict, section)

    if len(pruned) == 1:
        return pruned[0]
    if not pruned:
        raise ValueError(f"Candidate resolution removed all candidates (started with {cand}).")
    raise ValueError(f"Ambiguous candidates remain after resolution: {pruned} (started with {cand}).")


# -------------------------------
# Misc helpers used in TASK 3.e
# -------------------------------

def bfs_path(section: List[List[Any]], start: int, end: int) -> List[int]:
    visited = {start}
    queue = deque([(start, [start])])
    while queue:
        curr, path = queue.popleft()
        if curr == end:
            return path
        for (nbr, _bo) in section[curr][4]:
            if nbr not in visited:
                visited.add(nbr)
                queue.append((nbr, path + [nbr]))
    return []


def map_two_atoms_by_connectivity(section: List[List[Any]],
                                  local_a: int,
                                  local_b: int,
                                  martini_dict: Dict[str, List[Any]],
                                  final: List[str],
                                  generate_random_string) -> None:
    """
    Helper for the special 4-atom breakdown inside the 2-edge case.
    Maps exactly two atoms to a single bead based on their 2-atom path pattern.
    """
    ga = section[local_a][0]
    gb = section[local_b][0]

    # bond between them
    sym = ""
    for nbr, bo in section[local_a][4]:
        if nbr == local_b:
            sym = "=" if bo == 2 else "#" if bo == 3 else ""
            break

    path = (section[local_a][1] + sym + section[local_b][1]).upper()

    keys = find_linear_candidates(martini_dict, 7, path)
    if not keys:
        keys = find_linear_candidates(martini_dict, 11, path)

    if not keys:
        raise ValueError(
            f"Linear path of 2 from 4 cannot be mapped: no candidate bead found for 2-edge mapping (path: {path})."
        )

    bead = keys[0] + generate_random_string()
    final[ga] = bead
    final[gb] = bead


def build_and_assign_trio(section: List[List[Any]],
                          center_idx: int,
                          edge_pair: List[int],
                          martini_dict: Dict[str, List[Any]],
                          final: List[str],
                          generate_random_string) -> None:
    """
    Helper for the 4-edge case: construct a 3-atom path (edge-center-edge), pick a unique bead,
    and assign it to all three atoms.
    """
    eA, eB = edge_pair
    symA = section[eA][1]
    symC = section[center_idx][1]
    symB = section[eB][1]

    bond1 = next(bo for (nbr, bo) in section[eA][4] if nbr == center_idx)
    sym_bond1 = "=" if bond1 == 2 else "#" if bond1 == 3 else ""
    bond2 = next(bo for (nbr, bo) in section[center_idx][4] if nbr == eB)
    sym_bond2 = "=" if bond2 == 2 else "#" if bond2 == 3 else ""

    path_str = (symA + sym_bond1 + symC + sym_bond2 + symB).upper()

    candidate_keys = find_linear_candidates(martini_dict, 7, path_str)
    if not candidate_keys:
        candidate_keys = find_linear_candidates(martini_dict, 11, path_str)

    if len(candidate_keys) != 1:
        # attempt hydrogen-based disambiguation
        bead_key = resolve_unique_candidate(candidate_keys, martini_dict, [section[eA], section[center_idx], section[eB]])
    else:
        bead_key = candidate_keys[0]

    bead = bead_key + generate_random_string()
    for idx in (eA, center_idx, eB):
        final[section[idx][0]] = bead
        
# -------------------------------
# Graph / traversal helpers (TASK 3.d / 3.f)
# -------------------------------

def build_inner_graph_local(section: List[List[Any]]) -> Dict[int, List[int]]:
    """
    Build a local-index adjacency list from a section's inner connections.
    Nodes are local indices 0..len(section)-1.
    """
    return {i: [nbr for (nbr, _bo) in atom[4]] for i, atom in enumerate(section)}


def bfs_distance_graph(graph: Dict[int, List[int]], start: int, end: int) -> int:
    """
    BFS shortest-path distance in number of EDGES. Returns float('inf') if unreachable.
    """
    if start == end:
        return 0
    visited = {start}
    q = deque([(start, 0)])
    while q:
        cur, dist = q.popleft()
        for nbr in graph.get(cur, []):
            if nbr == end:
                return dist + 1
            if nbr not in visited:
                visited.add(nbr)
                q.append((nbr, dist + 1))
    return float('inf')


def bfs_path_graph(graph: Dict[int, List[int]], start: int, end: int) -> List[int]:
    """
    BFS path as a list of nodes (local indices). Empty list if unreachable.
    """
    if start == end:
        return [start]
    q = deque([(start, [start])])
    visited = {start}
    while q:
        cur, path = q.popleft()
        for nbr in graph.get(cur, []):
            if nbr == end:
                return path + [nbr]
            if nbr not in visited:
                visited.add(nbr)
                q.append((nbr, path + [nbr]))
    return []


def bfs_path_local(section: List[List[Any]], start: int, end: int) -> List[int]:
    """Convenience wrapper for bfs_path_graph() using a section."""
    return bfs_path_graph(build_inner_graph_local(section), start, end)


def induced_subgraph_is_connected(indices: List[int], section: List[List[Any]]) -> bool:
    """
    True if the subgraph induced by local indices is connected.
    """
    if not indices:
        return False
    graph = build_inner_graph_local(section)
    allowed = set(indices)
    start = indices[0]
    visited = {start}
    q = deque([start])
    while q:
        cur = q.popleft()
        for nbr in graph.get(cur, []):
            if nbr in allowed and nbr not in visited:
                visited.add(nbr)
                q.append(nbr)
    return visited == allowed


def induced_subgraph_components(indices: List[int], section: List[List[Any]]) -> List[List[int]]:
    """
    Connected components (lists of local indices) of the induced subgraph.
    """
    graph = build_inner_graph_local(section)
    remaining = set(indices)
    comps: List[List[int]] = []
    while remaining:
        start = next(iter(remaining))
        comp = []
        q = deque([start])
        seen = {start}
        while q:
            cur = q.popleft()
            comp.append(cur)
            for nbr in graph.get(cur, []):
                if nbr in remaining and nbr not in seen:
                    seen.add(nbr)
                    q.append(nbr)
        comps.append(comp)
        remaining -= set(comp)
    return comps


def dfs_trace_from_seed(section: List[List[Any]], seed: int) -> List[int]:
    """
    Replicates the Task 3.f nested dfs_trace(): trace a single chain starting at seed,
    always taking the first neighbor that is not visited and has degree < 3.
    """
    visited: set[int] = set()

    def _dfs(local_idx: int) -> List[int]:
        trace = [local_idx]
        visited.add(local_idx)
        for (nbr, _bond) in section[local_idx][4]:
            if nbr not in visited and len(section[nbr][4]) < 3:
                trace.extend(_dfs(nbr))
                break
        return trace

    return _dfs(seed)

