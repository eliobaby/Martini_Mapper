from typing import List, Dict, Any, Tuple, Optional
import warnings

# =================================================
# HELPER FUNCTIONS
# =================================================
from .non_aromatic_helpers import (
    build_cn_candidate_graph, 
    choose_pairs_for_nonbenzene_5_ring, 
    adjacent_in_allowed_set,
    nb6_assign_lone_carbon,
    nb6_try_assign_connected_pair_T,
    nb6_handle_count2,
    nb6_handle_count3,
)

from .non_ring_helpers import (
    trace_linear_path,
    trace_branch,
    get_final_edge,
    find_branch_candidates,
    resolve_unique_candidate,
    bfs_path,
    map_two_atoms_by_connectivity,
    build_and_assign_trio,
    build_inner_graph_local,
    bfs_distance_graph,
    bfs_path_graph,
    induced_subgraph_is_connected,
    induced_subgraph_components,
    dfs_trace_from_seed,
    dist_to_nearest_center,
)

from .ring_helpers import (
     outer_info,
     assign_bead,
     benzene_process_array3_two_passes,
     pick_bead_key,
     generate_random_string,
     find_bead,
     merge_phenol_to_diol,
     benzene_step8_pair_remaining_array1,
)

def build_section_index(section: List[List[Any]]):
    """Fast lookups for this section."""
    atoms_by_gid = {a[0]: a for a in section}
    gid_to_local = {a[0]: i for i, a in enumerate(section)}
    return atoms_by_gid, gid_to_local

def foreign_is_size1_nonring(full_mapping, outer_tup) -> bool:
    """True if outer points to section type 0, size 1."""
    foreign_sec, foreign_atom, _ = outer_info(full_mapping, outer_tup)
    return foreign_atom[2] == 0 and len(foreign_sec) == 1

def find_bead_by_val1(martini_dict: Dict[str, List[Any]], sect: int, val1: int) -> Optional[str]:
    for k, v in martini_dict.items():
        if len(v) > 1 and v[0] == sect and v[1] == val1:
            return k
    return None

def map_adjacent_SS_pairs(section, final, martini_dict) -> None:
    """If two unmapped sulfurs are adjacent, map them as SS immediately."""
    atoms_by_gid, _ = build_section_index(section)
    bead_key = find_bead(martini_dict, 5, "SS")
    if bead_key is None:
        return

    for atom in section:
        g = atom[0]
        if atom[1].upper() != "S" or final[g] != "":
            continue
        for nbr_local, _ in atom[4]:
            if 0 <= nbr_local < len(section):
                nbr_gid = section[nbr_local][0]
                if atoms_by_gid[nbr_gid][1].upper() == "S" and final[nbr_gid] == "":
                    assign_bead(final, bead_key, g, nbr_gid)
                    break

def _retry_with_cleared_ring(section, final, fn, *args, **kwargs):
    """
    Try fn(section, final, ...). If it fails, clear this ring's atoms in final and retry once.
    Returns updated final, or re-raises the second exception.
    """
    try:
        if fn == map_benzene_ring_section:
            return fn(section, final, *args, **kwargs, passes = 2)
        else:
            return fn(section, final, *args, **kwargs)
    except Exception as e1:
        for a in section:
            final[a[0]] = ""
        try:
            if fn == map_benzene_ring_section:
                return fn(section, final, *args, **kwargs, passes = 2)
            else:
                return fn(section, final, *args, **kwargs)
        except Exception as e2:
            # raise the second failure (more relevant after clearing), but keep context
            raise RuntimeError(f"Retry after clearing ring atoms failed: {e2}") from e1
                
# =============================================================================
# TASK 3.a: Mapping Benzene Ring Section
# =============================================================================
def map_benzene_ring_section(section: List[List[Any]],
                             final: List[str],
                             martini_dict: Dict[str, List[Any]],
                             full_mapping: List[List[List[Any]]],
                             passes: int = 1) -> List[str]:
    """
    Map a benzene ring section (section type == 2) using a tree‐like algorithm.
    Implements the following steps:
      1. Build array0: candidate atoms with type 'c/C' that have outer connections where the foreign atom's
         section type != 0, and that are unmapped.
      2. Pair atoms in array0 (using inner connections) and assign bead "TC5e" + random string.
      3. Build array1: atoms with no outer connections or outer connections of size 3+ or 2 that aren't CO
         that are unmapped;
     	if the amount of atoms in this array is odd: 
            first find any atom that is not connected to any other atoms in this array and remove it. If found then go straight to the pairing 
            if not found, find any atom that is only connecting to exactly 1 other atom in this array and remove it, if found then go to the pairing 
            if not found, throw an error saying benzene ring odd parity is bugged
            pair them (leftovers → array2) and assign beads ("TC5" if both C, "TN6a" if one N).
      4. Build array3: atoms that have outer connections where for at least one connection the foreign
         atom’s section is type 0 and of length 1.
         Then for each such atom:
            - If it has exactly one neighbor (via inner connections) found in array2, remove that neighbor
              from array2 and assign a bead based on the foreign atom type (using the first outer connection):
                  O: "SN6", N: "SN6d", S: "SC6", Cl: "SX3", I: "X1", C: "SC4"
            - If it has two neighbors, choose one (preferring a neighbor adjacent to an already mapped atom)
              and assign similarly.
            - If it has no neighbor in array2, use the bond order of the outer connection to decide:
                  if O with bond==2: "TN6a", if O with bond==1: "TN6", if C: "TC4" and more 
      5. Build array4: atoms with outer connections that are unmapped and for which the foreign section is
         of length 2; for each, if the foreign atom is O and its connected atom in that section is C,
         assign bead "SN2a" + random string.
      6. Build array5: every remaining unmapped atom; pair these using inner connections and assign beads
         ("TC5" if both C, "TN6a" if one N).
      7. Final check: ensure every atom is mapped.
    """
    # Fast lookups: global atom id -> atom record (keeps original list-based atom structure)
    atoms_by_gid, gid_to_local = build_section_index(section)
    # --- nH HANDLER (must run before STEP A) ---
    # If an unmapped nitrogen has exactly 1 attached H AND it has at least one INNER bond of order 1.5,
    # assign it immediately as TN6d.
    for atom in section:
        g = atom[0]
        if final[g] != "":
            continue
    
        if atom[1].upper() != "N":
            continue
    
        num_H = atom[6]  # [global_index, element, ring_status, outer_connection, inner_connection, isedge, num_H]
        if num_H != 1:
            continue
    
        # NEW requirement: N must have an inner aromatic bond (bond order 1.5) to qualify
        has_inner_15 = any(bo == 1.5 for _, bo in atom[4])
        if not has_inner_15:
            continue
    
        bead_key = find_bead(martini_dict, 5, "NH")
        assign_bead(final, bead_key, g)
        
    # ---------------------- Step 1: Build array0 ----------------------
    if passes == 1:
        array0 = []
        seen0 = set()
        for atom in section:
            idx = atom[0]
            # base criteria: carbon, unmapped, has any outer (foreign) connection
            if atom[1].lower() == 'c' and final[idx] == "" and atom[3]:
                valid_foreign = any(full_mapping[tup[0]][tup[1]][2] != 0 for tup in atom[3])
                if not valid_foreign:
                    continue
        
                # look for an inner neighbor that also meets those criteria
                for inner_tup in atom[4]:
                    nbr_global = section[inner_tup[0]][0]
                    if nbr_global in seen0 or nbr_global == idx:
                        continue
        
                    nbr_obj = atoms_by_gid[nbr_global]
                    if (nbr_obj[1].lower() == 'c'
                        and final[nbr_global] == ""
                        and nbr_obj[3]
                        and any(full_mapping[t[0]][t[1]][2] != 0 for t in nbr_obj[3])):
                        # both carbons form a valid pair—add them once
                        array0.extend([idx, nbr_global])
                        seen0.update([idx, nbr_global])
                        break
    
    # array0 should now contain paired indices: length 2, 4, or 6
    
    # ---------------------- Step 2: Pair atoms in array0 ----------------------
        pairs0 = []
        used0 = set()
        for atom in section:
            if atom[0] in array0 and atom[0] not in used0:
                for tup in atom[4]:
                    nbr_local = tup[0]
                    if 0 <= nbr_local < len(section):
                        neighbor = section[nbr_local]
                        if neighbor[0] in array0 and neighbor[0] not in used0:
                            pairs0.append((atom[0], neighbor[0]))
                            used0.add(atom[0])
                            used0.add(neighbor[0])
                            break
        for p in pairs0:
            atom_candidate = atoms_by_gid[p[0]]
        
            # pick bead key from martini_dict instead of hard-coding
            bead_key = None
            for tup in atom_candidate[3]:
                foreign_atom = full_mapping[tup[0]][tup[1]]
                order = foreign_atom[2]
                if order in (1, 2):
                    # find the dict key where val[0] == 8 and val[1] == bond order
                    bead_key = find_bead_by_val1(martini_dict, 8, order)
                    break
            assign_bead(final, bead_key, p[0], p[1])
        
    # ---------------------- Step 3: Build array3 ----------------------
    # For every atom with an outer connection (index3 non-empty) that is still unmapped,
    # and for at least one outer connection the foreign atom (in full_mapping) has section type 0
    # and the foreign section has length 1, and that foreign atom is unmapped.
    array3 = []
    for atom in section:
        if len(atom[3]) != 0 and final[atom[0]] == "":
            for tup in atom[3]:
                foreign_sec = full_mapping[tup[0]]
                foreign_atom = foreign_sec[tup[1]]
                if foreign_atom[2] == 0 and len(foreign_sec) == 1 and final[foreign_atom[0]] == "":
                    array3.append(atom[0])
                    break
    
    # ---------------------- Step 4: Build array4 ----------------------
    # Atoms with outer connections still unmapped where for at least one outer connection:
    #   - foreign section has length 2
    #   - the foreign atom in the tuple is section type 0 and unmapped
    #   - the foreign section contains exactly C and O, and the C–O bond is order 1
    #   - the inner atom connects to the foreign atom O (the tuple must point to O)
    
    array4 = []
    
    for atom in section:
        if len(atom[3]) != 0 and final[atom[0]] == "":
            for tup in atom[3]:
                si, li, bo = tup  # bo is the bond order between inner atom and this foreign atom
                foreign_sec = full_mapping[si]
                if len(foreign_sec) != 2:
                    continue
    
                foreign_atom = foreign_sec[li]
    
                # tuple foreign atom must be section type 0 and unmapped
                if foreign_atom[2] != 0 or final[foreign_atom[0]] != "":
                    continue
    
                # inner atom must connect to the foreign atom O
                if foreign_atom[1].upper() != "O":
                    continue
    
                # foreign section must contain exactly C and O
                elems = {a[1].upper() for a in foreign_sec}
                if elems != {"C", "O"}:
                    continue
    
                # C–O bond inside the foreign section must be order 1
                # Find the local indices of C and O in the foreign section
                c_local = next(i for i, a in enumerate(foreign_sec) if a[1].upper() == "C")
                o_local = next(i for i, a in enumerate(foreign_sec) if a[1].upper() == "O")
    
                # Check if C is bonded to O with order 1 (either direction)
                co_ok = any(nbr == o_local and bond == 1 for (nbr, bond) in foreign_sec[c_local][4]) or \
                        any(nbr == c_local and bond == 1 for (nbr, bond) in foreign_sec[o_local][4])
    
                if not co_ok:
                    continue
    
                array4.append(atom[0])
                break
    
    # ---------------------- Step 5: Process array4 (SN2a) ----------------------
    # NEW RULE:
    # Only map SN2a if array3 has size 1+ (instead of the old parity check).
    if len(array3) >= 1:
        for a_idx in array4:
            atom = atoms_by_gid[a_idx]
            foreign_tup = atom[3][0]  # take the first outer connection
            foreign_sec = full_mapping[foreign_tup[0]]
            foreign_atom = foreign_sec[foreign_tup[1]]
    
            # find the bead key for C(OC)
            bead_key = find_bead(martini_dict, 2, "C(OC)")
            targets = [atom[0], foreign_atom[0]]
            if foreign_atom[4]:
                nbr_local = foreign_atom[4][0][0]
                nbr_in_foreign = foreign_sec[nbr_local]
                targets.append(nbr_in_foreign[0])
            
            assign_bead(final, bead_key, *targets)
    # else: skip SN2a completely
    
    # ---------------------- Step 6: Process array3 (TWO PASSES) ----------------------
    # NEW:
    # - We run the loop twice.
    # - Candidate neighbors must:
    #   (c1) be Carbon (candidate atom type == 'C')
    #   (c2) be unmapped
    #   (c3) and satisfy:
    #        - has 0 outer connections
    #          OR
    #        - has an outer connection to a foreign section:
    #            (a) len >= 3, OR
    #            (b) len == 2 and not (C,O)
    #
    # - If len(nbr_candidates) == 2:
    #     * PASS 1: skip
    #     * PASS 2: handle (choose candidate adjacent to an already mapped atom; else first)
    #
    # - Every time we successfully assign a bead for an array3 atom, remove it from array3
    #   so PASS 2 only sees what's left.
    benzene_process_array3_two_passes(
        section=section,
        atoms_by_gid=atoms_by_gid,
        final=final,
        martini_dict=martini_dict,
        full_mapping=full_mapping,
        array3=array3,
    )
    # ---------------------- Step 7: Build array1 (NO SEARCH) ----------------------
    # NEW:
    # "the rest of the atoms can be put into array1 ... without any search whatsoever"
    array1 = [atom[0] for atom in section if final[atom[0]] == ""]
    
    # ---------------------- Step 8: Pair atoms from array1 (done LAST) ----------------------
    try:
        benzene_step8_pair_remaining_array1(section, atoms_by_gid, gid_to_local, array1, final, martini_dict)
    
        # Final Check
        for atom in section:
            if final[atom[0]] == "":
                raise ValueError("Benzene ring section not fully mappable!")
                
    except ValueError:
        # Fallback: if any ring atom has SN2a* bead, remove that entire bead-group and retry.
        sn2a_bead = None
        for atom in section:
            g = atom[0]
            if isinstance(final[g], str) and final[g].startswith("SN2a"):
                sn2a_bead = final[g]
                break
    
        if sn2a_bead is None:
            # no SN2a to remove; re-raise the original failure
            raise
    
        # remove every atom in the entire mapping that has that exact bead string
        for i, v in enumerate(final):
            if v == sn2a_bead:
                final[i] = ""
    
        # IMPORTANT: recompute array1 because final[] changed
        array1 = [atom[0] for atom in section if final[atom[0]] == ""]
        
        benzene_step8_pair_remaining_array1(section, atoms_by_gid, gid_to_local, array1, final, martini_dict)
    
    # step 9
    merge_phenol_to_diol(section, final, martini_dict, full_mapping)
    
    for atom in section:
        if final[atom[0]] == "":
            raise ValueError("Benzene ring section not fully mappable!")
    
    return final

# =============================================================================
# TASK 3.b: Mapping Non-Benzene 6-Ring Section
# =============================================================================
def map_nonbenzene_6_ring_section(
    section: List[List[Any]],
    final: List[str],
    martini_dict: Dict[str, List[Any]],
    full_mapping: List[List[List[Any]]]
) -> List[str]:
    """
    Revised mapping function for non-benzene 6-ring sections incorporating outer
    (foreign) connections.

    The algorithm proceeds in four major steps:
    
    1. Process ring–ring border nodes: For each unmapped node having an outer connection
       to a group with section != 0, check its inner neighbors:
         - Case B1: If exactly one inner neighbor qualifies and one of the nodes
           has exactly 2 outer connections with one of the non-shared foreign connections
           (foreign section == 0 and length == 1), assign the candidate, its inner neighbor,
           and that foreign atom to a bead (using the provided rules snippet).
         - Case B2: If both nodes have 2 such foreign connections, raise an error.
         - Case B3: Otherwise, simply assign the two atoms a "TC5" bead.
    
    2. Process double–bond connections: For each unmapped pair connected by an inner bond of
       order 2, check their outer connections:
         - Case D1: If one atom has a qualifying foreign connection (section == 0, size == 1),
           assign all three atoms to a bead using the provided snippet.
         - Case D2: If both atoms have such foreign connections, map each with its foreign
           partner using alternate rules.
         - Case D3: Else assign the pair a "TC5" bead.
    
    3. Process remaining unmapped nodes that have an outer connection to a non–ring (section 0)
       group of size 1:
         - Build an auxiliary array from eligible inner neighbors (based on two steps ahead).
         - If no eligible neighbor is found, map the node and its chosen foreign neighbor as a pair.
         - If exactly one is found, map the node, that neighbor, and the foreign neighbor together.
         - If two are found, raise an error.
    
    4. Final fallback: Count remaining unmapped atoms.
         - If exactly 6 remain, use the original (oxygen-based) mapping.
         - Otherwise, for counts 5, 4, 3, 2, or 1, fixed
    """
    # --- SPECIAL S–S HANDLER  ---
    # If two unmapped sulfurs sit next to each other in the ring, pair them immediately.
    map_adjacent_SS_pairs(section, final, martini_dict)
    atoms_by_gid, gid_to_local = build_section_index(section)
    
    # --- Step A: Process ring–ring border nodes
    for atom in section:
        has_outer = False
        for t in atom[3]:
            si, li = t[0], t[1]
            if 0 <= si < len(full_mapping) and 0 <= li < len(full_mapping[si]):
                fa = full_mapping[si][li]
                if len(fa) > 2 and fa[2] != 0:
                    has_outer = True
                    break
        
        if final[atom[0]] == "" and has_outer:
            candidate = atom
            matching_neighbors = []
            # Check candidate's inner neighbors for one that has an outer connection to the same ring
            for t in candidate[4]:
                nbr_local = t[0]
                if 0 <= nbr_local < len(section):
                    neighbor = section[nbr_local]
                    if final[neighbor[0]] == "":
                        # Does neighbor have an outer connection with the same section type as candidate?
                        for ft in neighbor[3]:
                            foreign_sec = full_mapping[ft[0]]
                            foreign_atom = foreign_sec[ft[1]]
                            if foreign_atom[2] == candidate[2]:
                                matching_neighbors.append(neighbor)
                                break
                            
            if len(matching_neighbors) == 0:
                print("ring-ring neighbors connected through only 1 atom, is ignore")
            elif len(matching_neighbors) == 1:
                inner_neighbor = matching_neighbors[0]
                cand_qual = [t for t in candidate[3] if foreign_is_size1_nonring(full_mapping, t)]
                neigh_qual = [t for t in inner_neighbor[3] if foreign_is_size1_nonring(full_mapping, t)]
                if candidate[1].upper() != 'C' or inner_neighbor[1].upper() != 'C':
                    print(f"ring–ring border fallback supports only C–C; found {candidate[1].upper()}-{inner_neighbor[1].upper()} skipping")
                    continue
                if (len(cand_qual) == 1 and len(neigh_qual) == 1):
                    # Case B2: Too complex intermediate section
                    # copy from Case D2: Both nodes have a qualifying foreign connection.
                    # Use each node's first qualifying connection
                    fc1 = cand_qual[0]
                    fc2 = neigh_qual[0]
                    _, foreign_atom1, _ = outer_info(full_mapping, fc1)
                    _, foreign_atom2, _ = outer_info(full_mapping, fc2)
                    # Use the same bond order from one of the nodes
                    bond_order = atom[3][0][2]
                    key1 = pick_bead_key(martini_dict, foreign_atom1[1].upper(), bond_order, kind='T')
                    key2 = pick_bead_key(martini_dict, foreign_atom2[1].upper(), bond_order, kind='T')
                    assign_bead(final, key1, candidate[0], foreign_atom1[0])
                    assign_bead(final, key2, inner_neighbor[0], foreign_atom2[0])
                elif (len(cand_qual) == 1) or (len(neigh_qual) == 1):
                    # Case B1: Use the one (candidate or neighbor) that has 2 outer connections.
                    fc = cand_qual[0] if len(cand_qual) == 1 else neigh_qual[0]
                    _, foreign_atom, _ = outer_info(full_mapping, fc)
                    bond_order = candidate[3][0][2]  # using candidate's first outer connection's bond order
                    key = pick_bead_key(martini_dict,
                                        foreign_atom[1].upper(),
                                        bond_order,
                                        kind='S')
                    assign_bead(final, key, candidate[0], inner_neighbor[0], foreign_atom[0])
                else:
                    # Case B3: only C–C allowed, lookup the bead with val[0]==1 and val[2]=="CC"
                    bead_key = find_bead(martini_dict, 1, "CC")
                    assign_bead(final, bead_key, candidate[0], inner_neighbor[0])
            else:
                # Ambiguous: more than one matching inner neighbor → skip this atom
                print(f"Ambiguous ring–ring border for atom {atom[0]}; skipping")
                continue
    # --- Step B: Process double–bond connections (inner bonds with order 2) ---
    processed_double = set()
    for atom in section:
        if final[atom[0]] != "":
            continue
    
        for t in atom[4]:
            nbr_local, bo = t
            if bo != 2 and bo != 1.5:
                continue
    
            if not (0 <= nbr_local < len(section)):
                continue
    
            neighbor = section[nbr_local]
            a_idx = atom[0]
            b_idx = neighbor[0]
    
            # skip if already mapped or already processed
            if final[b_idx] != "" or (a_idx, b_idx) in processed_double or (b_idx, a_idx) in processed_double:
                continue
    
            # pull the bond order directly from the inner‐bond tuple
            bond_order = bo
    
            # find qualifying foreign connections on each
            cand_qual  = [x for x in atom[3]     if foreign_is_size1_nonring(full_mapping, x)]
            neigh_qual = [x for x in neighbor[3] if foreign_is_size1_nonring(full_mapping, x)]
            # Case D1: both sides qualify
            if len(cand_qual) == 1 and len(neigh_qual) == 1:
                # FIX: ensure both ring atoms are carbon
                if atom[1].upper() != 'C' or neighbor[1].upper() != 'C':
                    raise ValueError("Double‐bond not mappable: both atoms must be C–C for D1")
                _, fa, _ = outer_info(full_mapping, cand_qual[0])
                _, fb, _ = outer_info(full_mapping, neigh_qual[0])
                key1 = pick_bead_key(martini_dict, fa[1].upper(), bond_order, kind='T')
                key2 = pick_bead_key(martini_dict, fb[1].upper(), bond_order, kind='T')
                assign_bead(final, key1, a_idx, fa[0])
                assign_bead(final, key2, b_idx, fb[0])
    
            # Case D2: exactly one side qualifies
            elif len(cand_qual) == 1 or len(neigh_qual) == 1:
                # FIX: ensure both ring atoms are carbon
                if atom[1].upper() != 'C' or neighbor[1].upper() != 'C':
                    raise ValueError("Double‐bond not mappable: both atoms must be C–C for D2")
                if len(cand_qual) == 1:
                    fc = cand_qual[0]
                    base_idx, other_idx = a_idx, b_idx
                else:
                    fc = neigh_qual[0]
                    base_idx, other_idx = b_idx, a_idx
                _, foreign, _ = outer_info(full_mapping, fc)
                key = pick_bead_key(martini_dict,
                                    foreign[1].upper(),
                                    bond_order,
                                    kind='S')
                assign_bead(final, key, base_idx, other_idx, foreign[0])
            # Case D3: neither side qualifies – choose bead type based on atom types
            else:
                # Fallback: look up a section-1 bead for the pair via the dictionary
                elem_a = atom[1].upper()
                elem_b = neighbor[1].upper()
                path = elem_a + elem_b

                # find a bead key in martini_dict where v[0]==1 and v[2] matches "CC", "CN"/"NC", or "NN"
                bead_key = find_bead(martini_dict, 1, path)
                assign_bead(final, bead_key, a_idx, b_idx)
            processed_double.add((a_idx, b_idx))
            # once we’ve assigned for this atom, break out to next atom
            break
    # --- Step C: Process remaining unmapped nodes with 1 or 2 non-ring size-1 foreign connections ---
    for i in range(2):
        for atom in section:
            if final[atom[0]] == "":
                # collect all size-1 non-ring outer tuples whose foreign atom is still unmapped
                f_tups = [t for t in atom[3]
                          if foreign_is_size1_nonring(full_mapping, t) and final[full_mapping[t[0]][t[1]][0]] == ""]
                # case: two size-1 foreigns → check all-C and assign via martini_dict
                if len(f_tups) == 2:
                    ft1, ft2 = f_tups
                    _, fa1, _ = outer_info(full_mapping, ft1)
                    _, fa2, _ = outer_info(full_mapping, ft2)
                    elem = atom[1].upper()
                    e1, e2 = fa1[1].upper(), fa2[1].upper()

                    if elem not in ("C", "O"):
                        raise ValueError("2 size-1 neighbors but central atom not C/O")
                    
                    # convenience
                    pair = (e1, e2)
                    pair_set = frozenset(pair)
                    
                    # Case 1: central C and ONE foreign carbon -> use S-key based on the *non-carbon* neighbor
                    if elem == "C" and ("C" in pair_set):
                        # pick the non-carbon if available, otherwise just pick one (old code picked ft2 when both are C)
                        if e1 != "C":
                            chosen_ft = ft1
                        else:
                            chosen_ft = ft2
                    
                        _, fa, bond = outer_info(full_mapping, chosen_ft)
                    
                        key = pick_bead_key(martini_dict, fa[1].upper(), bond, kind="S")
                        if not key:
                            raise ValueError(f"No S-key for {fa[1]}")
                        assign_bead(final, key, atom[0], fa1[0], fa2[0])
                                            
                    else:
                        # Case 2: central C or O with matching pair rules -> pick the bead pattern
                        pattern_map = {
                            # central C rules
                            ("C", frozenset({"O"})): "OCO",
                            ("C", frozenset({"N"})): "NCN",
                            ("C", frozenset({"F"})): "FCF",
                    
                            # central O rules
                            ("O", frozenset({"C"})): "COC",
                        }
                    
                        pattern = pattern_map.get((elem, pair_set))
                        if pattern is None:
                            raise ValueError(f"2 size-1 neighbors incompatible for central {elem}")
                    
                        bead_key = find_bead(martini_dict, 4, pattern)
                        assign_bead(final, bead_key, atom[0], fa1[0], fa2[0])
                
                    continue  # move on after handling this case
                # case: single size-1 foreign → proceed with original single-tuple logic
                if len(f_tups) == 1:
                    primary = f_tups[0]
                    _, foreign_atom, _ = outer_info(full_mapping, primary)
                    bond_order = primary[2]
                    # build array1 of candidate inner neighbors
                    array1 = []
                    for t in atom[4]:
                        nbr = section[t[0]]
                        if final[nbr[0]] == "":
                            # check if all of neighbor's inner connections remain unmapped
                            inner_unmapped = all(
                                final[section[t2[0]][0]] == "" for t2 in nbr[4]
                            )
                            # gather neighbor's size-1 non-ring foreign tuples unmapped
                            nf_tups = [x for x in nbr[3]
                                       if foreign_is_size1_nonring(full_mapping, x) and final[full_mapping[x[0]][x[1]][0]] == ""]
                            # include if neighbor is C and not inner_unmapped, with tuple rules
                            if nbr[1].upper() == 'C' and not inner_unmapped and not nf_tups:
                                array1.append(nbr)

                    inner_el = atom[1].upper()
                    outer_el = foreign_atom[1].upper()
                    
                    # 1) Handle N first (outside the len(array1) branching)
                    if inner_el == "N":
                        if outer_el == "C":
                            key = find_bead(martini_dict, 2, "N(C)")
                            if key is None:
                                raise ValueError("No bead found for section 2 pattern N(C)")
                            assign_bead(final, key, atom[0], foreign_atom[0])
                        else:
                            raise ValueError(
                                f"Unhandled N case: inner={inner_el}, outer={outer_el}"
                            )
                    
                    # 2) If C, then do the len(array1) == 0/1/2 logic
                    elif inner_el == "C":
                        if len(array1) == 0:
                            key = pick_bead_key(martini_dict, outer_el, bond_order, kind="T")
                            if key is None:
                                raise ValueError(f"No bead key found for inner={inner_el} outer={outer_el} bond={bond_order}")
                            assign_bead(final, key, atom[0], foreign_atom[0])
                    
                        elif len(array1) == 1:
                            key = pick_bead_key(martini_dict, outer_el, bond_order, kind="S")
                            if key is None:
                                raise ValueError(f"No bead key found for inner={inner_el} outer={outer_el} bond={bond_order}")
                            assign_bead(final, key, atom[0], array1[0][0], foreign_atom[0])
                    
                        elif len(array1) == 2:
                            if i == 1:
                                # Pick inner neighbor whose inner-neighbors do NOT already carry a “T”-bead in final[]
                                chosen = None
                                for nbr in array1:
                                    inner_globals = [section[t2[0]][0] for t2 in nbr[4]]
                                    if not any(final[g].startswith("T") for g in inner_globals if final[g] != ""):
                                        chosen = nbr
                                        break
                                if chosen is None:
                                    chosen = array1[0]
                    
                                bond_order = primary[2]
                                key = pick_bead_key(martini_dict, outer_el, bond_order, kind="S")
                                if key is None:
                                    raise ValueError(f"No bead key found for inner={inner_el} outer={outer_el} bond={bond_order}")
                                assign_bead(final, key, atom[0], chosen[0], foreign_atom[0])
                            else:
                                raise ValueError(f"Unhandled C case: len(array1)==2 but i={i}")
                    
                        else:
                            raise ValueError(f"Unhandled C case: len(array1)={len(array1)}")
                    
                    # 3) Neither N nor C
                    else:
                        raise ValueError(f"Unhandled inner element: inner={inner_el}, outer={outer_el}, len(array1)={len(array1)}")
                    
    # --- Step C.5: Special catch‐all for any unmapped inner atom that is not C or O,
    #                 but has a single‐sized non-ring (section 0) foreign neighbor which is C.
    #                 Map both atoms to a T‐type bead (bond order = 1).
    for atom in section:
        if final[atom[0]] != "":
            continue
        elem = atom[1].upper()
        if elem in ('C', 'O'):
            continue
        for ft in atom[3]:
            # ft = (foreign_section_index, foreign_atom_local_index, bond_order)
            _, foreign_atom, _ = outer_info(full_mapping, ft)
            if foreign_is_size1_nonring(full_mapping, ft) and foreign_atom[1].upper() == 'C':
                # Use bond_order = 1 and kind='T'
                key = pick_bead_key(martini_dict, elem, 1, kind='T')
                assign_bead(final, key, atom[0], foreign_atom[0])
                break   # move on to next atom
    # --- Step D: Final fallback mapping for remaining unmapped atoms
    remaining = [atom for atom in section if final[atom[0]] == ""]
    count = len(remaining)
    # ensure only C or O remain
    if any(atom[1].upper() not in ('C', 'O', 'N') for atom in remaining):
        raise ValueError("non-benzene 6-ring section has non C/O/N atom")
    if count == 6:
        # Fall back to the original oxygen-based mapping.
        array_O = [atom[0] for atom in section if atom[1].upper() == "O" and final[atom[0]] == ""]
    
        if len(array_O) == 1:
            # SN4a case: val[0]=3, val[2]=="COC", val[4]=="SC3"
            bead_key = find_bead(martini_dict, 3, "COC", tag="SC3")
            target_idx = array_O[0]
    
            # map oxygen + its two neighbors with one bead
            nbr_globals = [
                section[nbr_local][0]
                for nbr_local, _ in next(a for a in section if a[0] == target_idx)[4]
            ]
            assign_bead(final, bead_key, target_idx, *nbr_globals)
    
        elif len(array_O) == 2:
            # SN3a case: val[0]=3, val[2]=="COC", val[4]=="SN3a"
            bead_key = find_bead(martini_dict, 3, "COC", tag="SN3a")
            for idx in array_O:
                nbr_globals = [
                    section[nbr_local][0]
                    for nbr_local, _ in atoms_by_gid[idx][4]
                ]
                assign_bead(final, bead_key, idx, *nbr_globals)
    
        # now whatever remains
        array1 = [atom[0] for atom in section if final[atom[0]] == ""]
    
        # --- NEW CASE: 1 N and 5 C -> group CCN and CCC ---
        # i.e., make one 3-atom bead with pattern "CCN" and one 3-atom bead with pattern "CCC"
        if len(array1) == 6:
            n_nodes = [g for g in array1 if atoms_by_gid[g][1].upper() == "N"]
            c_nodes = [g for g in array1 if atoms_by_gid[g][1].upper() == "C"]
    
            if len(n_nodes) == 1 and len(c_nodes) == 5:
                n_g = n_nodes[0]
    
                n_adj_c = adjacent_in_allowed_set(section, gid_to_local, global_idx=n_g, allowed=set(c_nodes))
                if len(n_adj_c) < 2:
                    raise ValueError("1N+5C case: N does not have 2 adjacent carbons in this section")
    
                cA, cB = n_adj_c[0], n_adj_c[1]
                ccn_group = [cA, cB, n_g]
    
                ccc_group = [g for g in array1 if g not in ccn_group]
                if len(ccc_group) != 3 or any(atoms_by_gid[g][1].upper() != "C" for g in ccc_group):
                    raise ValueError("1N+5C case: leftover after CCN is not exactly 3 carbons")
    
                bead_key_ccn = find_bead(martini_dict, 3, "CCN")
                if bead_key_ccn is None:
                    raise ValueError("No section-3 bead found for pattern 'CCN'")
    
                bead_key_ccc = find_bead(martini_dict, 3, "CCC")
                if bead_key_ccc is None:
                    raise ValueError("No section-3 bead found for pattern 'CCC'")
    
                assign_bead(final, bead_key_ccn, *ccn_group)
                assign_bead(final, bead_key_ccc, *ccc_group)
    
            else:
                # Existing behavior for 6 remaining (all carbon etc.)
                # split into two SC3 groups of three
                bead_key = find_bead(martini_dict, 3, "CCC")
                g1, g2 = array1[:3], array1[3:]
                assign_bead(final, bead_key, *g1)
                assign_bead(final, bead_key, *g2)
    
        if len(array1) == 3:
            # SC3 case: val[0]=3, val[2]=="CCC"
            bead_key = find_bead(martini_dict, 3, "CCC")
            assign_bead(final, bead_key, *array1)
    
        elif len(array1) != 0:
            raise ValueError(f"Fallback for count=6 has {len(array1)} remaining — not implemented.")
    
    elif count == 5:
        # find the C-atom with exactly one mapped neighbor
        cand = None
        for atom in remaining:
            if final[atom[0]] != "":
                continue
            if atom[1].upper() != "C":
                continue
    
            # mapped neighbors
            mapped_nbrs = [section[t[0]] for t in atom[4]
                             if final[section[t[0]][0]] != ""]
    
            if len(mapped_nbrs) == 1:
                cand = (atom, mapped_nbrs[0])
                break

        if not cand:
            raise ValueError("Non-benzene 6-ring: no singleton C neighbor for count=5")
    
        atom, nbr = cand
        old_bead = final[nbr[0]]
    
        # find all indices in final with old_bead
        targets = [i for i, v in enumerate(final) if v == old_bead]
        new_bead = old_bead
        if new_bead.startswith('T'):
            new_bead = 'S' + new_bead[1:]
        else:
            raise ValueError("Non-benzene 6-ring: C neighbor next to no T beads")
    
        for i in targets:
            final[i] = new_bead
        final[atom[0]] = new_bead
    
        # now fall through to count=4 logic
        count = 4
        remaining = [atom for atom in section if final[atom[0]] == ""]
    
    if count == 4:
        # --- Special case: one atom is disconnected from the other three ---
        remaining_gids = {a[0] for a in remaining}
    
        # For each atom, see if it has at least one neighbor that is also in remaining
        connected_flags = {}
        for a in remaining:
            gid = a[0]
            has_neighbor_in_remaining = any(section[t[0]][0] in remaining_gids for t in a[4])
            connected_flags[gid] = has_neighbor_in_remaining
    
        isolated = [a for a in remaining if not connected_flags[a[0]]]
    
        if len(isolated) == 1:
            lone = isolated[0]
            triple = [a for a in remaining if a[0] != lone[0]]
    
            # Handle the connected triple using the same logic as count==3 (center + 2 neighbors)
            nb6_handle_count3(final, section, martini_dict, full_mapping, triple, assign_bead)
    
            # Then handle the isolated atom as count==1
            if final[lone[0]] == "":
                nb6_assign_lone_carbon(final, section, lone)
    
        else:
            # Original count==4 behavior: map any connected (atom, single-unmapped-neighbor) as a T-pair
            for atom in remaining:
                idx = atom[0]
                if final[idx] != "":
                    continue
    
                nbrs = [section[t[0]] for t in atom[4] if final[section[t[0]][0]] == ""]
                if len(nbrs) == 1:
                    # Only assign if the pair is actually connected (mutual), otherwise leave it for lower tiers.
                    _ = nb6_try_assign_connected_pair_T(section, final, martini_dict, atom, nbrs[0], assign_bead)
    
        remaining = [atom for atom in section if final[atom[0]] == ""]
        count = len(remaining)
    
    if count == 3:
        nb6_handle_count3(final, section, martini_dict, full_mapping, remaining, assign_bead)
        remaining = [atom for atom in section if final[atom[0]] == ""]
        count = len(remaining)
    
    if count == 2:
        nb6_handle_count2(final, section, martini_dict, remaining, assign_bead)
        remaining = [atom for atom in section if final[atom[0]] == ""]
        count = len(remaining)
    
    if count == 1:
        nb6_assign_lone_carbon(final, section, remaining[0])
        
    merge_phenol_to_diol(section, final, martini_dict, full_mapping)
    return final

# =============================================================================
# TASK 3.c: Mapping Non-Benzene 5-Ring Section
# =============================================================================
def map_nonbenzene_5_ring_section(section: List[List[Any]],
                                  final: List[str],
                                  martini_dict: Dict[str, List[Any]],
                                  full_mapping: List[List[List[Any]]]) -> List[str]:
    """
    Map a non-benzene 5-membered ring section.
    
    Pseudocode:
      (A) Build array0:
          For every atom that has:
            - element (index 1) is 'C' (or 'c')
            - final[atom[0]] is ""
            - AND it has at least one inner neighbor (from its inner connections at index 4)
              such that the neighbor's element is 'C' or 'N' AND the bond equals 2,
          add the pair (atom[0], neighbor[0]) as a tuple to array0.
      For each tuple in array0:
          - If both atoms are C: 
              if both atoms each have a foreign neighbor with the foreign sections of size 1: 
                  return an error saying 5-ring section is too difficult to be mapped and add a todo comment for later 
              if one atom has a foreign neighbor with the foreign section of size 1: 
                  assign a bead based on the foreign atom type (using the first outer connection): 
                      O: "SN6", N: "SN6d", S: "SC6", Cl: "SX3", I: "X1", C: "SC4", Br: SX2 
                      assign it + random string to all 3 atoms  
              else:
                assign bead "TC5" + random string. assign bead "TC5" + random string.
          - If one atom is N: assign bead "TN6a" + random string.
      
      (B) Build array1:
          For every atom that is unmapped (final[atom[0]] == "") and its element is not 'C', add its global index.
          Then:
           - If length of array1 is 1:
                * For the single atom:
                    - If element is S:
                        * If both of its first two inner neighbors are unmapped, generate a random string and assign bead "SC6" + that string to the atom and both neighbors.
                        * Else, assign "TC6" + random string only to the atom.
                    - If element is NH: 
                        * If both of its first two inner neighbors are unmapped, generate a random string and assign bead "TN6d" + that string to the atom and only 1 neighbor.
                        * Else, assign "TN6d" + random string to the only atom.
                    - If element is O: assign "TN4a" + random string to the atom and (only) to its first inner neighbor.
                    - If element is N and it has a foreign connection (len(atom[3]) > 0): 
                        Use the first outer connection to look up the foreign atom in full_mapping and assign "TN1" + random string to both.
           - If length of array1 is 2:
                * If both atoms have element O:
                    - Find a common neighbor (via inner connections) and assign "SN5a" + random string to both atoms and that common neighbor.
      
      (C) Build array2:
          For every atom that is unmapped and has element 'C', add its global index.
          We assume array2 will have either 5, 3, or 2 elements.
           - If length is 5:
                * Loop through array2, and find an atom that has:
                      - element 'C' (or 'c')
                      - at least one outer connection (len(atom[3]) != 0)
                  If found:
                      - If the foreign neighbor (from the first outer connection) has element 'C' and the size of its section is 1,
                        assign "SC3" + random string to that atom, the foreign neighbor, and its first inner neighbor.
                      - Else if the foreign neighbor has element 'O' and its section size is 1,
                        assign "SN6" + random string similarly.
                      - Else if the foreign neighbor’s section size is > 1, assign "TC3" + random string to that atom and its first inner neighbor.
                  If not found:
                      - take any connecting 2 atoms and assign TC3 + random string to those two atoms
                * After this, loop through array2 and remove atoms that are now mapped.
           - If the remaining size is 2 and all atoms have element 'C', assign "TC3" + random string.
           - If the remaining size is 3 and all atoms have element 'C', assign "SC3" + random string.
      
      (D) Final check: if any atom in the section remains unmapped, raise an error.
      
    Returns final
    """
    map_adjacent_SS_pairs(section, final, martini_dict)
    atoms_by_gid, gid_to_local = build_section_index(section)
    
    # --- nH HANDLER (must run before STEP A) ---
    # If an unmapped nitrogen has exactly 1 attached H AND it has at least one INNER bond of order 1.5,
    # assign it immediately as TN6d.
    for atom in section:
        g = atom[0]
        if final[g] != "":
            continue
    
        if atom[1].upper() != "N":
            continue
    
        num_H = atom[6]  # [global_index, element, ring_status, outer_connection, inner_connection, isedge, num_H]
        if num_H != 1:
            continue
    
        # NEW requirement: N must have an inner aromatic bond (bond order 1.5) to qualify
        has_inner_15 = any(bo == 1.5 for _, bo in atom[4])
        if not has_inner_15:
            continue
    
        bead_key = find_bead(martini_dict, 5, "NH")
        assign_bead(final, bead_key, g)
    # -----------------------
    # (A) Process double bonds among C and N atoms.
    # --- STEP A (revised): graph-based matching + required C/N pattern handling ---
    # Build a candidate graph on unmapped C/N atoms involved in (2 or 1.5) inner bonds,
    # then pick an edge-matching that (1) maximizes pairing, and (2) prefers CC pairs when possible.
    cand, elem_of, neighbors = build_cn_candidate_graph(section, final)
    pairs = choose_pairs_for_nonbenzene_5_ring(
        cand=cand,
        elem_of=elem_of,
        neighbors=neighbors,
        final=final,
    )

    # 2) assign beads for each pair
    for a, b in pairs:
        atomA = atoms_by_gid[a]
        atomB = atoms_by_gid[b]
        eA, eB = atomA[1].upper(), atomB[1].upper()
        
        if eA == "C" and eB == "C":
            # look for size-1 foreign neighbors
            fA = next((t for t in atomA[3] if len(full_mapping[t[0]]) == 1), None)
            fB = next((t for t in atomB[3] if len(full_mapping[t[0]]) == 1), None)
    
            if fA and fB:
                # two separate T-beads: one for each (inner↔foreign)
                for atom_inner, ft in ((atomA, fA), (atomB, fB)):
                    fa = full_mapping[ft[0]][ft[1]]
                    bo = ft[2]
                    key = pick_bead_key(martini_dict, fa[1].upper(), bo, kind="T")
                    if key is None:
                        raise ValueError(f"No T-key for foreign {fa[1].upper()} on double bond")
                    assign_bead(final, key, atom_inner[0], fa[0])
    
            elif fA or fB:
                # exactly one small foreign → S-bead trio
                ft = fA or fB
                fa = full_mapping[ft[0]][ft[1]]
                bo = ft[2]
                key = pick_bead_key(martini_dict, fa[1].upper(), bo, kind="S")
                if key is None:
                    raise ValueError(f"No S-key for foreign {fa[1].upper()} on double bond")
                assign_bead(final, key, a, b, fa[0])
            else:
                # no small foreign → default CC-double bead from dict
                key = find_bead(martini_dict, 1, "CC")
                if key is None:
                    raise ValueError("5-ring CC double not mappable (no bead for CC)")
                assign_bead(final, key, a, b)
    
        # ---- C–N ----
        elif ("C" in (eA,eB)) and ("N" in (eA,eB)):
            key = find_bead(martini_dict, 1, "CN")
            if key is None:
                raise ValueError("5-ring CN double not mappable (no bead for CN)")
            assign_bead(final, key, a, b)
    
        # ---- N–N ----
        else:
            # both must be N
            key = find_bead(martini_dict, 1, "NN")
            if key is None:
                raise ValueError("5-ring NN double not mappable (no bead for NN)")
            assign_bead(final, key, a, b)
            
    # -----------------------
    # (B) Build array1 for atoms not of type C.
    array1 = [atom[0] for atom in section if final[atom[0]] == "" and atom[1].upper() != "C"]
    if len(array1) == 1:
        idx = array1[0]
        atom = atoms_by_gid[idx]
    
        el = atom[1].upper()
        # S-center
        if el == "S":
            # if both first two inner neighbors unmapped → SC6, else TC6
            key = None
            if len(atom[4]) >= 2:
                inn1 = section[atom[4][0][0]]
                inn2 = section[atom[4][1][0]]
                if final[inn1[0]] == "" and final[inn2[0]] == "":
                    key = find_bead(martini_dict, 5, "CSC")   # SC6 → [0]=5, [2]="CSC"
            if key is None:
                key = find_bead(martini_dict, 5, "S")        # TC6  → [0]=5, [2]="S"
            if key is None:
                raise ValueError("No bead found for sulfur in 5-ring special case")
                # if SC6, also tag those two
            if key != find_bead(martini_dict, 5, "S"):
                assign_bead(final, key, atom[0], inn1[0], inn2[0])
            else:
                assign_bead(final, key, atom[0])
    
        # double-bond N(C) center
        elif el == "N" and len(atom[3]) > 0:
            # TN3 → [0]=6, [2]="N(C)"
            key = find_bead(martini_dict, 6, "N(C)")
            if key is None:
                raise ValueError("No bead found for N(C) in 6-ring special case")
            foreign = full_mapping[atom[3][0][0]][atom[3][0][1]]
            assign_bead(final, key, atom[0], foreign[0])
            
        # NH center
        elif el == "NH" or el == "N":
            # TN6d → [0]=5, [2]="NH"
            key = find_bead(martini_dict, 5, "NH")
            if key is None:
                raise ValueError("No bead found for NH in 5-ring special case")
            # if two inners available and both unmapped → tag one, else just center
            assign_bead(final, key, atom[0])
            
            if len(atom[4]) >= 2:
                inn1 = section[atom[4][0][0]]
                inn2 = section[atom[4][1][0]]
                if final[inn1[0]] == "" and final[inn2[0]] == "":
                    final[inn1[0]] = final[atom[0]]
    
        # O-center
        elif el == "O":
            # TP6a → [0]=5, [2]="OC"  for two unmapped inners → tag left
            # TN3a → [0]=5, [2]="O"   otherwise
            bead_key = None
            if len(atom[4]) >= 2:
                inn1 = section[atom[4][0][0]]
                inn2 = section[atom[4][1][0]]
                if final[inn1[0]] == "" and final[inn2[0]] == "":
                    bead_key = find_bead(martini_dict, 5, "OC")   # TP6a
            if bead_key is None:
                bead_key = find_bead(martini_dict, 5, "O")      # TN3a
            assign_bead(final, bead_key, atom[0])
            if bead_key == find_bead(martini_dict, 5, "OC"):
                # only the left neighbor
                final[inn1[0]] = final[atom[0]]  
    
        else:
            raise ValueError(f"Unrecognized single-atom case for element {el}")

    elif len(array1) == 2:
        # If both atoms are O.
        atom1 = next(a for a in section if a[0] == array1[0])
        atom2 = next(a for a in section if a[0] == array1[1])
        if atom1[1].upper() == "O" and atom2[1].upper() == "O":
            # find their common inner neighbor
            common = None
            for tup in atom1[4]:
                n1 = section[tup[0]]
                for tup2 in atom2[4]:
                    n2 = section[tup2[0]]
                    if n1[0] == n2[0]:
                        common = n1[0]
                        break
                if common is not None:
                    break
    
            if common is not None:
                # lookup SN5a in martini_dict: v[0]==4, v[2]=="OCO" (order-insensitive)
                bead_key = find_bead(martini_dict, 4, "OCO")
                assign_bead(final, bead_key, atom1[0], atom2[0], common)
                
    # -----------------------
    # (C) Build array2: for unmapped Cs.
    array2 = [atom[0] for atom in section if final[atom[0]] == "" and atom[1].upper() == "C"]
    if len(array2) == 5:
        # Try to handle the “one C with an outer” case exactly once.
        assigned = False
        for idx in array2:
            atom = atoms_by_gid[idx]
            if atom[3]:  # has at least one outer connection
                foreign_tup  = atom[3][0]
                foreign_sec  = full_mapping[foreign_tup[0]]
                foreign_atom = foreign_sec[foreign_tup[1]]
                
                if len(foreign_sec) == 1:
                    f_el = foreign_atom[1].upper()
                
                    if f_el == "O":
                        bo = foreign_tup[2]  # bond order to the foreign atom
                        if bo == 1:
                            pattern = "CC(O)"
                        elif bo == 2:
                            pattern = "CC(=O)"
                        else:
                            raise ValueError(f"Unsupported C–O bond order for 5-ring case: {bo}")
                    elif f_el == "C":
                        pattern = "CC(C)"
                    elif f_el == "N":
                        pattern = "CC(N)"
                    else:
                        raise ValueError("No section-6 bead found for pattern 'CC(?)'")
                
                    bead_key = find_bead(martini_dict, 6, pattern)
                else:
                    bead_key = find_bead(martini_dict, 5, "CC")
                
                # assign to the trio: this C, its foreign neighbour, and its first inner neighbour
                targets = [atom[0], foreign_atom[0]]
                if atom[4]:
                    inn = section[atom[4][0][0]]
                    targets.append(inn[0])
                assign_bead(final, bead_key, *targets)
                assigned = True
                break  # << don’t assign again!
                
        for idx in array2:
            if not assigned:
                # fallback: pair any two connected Cs once, with TC3 looked up from martini_dict
                for i, idx1 in enumerate(array2):
                    atom1 = next(a for a in section if a[0] == idx1)
                    for idx2 in array2[i+1:]:
                        if any(section[t][0] == idx2 for t,_ in atom1[4]):
                            # lookup the TC3 bead key: v[0]==5 and v[2]=="CC"
                            bead_key = find_bead(martini_dict, 5, "CC")
                            assign_bead(final, bead_key, idx1, idx2)
                            assigned = True
                            break
                    if assigned:
                        break
                if not assigned:
                    raise ValueError(
                       "Could not find any two connecting atoms in 5-ring array2 for TC3 assignment."
                    )
    
        # now remove those atoms from array2
        array2 = [idx for idx in array2 if final[idx] == ""]
    # At this point, array2 contains only the *unmapped* Cs after above.
    # If exactly two remain, pair them using TC3 from martini_dict; if three remain, use SC3.
    if len(array2) == 2:
        bead_key = find_bead(martini_dict, 5, "CC")
        assign_bead(final, bead_key, *array2)
        array2.clear()
    
    elif len(array2) == 3:
        bead_key = find_bead(martini_dict, 3, "CCC")
        assign_bead(final, bead_key, *array2)
        array2.clear()
        
    # -----------------------
    # (D) Final check with merging fallback
    for atom in section:
        idx = atom[0]
        if final[idx]=='' and atom[1].upper() == "C":
            # gather neighbors: inner + outer
            neigh = [section[t[0]][0] for t in atom[4]] + [full_mapping[t[0]][t[1]][0] for t in atom[3]]
            # try TC5 merge
            merged=False
            for prefix in ('TC5',):
                for n in neigh:
                    bead_str = final[n]
                    if bead_str.startswith(prefix):
                        # clear old
                        old = bead_str
                        ids = [i for i,v in enumerate(final) if v==old]
                        for i in ids: final[i] = ''
                        # pick new S-key
                        bond = atom[4][0][1] if atom[4] else None
                        key = pick_bead_key(martini_dict, atom[1].upper(), bond, kind='S')
                        if not key:
                            raise ValueError("Fallback mapping failed: no S-key for element")
                        bead = key + generate_random_string()
                        for i in ids: final[i]=bead
                        # apply to atom & neighbor
                        final[idx]=bead; final[n]=bead
                        merged=True; break
                if merged: break
            if merged: continue
            # try T->S convert
            for n in neigh:
                bead_str = final[n]
                if bead_str.startswith('TN6'): 
                    new = 'SN4' + bead_str[3:]
                    ids=[i for i,v in enumerate(final) if v==bead_str]
                    for i in ids: final[i]=new
                    final[idx]=new; final[n]=new
                    merged=True; break
                elif bead_str.startswith('T'):
                    new = 'S' + bead_str[1:]
                    ids=[i for i,v in enumerate(final) if v==bead_str]
                    for i in ids: final[i]=new
                    final[idx]=new; final[n]=new
                    merged=True; break
            if merged: continue
            # try S propagate
            for n in neigh:
                bead_str = final[n]
                if bead_str.startswith('S'):
                    ids=[i for i,v in enumerate(final) if v==bead_str]
                    for i in ids: final[i]=bead_str
                    final[idx]=bead_str; final[n]=bead_str
                    merged=True; break
            if merged: continue
            raise ValueError("Non-benzene 5-ring section not fully mappable!")
        elif final[idx]=='':
            raise ValueError("Non-benzene 5-ring section not fully mappable!")
    merge_phenol_to_diol(section, final, martini_dict, full_mapping)
    return final

# =============================================================================
# TASK 3.d: Determine if a Non-Ring Section Is 1-Bead Mappable
# =============================================================================
def is_non_ring_section_1bead_mappable(section: List[List[Any]]) -> bool:
    """
    Determine whether a non-ring section (section type == 0) is mappable as a single bead.

    Pseudocode:
      - Collect all edge atoms (atom[5] == True) into an array.
      - If exactly one edge atom exists, return True.
      - If none, throw an error.
      - If more than one, for each pair of edge atoms compute the shortest path (using BFS on the inner connections).
        If any such path is longer than 4 atoms, return False; otherwise, return True.
    """
    edge_locals = [i for i, atom in enumerate(section) if atom[5]]
    if len(edge_locals) == 1:
        return True
    if len(edge_locals) == 0:
        raise ValueError("Cannot map non-ring section: no edge atoms found!")

    graph = build_inner_graph_local(section)

    # Original behavior measured path length in *atoms* (nodes), starting at 1.
    # That is equivalent to (edge-distance + 1).
    for i in range(len(edge_locals)):
        for j in range(i + 1, len(edge_locals)):
            edge_dist = bfs_distance_graph(graph, edge_locals[i], edge_locals[j])
            if edge_dist + 1 > 4:
                return False
    return True

# =============================================================================
# TASK 3.e: Map Non-Ring Section Using 1–Bead Mapping
# =============================================================================
def map_non_ring_section_1bead(section: List[List[Any]],
                               final: List[str],
                               martini_dict: Dict[str, List[Any]],
                               full_mapping: List[List[List[Any]]]) -> List[str]:
    """
    Map a non‐ring section (section type == 0) that is 1‐bead mappable.
    
    This implementation follows the pseudocode:
      - Build array0 as the list of global indices of edge atoms (atom[5] == True).
      - Depending on whether there are 2, 3, or 4 edge atoms, trace a SMILES‐like path
        (using inner connections from index 4) and then look for candidate bead keys in martini_dict.
      - For a 2-edge case, we trace a single linear path.
      - For a 3-edge case, we locate the atom with 3 inner connections and trace each branch.
      - For a 4-edge case, we trace each branch from the atom with 4 inner connections and then
        assign special bead names (“SX4e” to the branch edges and “TP1d” to the remaining atoms)
        if the unmapped atoms are a C–O pair.
        
    If exactly one candidate bead is found (or a special case is met) the bead name (with a random tag)
    is assigned to all atoms in the section. Otherwise, an error is raised.
    """
    # ----------------------------------------------------------------
    # Build cache: global index -> (section index, local index)
    global_to_sec_loc: Dict[int, (int,int)] = {}
    for sec_idx, sec in enumerate(full_mapping):
        for loc_idx, atom_info in enumerate(sec):
            gi = atom_info[0]
            global_to_sec_loc[gi] = (sec_idx, loc_idx)
    # ----------------------------------------------------------------
    # --- Helper variables ---
    edge_globals = [a[0] for a in section if a[5]]
    num_edges = len(edge_globals)
    # --- fallback for single-edge sections ---
    # … earlier in map_non_ring_section_1bead …
    # global_to_sec_loc: Dict[int, (int,int)] = { … }
    if len(section) == 1:
        idx  = section[0][0]
        gi   = section[0][0]
        atom = section[0]

        # gather global neighbor indices (inner + outer)
        neigh = [ section[t[0]][0]           for t in atom[4] ] + \
                [ full_mapping[t[0]][t[1]][0] for t in atom[3] ]

        # Try merging into existing TC5/TC3 beads first
        for prefix in ("TC5","TC3"):
            for n in neigh:
                b = final[n]
                if b.startswith(prefix):
                    # clear the old group
                    old_tag = b
                    ids = [i for i,v in enumerate(final) if v == old_tag]
                    for i in ids:
                        final[i] = ""
                    # pick new S‐type key
                    bond = atom[4][0][1] if atom[4] else None
                    key  = pick_bead_key(martini_dict,
                                         atom[1].upper(),
                                         bond,
                                         kind='S')
                    if not key:
                        raise ValueError(f"Fallback: no S-key for {atom[1]} with bond {bond}")
                    bead = key + generate_random_string()
                    for i in ids:
                        final[i] = bead
                    final[gi] = final[n] = bead
                    return final

        # If it’s not C, we can’t handle it and therefore we will assign a provisional bead.
        if atom[1].upper() != "C":
            raise ValueError("1-edge non-ring not mappable (non-C)")
        else:
            # convert T* to S* (or borrow S*)
            for n in neigh:
                b = final[n]
                if b.startswith("T"):
                    new = "S" + b[1:]
                    ids = [i for i, v in enumerate(final) if v == b]
                    for i in ids:
                        final[i] = new
                    final[idx] = final[n] = new
                    return final
        
                elif b.startswith("S"):
                    ids = [i for i, v in enumerate(final) if v == b]
                    for i in ids:
                        final[i] = b
                    final[idx] = final[n] = b
                    return final
        
            raise ValueError("1-edge non-ring still not mappable")
      
    # CASE 1: 2 edge atoms.
    if num_edges == 2:
        start_global = edge_globals[0]
        start_local = next(i for i, atom in enumerate(section) if atom[0] == start_global)
        path_str = trace_linear_path(section, start_local).upper()

        # 1) Try section 7
        candidate_keys = [
            key for key, val in martini_dict.items()
            if val[0] == 7 and val[1] == 0
               and (val[2].upper() == path_str or val[2].upper() == path_str[::-1])
        ]
        candidate_keys = list(set(candidate_keys))

        # 2) Fallback to section 11 if needed
        if not candidate_keys:
            # ---- new: special split logic for 4‐atom sections ----
            # 2‑edge case
            if len(section) == 4:
                print(f"Linear path of 4 cannot be mapped: no candidate bead found for 2-edge mapping (path: {path_str}). Beginning to break down further.")
                # Build a little helper to map exactly two atoms by connectivity
                start_local = next(i for i,a in enumerate(section) if a[0] == edge_globals[0])
                nbr_local, _ = section[start_local][4][0]
                # Identify the other two locals
                used = { section[start_local][0], section[nbr_local][0] }
                others = [i for i,a in enumerate(section) if a[0] not in used]
                if len(others) != 2:
                    raise ValueError("Expected exactly 2 leftover atoms when splitting 4‑atom 2‑edge case.")
                # Map each sub‑pair
                map_two_atoms_by_connectivity(section, start_local, nbr_local, martini_dict, final, generate_random_string)
                map_two_atoms_by_connectivity(section, others[0], others[1], martini_dict, final, generate_random_string)
                return final

            else:
                candidate_keys = [
                    key for key, val in martini_dict.items()
                    if val[0] == 11 and val[1] == 0
                       and (val[2].upper() == path_str or val[2].upper() == path_str[::-1])
                ]
                candidate_keys = list(set(candidate_keys))
                if candidate_keys:
                    warnings.warn(
                        f"Warning: falling back to section 11 for 2-edge mapping (path: {path_str})",
                        UserWarning
                    )

        # 3) Commit or error
        if len(candidate_keys) == 1:
            rstr = generate_random_string()
            bead = candidate_keys[0] + rstr
            for atom in section:
                final[atom[0]] = bead
            return final

        # no candidates even after fallback
        if not candidate_keys:
            raise ValueError(
                f"Non-ring section cannot be mapped: no candidate bead found for 2-edge mapping (path: {path_str})."
            )
        else:
            # If more than one candidate remains, resolve it.
            # For CO/OC we use the foreign bead prefix rule first (SX4e -> TP1d, TC5 -> TN2a),
            # otherwise we disambiguate using the val[5] hydrogen-count rule.
            foreign_bead = None
            if path_str in ("CO", "OC"):
                for atom in section:
                    if atom[0] in edge_globals and atom[3]:
                        foreign_sec_index = atom[3][0][0]
                        index = atom[3][0][1]
                        foreign_atom = full_mapping[foreign_sec_index][index]
                        foreign_bead = final[foreign_atom[0]]
                        break

            bead_key = resolve_unique_candidate(
                candidate_keys, martini_dict, section, path_str=path_str, foreign_bead=foreign_bead
            )
            rstr = generate_random_string()
            for a in section:
                final[a[0]] = bead_key + rstr
            return final

    # CASE 2: 3 edge atoms.
    elif num_edges == 3:
        # Find the atom with 3 inner connections.
        center_local = None
        for i, atom in enumerate(section):
            if len(atom[4]) == 3:
                center_local = i
                break
        if center_local is None:
            raise ValueError("Non‐ring section with 3 edges does not have an atom with 3 inner connections.")
        # Trace each branch from the center atom.
        branch_paths = []
        for (nbr, _) in section[center_local][4]:
            branch_str = trace_branch(section, center_local, nbr)
            branch_paths.append(branch_str)

        # 1) Try section 7 matching
        candidate_keys = find_branch_candidates(
            martini_dict, 7, section[center_local][1], branch_paths
        )

        # 2) Fallback to section 11 if needed
        if not candidate_keys:
            candidate_keys = find_branch_candidates(
                martini_dict, 11, section[center_local][1], branch_paths
            )

        # 3) Commit or error (resolve to exactly one candidate)
        if not candidate_keys:
            raise ValueError(
                f"Non‐ring section cannot be mapped: no candidate bead found for 3-edge mapping "
                f"(center: {section[center_local][1]} and branches: {branch_paths})."
            )

        bead_key = resolve_unique_candidate(candidate_keys, martini_dict, section)

        rstr = generate_random_string()
        for atom in section:
            final[atom[0]] = bead_key + rstr
        return final

    # CASE 3: 4 edge atoms.
    elif num_edges == 4:
        center_local = None
        for i, atom in enumerate(section):
            if len(atom[4]) == 4:
                center_local = i
                break
        # 2) If no 4-connected center, look for exactly two atoms each with 3 inner bonds,
        #    of which exactly 2 neighbors are currently flagged as edges. Split those off.
        if center_local is None:
            three_centers = []
            for i, atom in enumerate(section):
                if len(atom[4]) == 3:
                    edge_nbrs = [nbr for (nbr, _) in atom[4] if section[nbr][5]]
                    if len(edge_nbrs) == 2:
                        three_centers.append((i, edge_nbrs))
            if len(three_centers) == 2:
                # Unpack the two 3-connected centers
                (c1, edges1), (c2, edges2) = three_centers

                build_and_assign_trio(section, c1, edges1, martini_dict, final, generate_random_string)
                # Map second trio (c2 with edges2)
                build_and_assign_trio(section, c2, edges2, martini_dict, final, generate_random_string)

                return final
            else:
                raise ValueError("Non-ring section with 4 edges has no 4-connected center or valid pair of 3-connected centers.")
        total_atoms = len(section)
        if total_atoms == 6:
            two_branch = None
            one_branches = []
            for (nbr_local, _) in section[center_local][4]:
                final_edge = get_final_edge(section, center_local, nbr_local)
                if final_edge is not None:
                    path = bfs_path(section, center_local, final_edge)
                    if len(path) == 3:
                        two_branch = path   # [center_local, intermediate_local, final_edge]
                    else:
                        one_branches.append((center_local, final_edge))
        
            if two_branch is None or len(one_branches) != 3:
                raise ValueError("4-edge 6-atom split: could not identify branches correctly.")
        
            _, int_local, edge_local = two_branch
        
            # 1) Remove the inner-bond entry from center<->intermediate, but
            #    “move” it as an outer-connection on `int_local` pointing back to center.
            #    First, find the bond order between them:
            bond_to_center = None
            for (nbr, bo) in section[int_local][4]:
                if nbr == center_local:
                    bond_to_center = bo
                    break
            # Remove from int_local’s inner list:
            section[int_local][4] = [
                (nbr, bo) for (nbr, bo) in section[int_local][4]
                if nbr != center_local
            ]
            # Remove from center_local’s inner list:
            section[center_local][4] = [
                (nbr, bo) for (nbr, bo) in section[center_local][4]
                if nbr != int_local
            ]
        
            # Now “move” that edge into int_local’s outer list. We need the index of this section in full_mapping:
            sec_idx = full_mapping.index(section)
            section[int_local][3].append((sec_idx, center_local, bond_to_center))
        
            # 2) Build sub2 = [center_local] + three single leaves, recurse on sub2 first.
            sub2_locals = [center_local] + [loc for (_, loc) in one_branches]
            sub2 = [section[i] for i in sub2_locals]
            # Mark only the three leaf nodes as edges; center_local stays isedge=False.
            for (_, leaf_loc) in one_branches:
                section[leaf_loc][5] = True
            sub2_reindexed = reindex_section(sub2, sub2_locals)
            final = map_non_ring_section_1bead(sub2_reindexed, final, martini_dict, full_mapping)
        
            # 3) Now build sub1 = [int_local, edge_local], recurse on sub1.
            #    int_local already has [3] updated to list its outer = center.
            section[int_local][5] = True
            section[edge_local][5] = True
            sub1_locals = [int_local, edge_local]
            sub1 = [section[i] for i in sub1_locals]
            sub1_reindexed = reindex_section(sub1, sub1_locals)
            final = map_non_ring_section_1bead(sub1_reindexed, final, martini_dict, full_mapping)
        
            return final

        # Part 2: if section has 5 atoms, treat like 3-edge but with val[1]==2
        elif total_atoms == 5:
            # All four neighbors of center are direct edges
            branch_paths = []
            for (nbr, _) in section[center_local][4]:
                if section[nbr][5]:
                    branch_str = trace_branch(section, center_local, nbr)
                    branch_paths.append(branch_str)
            if len(branch_paths) != 4:
                raise ValueError(f"Non‐ring 4-edge 5-atom mapping: expected 4 branches, found {len(branch_paths)}.")
            import re
            candidate_keys = []
            for key, val in martini_dict.items():
                if val[0] == 7 and val[1] == 2 and "(" in val[2]:
                    prefix, remainder = val[2].split("(", 1)
                    if prefix.strip() != section[center_local][1]:
                        continue
                    branch_segments = re.findall(r'\(([^)]+)\)', val[2])
                    if len(branch_segments) != len(branch_paths):
                        continue
                    unmatched = branch_segments.copy()
                    matched_all = True
                    for bp in branch_paths:
                        for bs in list(unmatched):
                            if bp == bs:
                                unmatched.remove(bs)
                                break
                        else:
                            matched_all = False
                            continue
                    if val[4] == "." and len(unmatched) <= 2:
                        candidate_keys.append(key)
                        break
                    if matched_all:
                        candidate_keys.append(key)
            candidate_keys = list(set(candidate_keys))
            if len(candidate_keys) == 1:
                rstr = generate_random_string()
                for atom in section:
                    final[atom[0]] = candidate_keys[0] + rstr
                return final
            elif not candidate_keys:
                warnings.warn(
                    f"Warning: falling back to section 11 for 4-edge mapping (center: {section[center_local][1]} and branches: {branch_paths})",
                    UserWarning
                )
                for key, val in martini_dict.items():
                    if val[0] == 11 and val[1] == 2 and "(" in val[2]:
                        prefix, remainder = val[2].split("(", 1)
                        if prefix.strip() != section[center_local][1]:
                            continue
                        branch_segments = re.findall(r'\(([^)]+)\)', val[2])
                        if len(branch_segments) != len(branch_paths):
                            continue
                        unmatched = branch_segments.copy()
                        matched_all = True
                        for bp in branch_paths:
                            for bs in list(unmatched):
                                if bp == bs:
                                    unmatched.remove(bs)
                                    break
                            else:
                                matched_all = False
                                continue
                        if val[4] == "." and len(unmatched) <= 2:
                            candidate_keys.append(key)
                            break
                        if matched_all:
                            candidate_keys.append(key)
                candidate_keys = list(set(candidate_keys))
                if len(candidate_keys) == 1:
                    rstr = generate_random_string()
                    for atom in section:
                        final[atom[0]] = candidate_keys[0] + rstr
                    return final
                elif not candidate_keys:
                    raise ValueError(
                        f"Non‐ring section cannot be mapped: no candidate bead found for 4-edge mapping (center: {section[center_local][1]} and branches: {branch_paths})."
                    )
                else:
                    raise ValueError(
                        f"Non‐ring section cannot be mapped: ambiguous candidate keys for 4-edge mapping (center: {section[center_local][1]} and branches: {branch_paths})."
                    )
            else:
                raise ValueError(
                    f"Non‐ring section cannot be mapped: ambiguous candidate keys for 4-edge mapping (center: {section[center_local][1]} and branches: {branch_paths})."
                )
        else:
            raise ValueError(f"Non‐ring 4-edge mapping only supported for 5 or 6 atoms; found {total_atoms}.")

    # CASE 4: 5 or 6 edge atoms.
    elif num_edges > 4:
        # find the two centers (atoms with >2 inner connections)
        centers = [i for i, atom in enumerate(section) if len(atom[4]) > 2]
        if len(centers) != 2:
            raise ValueError(f"Need exactly 2 centers to split; found {len(centers)}")
    
        c0, c1 = centers
    
        # disconnect the two centers in this section
        section[c0][4] = [(nbr, bo) for nbr, bo in section[c0][4] if nbr != c1]
        section[c1][4] = [(nbr, bo) for nbr, bo in section[c1][4] if nbr != c0]
    
        # make a local copy of the mapping so we don't clobber the caller's
        local_mapping = full_mapping[:]
        try:
            sec_idx = local_mapping.index(section)
        except ValueError:
            sec_idx = None
    
        # trace out each sub‑section by DFS (they’re now disconnected)
        sub_locals_list = []
        for center in centers:
            visited = {center}
            stack   = [center]
            sub_locs = []
            while stack:
                cur = stack.pop()
                sub_locs.append(cur)
                for nbr, _ in section[cur][4]:
                    if nbr not in visited:
                        visited.add(nbr)
                        stack.append(nbr)
            sub_locals_list.append(sorted(sub_locs))
    
        # rebuild each sub‑section with correct local indices
        subs = []
        for locs in sub_locals_list:
            # extract the raw atoms
            raw = [section[i] for i in locs]
            # reindex so raw[*][4] refers to positions within this small list
            sub = reindex_section(raw, locs)
            # recompute isedge flags
            fix_isedge(sub)
            subs.append(sub)
    
        # splice them into our local mapping (if we know where this section lived)
        if sec_idx is not None:
            local_mapping.pop(sec_idx)
            # insert in reverse so subs[0] ends up at sec_idx
            local_mapping.insert(sec_idx, subs[1])
            local_mapping.insert(sec_idx, subs[0])
    
        # recurse on each piece using our updated local_mapping
        for sub_section in subs:
            final = map_non_ring_section_1bead(
                sub_section,
                final,
                martini_dict,
                local_mapping
            )
        return final
    # Should not be reached.
    raise ValueError("Non‐ring section cannot be mapped: unhandled edge count case.")
    return final


# =============================================================================
# TASK 3.f: Map Non-Ring Sections That Are Too Long (Subdivision)
# =============================================================================
def fix_isedge(section: List[List[Any]]) -> None:
    """
    Recalculate the isedge flag for each atom in a section.
    An atom is considered an edge if it has 0 or 1 inner connections (index 4) within the section.
    """
    for atom in section:
        if len(atom[4]) <= 1:
            atom[5] = True
        else:
            atom[5] = False

def rebuild_all_outer_connections(full_mapping: List[List[List[Any]]],
                                  old_pair_to_gid: Dict[Tuple[int,int], int]) -> None:
    """
    Rewrites every atom[3] = [(sec_idx, local_idx, bond_order), ...] so that
    sec_idx/local_idx are correct after full_mapping has been structurally changed.

    Uses global atom id as the stable key.
    """
    # Build new gid -> (sec_idx, local_idx)
    gid_to_new_pair: Dict[int, Tuple[int,int]] = {}
    for si, sec in enumerate(full_mapping):
        for li, atom in enumerate(sec):
            gid_to_new_pair[atom[0]] = (si, li)

    # Rewrite all outer tuples everywhere
    for sec in full_mapping:
        for atom in sec:
            new_outer = []
            for (fsi, fli, bo) in atom[3]:
                # Translate old (sec,local) -> foreign gid (from before the split)
                fgid = old_pair_to_gid.get((fsi, fli))
                if fgid is None:
                    # Defensive: if something is missing, keep original tuple
                    # (or raise if you prefer strictness)
                    new_outer.append((fsi, fli, bo))
                    continue

                # Translate foreign gid -> new (sec,local)
                if fgid not in gid_to_new_pair:
                    raise ValueError(f"Outer-connection target global id {fgid} vanished after remap.")
                nsi, nli = gid_to_new_pair[fgid]
                new_outer.append((nsi, nli, bo))

            atom[3] = new_outer

def subdivide_non_ring_section_normal(section: List[List[Any]], 
                                      final: List[str],
                                      martini_dict: Dict[str, List[Any]],
                                      full_mapping: List[List[List[Any]]]) -> Tuple[List[List[Any]], List[List[Any]], List[int], List[int]]:
    """
    Normal subdivision using DFS from a single candidate edge.
    Returns:
      sub_section, remainder, sub_old_indices, rem_old_indices
    """
    # Build candidate_edges using local indices (atoms flagged as edge).
    candidate_edges = [i for i, atom in enumerate(section) if atom[5]]
    if not candidate_edges:
        candidate_edges = [0, len(section) - 1]
        
    # SPECIAL CASE for linear sections (exactly 2 edge nodes)
    if len(candidate_edges) == 2 and len(section) in (6, 9):
    
        start = candidate_edges[0]  # pick one edge to grow from
    
        # Find a 3-node path: start -> n1 -> n2
        # Because it's "linear", each internal node should have degree ~2 within this subgraph.
        sub_old_indices = [start]
    
        # Step 1: pick a neighbor (prefer one that is NOT the other edge if possible)
        nbrs0 = [nbr for (nbr, _bo) in section[start][4]]
        if not nbrs0:
            raise ValueError("Linear special-case: start edge has no inner neighbors")
    
        other_edge = candidate_edges[1]
        n1 = None
        for n in nbrs0:
            if n != other_edge:
                n1 = n
                break
        if n1 is None:
            n1 = nbrs0[0]  # fallback (only neighbor is the other edge)
    
        sub_old_indices.append(n1)
    
        # Step 2: from n1, pick the next neighbor that isn't start
        nbrs1 = [nbr for (nbr, _bo) in section[n1][4]]
        n2 = None
        for n in nbrs1:
            if n != start:
                n2 = n
                break
        if n2 is None:
            raise ValueError("Linear special-case: cannot extend path to 3 atoms")
    
        sub_old_indices.append(n2)
    
        # remainder = everything else
        rem_old_indices = [j for j in range(len(section)) if j not in sub_old_indices]
    
        sub_section = [section[i] for i in sub_old_indices]
        remainder   = [section[i] for i in rem_old_indices]
        return sub_section, remainder, sub_old_indices, rem_old_indices

    # NEW SPECIAL CASE: look for any “center” atom with ≥ 4 inner connections and ≥ 3 of those neighbors flagged as edges
    for i, atom in enumerate(section):
        if len(atom[4]) >= 4:
            # gather inner‐neighbors that are currently edges (atom[5] == True)
            edge_nbrs = [nbr for (nbr, _) in atom[4] if section[nbr][5]]
            if len(edge_nbrs) >= 3:
                center_local   = i
                chosen_edges   = edge_nbrs[:3]   # pick the first 3 edge‐neighbors
                # build sub_old_indices = [center] + those 3 edges
                sub_old_indices = [center_local] + chosen_edges
                # remainder = “all other locals”
                rem_old_indices = [j for j in range(len(section)) if j not in sub_old_indices]
                sub_section = [section[idx] for idx in sub_old_indices]
                remainder   = [section[idx] for idx in rem_old_indices]
                return sub_section, remainder, sub_old_indices, rem_old_indices

    # If there are many candidate edges, attempt a BFS-based split that yields mappable pieces.
    if len(candidate_edges) > 2:
        graph = build_inner_graph_local(section)

        # Make a copy of candidate_edges before removal.
        candidate_edges_copy = candidate_edges[:]

        # Remove too-close candidate edges (distance <= 4 edges).
        to_remove = set()
        for i in range(len(candidate_edges)):
            for j in range(i + 1, len(candidate_edges)):
                if bfs_distance_graph(graph, candidate_edges[i], candidate_edges[j]) <= 4:
                    to_remove.add(candidate_edges[i])
                    to_remove.add(candidate_edges[j])
        candidate_edges = [i for i in candidate_edges if i not in to_remove]
        
        # NEW: Filter out candidate edges that are too close to a center node (score must be > 1)
        center_nodes = {i for i, atom in enumerate(section) if len(atom[4]) >= 3}
        
        # Only apply the filter if we actually have center nodes to score against
        if center_nodes:
            scored_all = [(e, dist_to_nearest_center(section, center_nodes, e)) for e in candidate_edges]
            candidate_edges = [e for (e, d) in scored_all if d > 1]
        # else: no center nodes -> leave candidate_edges as-is
        
        # If pruning left multiple candidate edges, keep only the one "furthest" from a center node.
        # A "center node" is any node with >= 3 inner connections.
        if len(candidate_edges) > 1:
            center_nodes = {i for i, atom in enumerate(section) if len(atom[4]) >= 3}

            scored = [(e, dist_to_nearest_center(section, center_nodes, e)) for e in candidate_edges]
            maxd = max(d for (_e, d) in scored)

            # Keep only one edge (stable tie-break: smallest local index)
            best_edges = [e for (e, d) in scored if d == maxd]
            candidate_edges = [min(best_edges)]
            
        if not candidate_edges:
            # --- New procedure starts here ---
            # Reset candidate_edges from the copy.
            candidate_edges = candidate_edges_copy

            # Build an array of candidate edge pairs with their BFS distances.
            minimal_pairs = []
            min_distance = float('inf')
            for i in range(len(candidate_edges)):
                for j in range(i + 1, len(candidate_edges)):
                    d = bfs_distance_graph(graph, candidate_edges[i], candidate_edges[j])
                    if d < min_distance:
                        min_distance = d
                        minimal_pairs = [[candidate_edges[i], d, candidate_edges[j]]]
                    elif d == min_distance:
                        pair_sorted = sorted([candidate_edges[i], candidate_edges[j]])
                        # Avoid duplicates (regardless of order)
                        if pair_sorted not in [sorted([p[0], p[2]]) for p in minimal_pairs]:
                            minimal_pairs.append([candidate_edges[i], d, candidate_edges[j]])

            # Now try each candidate minimal pair.
            success = False
            chosen_sub_old_indices = None
            chosen_rem_old_indices = None
            for pair_info in minimal_pairs:
                edge1, _d, edge2 = pair_info
                path = bfs_path_graph(graph, edge1, edge2)
                if not path:
                    continue
                candidate_sub_old_indices = sorted(path)
                candidate_rem_old_indices = [i for i in range(len(section)) if i not in candidate_sub_old_indices]

                # First, ensure connectivity on the remainder:
                if not induced_subgraph_is_connected(candidate_rem_old_indices, section):
                    # The remainder is split into multiple connected components.
                    components = induced_subgraph_components(candidate_rem_old_indices, section)
                    if len(components) == 2:
                        if len(components[0]) == len(components[1]):
                            # Try both merging options.
                            option1_sub = sorted(candidate_sub_old_indices + components[0])
                            option1_rem = [i for i in range(len(section)) if i not in option1_sub]
                            option2_sub = sorted(candidate_sub_old_indices + components[1])
                            option2_rem = [i for i in range(len(section)) if i not in option2_sub]
                            try:
                                dummy_final = final[:]  # Test option1
                                section_index = full_mapping.index(section)
                                dummy_full_mapping = full_mapping[:]  
                                dummy_full_mapping.pop(section_index)
                                dummy_full_mapping.insert(section_index, [section[i] for i in option1_sub])
                                dummy_full_mapping.insert(section_index + 1, [section[i] for i in option1_rem])
                                # --- New reindexing step:
                                reindexed_sub = reindex_section([section[i] for i in option1_sub], option1_sub)
                                reindexed_rem = reindex_section([section[i] for i in option1_rem], option1_rem)
                                # Now update the boundary mapping on the reindexed candidate subdivisions.
                                updated_sub_section, updated_remainder = update_divided_section_mapping(
                                    reindexed_sub,
                                    reindexed_rem
                                )
                                dummy_final = map_non_ring_section_1bead(updated_sub_section, dummy_final, martini_dict, dummy_full_mapping)
                                dummy_final = map_non_ring_section_1bead(updated_remainder, dummy_final, martini_dict, dummy_full_mapping)
                                candidate_sub_old_indices = option1_sub
                                candidate_rem_old_indices = option1_rem
                            except Exception:
                                try:
                                    dummy_final = final[:]  # Test option2
                                    section_index = full_mapping.index(section)
                                    dummy_full_mapping = full_mapping[:]
                                    dummy_full_mapping.pop(section_index)
                                    dummy_full_mapping.insert(section_index, [section[i] for i in option2_sub])
                                    dummy_full_mapping.insert(section_index + 1, [section[i] for i in option2_rem])
                                    # --- New reindexing step:
                                    reindexed_sub = reindex_section([section[i] for i in option2_sub], option2_sub)
                                    reindexed_rem = reindex_section([section[i] for i in option2_rem], option2_rem)
                                    updated_sub_section, updated_remainder = update_divided_section_mapping(
                                        reindexed_sub,
                                        reindexed_rem
                                    )
                                    dummy_final = map_non_ring_section_1bead(updated_sub_section, dummy_final, martini_dict, dummy_full_mapping)
                                    dummy_final = map_non_ring_section_1bead(updated_remainder, dummy_final, martini_dict, dummy_full_mapping)
                                    candidate_sub_old_indices = option2_sub
                                    candidate_rem_old_indices = option2_rem
                                except Exception:
                                    continue
                        else:
                            # When the two components have different sizes, choose the smaller one.
                            if len(components[0]) > len(components[1]):
                                merged_sub = sorted(candidate_sub_old_indices + components[1])
                            else:
                                merged_sub = sorted(candidate_sub_old_indices + components[0])
                            rem_after_merge = [i for i in range(len(section)) if i not in merged_sub]
                            if not induced_subgraph_is_connected(rem_after_merge, section):
                                continue
                            candidate_sub_old_indices = merged_sub
                            candidate_rem_old_indices = rem_after_merge
                    else:
                        continue

                # *** New test: update boundary mapping and try mapping the candidate parts.
                try:
                    dummy_final = final[:]  # Make a copy of current final mapping.
                    section_index = full_mapping.index(section)
                    dummy_full_mapping = full_mapping[:]  # Copy the full mapping.
                    dummy_full_mapping.pop(section_index)
                    dummy_full_mapping.insert(section_index, [section[i] for i in candidate_sub_old_indices])
                    dummy_full_mapping.insert(section_index + 1, [section[i] for i in candidate_rem_old_indices])
                    # --- New reindexing step:
                    reindexed_sub = reindex_section([section[i] for i in candidate_sub_old_indices], candidate_sub_old_indices)
                    reindexed_rem = reindex_section([section[i] for i in candidate_rem_old_indices], candidate_rem_old_indices)
                    updated_sub_section, updated_remainder = update_divided_section_mapping(
                        reindexed_sub,
                        reindexed_rem
                    )
                    dummy_final = map_non_ring_section_1bead(updated_sub_section, dummy_final, martini_dict, dummy_full_mapping)
                    dummy_final = map_non_ring_section_1bead(updated_remainder, dummy_final, martini_dict, dummy_full_mapping)
                except Exception:
                    # Mapping candidate failed; try next minimal pair.
                    continue

                # If we reached this point without error, candidate mapping succeeded.
                chosen_sub_old_indices = candidate_sub_old_indices
                chosen_rem_old_indices = candidate_rem_old_indices
                success = True
                break

            if success and chosen_sub_old_indices is not None and chosen_rem_old_indices is not None:
                # Instead of proceeding with DFS tracing, immediately return this subdivision.
                sub_section = [section[i] for i in chosen_sub_old_indices]
                remainder = [section[i] for i in chosen_rem_old_indices]
                return sub_section, remainder, chosen_sub_old_indices, chosen_rem_old_indices
            else:
                raise ValueError("Non‐ring section division leads to unassignable sections.")
        # --- End of new procedure ---

    # Use the first candidate as seed for DFS tracing if candidate_edges is nonempty.
    seed = candidate_edges[0]
    trace_nodes = dfs_trace_from_seed(section, seed)

    # NEW CONDITION 1: If the trace is long enough and the fourth node has ≥ 3 inner connections, limit trace to 3 nodes.
    if len(trace_nodes) >= 4 and len(section[trace_nodes[3]][4]) >= 3:
        trace_nodes = trace_nodes[:3]
    if len(trace_nodes) > 4:
        trace_nodes = trace_nodes[:4]
    # NEW CONDITION 2: Ensure remainder has at least 2 atoms.
    while len(section) - len(trace_nodes) < 2 and len(trace_nodes) > 2:
        trace_nodes = trace_nodes[:-1]
    if len(section) - len(trace_nodes) < 2:
        raise ValueError("Subdivision not possible: remainder too small (normal method).")
    sub_old_indices = sorted(trace_nodes)
    rem_old_indices = [i for i in range(len(section)) if i not in sub_old_indices]
    sub_section = [section[i] for i in sub_old_indices]
    remainder = [section[i] for i in rem_old_indices]
    return sub_section, remainder, sub_old_indices, rem_old_indices

def subdivide_non_ring_section_multi(section: List[List[Any]]) -> Tuple[List[List[Any]], List[List[Any]], List[int], List[int]]:
    """
    Special subdivision: when candidate_edges has 2 or more elements.
    Find the shortest path between any two edge nodes (without the <=4 constraint),
    and use that path as the sub_section.
    Returns:
      sub_section, remainder, sub_old_indices, rem_old_indices
    """
    edge_nodes = [i for i, atom in enumerate(section) if atom[5]]
    if len(edge_nodes) < 2:
        raise ValueError("Special subdivision not possible: fewer than 2 edge nodes.")

    graph = build_inner_graph_local(section)

    best_path: Optional[List[int]] = None
    best_length = float('inf')

    for i in range(len(edge_nodes)):
        for j in range(i + 1, len(edge_nodes)):
            path = bfs_path_graph(graph, edge_nodes[i], edge_nodes[j])
            if not path:
                continue
            length = len(path) - 1
            if length < best_length:
                best_length = length
                best_path = path

    if best_path is None or len(best_path) < 2:
        raise ValueError("Special subdivision not possible: no valid path found.")

    sub_old_indices = sorted(best_path)
    rem_old_indices = [i for i in range(len(section)) if i not in sub_old_indices]

    if len(rem_old_indices) < 2:
        # If remainder too small, try removing the last node from best_path.
        sub_old_indices = best_path[:-1]
        rem_old_indices = [i for i in range(len(section)) if i not in sub_old_indices]
        if len(rem_old_indices) < 2:
            raise ValueError("Special subdivision not possible: remainder too small.")

    sub_section = [section[i] for i in sub_old_indices]
    remainder = [section[i] for i in rem_old_indices]
    return sub_section, remainder, sub_old_indices, rem_old_indices

def subdivide_non_ring_section(section: List[List[Any]],
                               final: List[str],
                               full_mapping: List[List[List[Any]]],
                               martini_dict: Dict[str, List[Any]]) -> Tuple[List[List[Any]], List[List[Any]], List[int], List[int]]:
    """
    Chooses the subdivision method based on candidate edges.
    If candidate_edges has <=1 element, use normal; if >=2, try normal first and fall back to multi-edge.
    """
    candidate_edges = [i for i, atom in enumerate(section) if atom[5]]
    if len(candidate_edges) <= 1:
        return subdivide_non_ring_section_normal(section, final, martini_dict, full_mapping)
    else:
        try:
            return subdivide_non_ring_section_normal(section, final, martini_dict, full_mapping)
        except ValueError as e:
            print("Critical error:", e)
            return subdivide_non_ring_section_multi(section)


def reindex_section(section: List[List[Any]], old_indices: List[int]) -> List[List[Any]]:
    """
    Reindex a subdivided section so that inner connections (index 4) refer to new local indices.
    'old_indices' is a list of the original positions (local indices) for the atoms in this new section.
    """
    mapping = {old: new for new, old in enumerate(old_indices)}
    new_section = []
    for new_idx, atom in enumerate(section):
        new_atom = list(atom)
        new_inner = []
        for (nbr_old, bond) in new_atom[4]:
            # Convert neighbor's old index to new local index if it appears in the mapping.
            if nbr_old in mapping:
                new_inner.append((mapping[nbr_old], bond))
        new_atom[4] = new_inner
        new_section.append(new_atom)
    return new_section


def update_divided_section_mapping(sub_section: List[List[Any]], remainder: List[List[Any]]
                                  ) -> Tuple[List[List[Any]], List[List[Any]]]:
    """
    TASK 3.g (Normal Version, simplified):
    Given a subdivision into two parts (sub_section and remainder) along with their original positions
    (sub_old_indices and rem_old_indices, respectively),
    update the mapping so that the inner connections are reindexed to refer to the new local positions.
    Then update the isedge flag for both parts.
    
    Outer connections (index 3) are left unchanged.
    """
    fix_isedge(sub_section)
    fix_isedge(remainder)
    return sub_section, remainder

def map_non_ring_section_long(section: List[List[Any]],
                              final: List[str],
                              martini_dict: Dict[str, List[Any]],
                              full_mapping: List[List[List[Any]]],
                              section_index: int) -> Tuple[List[str], List[List[List[Any]]]]:
    """
    TASK 3.f/3.g (Revised Special):
    For a non-ring section (type == 0) that is too long to be mapped as a single bead,
    subdivide the section (using either normal or multi-edge splitting as needed) and update the mapping.
    
    Process:
      - Subdivide the section using subdivide_non_ring_section.
      - Reindex each new section.
      - Depending on which method was used, update the boundary.
        (For multi-edge, use update_divided_section_mapping_multi; otherwise, use normal.)
      - Replace the section in full_mapping with the two new sections.
      - Recursively map each new section.
    Returns updated final and full mapping.
    """
    old_pair_to_gid = {}
    for si, sec in enumerate(full_mapping):
        for li, a in enumerate(sec):
            old_pair_to_gid[(si, li)] = a[0]
    sub_sec, rem_sec, sub_old_indices, rem_old_indices = subdivide_non_ring_section(section, final, full_mapping, martini_dict)
    sub_sec = reindex_section(sub_sec, sub_old_indices)
    rem_sec = reindex_section(rem_sec, rem_old_indices)
    
    sub_sec, rem_sec = update_divided_section_mapping(sub_sec, rem_sec)
    
    full_mapping.pop(section_index)
    full_mapping.insert(section_index, sub_sec)
    full_mapping.insert(section_index + 1, rem_sec)
    
    rebuild_all_outer_connections(full_mapping, old_pair_to_gid)
    if not is_non_ring_section_1bead_mappable(sub_sec):
        final, full_mapping = map_non_ring_section_long(sub_sec, final, martini_dict, full_mapping, section_index)
    else:
        final = map_non_ring_section_1bead(sub_sec, final, martini_dict, full_mapping)
    if not is_non_ring_section_1bead_mappable(rem_sec):
        final, full_mapping = map_non_ring_section_long(rem_sec, final, martini_dict, full_mapping, section_index + 1)
    else:
        final = map_non_ring_section_1bead(rem_sec, final, martini_dict, full_mapping)
    return final, full_mapping

# =============================================================================
# TASK 3.h.a: Map Non-Benzene 3-Ring Section
# =============================================================================
def map_nonbenzene_3_ring_section(section: List[List[Any]],
                                  final: List[str],
                                  martini_dict: Dict[str, List[Any]]) -> List[str]:
    """
    TASK 3.h.a: Map a non-benzene 3-membered ring section.

    For now, trace the atoms in the section and if they are all carbon atoms
    with exactly two inner connections and every inner bond has order 1,
    assign the bead "SC3" + random string to each atom in the section.
    Otherwise, if exactly one atom is oxygen (and the other two are carbon,
    with the same inner‐bond criteria), assign "SN3a" + random string to all
    atoms in the section. Otherwise, throw an error indicating that the
    3-ring is not mappable.
    """

    # Count how many are oxygen vs. carbon
    o_count = sum(1 for atom in section if atom[1].upper() == "O")
    c_count = sum(1 for atom in section if atom[1].upper() == "C")

    # Verify inner‐bond structure: every atom must have exactly 2 inner connections,
    # and each bond must be order 1.
    for atom in section:
        if len(atom[4]) != 2:
            raise ValueError("3-ring not mappable: Expected exactly 2 inner connections per atom.")
        for (_, bond) in atom[4]:
            if bond != 1:
                raise ValueError("3-ring not mappable: Found an inner bond order not equal to 1.")

    # Case A: All three are carbon → look up the SC3 bead in martini_dict (v[0]==9, v[2]=="CCC")
    if o_count == 0 and c_count == 3:
        bead_key = find_bead(martini_dict, 9, "CCC")
        bead = bead_key + generate_random_string()
        for atom in section:
            final[atom[0]] = bead
        return final

    # Case B: Exactly one oxygen and two carbons → look up the SN4a bead (v[0]==9, v[2]=="COC")
    if o_count == 1 and c_count == 2:
        bead_key = find_bead(martini_dict, 9, "COC")
        bead = bead_key + generate_random_string()
        for atom in section:
            final[atom[0]] = bead
        return final


    # Otherwise, cannot map
    raise ValueError("3-ring not mappable: Element composition not supported.")

# =============================================================================
# TASK 3.h.b: Map Non-Benzene 4-Ring Section
# =============================================================================
def map_nonbenzene_4_ring_section(section: List[List[Any]],
                                  final: List[str],
                                  martini_dict: Dict[str, List[Any]]) -> List[str]:
    raise ValueError("Non-benzene ring section of size 4 is not supported.")


# =============================================================================
# Top-Level Mapping Function (Tree-Like Algorithm)
# =============================================================================
def map_martini_beads(mapping: List[List[List[Any]]],
                      final: List[str],
                      martini_dict: Dict[str, List[Any]]) -> List[str]:
    """
    Two-pass ring handling + then non-ring sections.
    """
    final2 = final
    try:
        skipped = []
    
        # --- PASS 1: Benzene rings (type 2) ---
        for idx, section in enumerate(mapping):
            if section and section[0][2] == 2 and any(final[a[0]] == "" for a in section):
                try:
                    new_final = map_benzene_ring_section(section, final.copy(), martini_dict, mapping)
                    final = new_final
                except Exception as e:
                    # record for second pass
                    print(f"mapping section {idx} failed:", e)
                    skipped.append(("benzene", idx))
                    
        # --- PASS 1: Non-benzene rings (type 1) ---
        for idx, section in enumerate(mapping):
            if section and section[0][2] == 1 and any(final[a[0]] == "" for a in section):
                size = len(section)
                try:
                    if size == 6:
                        new_final = map_nonbenzene_6_ring_section(section, final.copy(), martini_dict, mapping)
                    elif size == 5:
                        try:
                            new_final = map_nonbenzene_5_ring_section(section, final.copy(), martini_dict, mapping)
                        except Exception as e:
                            print(f"mapping section {idx} failed as a 5-ring section, attempt as if it is a 6-ring section:", e)
                            new_final = map_nonbenzene_6_ring_section(section, final.copy(), martini_dict, mapping)
                    elif size == 4:
                        new_final = map_nonbenzene_6_ring_section(section, final.copy(), martini_dict, mapping)
                    elif size == 3:
                        new_final = map_nonbenzene_3_ring_section(section, final.copy(), martini_dict)
                    else:
                        raise ValueError(f"Unexpected non-benzene ring size {size}")
                    final = new_final
                except Exception as e:
                    print(f"mapping section {idx} failed:", e)
                    skipped.append(("nonbenzene", idx))
    
        # --- PASS 2: Retry skipped ---
        for sect_type, idx in skipped:
            section = mapping[idx]
            if sect_type == "benzene":
                try:
                    final = map_benzene_ring_section(section, final, martini_dict, mapping)
                except Exception:
                    # PASS 2: if it fails, clear ring atoms in final and retry once
                    final = _retry_with_cleared_ring(
                        section, final,
                        map_benzene_ring_section,
                        martini_dict, mapping
                    )
            else:
                size = len(section)
                if size == 6:
                    try:
                        final = map_nonbenzene_6_ring_section(section, final, martini_dict, mapping)
                    except Exception:
                        # PASS 2: if it fails, clear ring atoms in final and retry once
                        final = _retry_with_cleared_ring(
                            section, final,
                            map_nonbenzene_6_ring_section,
                            martini_dict, mapping
                        )
                elif size == 5:
                    # PASS 2: try 5-ring first; if fails, clear and retry 6-ring; last resort try 5-ring (with cleared final)
                    try:
                        final = map_nonbenzene_5_ring_section(section, final, martini_dict, mapping)
                    except Exception:
                        # clear + retry 5-ring
                        try:
                            final = _retry_with_cleared_ring(
                                section, final,
                                map_nonbenzene_5_ring_section,
                                martini_dict, mapping
                            )
                        except Exception as e2:
                            # last resort: 6-ring mapping with final already cleared
                            try:
                                final = _retry_with_cleared_ring(
                                    section, final,
                                    map_nonbenzene_6_ring_section,
                                    martini_dict, mapping
                                )
                            except Exception as e3:
                                raise RuntimeError(
                                    f"PASS2 5-ring failed after (6-ring -> clear+retry 6-ring -> 5-ring): {e3}"
                                ) from e2
                elif size == 4:
                    final = map_nonbenzene_6_ring_section(section, final, martini_dict, mapping)
                elif size == 3:
                    final = map_nonbenzene_3_ring_section(section, final, martini_dict)
        # --- Then the non-ring sections as before ---
        i = 0
        while i < len(mapping):
            section = mapping[i]
            if section and section[0][2] == 0:
                # keep going until every atom in this section is filled
                while any(final[a[0]] == "" for a in section):
                    if is_non_ring_section_1bead_mappable(section):
                        final = map_non_ring_section_1bead(section, final, martini_dict, mapping)
                    else:
                        final, mapping = map_non_ring_section_long(section, final, martini_dict, mapping, i)
                        section = mapping[i]
            i += 1
    except Exception:
        final = final2
        skipped = []
    
        # --- PASS 1: Benzene rings (type 2) ---
        for idx, section in enumerate(mapping):
            if section and section[0][2] == 2 and any(final[a[0]] == "" for a in section):
                try:
                    new_final = map_benzene_ring_section(section, final.copy(), martini_dict, mapping)
                    final = new_final
                except Exception as e:
                    # record for second pass
                    print(f"mapping section {idx} failed:", e)
                    skipped.append(("benzene", idx))
    
        # --- PASS 1: Non-benzene rings (type 1) ---
        for idx, section in enumerate(mapping):
            if section and section[0][2] == 1 and any(final[a[0]] == "" for a in section):
                size = len(section)
                try:
                    if size == 6 or size == 5 or size == 4:
                        new_final = map_nonbenzene_6_ring_section(section, final.copy(), martini_dict, mapping)
                    elif size == 3:
                        new_final = map_nonbenzene_3_ring_section(section, final.copy(), martini_dict)
                    else:
                        raise ValueError(f"Unexpected non-benzene ring size {size}")
                    final = new_final
                except Exception as e:
                    print(f"mapping section {idx} failed:", e)
                    skipped.append(("nonbenzene", idx))
    
        # --- PASS 2: Retry skipped ---
        for sect_type, idx in skipped:
            section = mapping[idx]
            if sect_type == "benzene":
                final = map_benzene_ring_section(section, final, martini_dict, mapping)
            else:
                size = len(section)
                if size == 6 or size == 5 or size == 4:
                    final = map_nonbenzene_6_ring_section(section, final, martini_dict, mapping)
                elif size == 3:
                    final = map_nonbenzene_3_ring_section(section, final, martini_dict)
        # --- Then the non-ring sections as before ---
        i = 0
        while i < len(mapping):
            section = mapping[i]
            if section and section[0][2] == 0:
                # keep going until every atom in this section is filled
                while any(final[a[0]] == "" for a in section):
                    if is_non_ring_section_1bead_mappable(section):
                        final = map_non_ring_section_1bead(section, final, martini_dict, mapping)
                    else:
                        final, mapping = map_non_ring_section_long(section, final, martini_dict, mapping, i)
                        section = mapping[i]
            i += 1
    # sanity check
    if any(name == "" for name in final):
        raise ValueError("Final mapping incomplete; some atoms remain unmapped!")
    return final
