from typing import List, Tuple, Optional

# TASK 1a: Parse the SMILES string into tokens.
def parse_smiles(smiles: str) -> List[str]:
    """
    Tokenize the SMILES string into meaningful tokens:
    - Multi-character elements: Cl, Br, NH
    - Ring closures: '1'-'9' and '%dd' for 10+
    - Element symbols and lowercase aromatic symbols
    Ignores parentheses, bond symbols, stereochemistry markers, charges, and dot separators.

    Examples:
      c1c(CC)cccc1 -> ['c','1','c','C','C','c','c','c','c','1']
      C12CC      -> ['C','12','C','C']
    """
    tokens: List[str] = []
    i = 0
    n = len(smiles)
    while i < n:
        ch = smiles[i]
        # Skip ignored characters
        if smiles.startswith('Hg', i):
            tokens.append('Hg')
            i += 2
            continue
        if ch in '()=#$@H[]+-/\\':
            i += 1
            continue
        # Multi-letter element tokens
        if smiles.startswith('Cl', i):
            tokens.append('Cl')
            i += 2
            continue
        if smiles.startswith('Br', i):
            tokens.append('Br')
            i += 2
            continue
        if smiles.startswith('NH', i):
            tokens.append('NH')
            i += 2
            continue
        if smiles.startswith('Pb', i):
            tokens.append('Pb')
            i += 2
            continue
        if smiles.startswith('Si', i):
            tokens.append('Si')
            i += 2
            continue
        if smiles.startswith('Si', i):
            tokens.append('Si')
            i += 2
            continue
        if smiles.startswith('Ge', i):
            tokens.append('Ge')
            i += 2
            continue
        if smiles.startswith('Sn', i):
            tokens.append('Sn')
            i += 2
            continue
        # Ring closure: multi-digit
        if ch == '%':
            if i + 2 < n and smiles[i+1].isdigit() and smiles[i+2].isdigit():
                tokens.append(smiles[i+1:i+3])
                i += 3
                continue
            else:
                raise ValueError("% is not supposed to be here in SMILES")
        # Ring closure: single-digit
        if ch.isdigit():
            tokens.append(ch)
            i += 1
            continue
        # Otherwise, element or aromatic symbol
        tokens.append(ch)
        i += 1
    return tokens


# TASK 1b: Given a token index in the smiles_tokens array, return the corresponding index in the filtered (non-digit) sequence.
def get_property_index(token_index: int, smiles_tokens: List[str]) -> int:
    """
    Count non-digit tokens preceding the given token_index.
    For example, with tokens ['c','1','c','C',...]:
      get_property_index(2, ...) returns 1 because:
        token0 'c' (non-digit, count=1), token1 '1' (skip), token2 'c' becomes the 2nd non-digit => index 1.
    """
    count = 0
    for i in range(token_index):
        if not smiles_tokens[i].isdigit():
            count += 1
    return count

# TASK 1e: Given a property index (i.e. index in the non-digit sequence) return the corresponding token.
def get_atom_element(property_index: int, smiles_tokens: List[str]) -> str:
    """
    Skip digits and return the token that is the (property_index)-th non-digit token.
    For example, with tokens ['c','1','c','C','C',...], property index 2 returns 'C'.
    """
    count = 0
    for token in smiles_tokens:
        if not token.isdigit():
            if count == property_index:
                return token
            count += 1
    return ""

# TASK 1c: Given a starting atom index and a connectivity (adjacency) matrix, use DFS to find a cycle
# that goes through exactly 4 or 5 distinct nodes (i.e. a ring).
def find_ring(start_index: int, 
              connectivity: List[List[int]], 
              solution: List[List[List]], 
              allowed_lengths: List[int] = [3, 4, 5, 6]) -> List[int]:
    """
    Performs a DFS to find a cycle (ring) starting at start_index that goes through a number of
    distinct nodes matching one of the allowed_lengths. It will reject cycles for which all node
    indices are already present in one section of the solution.
    
    Returns a list of node indices forming the cycle (without repeating the start at the end).
    If no cycle is found, returns an empty list.
    """
    def is_duplicate(candidate: List[int], solution: List[List[List]]) -> bool:
        """Check if the candidate cycle (as a set of indices) is already present in any section."""
        if len(solution) == 0: return False
        section_indices = []
        candidate_set = set(candidate)
        for section in solution:
            # Each atom is represented as [global_index, ...]
            for atom in section:
                section_indices.append(atom[0])
            # If all candidate indices are found in this section, it's a duplicate.
            if candidate_set.issubset(section_indices):
                return True
        return False

    def dfs(current: int, path: List[int]) -> Optional[List[int]]:
        # Check for cycle: if we're back at start, and the path length is allowed (and not trivial)
        if current == start_index and len(path) in allowed_lengths and len(path) > 1:
            candidate = path[:]  # candidate cycle (without repeating the start at the end)
            if not is_duplicate(candidate, solution):
                return candidate
            # Else, continue searching for an alternative cycle.
        # Prevent going too deep:
        if len(path) > max(allowed_lengths):
            return None
        for neighbor, bond in enumerate(connectivity[current]):
            if bond == 0:
                continue
            # Check if closing the cycle is possible:
            if neighbor == start_index and len(path) in allowed_lengths:
                candidate = path[:]  # candidate cycle found
                if not is_duplicate(candidate, solution):
                    return candidate
                # Otherwise, do not return this duplicate; continue searching.
            # Avoid revisiting nodes in the current path.
            if neighbor in path:
                continue
            new_path = path + [neighbor]
            result = dfs(neighbor, new_path)
            if result is not None:
                return result
        return None

    cycle = dfs(start_index, [start_index])
    return cycle if cycle is not None else []

# TASK 1d: Append a ring section into the solution.
def append_ring_section(ring_indices: List[int],
                        connectivity: List[List[int]],
                        properties: List[List],
                        smiles_tokens: List[str],
                        is_benzene: bool,
                        solution: List[List[List]]) -> Tuple[List[List[List]], List[List]]:
    """
    Creates a new section (a ring) where each atom is represented as:
      [global_index, element, ring_status, outer_connection, inner_connection, isedge]
    ring_status: 2 if benzene, otherwise 1.
    inner_connection: For each neighbor (within ring_indices), record (local index, bond value).
    Also marks the atom in properties as "mapped" by setting its element to "X".
    """
    section = []
    ring_status = 2 if is_benzene else 1
    for idx in ring_indices:
        element = get_atom_element(idx, smiles_tokens)
        outer_connection = []  # to be updated later in TASK 1h
        inner_connection = []
        # For each neighbor in connectivity row that is within the ring
        for j, bond in enumerate(connectivity[idx]):
            if bond != 0 and j in ring_indices:
                local_j = ring_indices.index(j)
                inner_connection.append((local_j, bond))
        isedge = False  # default set to False; will update later
        num_H = properties[idx][3]
        atom_rep = [idx, element, ring_status, outer_connection, inner_connection, isedge, num_H]
        section.append(atom_rep)
        # Mark as mapped (set element to "X")
        properties[idx][0] = "X"
    solution.append(section)
    return solution, properties

# TASK 1f: Use BFS to find connected nodes (by connectivity) that are not yet mapped.
def find_connected_non_mapped(start_index: int,
                              properties: List[List],
                              connectivity: List[List[int]]) -> List[int]:
    """
    Starting at start_index, perform BFS on the connectivity graph and return all node indices that
    are connected and are not marked as "X" in the properties matrix.
    """
    visited = set()
    queue = [start_index]
    while queue:
        current = queue.pop(0)
        if current in visited:
            continue
        if properties[current][0] == "X":
            continue
        visited.add(current)
        for neighbor, bond in enumerate(connectivity[current]):
            if bond != 0 and properties[neighbor][0] != "X" and neighbor not in visited:
                queue.append(neighbor)
    return list(visited)

# TASK 1g: Append a non-ring section (subgraph) into the solution.
def append_non_ring_section(section_indices: List[int],
                            connectivity: List[List[int]],
                            properties: List[List],
                            smiles_tokens: List[str],
                            solution: List[List[List]]) -> Tuple[List[List[List]], List[List]]:
    """
    Similar to append_ring_section but for a non-ring section.
    Here ring_status is 0.
    """
    section = []
    ring_status = 0
    for idx in section_indices:
        element = get_atom_element(idx, smiles_tokens)
        outer_connection = []  # to be updated later
        inner_connection = []
        for j, bond in enumerate(connectivity[idx]):
            if bond != 0 and j in section_indices:
                local_j = section_indices.index(j)
                inner_connection.append((local_j, bond))
        isedge = False  # default set to False; will update later
        num_H = properties[idx][3]
        atom_rep = [idx, element, ring_status, outer_connection, inner_connection, isedge, num_H]
        section.append(atom_rep)
        properties[idx][0] = "X"
    solution.append(section)
    return solution, properties

# TASK 1h: Update outer (foreign) connections between sections.
def update_outer_connections(connectivity: List[List[int]],
                             solution: List[List[List]]) -> List[List[List]]:
    """
    For each atom in each section, inspect its connectivity row.
    For each neighbor (global index) that is not in the same section, find the section and local index
    where that neighbor is mapped. Append a tuple (section_index, local_index, bond) to the atom's
    outer_connection list.
    """
    for sec_idx, section in enumerate(solution):
        for atom_idx, atom in enumerate(section):
            global_index = atom[0]
            # Check each possible connection for this global index.
            for neighbor, bond in enumerate(connectivity[global_index]):
                if bond == 0:
                    continue
                # Skip if the neighbor is in the same section.
                if any(neighbor == other_atom[0] for other_atom in section):
                    continue
                # Search among all sections.
                found = False
                for other_sec_idx, other_section in enumerate(solution):
                    if other_sec_idx == sec_idx:
                        continue
                    for local_idx, other_atom in enumerate(other_section):
                        if other_atom[0] == neighbor:
                            atom[3].append((other_sec_idx, local_idx, bond))
                            found = True
                            break
                    if found:
                        break
                if not found:
                    raise ValueError(f"Foreign index {neighbor} not found in any section")
    return solution

# New TASK: Update isedge flags for each atom.
def update_isedge_flags(solution: List[List[List]]) -> List[List[List]]:
    """
    Loop through the molecule, and for each atom (in each section), set isedge to True
    if the atom's inner_connection length is 0 or 1, otherwise False.
    """
    for section in solution:
        for atom in section:
            # atom structure: [global_index, element, ring_status, outer_connection, inner_connection, isedge]
            inner_connections = atom[4]
            atom[5] = len(inner_connections) <= 1
    return solution

# Main function that maps the whole molecule.
def map_molecule(smiles_tokens: List[str],
                 connectivity: List[List[int]],
                 properties: List[List]) -> List[List[List]]:
    """
    Overall mapping function.
    
    1. Parse the SMILES string (TASK 1a).
    2. For each ring-digit (e.g. '1', '2', ...) found in the tokens, find the corresponding ring.
       (Uses TASK 1b, 1c, 1d – note that the algorithm finds the atom before the digit and then
       finds a ring cycle starting from its property index. Also, determine benzene status by checking
       if the token is lowercase 'c'.)
    3. Then, for every unmapped atom (not marked as "X" in properties) use BFS (TASK 1f) and map
       non-ring sections (TASK 1g).
    4. Finally, update outer connections (TASK 1h) and update isedge flags.
    """
    solution = []
    # Process rings by scanning for ring digits (starting at '1' upwards)
    ring_num = 1
    while str(ring_num) in smiles_tokens:
        # Get the first occurrence of this ring digit.
        indices = [i for i, token in enumerate(smiles_tokens) if token == str(ring_num)]
        if not indices:
            break
        pos = indices[0]
        # The atom token preceding the ring digit is assumed to be part of the ring.
        if pos > 0 and not smiles_tokens[pos-1].isdigit():
            # Get the corresponding property index from the token position.
            prop_idx = get_property_index(pos-1, smiles_tokens)
            # Find a cycle (ring) starting from this property index.
            ring_cycle = find_ring(prop_idx, connectivity, solution)
            if ring_cycle:
                if len(ring_cycle) == 6:
                    # Condition 2: Count atoms in the ring that have at least one inner connectivity with bond order >= 1.5.
                    count_high_connectivity = sum(
                        1 for atom in ring_cycle if any(
                            connectivity[atom][other] >= 1.5 for other in ring_cycle if other != atom
                        )
                    )
                    
                    is_benzene = (count_high_connectivity >= 4)
                else:
                    is_benzene = False
                solution, properties = append_ring_section(ring_cycle, connectivity, properties,
                                                             smiles_tokens, is_benzene, solution)
                ring_cycle = find_ring(prop_idx, connectivity, solution)
                if ring_cycle:
                    # Determine if benzene by checking if any bond in the ring cycle has a bond order of 1.5.
                    if len(ring_cycle) == 6:
                        # Condition 2: Count atoms in the ring that have at least one inner connectivity with bond order >= 1.5.
                        count_high_connectivity = sum(
                            1 for atom in ring_cycle if any(
                                connectivity[atom][other] >= 1.5 for other in ring_cycle if other != atom
                            )
                        )
                        
                        is_benzene = (count_high_connectivity >= 4)
                    else: 
                        is_benzene = False
                    solution, properties = append_ring_section(ring_cycle, connectivity, properties,
                                                                 smiles_tokens, is_benzene, solution)
                    if ring_cycle:
                        # Determine if benzene by checking if any bond in the ring cycle has a bond order of 1.5.
                        if len(ring_cycle) == 6:
                            # Condition 2: Count atoms in the ring that have at least one inner connectivity with bond order >= 1.5.
                            count_high_connectivity = sum(
                                1 for atom in ring_cycle if any(
                                    connectivity[atom][other] >= 1.5 for other in ring_cycle if other != atom
                                )
                            )
                            
                            is_benzene = (count_high_connectivity >= 4)
                        else: 
                            is_benzene = False
                        solution, properties = append_ring_section(ring_cycle, connectivity, properties,
                                                                     smiles_tokens, is_benzene, solution)
                        if ring_cycle:
                            # Determine if benzene by checking if any bond in the ring cycle has a bond order of 1.5.
                            if len(ring_cycle) == 6:
                                # Condition 2: Count atoms in the ring that have at least one inner connectivity with bond order >= 1.5.
                                count_high_connectivity = sum(
                                    1 for atom in ring_cycle if any(
                                        connectivity[atom][other] >= 1.5 for other in ring_cycle if other != atom
                                    )
                                )
                                
                                is_benzene = (count_high_connectivity >= 4)
                            else: 
                                is_benzene = False
                            solution, properties = append_ring_section(ring_cycle, connectivity, properties,
                                                                         smiles_tokens, is_benzene, solution)
                            if ring_cycle:
                                # Determine if benzene by checking if any bond in the ring cycle has a bond order of 1.5.
                                if len(ring_cycle) == 6:
                                    # Condition 2: Count atoms in the ring that have at least one inner connectivity with bond order >= 1.5.
                                    count_high_connectivity = sum(
                                        1 for atom in ring_cycle if any(
                                            connectivity[atom][other] >= 1.5 for other in ring_cycle if other != atom
                                        )
                                    )
                                    
                                    is_benzene = (count_high_connectivity >= 4)
                                else: 
                                    is_benzene = False
                                solution, properties = append_ring_section(ring_cycle, connectivity, properties,
                                                                             smiles_tokens, is_benzene, solution)
        ring_num += 1

    # Process non-ring sections: for every atom not yet mapped (i.e. not marked "X")
    n = len(properties)
    for i in range(n):
        if properties[i][0] != "X":
            section_indices = find_connected_non_mapped(i, properties, connectivity)
            solution, properties = append_non_ring_section(section_indices, connectivity, properties,
                                                             smiles_tokens, solution)
    # Update outer connections between sections.
    solution = update_outer_connections(connectivity, solution)
    # Update the isedge flag for each atom.
    solution = update_isedge_flags(solution)
    return solution
