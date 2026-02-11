from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Optional, Tuple, Dict, List, Set

from rdkit import Chem


# -------------------- Parse name + SMILES from .txt --------------------
COMPOUND_RE = re.compile(r"^\s*Compound\s*:\s*(.+?)\s*$", re.IGNORECASE)
SMILES_RE = re.compile(r"^\s*SMILES\s*:\s*(.+?)\s*$", re.IGNORECASE)


def parse_compound_and_smiles(txt_path: Path) -> Tuple[Optional[str], Optional[str]]:
    name: Optional[str] = None
    smiles: Optional[str] = None
    with txt_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if name is None:
                m = COMPOUND_RE.match(line)
                if m:
                    name = m.group(1).strip()
                    continue
            if smiles is None:
                m = SMILES_RE.match(line)
                if m:
                    smiles = m.group(1).strip()
                    continue
            if name is not None and smiles is not None:
                break
    return name, smiles


# -------------------- Parse bead table from .txt --------------------
ALL_BEADS_HDR = re.compile(r"^\s*All beads\s*$", re.IGNORECASE)
BEAD_CONNS_HDR = re.compile(r"^\s*Bead connections\s*$", re.IGNORECASE)
DASHES = re.compile(r"^-{3,}\s*$")


def _split_int_list(s: str) -> List[int]:
    s = s.strip()
    if not s:
        return []
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def parse_beads(txt_path: Path) -> List[dict]:
    """
    Reads only the 'All beads' table and returns:
      beads: [{bead_index:int, bead_type:str, atoms_idx:List[int]}]
    """
    beads: List[dict] = []
    in_beads = False
    header_seen = False

    with txt_path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")

            if ALL_BEADS_HDR.match(line):
                in_beads = True
                header_seen = False
                continue
            if BEAD_CONNS_HDR.match(line):
                in_beads = False
                continue
            if DASHES.match(line.strip()):
                continue
            if not in_beads:
                continue

            # header
            if not header_seen:
                if "bead_type" in line and "atoms_idx" in line:
                    header_seen = True
                continue

            if not line.strip():
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            bead_index = int(parts[0])
            bead_type = parts[1]
            atoms_idx = _split_int_list(parts[2])

            beads.append({"bead_index": bead_index, "bead_type": bead_type, "atoms_idx": atoms_idx})

    return beads


def parse_bead_connections(txt_path: Path) -> List[tuple]:
    """
    Reads only the 'Bead connections' table and returns a list of undirected edges:
      edges: [(i:int, j:int, bond:float)]
    Notes:
      - The 'bond' column is the bead-bead bond order indicator (e.g., 1.0 single, 1.5 aromatic-like).
      - This function does NOT infer multiplicity; each row is a single bead-bead edge.
    """
    edges: List[tuple] = []
    in_conns = False
    header_seen = False

    with txt_path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")

            if BEAD_CONNS_HDR.match(line):
                in_conns = True
                header_seen = False
                continue

            if ALL_BEADS_HDR.match(line):
                # if another section starts, stop reading connections
                if in_conns:
                    break

            if DASHES.match(line.strip()):
                continue
            if not in_conns:
                continue

            if not header_seen:
                if "to_idx" in line and "bond" in line:
                    header_seen = True
                continue

            if not line.strip():
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            i = int(parts[0])
            j = int(parts[1])
            bond = float(parts[2])
            a, b = (i, j) if i < j else (j, i)
            edges.append((a, b, bond))

    # de-duplicate (keep the first occurrence, but preserve bond if repeated identically)
    uniq = {}
    for a, b, bond in edges:
        key = (a, b)
        if key not in uniq:
            uniq[key] = bond
        else:
            # if conflicting, keep the larger bond order (e.g., 1 vs 2)
            if abs(uniq[key] - bond) > 1e-6:
                uniq[key] = max(uniq[key], bond)

    return [(a, b, uniq[(a, b)]) for (a, b) in sorted(uniq.keys())]


# -------------------- CGsmiles builders --------------------

def bond_order_symbol(n: int) -> str:
    """
    Map number of AA connections between beads (0..4) to CG bond-order symbols:
      0 -> '.' (virtual edge)
      1 -> ''  (implicit single, we omit symbol)
      2 -> '='
      3 -> '#'
      4 -> '$'
    """
    if n <= 0:
        return "."
    if n == 1:
        return ""
    if n == 2:
        return "="
    if n == 3:
        return "#"
    return "$"

def bond_value_symbol(bond: float) -> str:
    """
    Map bead-level bond order values from the 'Bead connections' table to CGsmiles bond symbols.

    We keep the CG graph SIMPLE (no parallel edges). Aromatic-like 1.5 is treated as a normal
    connection at the CG layer (no explicit ':'), because ring-closure numbering already
    captures the topology and many downstream uses don't require aromatic symbols here.

      ~1.0 -> '' (implicit single)
      ~1.5 -> '' (implicit, aromatic-like)
      ~2.0 -> '='
      ~3.0 -> '#'
      else -> '' (fallback)
    """
    if abs(bond - 2.0) < 1e-6:
        return "="
    if abs(bond - 3.0) < 1e-6:
        return "#"
    return ""



def ring_token(r: int) -> str:
    return str(r) if r <= 9 else f"%{r}"


def derive_bead_adjacency_from_smiles(
    mol: Chem.Mol, beads: List[dict]
) -> Tuple[Dict[Tuple[int, int], int], Dict[int, int]]:
    """
    Returns:
      edge_counts: (bead_i, bead_j) -> number of AA bonds crossing between those beads
      atom_to_bead: atom_idx -> bead_index
    """
    atom_to_bead: Dict[int, int] = {}
    for b in beads:
        bi = b["bead_index"]
        for aidx in b["atoms_idx"]:
            atom_to_bead[aidx] = bi

    edge_counts: Dict[Tuple[int, int], int] = {}
    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a not in atom_to_bead or b not in atom_to_bead:
            continue
        ba = atom_to_bead[a]
        bb = atom_to_bead[b]
        if ba == bb:
            continue
        i, j = (ba, bb) if ba < bb else (bb, ba)
        edge_counts[(i, j)] = edge_counts.get((i, j), 0) + 1

    return edge_counts, atom_to_bead


def _suffix_letter(i: int) -> str:
    # 0->A, 1->B ...
    return chr(ord("A") + i)


def assign_bead_labels(beads: List[dict], frag_by_bead_i: Dict[int, str]) -> Dict[int, str]:
    """
    User rule:
      - Default label is bead_type (e.g., TC5)
      - If the same bead_type maps to different fragment SMILES, label as TC5A, TC5B, ...
        (suffix based on unique fragment string)
    """
    by_type: Dict[str, Dict[str, List[int]]] = {}
    for b in beads:
        i = b["bead_index"]
        t = b["bead_type"]
        frag = frag_by_bead_i.get(i, "")
        by_type.setdefault(t, {}).setdefault(frag, []).append(i)

    out: Dict[int, str] = {}
    for t, frag_groups in by_type.items():
        frags = sorted(frag_groups.keys())
        if len(frags) == 1:
            for i in frag_groups[frags[0]]:
                out[i] = t
        else:
            for idx, frag in enumerate(frags):
                suf = _suffix_letter(idx)
                for i in frag_groups[frag]:
                    out[i] = f"{t}{suf}"
    return out


def build_graph_block(bead_indices: List[int], edges: List[tuple], bead_label_by_i: Dict[int, str]) -> str:
    """
    Build first-resolution CGsmiles graph.

    Nodes: [#<label>] where <label> is bead_type or bead_typeA/B...
    Edges: taken directly from the 'Bead connections' table (simple graph; no parallel edges).
    """
    if not bead_indices:
        return "{}"

    bead_set = set(bead_indices)

    adj: Dict[int, List[int]] = {i: [] for i in bead_indices}
    sym: Dict[Tuple[int, int], str] = {}

    for i, j, bond in edges:
        if i in bead_set and j in bead_set:
            adj[i].append(j)
            adj[j].append(i)
            a, b = (i, j) if i < j else (j, i)
            sym[(a, b)] = bond_value_symbol(bond)

    for k in adj:
        adj[k].sort()

    start = min(bead_indices)
    visited: Set[int] = set()
    parent: Dict[int, int] = {}
    token_pos: Dict[int, int] = {}

    parts: List[str] = []

    ring_ids: Dict[Tuple[int, int], int] = {}
    next_ring_id = 1

    def edge_symbol(u: int, v: int) -> str:
        a, b = (u, v) if u < v else (v, u)
        return sym.get((a, b), "")

    def emit_node(u: int):
        nonlocal next_ring_id
        visited.add(u)

        token = f"[#{bead_label_by_i[u]}]"
        parts.append(token)
        token_pos[u] = len(parts) - 1

        # ring closures: any visited neighbor that isn't the parent
        # IMPORTANT: emit each ring-closure digit exactly once per edge (otherwise you can
        # accidentally duplicate closures and appear to duplicate nodes in downstream parsers).
        for v in adj.get(u, []):
            if v in visited and parent.get(u) != v:
                a, b = (u, v) if u < v else (v, u)
                if (a, b) not in ring_ids:
                    ring_ids[(a, b)] = next_ring_id
                    next_ring_id += 1
                    rid = ring_ids[(a, b)]
                    rt = ring_token(rid)
                    parts[token_pos[u]] += rt
                    # v is already emitted earlier; we can safely edit its token in-place once.
                    if v in token_pos:
                        parts[token_pos[v]] += rt

        first_child = True
        for v in adj.get(u, []):
            if v in visited:
                continue
            parent[v] = u
            bsym = edge_symbol(u, v)
            if first_child:
                if bsym:
                    parts.append(bsym)
                emit_node(v)
                first_child = False
            else:
                parts.append("(")
                if bsym:
                    parts.append(bsym)
                emit_node(v)
                parts.append(")")


    emit_node(start)

    # disconnected components, if any
    for comp in sorted(i for i in bead_indices if i not in visited):
        parts.append(".")
        parent[comp] = -1
        emit_node(comp)

    return "{" + "".join(parts) + "}"


def _replace_mapped_dummy(smiles: str, atom_map: int, connector: str) -> str:
    """
    Replace ONE specific RDKit mapped dummy atom token [*:atom_map] into a CGsmiles bonding connector.

    Handles both:
      - branch dummy: (...[*:k]...) -> ...<bond><connector>...
      - inline dummy: X[*:k] or X=[*:k] -> X<bond><connector>

    We remove the dummy token entirely; connector becomes an annotation on the neighbor atom.
    """
    # 1) branch form: "( [*:k] )" optionally with bond symbol
    pat_branch = re.compile(rf"\(\s*([=\-#:\$\\/\.]?)\s*\[\*:{atom_map}\]\s*\)")
    new_smiles, n = pat_branch.subn(lambda m: (m.group(1) + connector), smiles, count=1)
    if n:
        return new_smiles

    # 2) inline form: "=[*:k]" or "[*:k]" after atom
    pat_inline = re.compile(rf"([=\-#:\$\\/\.]?)\[\*:{atom_map}\]")
    new_smiles, n = pat_inline.subn(lambda m: (m.group(1) + connector), smiles, count=1)
    if n:
        return new_smiles

    # if not found, return unchanged (caller can decide)
    return smiles


def build_bead_fragment_with_edge_connectors(
    mol: Chem.Mol,
    bead_atoms: List[int],
    atom_to_bead: Dict[int, int],
    bead_i: int,
) -> str:
    """
    Build an atomistic fragment for this bead and add CGsmiles bonding connectors [$]
    for each cut bond that leaves the bead.

    IMPORTANT:
      - Output must contain NO '*' dummy atoms.
      - Connectors are undirected '$' and uses the undirected bonding operator '[$]' for each cut bond.
    """
    subset = set(bead_atoms)
    if not subset:
        return ""

    # Collect cut-bonds in a deterministic order
    cut_info: List[Tuple[int, int]] = []  # (inside_atom_idx, neighbor_bead_idx)
    for ai in sorted(subset):
        a = mol.GetAtomWithIdx(ai)
        for nb in a.GetNeighbors():
            aj = nb.GetIdx()
            if aj not in subset and aj in atom_to_bead:
                cut_info.append((ai, atom_to_bead[aj]))

    # Build fragment: subset atoms + one mapped dummy per cut bond
    rw = Chem.RWMol()
    old_to_new: Dict[int, int] = {}
    for aidx in sorted(subset):
        a = mol.GetAtomWithIdx(aidx)
        na = Chem.Atom(a.GetAtomicNum())
        na.SetFormalCharge(a.GetFormalCharge())
        na.SetIsAromatic(False)
        old_to_new[aidx] = rw.AddAtom(na)

    # internal bonds
    for bond in mol.GetBonds():
        a0 = bond.GetBeginAtomIdx()
        a1 = bond.GetEndAtomIdx()
        if a0 in subset and a1 in subset:
            rw.AddBond(old_to_new[a0], old_to_new[a1], bond.GetBondType())
            # Do NOT propagate aromatic flags into fragments. If a bead contains only a portion
            # of an aromatic ring, those flags can become invalid after cutting and RDKit will
            # error with "non-ring atom marked aromatic".

    # dummy attachments with atom-map numbers 1..N
    for k, (ai, _nb_bead) in enumerate(cut_info, start=1):
        d = Chem.Atom(0)  # [*]
        d.SetAtomMapNum(k)  # [*:k]
        d_idx = rw.AddAtom(d)
        rw.AddBond(old_to_new[ai], d_idx, Chem.BondType.SINGLE)

    frag = rw.GetMol()
    Chem.SanitizeMol(frag, catchErrors=True)

    # Force Kekulé (non-aromatic) output so CGsmiles doesn't reject fragments with connectors
    try:
        Chem.Kekulize(frag, clearAromaticFlags=True)
    except Exception:
        pass

    smi = Chem.MolToSmiles(frag, isomericSmiles=True, kekuleSmiles=True)

    # Replace each mapped dummy [*:k] with its bonding connector
    for k, (_ai, nb_bead) in enumerate(cut_info, start=1):
        conn = "[$]"
        new_smi = _replace_mapped_dummy(smi, k, conn)
        smi = new_smi

    # Safety: no stray '*' should remain
    if "*" in smi:
        raise ValueError(f"Unconverted dummy atom(s) in fragment for bead {bead_i}: {smi}")

    return smi


def build_definition_block_unique_labels(
    beads: List[dict],
    bead_label_by_i: Dict[int, str],
    frag_by_bead_i: Dict[int, str],
) -> str:
    """
    Emit a single fragment definition per UNIQUE bead label:
      {#TC5=...,#TN6aA=...}
    """
    label_to_frag: Dict[str, str] = {}
    for b in beads:
        i = b["bead_index"]
        lab = bead_label_by_i[i]
        frag = frag_by_bead_i.get(i, "")
        if lab not in label_to_frag:
            label_to_frag[lab] = frag
        else:
            # Should match if the label assignment is correct
            if label_to_frag[lab] != frag:
                raise ValueError(f"Label collision: {lab} has multiple different fragments.")

    items = [f"#{lab}={label_to_frag[lab]}" for lab in sorted(label_to_frag.keys())]
    return "{" + ",".join(items) + "}"


def build_cgsmiles_from_txt(smiles: str, beads: List[dict], edges: List[tuple]) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit failed to parse SMILES: {smiles}")
    
    # IMPORTANT: kekulize whole molecule before fragmenting
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        # If RDKit can't kekulize, we'll still try, but some fragments may fail
        pass

    # atom -> bead mapping is determined by the bead table
    atom_to_bead: Dict[int, int] = {}
    for b in beads:
        bi = b["bead_index"]
        for aidx in b["atoms_idx"]:
            atom_to_bead[aidx] = bi


    # Build per-bead fragments (with NO '*' and with edge-labeled connectors)
    frag_by_bead_i: Dict[int, str] = {}
    bead_indices = sorted(b["bead_index"] for b in beads)

    for b in beads:
        i = b["bead_index"]
        frag_by_bead_i[i] = build_bead_fragment_with_edge_connectors(
            mol=mol,
            bead_atoms=b["atoms_idx"],
            atom_to_bead=atom_to_bead,
            bead_i=i,
        )

    # Assign labels per user rule (no instance labels unless same type has different fragments)
    bead_label_by_i = assign_bead_labels(beads, frag_by_bead_i)

    # Coarse graph + unique definitions
    graph = build_graph_block(bead_indices, edges, bead_label_by_i)
    defs = build_definition_block_unique_labels(beads, bead_label_by_i, frag_by_bead_i)
    return f"{graph}.{defs}"


# -------------------- Walk folders + write ONLY full_list.csv --------------------
def build_full_list_csv(root_dir: str, out_csv: str = "full_list.csv") -> Path:
    root = Path(root_dir)
    if not root.exists():
        raise FileNotFoundError(f"Root directory not found: {root}")

    rows: List[Dict[str, str]] = []

    xtb_dirs = sorted(p for p in root.rglob("xtb") if p.is_dir())
    if not xtb_dirs:
        out_path = root / out_csv
        with out_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["Name", "smiles", "cgsmiles"])
            w.writeheader()
        return out_path

    for xtb_dir in xtb_dirs:
        for sub in sorted(p for p in xtb_dir.iterdir() if p.is_dir()):
            mol_folder = sub.name
            txt_path = sub / f"{mol_folder}.txt"
            if not txt_path.exists():
                continue

            name, smiles = parse_compound_and_smiles(txt_path)
            name = name or ""
            smiles = smiles or ""

            beads = parse_beads(txt_path)
            edges = parse_bead_connections(txt_path)

            cgsmiles = ""
            if smiles and beads:
                try:
                    cgsmiles = build_cgsmiles_from_txt(smiles, beads, edges)
                except Exception as e:
                    print(f"[WARN] Failed CGsmiles for {txt_path}: {e}")
                    cgsmiles = ""

            # Fallback: single bead, no mapping info
            if not cgsmiles and smiles:
                cgsmiles = "{[#B0]}.{#B0=" + smiles + "}"

            rows.append({"Name": name, "smiles": smiles, "cgsmiles": cgsmiles})

    out_path = root / out_csv
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["Name", "smiles", "cgsmiles"])
        w.writeheader()
        w.writerows(rows)

    print(f"[OK] Wrote {out_path} with {len(rows)} rows.")
    return out_path


if __name__ == "__main__":
    ROOT_DIR = "Working_Molecules"
    build_full_list_csv(ROOT_DIR, out_csv="full_list.csv")
