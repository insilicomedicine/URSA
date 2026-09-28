from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable
from typing import Sequence


def _try_canonicalize_smiles(smiles: str) -> str:
    """
    Returns a canonicalized SMILES string when RDKit is available.
    Falls back to the input string if RDKit is not installed or parsing fails.
    """

    try:
        # Delayed import to keep this utility lightweight if RDKit isn't used
        from rdkit import Chem  # type: ignore
    except Exception:
        return smiles

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        return Chem.MolToSmiles(mol, canonical=True, ignoreAtomMapNumbers=True)
    except Exception:
        return smiles


def _split_reaction_smiles(reaction_smiles: str) -> tuple[list[str], list[str]]:
    """
    Splits a reaction SMILES of the form 'R1.R2...>>P1.P2...' into
    (reactants, products) lists of SMILES. Empty parts are ignored.
    """

    if ">>" not in reaction_smiles:
        raise ValueError(f"Reaction SMILES must contain '>>': {reaction_smiles!r}")

    left, right = reaction_smiles.split(">>", 1)
    reactants = [s for s in left.split(".") if s]
    products = [s for s in right.split(".") if s]
    return reactants, products


@dataclass
class _MoleculeRef:
    index: int
    smiles: str


def _assign_indices(
    reactions: Sequence[tuple[list[str], list[str]]],
) -> tuple[dict[str, int], str]:
    """
    Assigns integer indices to unique molecules across reactions and returns
    (smiles_to_index, final_product_smiles). The final product is inferred as
    a product that never appears as a reactant across the route.

    The final product is always assigned index 0; all other molecules get
    indices starting from 1.
    """

    reactant_set: set[str] = set()
    product_set: set[str] = set()

    for reactants, products in reactions:
        for smi in reactants:
            reactant_set.add(_try_canonicalize_smiles(smi))
        for smi in products:
            product_set.add(_try_canonicalize_smiles(smi))

    # Final product(s) are products not used as reactants
    candidate_targets = list(product_set - reactant_set)
    if not candidate_targets:
        raise ValueError(
            "Could not infer a final product: every product is also a reactant."
        )

    # Heuristic: if multiple candidates, choose the one that appears as a product
    # in the first step that produces a candidate (often the most downstream one
    # is referenced later). This keeps behavior deterministic.
    candidate_rank: dict[str, int] = {}
    for step_idx, (_, products) in enumerate(reactions):
        canon_products = [_try_canonicalize_smiles(s) for s in products]
        for smi in canon_products:
            if smi in candidate_targets and smi not in candidate_rank:
                candidate_rank[smi] = step_idx

    candidate_targets.sort(key=lambda s: candidate_rank.get(s, 1_000_000))
    final_product = (
        candidate_targets[-1] if len(candidate_targets) > 1 else candidate_targets[0]
    )

    # Assign indices: final product gets 0
    smiles_to_index: dict[str, int] = {final_product: 0}

    next_index = 1
    for reactants, products in reactions:
        for smi in list(reactants) + list(products):
            canon = _try_canonicalize_smiles(smi)
            if canon not in smiles_to_index:
                smiles_to_index[canon] = next_index
                next_index += 1

    # Ensure final product stays 0 (in case ordering assigned it otherwise)
    if smiles_to_index[final_product] != 0:
        # Find who has 0 and swap
        for smi, idx in list(smiles_to_index.items()):
            if idx == 0:
                smiles_to_index[smi] = smiles_to_index[final_product]
        smiles_to_index[final_product] = 0

    return smiles_to_index, final_product


def convert_reaction_smiles_to_route(
    reaction_smiles_list: Iterable[str],
    labels: Sequence[str] | None = None,
) -> dict:
    """
    Converts a list of reaction SMILES strings into a route dictionary that is
    compatible with build_cdxml_document(). The resulting dict has at least the
    keys required by the CDXML renderer:

    - building_blocks: list of {"molecule": {"index": int, "smiles": str}}
    - steps: list of steps; each step has
        - product: {"molecule": {"index": int, "smiles": str}}
        - reactants: list of {"molecule": {"index": int, "smiles": str}}

    The converter infers building blocks as molecules that never appear as a
    product in any step. The final product is inferred as a product that never
    reappears as a reactant and is assigned index 0.
    """

    # Normalize and parse all reactions once
    parsed: list[tuple[list[str], list[str]]] = []
    rxn_list = list(reaction_smiles_list)
    for rxn in rxn_list:
        reactants, products = _split_reaction_smiles(rxn)
        if not products:
            raise ValueError(f"Reaction has no products: {rxn!r}")
        parsed.append((reactants, products))

    if labels is not None and len(labels) != len(rxn_list):
        raise ValueError(
            f"labels length ({len(labels)}) must match reactions length "
            f"({len(rxn_list)})"
        )

    smiles_to_index, final_product = _assign_indices(parsed)

    # Determine building blocks (reactants that never appear as products)
    all_reactants: set[str] = set()
    all_products: set[str] = set()
    for reactants, products in parsed:
        all_reactants.update(_try_canonicalize_smiles(s) for s in reactants)
        all_products.update(_try_canonicalize_smiles(s) for s in products)

    building_block_smiles = sorted(
        all_reactants - all_products, key=lambda s: smiles_to_index[s]
    )

    def _molecule_entry(smiles: str) -> dict:
        canon = _try_canonicalize_smiles(smiles)
        return {
            "molecule": {
                "index": smiles_to_index[canon],
                "smiles": canon,
            }
        }

    steps: list[dict] = []
    for step_idx, (reactants, products) in enumerate(parsed, start=1):
        # If multiple products are present, keep the first one as the step's product
        # to satisfy the tree builder expectations. Others (if any) are ignored.
        product_canon = _try_canonicalize_smiles(products[0])
        step_entry: dict = {
            "product": _molecule_entry(product_canon),
            "reactants": [_molecule_entry(r) for r in reactants],
            "reaction_meta": {"name": f"step-{step_idx}"},
        }
        if labels is not None:
            step_entry["label"] = labels[step_idx - 1]
        steps.append(step_entry)

    route: dict = {
        "building_blocks": [_molecule_entry(s) for s in building_block_smiles],
        "steps": steps,
    }

    # Guarantee that a step exists that yields the final product with index 0
    # (this should already be true if the inputs are consistent)
    has_root = any(step["product"]["molecule"]["index"] == 0 for step in route["steps"])
    if not has_root:
        raise ValueError(
            "Inconsistent reactions: none of the steps produces "
            "the inferred final product."
        )

    return route


def reaction_smiles_to_cdxml(
    reaction_smiles_list: Iterable[str],
    labels: Sequence[str] | None = None,
) -> str:
    """
    Convenience helper to convert reaction SMILES to a CDXML document string.
    """

    from .cdxml_doc import build_cdxml_document

    route = convert_reaction_smiles_to_route(reaction_smiles_list, labels=labels)
    return build_cdxml_document([route])
