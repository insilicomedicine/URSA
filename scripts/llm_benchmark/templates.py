from __future__ import annotations

import random

QUERY_TEMPLATES = (
    (
        "Please propose a plausible synthetic route for the target molecule "
        "<smiles>{query}</smiles>. The route must converge to commercially "
        "available starting materials (purchasable building blocks)."
    ),
    (
        "Please propose a stepwise retrosynthetic plan for the target molecule "
        "<smiles>{query}</smiles> in which the final precursors are limited to "
        "commercially supplied reagents and building blocks."
    ),
    (
        "Please propose a convergent retrosynthetic pathway for the target "
        "molecule <smiles>{query}</smiles>. It is mandatory that the route "
        "ends on commercially available starting materials, not bespoke "
        "intermediates."
    ),
    (
        "Design the entire synthetic route to reach the target molecule "
        "<smiles>{query}</smiles> from catalog chemicals. Explicitly ensure "
        "the retrosynthesis bottoms out in commercially accessible building "
        "blocks."
    ),
    (
        "For the target molecule <smiles>{query}</smiles>, outline a "
        "retrosynthetic analysis and propose a synthetic route. The route "
        "must terminate in reagents and building blocks that are commercially "
        "available from standard suppliers. Please propose the full route."
    ),
    (
        "Please propose a sequence of retrosynthetic disconnections that "
        "constitutes a valid synthetic route for the target molecule "
        "<smiles>{query}</smiles>. The route must lead to commercially "
        "available starting materials."
    ),
    (
        "Please perform retrosynthetic analysis and propose a synthetic route "
        "for the molecule <smiles>{query}</smiles>. The route is acceptable "
        "only if every terminal precursor can be sourced as a commercially "
        "available building block."
    ),
    (
        "Please propose a practical multi-step retrosynthetic plan for the "
        "target molecule <smiles>{query}</smiles>, in which all reactions "
        "ultimately rest on commercially catalogued building blocks as the "
        "inputs."
    ),
    (
        "Please propose the retrosynthetic pathway from the target molecule "
        "<smiles>{query}</smiles> down to commercially available building "
        "blocks."
    ),
    (
        "Please propose a retrosynthetic scheme for the target molecule "
        "<smiles>{query}</smiles>. The route must lead to commercially "
        "purchasable reagents and building blocks."
    ),
)

FORMAT_INSTRUCTION = (
    "Return only <synthesis_step> blocks, one per reaction. Use the following "
    "route only as an output-format example for the target "
    "CC1CCN(Cc2ccccc2COC2(C(F)F)CCCC2)CC1:\n"
    "<synthesis_step><product><name>query</name><smiles>"
    "CC1CCN(Cc2ccccc2COC2(C(F)F)CCCC2)CC1</smiles></product>"
    "<reactant><smiles>FC(F)C1(O[Na])CCCC1</smiles>"
    "<name>intermediate 1</name></reactant>"
    "<reactant><smiles>CC1CCN(Cc2ccccc2CCl)CC1</smiles>"
    "<name>intermediate 2</name></reactant></synthesis_step>"
    "<synthesis_step><product><name>intermediate 1</name>"
    "<smiles>FC(F)C1(O[Na])CCCC1</smiles></product>"
    "<reactant><smiles>OC1(C(F)F)CCCC1</smiles>"
    "<name>building block 1</name></reactant></synthesis_step>"
    "<synthesis_step><product><name>intermediate 2</name>"
    "<smiles>CC1CCN(Cc2ccccc2CCl)CC1</smiles></product>"
    "<reactant><smiles>CC1CCN(Cc2ccccc2CO)CC1</smiles>"
    "<name>building block 2</name></reactant></synthesis_step>\n"
    "Use query for the target, intermediate N for a reactant that must still "
    "be prepared, building block N for a purchasable one. Each intermediate "
    "N is the product of some later block; every branch bottoms out in a "
    "building block. Add no other text."
)


def build_prompt(target_smiles: str) -> str:
    """Build one randomized benchmark prompt for a target."""
    query = random.choice(QUERY_TEMPLATES).replace("{query}", target_smiles)
    return f"{query}\n\n{FORMAT_INSTRUCTION}"
