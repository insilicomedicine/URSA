# based on https://chemapps.stolaf.edu/iupac/cdx/sdk/IntroCDXML.htm specification
# (severely outdated, but still working)
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator
from typing import Optional

from lxml import etree  # nosec: B410

from .cdxml_elements import CDXMLBoundingBox
from .cdxml_elements import CDXMLLine
from .cdxml_elements import CDXMLText
from .cdxml_molecule import CDXMLMolecule
from .molecule import Molecule

# only used for trusted data

BRACKET_WIDTH = 5
BRACKET_MARGIN = 10
ARROW_LENGTH = 200


@dataclass
class CDXMLBracket:
    """
    Describes a bracket connecting a molecule (product) to a group of
    molecules (reactants) with a label (reaction name).

    :param lines: The lines of the bracket;
    :param label: The label of the bracket;
    """

    lines: list[CDXMLLine]
    label: Optional[CDXMLText] = None

    @classmethod
    def group_molecules(
        cls,
        mol_from: CDXMLMolecule,
        mols_to: list[CDXMLMolecule],
        label_text: Optional[str] = None,
    ) -> CDXMLBracket:
        """
        Creates a CDXMLBracket from a molecule (product) and a group of
        molecules (reactants) with a label (reaction name).

        :param mol_from: The starting molecule (product);
        :param mols_to: The ending group of molecules (reactants);
        :param label_text: The label (reaction name);
        :return: The CDXMLBracket object;
        """

        mol_bboxes = [mol.get_bounding_box() for mol in mols_to]

        alpha = (
            min(mol_bbox.min_x for mol_bbox in mol_bboxes)
            - BRACKET_MARGIN
            - BRACKET_WIDTH
        )
        beta = min(mol_bbox.center_y for mol_bbox in mol_bboxes)
        gamma = max(mol_bbox.center_y for mol_bbox in mol_bboxes)

        lines = []
        label: Optional[CDXMLText] = None

        lines.append(
            CDXMLLine(
                beg_pos_x=mol_from.get_bounding_box().max_x + BRACKET_MARGIN,
                beg_pos_y=mol_from.get_bounding_box().center_y,
                end_pos_x=alpha,
                end_pos_y=mol_from.get_bounding_box().center_y,
            )
        )

        if label_text:
            label = CDXMLText(
                text=label_text,
                pos_x=mol_from.get_bounding_box().max_x + BRACKET_MARGIN,
                pos_y=mol_from.get_bounding_box().center_y - 20,
            )

        for mol_bbox in mol_bboxes:
            lines.append(
                CDXMLLine(
                    beg_pos_x=alpha + BRACKET_MARGIN,
                    beg_pos_y=mol_bbox.center_y,
                    end_pos_x=alpha,
                    end_pos_y=mol_bbox.center_y,
                )
            )

        lines.append(
            CDXMLLine(
                beg_pos_x=alpha,
                beg_pos_y=beta,
                end_pos_x=alpha,
                end_pos_y=gamma,
            )
        )

        return cls(
            lines=lines,
            label=label,
        )

    def as_cdxml_element(self) -> etree.Element:
        """
        Returns the CDXML element for the bracket.
        :return: The CDXML element for the bracket;
        """

        element = etree.Element("group")
        for line in self.lines:
            element.append(line.as_cdxml_element())

        if self.label:
            element.append(self.label.as_cdxml_element())

        return element


class CDXMLTree:
    """
    Describes a tree branch. Recursively describes a synthesis tree.

    :param root_molecule: The root molecule (product);
    :param reaction_name: The reaction name;
    :param children: The children of the root molecule (reactants);
    """

    root_molecule: CDXMLMolecule
    reaction_name: Optional[str]
    children: list[CDXMLTree]

    def __init__(
        self,
        root_molecule: CDXMLMolecule,
        children: list[CDXMLTree],
        reaction_name: Optional[str] = None,
    ):
        """
        Initializes a CDXMLTree.
        :param root_molecule: The root molecule (product);
        :param reaction_name: The reaction name;
        :param children: The children of the root molecule (reactants);
        """

        self.root_molecule = root_molecule
        self.reaction_name = reaction_name
        self.children = children
        self._arrange()

    def _arrange(self):
        """
        Arranges by stacking the children branches right of the root molecule.
        The root molecule is translated to the middle-left of the stack.
        """

        if not self.children:
            return

        stack_width = max(
            child.root_molecule.get_bounding_box().width for child in self.children
        )
        stack_height = max(
            self.root_molecule.get_bounding_box().height,
            sum(child.get_bounding_box().height for child in self.children),
        )

        self.root_molecule.translate(
            x=0,
            y=-self.root_molecule.get_bounding_box().height / 2 + stack_height / 2,
        )

        if len(self.children) == 1:
            self.children[0].translate(
                x=self.root_molecule.get_bounding_box().max_x + ARROW_LENGTH,
                y=-self.children[0].get_bounding_box().height / 2 + stack_height / 2,
            )
        else:
            curr_stack_height = 0
            for child in self.children:
                child.translate(
                    x=(
                        self.root_molecule.get_bounding_box().max_x
                        + (stack_width - child.root_molecule.get_bounding_box().width)
                        / 2
                        + ARROW_LENGTH
                    ),
                    y=curr_stack_height,
                )
                curr_stack_height += child.get_bounding_box().height

    def get_bounding_box(self) -> CDXMLBoundingBox:
        """
        Returns the bounding box of the tree.
        :return: The bounding box of the tree;
        """

        if not self.children:
            return self.root_molecule.get_bounding_box()

        bboxes = [self.root_molecule.get_bounding_box()] + [
            child.get_bounding_box() for child in self.children
        ]

        min_x = min(bbox.min_x for bbox in bboxes)
        min_y = min(bbox.min_y for bbox in bboxes)
        max_x = max(bbox.max_x for bbox in bboxes)
        max_y = max(bbox.max_y for bbox in bboxes)

        return CDXMLBoundingBox(min_x, min_y, max_x, max_y)

    def translate(self, x: float, y: float) -> CDXMLTree:
        """
        Translates the tree by the given x and y coordinates (in-place).
        :param x: The x coordinate to translate the tree by;
        :param y: The y coordinate to translate the tree by;
        :return: The translated tree (same object);
        """

        self.root_molecule.translate(x, y)
        for child in self.children:
            child.translate(x, y)
        return self

    def as_cdxml_element(self) -> etree.Element:
        """
        Returns the CDXML element for the tree.
        :return: The CDXML element for the tree;
        """

        element = etree.Element("group")
        element.append(self.root_molecule.as_cdxml_element())
        for child in self.children:
            element.append(child.as_cdxml_element())

        if not self.children:
            return element

        bracket = CDXMLBracket.group_molecules(
            mol_from=self.root_molecule,
            mols_to=[child.root_molecule for child in self.children],
            label_text=self.reaction_name,
        )

        element.append(bracket.as_cdxml_element())

        return element

    @classmethod
    def _build_tree(cls, route: dict) -> dict:
        """
        Builds a tree of molecules from a route dict as stored in the
        retrosynthesis database.

        :param route: The route;
        :return: A hierarchical dictionary of molecules with the following
            structure:
            ```{
                "smi": product_smiles,
                "label": product_label,
                "outline": whether the molecule should be outlined,
                "reaction_name": the reaction name if applicable,
                "children": a list of child molecules,
            }```
        """

        tree = {}

        building_block_indices = {
            bb["molecule"]["index"] for bb in route["building_blocks"]
        }

        for step in route["steps"]:
            product = step["product"]
            reactants = [r for r in step["reactants"]]

            for mol in [product] + reactants:
                if mol["molecule"]["index"] not in tree:
                    tree[mol["molecule"]["index"]] = {
                        "smi": mol["molecule"]["smiles"],
                        "label": str(mol["molecule"]["index"]),
                        "outline": mol["molecule"]["index"] in building_block_indices,
                        "reaction_name": None,
                        "children": [],
                    }

            product_idx = product["molecule"]["index"]
            for reactant in reactants:
                # Prefer a custom step label if present; otherwise no label
                tree[product_idx]["reaction_name"] = step.get("label")
                tree[product_idx]["children"].append(
                    tree[reactant["molecule"]["index"]]
                )

        return tree[0]

    @classmethod
    def _build_cdxml_tree(cls, index_tree: dict, obj_id_it: Iterator[int]) -> CDXMLTree:
        """
        Builds a CDXMLTree from a hierarchical dictionary of molecules.
        :param index_tree: The hierarchical dictionary of molecules as
            returned by `_build_tree`;
        :param obj_id_it: An iterator for generating object IDs;
        :return: The CDXMLTree;
        """

        return cls(
            root_molecule=CDXMLMolecule.from_molecule(
                Molecule(index_tree["smi"]),
                obj_id_it,
                label=index_tree["label"],
                outline=index_tree["outline"],
            ),
            children=[
                cls._build_cdxml_tree(child, obj_id_it)
                for child in index_tree["children"]
            ],
            reaction_name=index_tree["reaction_name"],
        )

    @classmethod
    def from_route(
        cls,
        route: dict,
        obj_id_it: Iterator[int],
    ) -> CDXMLTree:
        """
        Builds a CDXMLTree from a route dict representation.
        :param route: The route dict;
        :param obj_id_it: An iterator for generating object IDs;
        :return: The CDXMLTree;
        """

        return cls._build_cdxml_tree(cls._build_tree(route), obj_id_it)
