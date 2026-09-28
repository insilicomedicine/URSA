# based on https://chemapps.stolaf.edu/iupac/cdx/sdk/IntroCDXML.htm specification
# (severely outdated, but still working)
from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Iterator
from typing import Optional

import numpy as np
import rdkit.Chem.AllChem as Chem
from lxml import etree  # nosec: B410
from rdkit.Chem.Draw.rdMolDraw2D import PrepareMolForDrawing
from rdkit.Chem.rdchem import GetPeriodicTable
from rdkit.Chem.rdMolTransforms import GetBondLength

from .cdxml_elements import CDXMLBoundingBox
from .cdxml_elements import CDXMLRectangle
from .cdxml_elements import CDXMLText
from .molecule import Molecule

# only used for trusted data

MOLECULE_PADDING = 25


# constants for bond order
BOND_ORDER_DICT: dict[Chem.BondType, float] = {
    Chem.BondType.DOUBLE: 2,
    Chem.BondType.TRIPLE: 3,
    Chem.BondType.AROMATIC: 1.5,
}

# constants for bond stereo
BOND_STEREO_DICT: dict[Chem.BondStereo, str] = {
    Chem.BondStereo.STEREOCIS: "Z",
    Chem.BondStereo.STEREOZ: "Z",
    Chem.BondStereo.STEREOTRANS: "E",
    Chem.BondStereo.STEREOE: "E",
}

# constants for bond direction
BOND_DIR_DICT: dict[Chem.BondDir, str] = {
    Chem.BondDir.BEGINDASH: "WedgedHashBegin",
    Chem.BondDir.BEGINWEDGE: "WedgeBegin",
}


@dataclass
class CDXMLAtom:
    """
    Describes an atom in a molecule.

    :param obj_id: The object ID of the atom;
    :param pos_x: The x coordinate of the atom;
    :param pos_y: The y coordinate of the atom;
    :param atomic_num: The atomic number of the atom;
    :param num_hs: The number of hydrogen atoms attached to the atom;
    :param atom_stereo: The stereochemistry of the atom;
    :param charge: The charge of the atom;
    :param isotope: The isotope of the atom;
    """

    obj_id: int
    pos_x: float
    pos_y: float
    atomic_num: int
    num_hs: int
    atom_stereo: str
    charge: int
    isotope: int

    def _prop_dict(self) -> dict[str, str]:
        """
        Returns a dictionary of properties for the atom element.
        :return: A dictionary of properties for the atom element;
        """

        props = {
            "id": str(self.obj_id),
            "p": f"{self.pos_x:.2f} {self.pos_y:.2f}",
            "Element": str(self.atomic_num),
        }

        # not adding hydrogens for carbons
        if self.atomic_num == 6 and self.num_hs != 0:
            props["NumHydrogens"] = str(self.num_hs)
        if self.atom_stereo != "N":
            props["AS"] = self.atom_stereo
        if self.charge != 0:
            props["Charge"] = str(self.charge)
        if self.isotope != 0:
            props["Isotope"] = str(self.isotope)

        return props

    def _label(self) -> str:
        """
        Returns a string label for the atom describing its chemical properties
        (i.e. charge, number of hydrogens, isotope, etc.).
        :return: A string label for the atom;
        """

        # label is required for heteroatoms and charged carbons
        if self.atomic_num == 6 and self.charge == 0:
            return ""

        # label deuterium correctly
        symbol = (
            "D"
            if self.atomic_num == 1 and self.isotope == 2
            else GetPeriodicTable().GetElementSymbol(self.atomic_num)
        )

        # label hydrogens
        hydrogen_str = ""
        if self.num_hs > 0:
            hydrogen_str = "H" if self.num_hs == 1 else f"H{self.num_hs}"

        # label charge
        charge_str = ""
        if self.charge != 0:
            charge_sign = "+" if self.charge > 0 else "-"
            charge_value = abs(self.charge) if abs(self.charge) != 1 else ""
            charge_str = f"{charge_sign}{charge_value}"

        return f"{symbol}{hydrogen_str}{charge_str}"

    def as_cdxml_element(self) -> etree.Element:
        """
        Returns the CDXML element for the atom.
        :return: The CDXML element for the atom;
        """

        atom_element = etree.Element("n", **self._prop_dict())
        element_label = self._label()

        if not element_label:
            return atom_element

        t_element = etree.SubElement(atom_element, "t")
        s_element = etree.SubElement(
            t_element,
            "s",
            font="1",
            size="10.0",
            face="96",
            color="0",
        )
        s_element.text = element_label

        return atom_element


@dataclass
class CDXMLBond:
    """
    Describes a bond between two atoms.

    :param begin_atom_obj_id: The object ID of the beginning atom;
    :param end_atom_obj_id: The object ID of the ending atom;
    :param order: The order of the bond;
    :param bond_stereo: The stereochemistry of the bond;
    :param bond_dir: The direction of the bond;
    """

    begin_atom_obj_id: int
    end_atom_obj_id: int
    order: Chem.BondType
    bond_stereo: Chem.BondStereo
    bond_dir: Chem.BondDir

    def _prop_dict(self) -> dict[str, str]:
        """
        Returns a dictionary of properties for the bond element.
        :return: A dictionary of properties for the bond element;
        """

        props = {
            "B": str(self.begin_atom_obj_id),
            "E": str(self.end_atom_obj_id),
        }

        if self.order != Chem.BondType.SINGLE:
            props["Order"] = str(BOND_ORDER_DICT[self.order])

        if self.bond_stereo != Chem.BondStereo.STEREONONE:
            props["BS"] = BOND_STEREO_DICT[self.bond_stereo]

        if self.bond_dir in BOND_DIR_DICT:
            props["Display"] = BOND_DIR_DICT[self.bond_dir]

        return props

    def as_cdxml_element(self) -> etree.Element:
        """
        Returns the CDXML element for the bond.
        :return: The CDXML element for the bond;
        """

        return etree.Element("b", **self._prop_dict())


@dataclass
class CDXMLMolecule:
    """
    Describes a molecule.

    :param atoms: The atoms in the molecule;
    :param bonds: The bonds in the molecule;
    :param label: The label of the molecule;
    :param outline: Whether the molecule is should be outlined;
    """

    atoms: list[CDXMLAtom]
    bonds: list[CDXMLBond]
    label: Optional[str]
    outline: bool

    @classmethod
    def from_molecule(
        cls,
        mol: Molecule,
        obj_id_it: Iterator[int],
        label: Optional[str] = None,
        outline: bool = False,
    ) -> CDXMLMolecule:
        """
        Creates a CDXMLMolecule from a Molecule object.

        :param mol: The Molecule object;
        :param obj_id_it: An iterator for generating object IDs;
        :param label: The label of the molecule;
        :param outline: Whether the molecule is should be outlined;
        :return: The CDXMLMolecule object;
        """

        rdmol = copy.deepcopy(mol.canonical_rdmol)
        rdmol = PrepareMolForDrawing(
            rdmol, kekulize=True, addChiralHs=True, wedgeBonds=True
        )
        conformer = rdmol.GetConformer()
        coords = np.array(conformer.GetPositions())[:, :2]

        if rdmol.GetNumBonds() > 0:
            avg_bond_length = np.mean(
                [
                    GetBondLength(
                        conformer,
                        bond.GetBeginAtomIdx(),
                        bond.GetEndAtomIdx(),
                    )
                    for bond in rdmol.GetBonds()
                ]
            )
            coords = coords / avg_bond_length * 14.40

        coords = coords * [[1, -1]]
        coords = coords - coords.min(axis=0)

        atoms: list[CDXMLAtom] = []
        for rdatom in rdmol.GetAtoms():
            atom = CDXMLAtom(
                obj_id=next(obj_id_it),
                pos_x=coords[rdatom.GetIdx(), 0] + MOLECULE_PADDING,
                pos_y=coords[rdatom.GetIdx(), 1] + MOLECULE_PADDING,
                atomic_num=rdatom.GetAtomicNum(),
                num_hs=rdatom.GetTotalNumHs(),
                atom_stereo=rdatom.GetPropsAsDict().get("_CIPCode", "N"),
                charge=rdatom.GetFormalCharge(),
                isotope=rdatom.GetIsotope(),
            )
            atoms.append(atom)

        bonds: list[CDXMLBond] = []
        for rdbond in rdmol.GetBonds():
            bond = CDXMLBond(
                begin_atom_obj_id=atoms[rdbond.GetBeginAtomIdx()].obj_id,
                end_atom_obj_id=atoms[rdbond.GetEndAtomIdx()].obj_id,
                order=rdbond.GetBondType(),
                bond_stereo=rdbond.GetStereo(),
                bond_dir=rdbond.GetBondDir(),
            )
            bonds.append(bond)

        return cls(
            atoms=atoms,
            bonds=bonds,
            label=label,
            outline=outline,
        )

    def translate(self, x: float, y: float) -> CDXMLMolecule:
        """
        Translates the molecule by the given x and y coordinates (in-place).

        :param x: The x coordinate to translate the molecule by;
        :param y: The y coordinate to translate the molecule by;
        :return: The translated molecule (same object);
        """

        for atom in self.atoms:
            atom.pos_x += x
            atom.pos_y += y

        return self

    def get_bounding_box(self) -> CDXMLBoundingBox:
        """
        Returns the bounding box of the molecule.
        :return: The bounding box of the molecule;
        """

        label_padding = 20 if self.label else 0

        min_x = min(atom.pos_x for atom in self.atoms) - MOLECULE_PADDING
        min_y = min(atom.pos_y for atom in self.atoms) - MOLECULE_PADDING
        max_x = max(atom.pos_x for atom in self.atoms) + MOLECULE_PADDING
        max_y = (
            max(atom.pos_y for atom in self.atoms) + MOLECULE_PADDING + label_padding
        )

        return CDXMLBoundingBox(min_x, min_y, max_x, max_y)

    def as_cdxml_element(self) -> etree.Element:
        """
        Returns the CDXML element for the molecule.
        :return: The CDXML element for the molecule;
        """

        element = etree.Element("fragment")
        for atom in self.atoms:
            element.append(atom.as_cdxml_element())
        for bond in self.bonds:
            element.append(bond.as_cdxml_element())

        if self.label:
            element.append(
                CDXMLText(
                    text=self.label,
                    pos_x=self.get_bounding_box().min_x + MOLECULE_PADDING,
                    pos_y=self.get_bounding_box().max_y - 15,
                ).as_cdxml_element()
            )

        if self.outline:
            element.append(
                CDXMLRectangle.from_bbox(self.get_bounding_box()).as_cdxml_element()
            )

        return element
