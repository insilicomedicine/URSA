from dataclasses import dataclass


@dataclass(frozen=True)
class BuildingBlock:
    """A starting material checked against the building-block catalog.

    Represents a leaf node (:attr:`~ursa.RetrosyntheticNode.is_starting_material`
    ``== True``) together with the result of looking it up in the catalog
    provided to :class:`~ursa.BuildingBlockChecker`. A starting material
    that is ``found_in_catalog`` is considered a *building block* in the
    chemical sense (commercially available).

    :param smiles: SMILES of the starting-material molecule.
    :type smiles: str
    :param found_in_catalog: ``True`` if ``smiles`` was found in the
        building-block catalog, ``False`` otherwise.
    :type found_in_catalog: bool
    """

    smiles: str
    found_in_catalog: bool
