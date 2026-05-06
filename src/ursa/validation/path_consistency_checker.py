from ..basic import RetrosyntheticNode
from ..basic import RetrosyntheticPath
from ..configs import ValidationConfig


class PathConsistencyChecker:
    """Checks structural consistency of a retrosynthetic tree.

    A path is consistent when:

    * The tree has exactly one root (the target molecule).
    * All molecules are reachable from the root via product-to-reactant
      relationships.
    * For every internal node, the ``reaction_smiles`` product matches
      the node's ``smiles``, and the reactants correspond to
      the node's children.
    """

    def check(self, path: RetrosyntheticPath) -> bool:
        """Return ``True`` if ``path`` is structurally consistent.

        Traverses the full tree and verifies connectivity and
        product-reactant correspondence at every internal node.

        :param path: The retrosynthetic path to validate.
        :type path: RetrosyntheticPath

        :return: ``True`` if the tree passes all consistency checks.
        :rtype: bool
        """
        if path.depth < ValidationConfig.min_path_depth.value:
            return False
        return self._check_node(path.root)

    def _check_node(self, node: RetrosyntheticNode) -> bool:
        """Recursively check a single node and its subtree.

        For an internal node, verifies that its ``reaction_smiles``
        product is consistent with ``node.smiles`` and that
        all children are reachable.

        :param node: The node to validate.
        :type node: RetrosyntheticNode

        :return: ``True`` if the node and its entire subtree are valid.
        :rtype: bool
        """
        if not node.is_valid:
            return False
        return all(self._check_node(child) for child in node.children)
