from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

from .node import RetrosyntheticNode


@dataclass(frozen=True)
class RetrosyntheticPath:
    """A retrosynthetic synthesis tree rooted at the target molecule.

    The tree is represented by its root node (the target molecule).
    Each internal node stores the forward reaction that transforms its
    children into the node's molecule. Leaf nodes are building blocks.

    :param path_id: Unique identifier for this retrosynthetic path.
    :type path_id: str
    :param root: Root node of the retrosynthetic tree (target molecule).
    :type root: RetrosyntheticNode
    """

    path_id: str
    root: RetrosyntheticNode

    @staticmethod
    def _dfs(node: RetrosyntheticNode) -> list[RetrosyntheticNode]:
        """Return all nodes rooted at ``node`` in DFS pre-order.

        :param node: Subtree root to traverse.
        :type node: RetrosyntheticNode
        :return: Flat list of all nodes, root first.
        :rtype: list[RetrosyntheticNode]
        """
        result = [node]
        for child in node.children:
            result.extend(RetrosyntheticPath._dfs(child))
        return result

    @staticmethod
    def _node_depth(node: RetrosyntheticNode) -> int:
        """Return the maximum depth of the subtree rooted at ``node``.

        :param node: Subtree root.
        :type node: RetrosyntheticNode
        :return: Number of edges on the longest root-to-leaf path.
        :rtype: int
        """
        if not node.children:
            return 0
        return 1 + max(RetrosyntheticPath._node_depth(c) for c in node.children)

    def get_all_steps(self) -> tuple[RetrosyntheticNode, ...]:
        """Return all internal (non-leaf) nodes in depth-first order.

        Each returned node represents one retrosynthetic step: it has a
        non-empty ``reaction_smiles`` and at least one child.

        :return: Tuple of all internal nodes in DFS pre-order.
        :rtype: tuple[RetrosyntheticNode, ...]
        """
        return tuple(n for n in self._dfs(self.root) if n.children)

    def get_starting_materials(self) -> tuple[RetrosyntheticNode, ...]:
        """Return all leaf nodes (starting materials) in depth-first order.

        Leaf nodes have ``is_starting_material == True`` (no children).
        Whether each starting material is in the building-block catalog is
        determined separately by :class:`~ursa.BuildingBlockChecker`.

        :return: Tuple of all leaf nodes in DFS pre-order.
        :rtype: tuple[RetrosyntheticNode, ...]
        """
        return tuple(n for n in self._dfs(self.root) if not n.children)

    @cached_property
    def depth(self) -> int:
        """Return the maximum depth of the retrosynthetic tree.

        Depth is defined as the number of edges on the longest path from
        the root to any leaf node.

        :return: Maximum tree depth.
        :rtype: int
        """
        return self._node_depth(self.root)

    @cached_property
    def num_steps(self) -> int:
        """Return the total number of retrosynthetic steps (internal nodes).

        :return: Number of internal nodes in the tree.
        :rtype: int
        """
        return len(self.get_all_steps())

    @classmethod
    def from_dict(cls, data: dict, path_id: str) -> RetrosyntheticPath:
        """Construct a :class:`RetrosyntheticPath` from a RetroCast JSON route.

        The top-level object contains a ``target`` molecule with a recursive
        ``synthesis_step`` / ``reactants`` tree. Example::

            {
                "target": {
                    "smiles": "<target SMILES>",
                    "inchikey": "<InChIKey>",
                    "synthesis_step": {
                        "reactants": [ { ... }, ... ]
                    },
                    "is_leaf": false
                },
                "rank": 1,
                "length": 3
            }

        Extra keys (``rank``, ``length``, ``leaves``, etc.) are ignored.

        :param data: Parsed RetroCast route dict containing a ``target`` key.
        :type data: dict
        :param path_id: Unique identifier to assign to the resulting path.
        :type path_id: str

        :raises KeyError: If ``target`` is missing from ``data``.

        :return: A fully populated :class:`RetrosyntheticPath`.
        :rtype: RetrosyntheticPath
        """
        try:
            target = data["target"]
        except KeyError as exc:
            raise KeyError("RetroCast route dict must contain a 'target' key.") from exc
        root = RetrosyntheticNode._from_mol_dict(target)
        return cls(path_id=path_id, root=root)
