from typing import Protocol

from ..basic import RetrosyntheticPath


class Checker(Protocol):
    """Protocol for retrosynthetic path checkers.

    All checkers accept a :class:`~ursa.RetrosyntheticPath` and return
    a boolean indicating whether the path passes the check.
    """

    def check(self, path: RetrosyntheticPath) -> bool:
        """Check a retrosynthetic path.

        :param path: The retrosynthetic path to check.
        :type path: RetrosyntheticPath

        :return: ``True`` if the path passes the check, ``False`` otherwise.
        :rtype: bool
        """
        ...
