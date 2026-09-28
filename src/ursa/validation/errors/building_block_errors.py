from .validation_errors import ValidationError


class BuildingBlockError(ValidationError):
    """Base class for building-block catalog errors."""


class CatalogLoadError(BuildingBlockError):
    """Raised when the building-block catalog file cannot be loaded.

    :param catalog_path: Path to the catalog file that failed to load.
    :type catalog_path: str
    :param msg: Human-readable error message.
    :type msg: str | None
    """

    def __init__(self, catalog_path: str, msg: str | None = None) -> None:
        """Initialize CatalogLoadError.

        :param catalog_path: Path to the catalog file that failed to load.
        :type catalog_path: str
        :param msg: Human-readable error message.
        :type msg: str | None
        """
        self.catalog_path = catalog_path
        self.msg = msg or f"Failed to load building-block catalog from {catalog_path!r}"
        super().__init__(self.msg)

    def __str__(self) -> str:
        return self.msg

    def __repr__(self) -> str:
        return (
            f"CatalogLoadError(catalog_path={self.catalog_path!r}, "
            f"msg={self.msg!r})"
        )


class UnsupportedBuildingBlockMatchPolicyError(BuildingBlockError):
    """Raised when :class:`~ursa.BuildingBlockChecker` gets an unknown policy.

    :param policy: The unsupported match-policy value.
    :type policy: object
    """

    def __init__(self, policy: object) -> None:
        """Initialize UnsupportedBuildingBlockMatchPolicyError.

        :param policy: The unsupported match-policy value.
        :type policy: object
        """
        self.policy = policy
        super().__init__(f"unsupported BB match policy: {policy!r}")

    def __repr__(self) -> str:
        return f"UnsupportedBuildingBlockMatchPolicyError(policy={self.policy!r})"
