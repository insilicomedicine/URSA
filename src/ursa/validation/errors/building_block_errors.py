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
