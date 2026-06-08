class UrsaError(Exception):
    """Base class for all Ursa errors."""


class CliError(UrsaError):
    """Base class for CLI errors."""


class InvalidTargetKeyError(CliError):
    """Raised when a predictions-file key cannot be parsed as valid SMILES."""

    def __init__(self, key: str) -> None:
        """Initialize InvalidTargetKeyError.

        :param key: Target identifier / SMILES from the raw predictions file.
        :type key: str
        """
        self.key = key
        super().__init__(f"Cannot build target from key {key!r}")

    def __repr__(self) -> str:
        return f"InvalidTargetKeyError(key={self.key!r})"
