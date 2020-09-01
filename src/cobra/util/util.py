"""General utilities used across the package."""

from typing import Any

from depinfo import print_dependencies


def format_long_string(string: str, max_length: int = 50) -> str:
    """Shorten long string into a small string with ellipsis.

    Parameters
    ----------
    string: str
        The long string to shorten.
    max_length: int, optional
        The maximum length after which to append ellipsis (default 50).

    Returns
    -------
    str
        The shortened string.

    """
    if len(string) > max_length:
        string = f"{string[: max_length - 3]}..."
    return string


class AutoVivification(dict):
    """
    Implementation of Perl's autovivification feature.

    Notes
    -----
    For more information, check https://stackoverflow.com/a/652284/280182 .

    """

    def __getitem__(self, item: Any) -> Any:
        """Retrieve if item is found, else add it.

        Parameters
        ----------
        item: Any
            The object to look for.

        Returns
        -------
        Any
            The retrieved object.

        """
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def show_versions() -> None:
    """Print dependency information."""
    print_dependencies("cobra")
