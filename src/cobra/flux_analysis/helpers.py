"""Helper functions for all flux analysis methods."""

from typing import TYPE_CHECKING, Optional


if TYPE_CHECKING:
    from cobra import Model


def normalize_cutoff(model: "Model", zero_cutoff: Optional[float] = None) -> float:
    """Return a valid zero cutoff value.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on.
    zero_cutoff : positive float, optional
        The zero cutoff value. If not specified, defaults to
        `model.tolerance` (default None).

    Returns
    -------
    float
        The normalized zero cutoff value.

    Raises
    ------
    ValueError
        If the specified `zero_cutoff` is lesser than `model.tolerance`.

    """
    if zero_cutoff is None:
        return model.tolerance
    else:
        if zero_cutoff < model.tolerance:
            raise ValueError(
                "The chosen zero cutoff cannot be less than the model's "
                "tolerance value."
            )
        else:
            return zero_cutoff
