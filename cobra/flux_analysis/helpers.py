# -*- coding: utf-8 -*-

"""Helper functions for all flux analysis methods."""

from __future__ import absolute_import

import logging


LOGGER = logging.getLogger(__name__)


def normalize_cutoff(model, zero_cutoff=None):
    """Return a valid zero cutoff value."""
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


def relax_model_bounds(model, bigM=1e4):
    """
    Relax all upper and lower bounds in the model.
    All positive upper bounds will become bigM.
    All negative lower bounds will become -bigM.
    All positive lower bounds and negative upper bounds will become zero.

    Parameters
    ----------
    model: cobra.Model
        cobra model. It *will* be modified.
    bigM: float, optional
        a large constant for relaxing the model bounds, default 1e4.

    Returns
    -------
    Nothing

    """

    for r in model.reactions:
        r.upper_bound = bigM if r.upper_bound > 0 else 0
        r.lower_bound = -bigM if r.lower_bound < 0 else 0
