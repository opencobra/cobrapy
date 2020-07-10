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
