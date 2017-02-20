# -*- coding: utf-8 -*-
from __future__ import absolute_import

from cobra.exceptions import DefunctError


def set_up_optknock(*args, **kwargs):
    raise DefunctError('set_up_optknock',
                       'cameo.strain_design.OptKnock',
                       'https://github.com/biosustain/cameo')


def dual_problem(*args, **kwargs):
    raise DefunctError('dual_problem')
