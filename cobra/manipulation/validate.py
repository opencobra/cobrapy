# -*- coding: utf-8 -*-

from __future__ import absolute_import

from math import isinf, isnan
from warnings import warn

NOT_MASS_BALANCED_TERMS = {"SBO:0000627",  # EXCHANGE
                           "SBO:0000628",  # DEMAND
                           "SBO:0000629",  # BIOMASS
                           "SBO:0000631",  # PSEUDOREACTION
                           "SBO:0000632",  # SINK
                           }


def check_mass_balance(model):
    unbalanced = {}
    for reaction in model.reactions:
        if reaction.annotation.get("SBO") not in NOT_MASS_BALANCED_TERMS:
            balance = reaction.check_mass_balance()
            if balance:
                unbalanced[reaction] = balance
    return unbalanced


# no longer strictly necessary, done by optlang solver interfaces
def check_reaction_bounds(model):
    warn("no longer necessary, done by optlang solver interfaces",
         DeprecationWarning)
    errors = []
    for reaction in model.reactions:
        if reaction.lower_bound > reaction.upper_bound:
            errors.append("Reaction '%s' has lower bound > upper bound" %
                          reaction.id)
        if isinf(reaction.lower_bound):
            errors.append("Reaction '%s' has infinite lower_bound" %
                          reaction.id)
        elif isnan(reaction.lower_bound):
            errors.append("Reaction '%s' has NaN for lower_bound" %
                          reaction.id)
        if isinf(reaction.upper_bound):
            errors.append("Reaction '%s' has infinite upper_bound" %
                          reaction.id)
        elif isnan(reaction.upper_bound):
            errors.append("Reaction '%s' has NaN for upper_bound" %
                          reaction.id)
    return errors


def check_metabolite_compartment_formula(model):
    errors = []
    for met in model.metabolites:
        if met.compartment is not None and \
                met.compartment not in model.compartments:
            errors.append("Metabolite '%s' compartment '%s' not found" %
                          (met.id, met.compartment))
        if met.formula is not None and len(met.formula) > 0:
            if not met.formula.isalnum():
                errors.append("Metabolite '%s' formula '%s' not alphanumeric" %
                              (met.id, met.formula))
    return errors
