# -*- coding: utf-8 -*-

from __future__ import absolute_import


NOT_MASS_BALANCED_TERMS = {"SBO:0000627",  # EXCHANGE
                           "SBO:0000628",  # DEMAND
                           "SBO:0000629",  # BIOMASS
                           "SBO:0000631",  # PSEUDOREACTION
                           "SBO:0000632",  # SINK
                           }


def check_mass_balance(model):
    unbalanced = {}
    for reaction in model.reactions:
        if reaction.annotation.get("sbo") not in NOT_MASS_BALANCED_TERMS:
            balance = reaction.check_mass_balance()
            if balance:
                unbalanced[reaction] = balance
    return unbalanced


def check_metabolite_compartment_formula(model):
    errors = []
    for met in model.metabolites:
        if met.formula is not None and len(met.formula) > 0:
            if not met.formula.isalnum():
                errors.append("Metabolite '%s' formula '%s' not alphanumeric" %
                              (met.id, met.formula))
    return errors
