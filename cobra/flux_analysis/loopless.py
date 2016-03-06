from ..core import Reaction, Metabolite
from ..manipulation.modify import convert_to_irreversible
from six import iteritems


def construct_loopless_model(cobra_model):
    """construct a loopless model

    This adds MILP constraints to prevent flux from proceeding in a loop, as
    done in http://dx.doi.org/10.1016/j.bpj.2010.12.3707
    Please see the documentation for an explanation of the algorithm.

    This must be solved with an MILP capable solver.

    """
    # copy the model and make it irreversible
    model = cobra_model.copy()
    convert_to_irreversible(model)
    max_ub = max(model.reactions.list_attr("upper_bound"))
    # a dict for storing S^T
    thermo_stoic = {"thermo_var_" + metabolite.id: {}
                    for metabolite in model.metabolites}
    # Slice operator is so that we don't get newly added metabolites
    original_metabolites = model.metabolites[:]
    for reaction in model.reactions[:]:
        # Boundary reactions are not subjected to these constraints
        if len(reaction._metabolites) == 1:
            continue
        # populate the S^T dict
        bound_id = "thermo_bound_" + reaction.id
        for met, stoic in iteritems(reaction._metabolites):
            thermo_stoic["thermo_var_" + met.id][bound_id] = stoic
        # I * 1000 > v --> I * 1000 - v > 0
        reaction_ind = Reaction(reaction.id + "_indicator")
        reaction_ind.variable_kind = "integer"
        reaction_ind.upper_bound = 1
        reaction_ub = Metabolite(reaction.id + "_ind_ub")
        reaction_ub._constraint_sense = "G"
        reaction.add_metabolites({reaction_ub: -1})
        reaction_ind.add_metabolites({reaction_ub: max_ub})
        # This adds a compensating term for 0 flux reactions, so we get
        # S^T x - (1 - I) * 1001 < -1 which becomes
        # S^T x < 1000 for 0 flux reactions and
        # S^T x < -1 for reactions with nonzero flux.
        reaction_bound = Metabolite(bound_id)
        reaction_bound._constraint_sense = "L"
        reaction_bound._bound = max_ub
        reaction_ind.add_metabolites({reaction_bound: max_ub + 1})
        model.add_reaction(reaction_ind)
    for metabolite in original_metabolites:
        metabolite_var = Reaction("thermo_var_" + metabolite.id)
        metabolite_var.lower_bound = -max_ub
        model.add_reaction(metabolite_var)
        metabolite_var.add_metabolites(
            {model.metabolites.get_by_id(k): v
             for k, v in iteritems(thermo_stoic[metabolite_var.id])})
    return model
