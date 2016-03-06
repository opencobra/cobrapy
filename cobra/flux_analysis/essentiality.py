from warnings import warn
try:
    from cobra.flux_analysis.moma import moma
except:
    warn("moma does not appear to be functional on your system")
from cobra.manipulation import initialize_growth_medium


def assess_medium_component_essentiality(cobra_model, the_components=None,
                                         the_medium=None,
                                         medium_compartment='e', solver='glpk',
                                         the_condition=None, method='fba'):
    """Determines which components in an in silico medium are essential for
    growth in the context of the remaining components.

    cobra_model: A Model object.

    the_components: None or a list of external boundary reactions that will be
    sequentially disabled.

    the_medium: Is None, a string, or a dictionary.  If a string then the
    initialize_growth_medium function expects that the_model has an attribute
    dictionary called media_compositions, which is a dictionary of dictionaries
    for various medium compositions.  Where a medium composition is a
    dictionary of external boundary reaction ids for the medium components and
    the external boundary fluxes for each medium component.

    medium_compartment: the compartment in which the boundary reactions
    supplying the medium components exist

    NOTE: that these fluxes must be negative because the convention is
    backwards means something is feed into the system.

    solver: 'glpk', 'gurobi', or 'cplex'

    returns: essentiality_dict:  A dictionary providing the maximum growth rate
    accessible when the respective component is removed from the medium.

    """
    if method.lower() == 'moma':
        wt_model = cobra_model.copy()
    cobra_model = cobra_model.copy()

    if isinstance(the_medium, str):
        try:
            the_medium = cobra_model.media_compositions[the_medium]
        except:
            raise Exception(
                the_medium + " is not in cobra_model.media_compositions")
    if the_medium is not None:
        initialize_growth_medium(cobra_model, the_medium, medium_compartment)
        if the_components is None:
            the_components = the_medium.keys()
    if not the_components:
        raise Exception("You need to specify the_components or the_medium")
    essentiality_dict = {}
    for the_component in the_components:
        the_reaction = cobra_model.reactions.get_by_id(the_component)
        original_lower_bound = float(the_reaction.lower_bound)
        the_reaction.lower_bound = 0.
        if method.lower() == 'fba':
            cobra_model.optimize(solver=solver)
            objective_value = cobra_model.solution.f
        elif method.lower() == 'moma':
            objective_value = moma(wt_model, cobra_model, solver=solver)[
                'objective_value']
        essentiality_dict[the_component] = objective_value
        the_reaction.lower_bound = original_lower_bound

    return(essentiality_dict)
