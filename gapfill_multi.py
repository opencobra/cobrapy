def gapfill_multi(model, missing_genes, **kwargs):
    
    """
    
    Generate a list of gapfilling reactions from a list of missing genes for a strain-specific model.
    
    :param model: COBRA model for the base strain with the objective coefficient for the reaction of interest (e.g. biomass reaction) set to 1.
    
    :param missing_genes: list of genes with no homologs in the strain of interest.
    
    :param lower_bound: minimum allowable yield of gapfilled model.
    
    :param biomass: override the current model settings and temporarily assign the objective coefficient for a function of interest to 1.
    
    :return: a list of gapfilling reactions.
    
    """
    
    #     model.solver = 'gurobi'

    if 'lower_bound' in kwargs.keys():
        lower_bound = kwargs['lower_bound']
    else:
        lower_bound = model.optimize().objective_value*0.5
        
    biomass_reactions = [rx.id for rx in model.reactions if rx.objective_coefficient == 1]
    if 'biomass' in kwargs.keys():
        biomass = kwargs['biomass']
        if len(biomass_reactions) > 1:
            for rx in set(biomass_reactions) - {biomass}:
                model.reactions.get_by_id(rx).objective_coefficient = 0
                 
    else:
        if len(biomass_reactions) > 1:
            raise Exception("This model has more than one objective. \n Please adjust the objective coefficient to 1 for the chosen objective reaction (e.g. biomass or ATP) and 0 for the rest of the reactions, \n or specify the reaction ID to use as an objective.")
        if len(biomass_reactions) > 1:
            raise Exception("The model doesn't have an objective function. Please set the appropriate objective coefficient to 1, or specify the reaction ID to use as an objective.")
        biomass = biomass_reactions[0]
        
        
    model.solver.configuration.tolerances.feasibility = 1e-9
    constraints = []
    indicators = []

    for rx in cobra.manipulation.find_gene_knockout_reactions(model, missing_genes):

        indicator = model.problem.Variable('%s_i'%rx.id , type = 'binary')
        indicators.append(indicator)

        new_cstr1 = model.problem.Constraint( rx.flux_expression - rx.upper_bound*indicator ,ub = 0)
        new_cstr2 = model.problem.Constraint(-rx.flux_expression + rx.lower_bound*indicator ,ub = 0)
        constraints += [new_cstr1, new_cstr2]
        model.add_cons_vars([new_cstr1, new_cstr2, indicator])

    model.reactions.get_by_id(biomass).lower_bound = lower_bound
    model.objective = model.problem.Objective(-sum(indicators))
    sol = model.optimize()
    indicator_results = [ind.name[:-2] for ind in indicators if ind.primal != 0.0]
    
    # removing changes to model
    model.remove_cons_vars(constraints+indicators)
    for rx in set(biomass_reactions):
        model.reactions.get_by_id(rx).objective_coefficient = 1   
        
    return indicator_results