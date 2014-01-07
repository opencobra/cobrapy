import cobra
from keggIO import import_kegg_reactions
from cobra import Model, Reaction, Metabolite

class SUXModelMILP(cobra.Model):
    """Model with additional Universal and Exchange reactions.
    Adds corresponding dummy reactions and dummy metabolites for each added
    reaction which are used to impose MILP constraints to minimize the
    total number of added reactions. See the figure for more
    information on the structure of the matrix.
    """
    def __init__(self, model, Universal=None,
            threshold=0.05, penalties={"Universal": 1, "Exchange": 1, "Demand": 1}, dm_rxns=0, ex_rxns=1):
        cobra.Model.__init__(self, "")
        # store parameters
        self.threshold = threshold
        self.penalties = penalties
        # want to only operate on a copy of Universal
        
        if Universal is None:
            Universal = cobra.Model("Universal_Reactions")
        else:
            Universal = Universal.copy()  #did this before - so as not to mess up the Universal Model
            
        
        
        
        
        # SUX += Exchange (when exchange generator has been written)
        # For now, adding exchange reactions to Universal - could add to a new model called exchange and allow their addition or not....
        if ex_rxns == 1:
            ex_reactions = [x for x in model.reactions if x.startswith('EX_')]
        else: 
            ex_reactions = [] 
            
        # ADD ALL EXCHANGE REACTIONS TO UNIVERSAL MODEL
        for r in ex_reactions:
            if r.lower_bound >= 0:
                rxn = r.copy()
                r.remove_from_model(model)
                rxn.lower_bound = -1000
                
                Universal.add_reaction(rxn)
        
        if dm_rxns == 1:
            # ADD DEMAND REACTIONS FOR ALL METABOLITES TO UNIVERSAL MODEL
            for m in model.metabolites:
                rxn = Reaction('DM_' + m.id)
                rxn.lower_bound = -1000
                rxn.upper_bound = 1000
                rxn.add_metabolites({m: -1.0})
                Universal.add_reaction(rxn)
    
        cobra.manipulation.modify.convert_to_irreversible(Universal)
        
        for rxn in Universal.reactions:
            if rxn.startswith('EX_'):
                rxn.notes["gapfilling_type"] = "Exchange"
            elif rxn.startswith('DM_'):
                rxn.notes["gapfilling_type"] = "Demand"
            else:
                rxn.notes["gapfilling_type"] = "Universal"    
        
        self += model
        self += Universal

        # Add MILP dummy reactions
        v = 1000  # maximum flux in a reaction
        # threshold = 0.05
        dummy_reactions = []
        # all reactions with an index < len(model.reactions) were original
        self.original_reactions = self.reactions[:len(model.reactions)]
        self.added_reactions = self.reactions[len(model.reactions):]

        # add in the dummy reactions for each added reaction
        # a dict will map from each added reaction (the key) to
        # the dummy reaction (the value)
        self._dummy_reaction_map = {}
        for reaction in self.added_reactions:
            dummy_metabolite = cobra.Metabolite("dummy_met_" + reaction.id)
            dummy_metabolite._constraint_sense = "L"
            reaction.add_metabolites({dummy_metabolite: 1})
            the_dummy_reaction = cobra.Reaction("dummy_rxn_" + reaction.id)
            the_dummy_reaction.add_metabolites({dummy_metabolite: -1 * v})
            the_dummy_reaction.lower_bound = 0
            the_dummy_reaction.upper_bound = 1
            the_dummy_reaction.variable_kind = "integer"
            dummy_reactions.append(the_dummy_reaction)
            self._dummy_reaction_map[reaction] = the_dummy_reaction
        self.add_reactions(dummy_reactions)
        # add in the dummy metabolite for the actual objective function
        self.objective_metabolite = cobra.Metabolite(
            "dummy_metabolite_objective_function")
        self.objective_metabolite._constraint_sense = "G"
        self.objective_metabolite._bound = self.threshold
        self._update_objectives()
        # make .add_reaction(s) call the ._add_reaction(s) functions
        self.add_reaction = self._add_reaction # CHANGED THIS FROM self.add_reation on 9-24-2012
        self.add_reactions = self._add_reactions # CHANGED THIS FROM self.add_reations on 9-24-2012
        # embed()
        
    def _update_objectives(self):
        """Update the metabolite which encodes the objective function
        with the objective coefficients for the reaction, and impose
        penalties for added reactions.
        """
        for reaction in self.original_reactions:
            reaction.add_metabolites({self.objective_metabolite: \
                                    reaction.objective_coefficient})
            reaction.objective_coefficient = 0
        # now make the objective coefficient the penalty
        for reaction in self.added_reactions:
            reaction.objective_coefficient += \
                self.penalties[reaction.notes["gapfilling_type"]]

    def _add_reaction(self, reaction):
        cobra.Model.add_reaction(self, reaction)
        self.original_reactions.append(reaction)
        self._update_objectives()

    def _add_reactions(self, reactions):
        cobra.Model.add_reactions(self, reactions)
        self.original_reactions.extend(reaction)
        self._update_objectives()

    def solve(self, iterations=1, debug=0):
        """solve the MILP problem"""
        used_reactions = {}
        numeric_error_cutoff = 0.0001
        self._update_objectives()
        # TODO implement iterations part of the algorithm
        for i in range(iterations):
            used_reactions[i] = []
            self.optimize(objective_sense="minimize")

            if debug: print 'iteration', i
            for reaction in self.added_reactions:
                # The dummy reaction should have a flux of either 0 or 1.
                # If it is 1 (nonzero), then the reaction was used in
                # the solution.
                
                if self.solution.x_dict[self._dummy_reaction_map[
                        reaction].id] > numeric_error_cutoff:
                    used_reactions[i].append(reaction)
                    reaction.objective_coefficient += self.penalties[reaction.notes["gapfilling_type"]]
                    if debug: print '\t', reaction, reaction.objective_coefficient

        return used_reactions
        


def growMatch(model, Universal=None, iterations=1, debug=0, dm_rxns=0, ex_rxns=1):
    """runs growMatch"""
    # if Universal is None:
        # Universal = import_kegg_reactions()
    model = model.copy()
    SUX = SUXModelMILP(model, Universal, dm_rxns=dm_rxns, ex_rxns=ex_rxns)
    used_reactions = SUX.solve(iterations=iterations, debug=debug)
    return used_reactions


def SMILEY(model, metabolite_id, Universal=None):
    """
    runs the SMILEY algorithm to determine which gaps should be
    filled in order for the model to create the metabolite with the
    given metabolite_id.

    This function is good for running the algorithm once. For more fine-
    grained control, create a SUXModelMILP object, add a demand reaction
    for the given metabolite_id, and call the solve function on the
    SUXModelMILP object.
    """
    if Universal is None:
        Universal = import_kegg_reactions()
    SUX = SUXModelMILP(model, Universal)
    # change the objective to be the metabolite
    for reaction in SUX.original_reactions:
        reaction.objective_coefficient = 0
    demand_reaction = cobra.Reaction("SMILEY_DEMAND_RXN_%s" % metabolite_id)
    demand_reaction.objective_coefficient = 1
    demand_reaction.add_metabolites(
        {SUX.metabolites[SUX.metabolites.index(metabolite_id)]: -1})
    SUX.add_reaction(demand_reaction)
    used_reactions = SUX.solve()
    for reaction in used_reactions:
        print reaction, SUX.solution.x_dict[reaction.id]
    return used_reactions
    
def add_gap_reactions(model, results):
    for r in results:
        reverse = 0
        #determine if reaction should be reversible or not...
        index = r.id.find('reverse') # =-1 if not found
        if index >= 0:
            rxn = r.id[:index-1]
            reverse = 1 # need to add functionality here...
        else:
            rxn = r.id
        
        # TO DO determine if reaction is already in model
                      
        #add reaction or adjust bounds
        model.reactions[model.reactions.index(rxn)].lower_bound = -1
        
    model.optimize()
        
    try:
        if model.solution.f < 0.00001:
            model.solution.f = 0
        print '\nNew model simulated growth rate is %1.3f'%model.solution.f
        
        # fout.write(org_group + ',' + org_name + ',' + org_id + ',' + str(rxn) + ',' + str(model.solution.f) + '\n')
    except:
        print 'no growth'
        
    return model
    


if __name__ == "__main__":
    from cPickle import load
    from os.path import join, abspath, dirname
    from time import time

    import cobra

    test_file_path = join(dirname(cobra.__file__), "test", "data", \
                          "salmonella.pickle")
    test_file = open(test_file_path, "rb")
    model = load(test_file)
    test_file.close()
    
    model.optimize(solver='glpk')
    
  
    # from cobra.flux_analysis.single_deletion import *
    # deletions = single_reaction_deletion(model)
    # deletions show that SDPTA is essential to growth:
    r = model.reactions[model.reactions.index('SDPTA')]
    test = r.copy()
    r.remove_from_model(model)
    
    
    tic = time()
    # Universal = import_kegg_reactions()
    Universal = cobra.Model("Kegg_Universal_Reactions")
    # embed()
    
    
    Universal.add_reaction(test)
    
    
    toc = time()
    print "%.2f sec to import kegg reactions" % (toc - tic)
    tic = toc
    print "growMatch: "
    results = growMatch(model, Universal)
    # toc = time()
    # print "%.2f sec for growmatch" % (toc - tic)
    # tic = toc
    # print "SMILEY results"
    # SMILEY(model, "atp_c", Universal)
    # toc = time()
    # print "%.2f sec for smiley" % (toc - tic)
