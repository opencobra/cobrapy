from warnings import warn
from copy import deepcopy, copy
from ..external.six import iteritems, string_types
from ..solvers import optimize
from .Object import Object
from .Formula import Formula
from .Solution import Solution
from .DictList import DictList



# Note, when a reaction is added to the Model it will no longer keep personal
#instances of it's Metabolites, it will reference Model.metabolites to improve
#performance.  When doing this, take care to monitor the metabolite coefficients.
#Do the same for Model.reactions[:].genes and Model.genes

class Model(Object):
    """Metabolic Model

    Refers to Metabolite, Reaction, and Gene Objects.
    """
    def __setstate__(self, state):
        """Make sure all cobra.Objects in the model point to the model"""
        self.__dict__.update(state)
        [[setattr(x, '_model', self)
          for x in self.__dict__[y]]
         for y in ['reactions', 'genes', 'metabolites']]

    def __init__(self, description=None):
        if isinstance(description, Model):
            self.__dict__ = description.__dict__
        else:
            Object.__init__(self, description)
            self.description = self.id
            self._trimmed = False #This might get changed to a dict of 
            #gene:[reactions in which the gene participates]
            self._trimmed_genes = None #This will be integrated with _trimmed
            self._trimmed_reactions = None #as will this
            self.legacy_format = False #DEPRECATED
            #Allow the creation of an empty object which will facilitate development
            #of SBML parsers and other development issues.
            self.genes = DictList()
            self.reactions = DictList() #A list of cobra.Reactions
            self.metabolites = DictList() #A list of cobra.Metabolites
            #genes based on their ids {Gene.id: Gene}
            self.compartments = {}
            self.solution = Solution(None)


    def __add__(self, other_model):
        """Adds two models. +

        The issue of reactions being able to exists in multiple Models now arises, the same
        for metabolites and such.  This might be a little difficult as a reaction with the
        same name / id in two models might have different coefficients for their metabolites
        due to pH and whatnot making them different reactions.

        """
        new_model = self.copy()
        new_reactions = deepcopy(other_model.reactions)
        new_model.add_reactions(new_reactions)
        new_model.id = self.id + '_' + other_model.id
        return new_model

    def __iadd__(self, other_model):
        """Adds a Model to this model +=

        The issue of reactions being able to exists in multiple Models now arises, the same
        for metabolites and such.  This might be a little difficult as a reaction with the
        same name / id in two models might have different coefficients for their metabolites
        due to pH and whatnot making them different reactions.

        """
        new_reactions = deepcopy(other_model.reactions)
        self.add_reactions(new_reactions)
        self.id = self.id + '_' + other_model.id
        return self

    def guided_copy(self):
        warn("deprecated")
        return self.copy()

    def copy(self, print_time=False):
        """Provides a partial 'deepcopy' of the Model.  All of the Metabolite, Gene,
        and Reaction objects are created anew but in a faster fashion than deepcopy
        """
        if print_time is not False:
            warn("print_time is a deprecated option")
        new = self.__class__()
        do_not_copy = {"metabolites", "reactions", "genes"}
        for attr in self.__dict__:
            if attr not in do_not_copy:
                new.__dict__[attr] = self.__dict__[attr]

        new.metabolites = DictList()
        do_not_copy = {"_reaction", "_model"}
        for metabolite in self.metabolites:
            new_met = metabolite.__class__()
            for attr, value in iteritems(metabolite.__dict__):
                if attr not in do_not_copy:
                    new_met.__dict__[attr] = copy(value) if attr == "formula" else value
            new_met._model = new
            new.metabolites.append(new_met)

        new.genes = DictList()
        for gene in self.genes:
            new_gene = gene.__class__(None)
            for attr, value in iteritems(gene.__dict__):
                if attr not in do_not_copy:
                    new_gene.__dict__[attr] = copy(value) if attr == "formula" else value
            new_gene._model = new
            new.genes.append(new_gene)

        new.reactions = DictList()
        do_not_copy = {"_model", "_metabolites", "_genes"}
        for reaction in self.reactions:
            new_reaction = reaction.__class__()
            for attr, value in iteritems(reaction.__dict__):
                if attr not in do_not_copy:
                    new_reaction.__dict__[attr] = value
            new_reaction._model = new
            new.reactions.append(new_reaction)
            # update awareness
            for metabolite, stoic in iteritems(reaction._metabolites):
                new_met = new.metabolites.get_by_id(metabolite.id)
                new_reaction._metabolites[new_met] = stoic
                new_met._reaction.add(new_reaction)
            for gene in reaction._genes:
                new_gene = new.genes.get_by_id(gene.id)
                new_reaction._genes.add(new_gene)
                new_gene._reaction.add(new_reaction)
        return new


    def add_metabolites(self, metabolite_list):
        """Will add a list of metabolites to the the object, if they do not
        exist and then expand the stochiometric matrix

        metabolite_list: A list of :class:`~cobra.core.Metabolite` objects

        """
        if not hasattr(metabolite_list, '__iter__'):
            metabolite_list = [metabolite_list]
        #First check whether the metabolites exist in the model
        metabolite_list = [x for x in metabolite_list
                           if x.id not in self.metabolites]
        [setattr(x, '_model', self) for x in metabolite_list]
        self.metabolites += metabolite_list


    def _update_reaction(self, reaction):
        """Updates everything associated with the reaction.id of reaction.

        reaction: A cobra.Reaction object, or a list of these objects.

        """
        warn("WARNING: To be modified.  It will be faster to use properties of the DictList " +\
             "self.reactions to find a reaction and update the matrices " +\
             "This function is only used after the Model has been " +\
             "converted to matrices.  It is typically faster to access the objects" +\
             "in the Model directly.  This function will eventually moved to another" +\
             "module for advanced users due to the potential for mistakes.")

        if not hasattr(reaction, '__iter__'):
            reaction = [reaction]
        for the_reaction in reaction:
            if the_reaction.id not in self.reactions:
                warn(the_reaction.id + ' is not in the model')
                continue
            reaction_index = self.reactions.index(the_reaction.id)
            self.reactions[reaction_index] = the_reaction

    def update(self):
        """Non functional.  Model.update is moved to ArrayBasedModel.  Please use
        the to_array_based_model property to create an ArrayBasedModel.
        
        """
        raise Exception("Model.update is moved to ArrayBasedModel.  Please use \n"
                        "the to_array_based_model property to create an ArrayBasedModel.")


    def add_reaction(self, reaction):
        """Will add a cobra.Reaction object to the model, if
        reaction.id is not in self.reactions.

        reaction: A :class:`~cobra.core.Reaction` object

        """
        self.add_reactions([reaction])


    def add_reactions(self, reaction_list):
        """Will add a cobra.Reaction object to the model, if
        reaction.id is not in self.reactions.

        reaction_list: A list of :class:`~cobra.core.Reaction` objects
      
        """
        #Only add the reaction if one with the same ID is not already
        #present in the model.
        
        # This function really should not used for single reactions
        if not hasattr(reaction_list, "__len__"):
            reaction_list = [reaction_list]
            warn("Use add_reaction for single reactions")

        reaction_list = DictList(reaction_list)
        reactions_in_model = [i.id for i in reaction_list if self.reactions.has_id(i.id)]
        
        if len(reactions_in_model) > 0:
            raise Exception("Reactions already in the model: " + \
                ", ".join(reactions_in_model))

        # Add reactions. Also take care of genes and metabolites in the loop
        for reaction in reaction_list:
            reaction._model = self  # the reaction now points to the model
            # keys() is necessary because the dict will be modified during
            # the loop
            for metabolite in list(reaction._metabolites.keys()):
                # if the metabolite is not in the model, add it
                # should we be adding a copy instead.
                if not self.metabolites.has_id(metabolite.id):
                    self.metabolites.append(metabolite)
                    metabolite._model = self
                    # this should already be the case. Is it necessary?
                    metabolite._reaction = set([reaction])
                # A copy of the metabolite exists in the model, the reaction
                # needs to point to the metabolite in the model.
                else:
                    stoichiometry = reaction._metabolites.pop(metabolite)
                    model_metabolite = self.metabolites.get_by_id(metabolite.id)
                    reaction._metabolites[model_metabolite] = stoichiometry
                    model_metabolite._reaction.add(reaction)

            for gene in list(reaction._genes):
                # If the gene is not in the model, add it
                if not self.genes.has_id(gene.id):
                    self.genes.append(gene)
                    gene._model = self
                    # this should already be the case. Is it necessary?
                    gene._reaction = set([reaction])
                # Otherwise, make the gene point to the one in the model
                else:
                    model_gene = self.genes.get_by_id(gene.id)
                    if model_gene is not gene:
                        reaction._dissociate_gene(gene)
                        reaction._associate_gene(model_gene)

        self.reactions += reaction_list


    def to_array_based_model(self, deepcopy_model=False, **kwargs):
        """Makes a :class:`~cobra.core.ArrayBasedModel` from a cobra.Model which
        may be used to perform linear algebra operations with the
        stoichiomatric matrix.

        deepcopy_model: Boolean.  If False then the ArrayBasedModel points
        to the Model
        
        """
        from .ArrayBasedModel import ArrayBasedModel
        return ArrayBasedModel(self, deepcopy_model=deepcopy_model, **kwargs)


    def optimize(self, objective_sense='maximize',
                 solver=None, 
                 error_reporting=None, quadratic_component=None,
                 tolerance_optimality=1e-6, tolerance_feasibility=1e-6,
                 tolerance_barrier=1e-10,  **kwargs):
        """Optimize self for self._objective_coefficients or new_objective.

        NOTE: Only the most commonly used parameters are presented here.  Additional
        parameters for cobra.solvers may be available and specified with the
        appropriate keyword=value.

        new_objective: Reaction, String, or Integer referring to a reaction in
        cobra_model.reactions to set as the objective.  In the case where the new_objective
        is a linear combination of Reactions then new_objective should be a dictionary where
        the key is the Reaction and the value is the objective coefficient
        (e.g., 0.1 reaction_1 + 0.5 reaction_2 would be new_objective = {reaction_1: 0.1, reaction_2: 0.5})
        
        objective_sense: 'maximize' or 'minimize'
        
        the_problem: None or a problem object for the specific solver that can be used to hot
        start the next solution.

        solver: 'glpk', 'gurobi', or 'cplex'


        quadratic_component: None or 
        scipy.sparse.dok of dim(len(cobra_model.reactions),len(cobra_model.reactions))
        If not None:
          Solves quadratic programming problems for cobra_models of the form:
          cobra_model = ArrayBasedModel(cobra_model)
          minimize: 0.5 * x' * quadratic_component * x + cobra_model._objective_coefficients' * x
          such that,
            cobra_model._lower_bounds <= x <= cobra_model._upper_bounds
            cobra_model._S * x (cobra_model._constraint_sense) cobra_model._b


        #See cobra.flux_analysis.solvers for more info on the following parameters.  Also,
        refer to your solver's manual
        
        tolerance_optimality: Solver tolerance for optimality.
            
        tolerance_feasibility: Solver tolerance for feasibility.

        tolerance_barrier: Solver tolerance for barrier method

        lp_method: Solver method to solve the problem

        #End solver parameters
        
        **kwargs: See additional parameters for your specific solver module in
        cobra.solvers

        
        """
        if "new_objective" in kwargs:
            warn("new_objective is deprecated. Use Model.change_objective")
            self.change_objective(kwargs.pop("new_objective"))
        the_solution = optimize(self, solver=solver,
                                objective_sense=objective_sense,
                                quadratic_component=quadratic_component,
                                tolerance_optimality=tolerance_optimality,
                                tolerance_feasibility=tolerance_feasibility,
                                tolerance_barrier=tolerance_barrier,
                                **kwargs)
        self.solution = the_solution
        return the_solution


    def remove_reactions(self, the_reactions):
        """
        the_reactions: instance or list of cobra.Reactions or strings of
        self.reactions[:].id.

        """
        if not hasattr(the_reactions, '__iter__') or \
               hasattr(the_reactions, 'id'):
            the_reactions = [the_reactions]
        if hasattr(the_reactions[0], 'id'):
            the_reactions = [x.id for x in the_reactions]
        reactions_to_delete = []
        for the_reaction in the_reactions:
            try:
                the_reaction = self.reactions[self.reactions.index(the_reaction)]
                the_reaction.remove_from_model()
            except:
                warn('%s not in %s'%(the_reaction, self))

    def repair(self):
        """Update all indexes and pointers in a model"""
        # DictList indexes
        self.reactions._generate_index()
        self.metabolites._generate_index()
        self.genes._generate_index()
        return  # TODO update the pointers as well

    def change_objective(self, objectives):
        """Change the objective in the cobrapy model.
        
        objectives: A list or a dictionary.  If a list then
        a list of reactions for which the coefficient in the
        linear objective is set as 1.  If a dictionary then the
        key is the reaction and the value is the linear coefficient
        for the respective reaction.

        """
        # I did not want to refactor code just to rename the variable, but this
        # way the API uses the variable "objectives"
        the_objectives = objectives
        # set all objective coefficients to 0 initially
        for x in self.reactions:
            x.objective_coefficient = 0.
        # update the objective coefficients if a dict is passed in
        if hasattr(the_objectives, "items"):
            for the_reaction, the_coefficient in iteritems(the_objectives):
                if isinstance(the_reaction, int):
                    the_reaction = self.reactions[the_reaction]
                else:
                    if hasattr(the_reaction, 'id'):
                        the_reaction = the_reaction.id
                    the_reaction = self.reactions.get_by_id(the_reaction)
                the_reaction.objective_coefficient = the_coefficient
        # If a list (or a single reaction is passed in), each reaction gets
        # 1 for the objective coefficent.
        else:
            # Allow for objectives to be constructed from multiple reactions
            if not hasattr(the_objectives, "__iter__") or \
                    isinstance(the_objectives, string_types):
                the_objectives = [the_objectives]
            for the_reaction in the_objectives:
                if isinstance(the_reaction, int):
                    the_reaction = self.reactions[the_reaction]
                else:
                    if hasattr(the_reaction, 'id'):
                        the_reaction = the_reaction.id
                    the_reaction = self.reactions.get_by_id(the_reaction)
                the_reaction.objective_coefficient = 1.
