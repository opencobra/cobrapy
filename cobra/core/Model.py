#cobra.core.Model.py
#TODO: Improve copy time.  Perhaps use references for metabolites in
#reactions.
#Here we're dealing with the fact that numpy isn't supported in jython
#and we've made some numjy modules that use the cern.colt matrices to
#mimic the used numpy behavior.
from warnings import warn
from copy import deepcopy
from ..solvers import optimize
from .Object import Object
from .Formula import Formula
from .DictList import DictList
#*********************************************************************************
#
#TODO: Implement self.reactions[:]._boundary and use

# Note, there are some issues with using a gene as opposed to a protein in the
#nomenclature; this will be addressed after the first iteration of cleaning.
#
# Note, when a reaction is added to the Model it will no longer keep personal
#instances of it's Metabolites, it will reference Model.metabolites to improve
#performance.  When doing this, take care to monitor the metabolite coefficients.
#Do the same for Model.reactions[:].genes and Model.genes
#
#*********************************************************************************

#######################
#BEGIN Class Model
#
class Model(Object):
    """Model is a class for analyzing metabolic models with
    the COBRA toolbox developed in the Palsson Lab at UCSD. Make
    all of the objects (Reaction, Metabolite, ...) to make OOP easier.
    
    """
    def __setstate__(self, state):
        """Make sure all cobra.Objects in the model point to the model
        
        TODO: Make sure that the genes and metabolites referenced by
        the reactions.(genes|metabolites) point to the model's genes
        and metabolites.
        
        """
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

    def copy(self, print_time=False):
        """Provides a partial 'deepcopy' of the Model.  All of the Metabolite, Gene,
        and Reaction objects are created anew but in a faster fashion than deepcopy

        print_time: Boolean used for debugging

        """
        the_copy = Object.guided_copy(self)
        the_copy.metabolites = None
        the_copy.reactions = None
        the_copy.compartments = deepcopy(self.compartments)
        if print_time:
            from time import time
            start_time = time()
        the_metabolites = DictList([x.guided_copy(the_copy)
                                    for x in self.metabolites])
        if print_time:
            print 'Metabolite guided copy: %1.4f'%(time() - start_time)
            start_time = time()
        the_genes = DictList([x.guided_copy(the_copy)
                              for x in self.genes])
        if print_time:
            print 'Gene guided copy: %1.4f'%(time() - start_time)
            start_time = time()
        #TODO: See if we can use the DictList objects instead
        ## metabolite_dict = dict([(k.id, k)
        ##                          for k in the_metabolites])
        ## gene_dict = dict([(k.id, k)
        ##                          for k in the_genes])
        the_reactions = DictList([x.guided_copy(the_copy, the_metabolites._object_dict,
                                                the_genes._object_dict)
                                  for x in self.reactions])

        if print_time:
            print 'Reaction guided copy: %1.4f'%(time() - start_time)
        the_copy.reactions = the_reactions
        the_copy.genes = the_genes
        the_copy.metabolites = the_metabolites
        return the_copy
        

        
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
                print the_reaction.id + ' is not in the model\n'
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
        self.add_reactions(reaction)

        
    def add_reactions(self, reaction_list):
        """Will add a cobra.Reaction object to the model, if
        reaction.id is not in self.reactions.

        reaction_list: A :class:`~cobra.core.Reaction` object or a list of them
      
        """
        #Only add the reaction if one with the same ID is not already
        #present in the model.
        if type(reaction_list) not in [tuple, list, set, DictList]:
            reaction_list = [reaction_list]
        #TODO: Use the DictList properties
        reactions_in_model = set([x.id
                                  for x in reaction_list]).intersection([x.id
                                                                       for x in self.reactions])
        if len(reactions_in_model) > 0:
            print '%i of %i reaction(s) %s already in the model'%(len(reactions_in_model),
                                                          len(reaction_list), repr(reactions_in_model))
            return
        #TODO: Consider using DictList's here, just make sure that the items get appended
        #to self.metabolites, self.genes
        metabolite_dict = {}
        gene_dict = {}
        [metabolite_dict.update(dict([(y.id, y) for y in x._metabolites]))
         for x in reaction_list]
        new_metabolites = [metabolite_dict[x]
                           for x in set(metabolite_dict).difference(self.metabolites._dict)]
        if new_metabolites:
            self.add_metabolites(new_metabolites)

        [gene_dict.update(dict([(y.id, y) for y in x._genes]))
         for x in reaction_list]
        new_genes = [gene_dict[x]
                           for x in set(gene_dict).difference(self.genes._dict)]
        if new_genes:
            self.genes += DictList(new_genes)
            [setattr(x, '_model', self)
             for x in new_genes]

        #This might slow down performance
        #Make sure each reaction knows that it is now part of a Model and uses
        #metabolites in the Model and genes in the Model
        for the_reaction in reaction_list:
            the_reaction._model = self
            the_reaction._metabolites = dict([(self.metabolites.get_by_id(k.id), v)
                                             for k, v in the_reaction._metabolites.iteritems()])
            the_reaction._genes = dict([(self.genes.get_by_id(k.id), v)
                                             for k, v in the_reaction._genes.iteritems()])
            #Make sure the metabolites and genes are aware of the reaction
            the_reaction._update_awareness()
            
        #Add the reactions to the Model
        self.reactions += reaction_list


    def to_array_based_model(self, deepcopy_model=False):
        """Makes a :class:`~cobra.core.ArrayBasedModel` from a cobra.Model which
        may be used to perform linear algebra operations with the
        stoichiomatric matrix.

        deepcopy_model: Boolean.  If False then the ArrayBasedModel points
        to the Model
        
        """
        from .ArrayBasedModel import ArrayBasedModel
        return ArrayBasedModel(self, deepcopy_model=deepcopy_model)


    def optimize(self, new_objective=None, objective_sense='maximize',
                 the_problem=None, solver='glpk', 
                 error_reporting=None, quadratic_component=None,
                 tolerance_optimality=1e-6, tolerance_feasibility=1e-6,
                 tolerance_barrier=1e-10,  **kwargs):
        """Optimize self for self._objective_coefficients or new_objective.

        NOTE: Only the most commonly used parameters are presented here.  Additional
        parameters for cobra.solvers may be available and specified with the
        appropriate keyword=value.

        new_objective: Reaction, String, or Integer referring to a reaction in
        cobra_model.reactions to set as the objective.  Currently, only supports single
        objective coeffients.  Will expand to include mixed objectives.
        
        objective_sense: 'maximize' or 'minimize'
        
        the_problem: None or a problem object for the specific solver that can be used to hot
        start the next solution.

        solver: 'glpk', 'gurobi', or 'cplex'
        
        error_reporting: None or True to disable or enable printing errors encountered
        when trying to find the optimal solution.

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
        the_solution = optimize(self, solver=solver, new_objective=new_objective,
                                objective_sense=objective_sense,
                                the_problem=the_problem,
                                error_reporting=error_reporting,
                                quadratic_component=quadratic_component,
                                tolerance_optimality=tolerance_optimality,
                                tolerance_feasibility=tolerance_feasibility,
                                tolerance_barrier=tolerance_barrier,
                                **kwargs)
        return the_solution

    

    def _update_metabolite_formula(self, metabolite_name, metabolite_formula):
        """Associate metabolite_formula with all self.metabolite_names that
        match metabolite_name.

        metabolite_name: A string from self.metabolite_names

        metabolite_formula: A string specifying the chemical formula for
        metabolite_name.
        
        TODO:  This should be moved to a separate module
        
        """
        if not isinstance(metabolite_formula, Formula):
            metabolite_formula = Formula(metabolite_formula)
        for the_metabolite in self.metabolites:
            if the_metabolite.name == metabolite_name:
                the_metabolite.formula = metabolite_formula


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
                the_reaction.remove_from_model(self)
            except:
                print '%s not in %s'%(the_reaction, self)
        
#
#END Class Model
#####################

