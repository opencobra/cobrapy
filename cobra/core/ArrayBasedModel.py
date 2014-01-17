#Dresses a cobra.Model with arrays and vectors so that linear algebra operations
#can be carried out on 
from sys import maxint
from copy import deepcopy
from warnings import warn

from numpy import array, hstack, ndarray
from scipy.sparse import lil_matrix, dok_matrix

from .Model import Model


class ArrayBasedModel(Model):
    """ArrayBasedModel is a class that adds arrays and vectors to
    a cobra.Model to make it easier to perform linear algebra operations.
    
    S : :class:`scipy.sparse.lil_matrix` (in CPython)
        Stoichiometric matrix of the model

    lower_bounds

    upper_bounds

    objective_coefficients

    b

    constraint_sense

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

    def __init__(self, description=None, deepcopy_model=False, matrix_type='scipy.lil_matrix'):
        """
        description: None | String | cobra.Model

        deepcopy_model: Boolean.  If True and description is
        a cobra.Model then make a deepcopy of the Model before
        creating the ArrayBasedModel.

        matrix_type: String. Specifies which type of backend matrix to use for self.S.
        Currently, 'scipy.lil_matrix' or 'scipy.dok_matrix'
        
        """
        if deepcopy_model and isinstance(description, Model):
            description = description.copy()
        Model.__init__(self, description)
        self.S = None
        self.matrix_type = matrix_type
        self.b = None
        self.constraint_sense = None
        self.update()

    @property
    def lower_bounds(self):
        return self._lower_bounds

    @property
    def upper_bounds(self):
        return self._upper_bounds

    @property
    def objective_coefficients(self):
        return self._objective_coefficients


    def copy(self, print_time=False):
        """Provides a partial 'deepcopy' of the Model.  All of the Metabolite, Gene,
        and Reaction objects are created anew but in a faster fashion than deepcopy

        print_time: Boolean used for debugging

        """
        the_copy = Model.copy(self, print_time)
        the_copy.update()
        return the_copy
        

    def add_metabolites(self, metabolite_list,
                        expand_stoichiometric_matrix=True):
        """Will add a list of metabolites to the the object, if they do not
        exist and then expand the stochiometric matrix

        metabolite_list: A list of :class:`~cobra.core.Metabolite` objects

        expand_stoichimetric_matrix: Boolean.  If True and self.S is
        not None then it will add rows to self.S.  self.S must be
        created after adding reactions and metabolites to self before
        it can be expanded.  Trying to expand self.S when self only
        contains metabolites is ludacris.

        """
        Model.add_metabolites(metabolite_list)
        self.constraint_sense = [x._constraint_sense for x in self.metabolites]
        if self.S is not None and expand_stoichiometric_matrix:
            s_expansion = len(self.metabolites) - self.S.shape[0]
            if s_expansion > 0:
                self.S = self.S.todok()
                self.S.resize((self.S.shape[0] + s_expansion,
                               self.S.shape[1]))
                self.S = self.S.tolil()


    def _update_reaction(self, reaction):
        """Updates everything associated with the reaction.id of reaction.

        reaction: A cobra.Reaction object, or a list of these objects.

        WARNING: This function is only used after the Model has been
        converted to matrices.  It is typically faster to access the objects
        in the Model directly.  This function will eventually moved to another
        module for advanced users due to the potential for mistakes.

    
        """
        if not hasattr(reaction, '__iter__'):
            reaction = [reaction]
        Model._update_reaction(reaction)
        for the_reaction in reaction:
            try:
                reaction_index = self.reactions.index(the_reaction.id)
            except KeyError:
                print the_reaction.id + ' is not in the model'
                continue

            #zero reaction stoichiometry column
            the_column = self.S[:, reaction_index]
            for nonzero_index in the_column.nonzero()[0]:
                the_column[nonzero_index, 0] = 0
            self._lower_bounds[reaction_index] = the_reaction.lower_bound
            self._upper_bounds[reaction_index] = the_reaction.upper_bound
            self.objective_coefficients[reaction_index] = the_reaction.objective_coefficient
            self.add_metabolites(the_reaction._metabolites)
            #Make sure that the metabolites are the ones contained in the model
            the_reaction._metabolites = [self.metabolite.get_by_id(x.id)
                                        for x in the_reaction._metabolites]
            #Update the stoichiometric matrix
            metabolite_indices = map(self.metabolites.index, the_reaction._metabolites)
            for (index, metabolite_index) in enumerate(metabolite_indices):
                self.S[metabolite_index, reaction_index] = the_reaction.stoichiometric_coefficients[index]


        
    def add_reactions(self, reaction_list, update_matrices=False):
        """Will add a cobra.Reaction object to the model, if
        reaction.id is not in self.reactions.

        reaction_list: A :class:`~cobra.core.Reaction` object or a list of them

        update_matrices:  Boolean.  If true populate / update matrices
        S, lower_bounds, upper_bounds, .... Note this is slow to run
        for very large models and using this option with repeated calls
        will degrade performance.  Better to call self.update() after
        adding all reactions.

        
         If the stoichiometric matrix is initially empty then initialize a 1x1
         sparse matrix and add more rows as needed in the self.add_metabolites
         function

        """
        Model.add_reactions(reaction_list)
        if update_matrices:
            self._update_matrices(reaction_list)

    def _construct_matrices(self):
        """Large sparse matrices take time to construct and to read / write.
        This function allows one to let the model exists without cobra_model.S
        and then generate it at needed.
        
        """
        self._update_matrices() #This does basic construction as well.

    def _update_reaction_vectors(self):
        """regenerates the lower_bounds, upper_bounds,
        and objective_coefficients vectors.

        WARNING: This function is only used after the Model has been
        converted to matrices.  It is typically faster to access the objects
        in the Model directly.  This function will eventually moved to another
        module for advanced users due to the potential for mistakes.


        """
        self._lower_bounds = LinkedArray(self.reactions, "lower_bound")
        self._upper_bounds = LinkedArray(self.reactions, "upper_bound")
        self._objective_coefficients = LinkedArray(self.reactions,
                                                   "objective_coefficient")


    # TODO deprecate and use @property for _b and constraint sense
    def _update_metabolite_vectors(self):
        """regenerates _b and _constraint_sense

        WARNING: This function is only used after the Model has been
        converted to matrices.  It is typically faster to access the objects
        in the Model directly.  This function will eventually moved to another
        module for advanced users due to the potential for mistakes.

        """
        _b, _constraint_sense = [], []
        [(_b.append(x._bound),
          _constraint_sense.append(x._constraint_sense))
         for x in self.metabolites]
        self.b = array(_b)
        self.constraint_sense = _constraint_sense
 
    # TODO deprecate and use @property
    def _update_matrices(self, reaction_list=None):
        """
        reaction_list: None or a list of cobra.Reaction objects that are in
        self.reactions.  If None then reconstruct the whole matrix.

        NOTE: reaction_list is assumed to be at the end of self.reactions.

        In the future, we'll be able to use reactions from anywhere in the
        list

        
        WARNING: This function is only used after the Model has been
        converted to matrices.  It is typically faster to access the objects
        in the Model directly.  This function will eventually moved to another
        module for advanced users due to the potential for mistakes.

        """
        #Pretty much all of these things are unnecessary to use the objects and
        #interact with the optimization solvers.  It might be best to move them
        #to linear algebra modules.
        #If no reactions are present in the Model, initialize the arrays
        if not self.S or reaction_list is None:
            reaction_list = self.reactions
            SMatrix = SMatrix_classes[self.matrix_type]
            self.S = SMatrix((len(self.metabolites),
                              len(self.reactions)), model=self) 
            self._update_reaction_vectors()
        else: #Expand the arrays to accomodate the new reaction
            self.S.resize((len(self.metabolites),
                           len(self.reactions)))
            lower_bounds = array([x.lower_bound
                                  for x in reaction_list])
            upper_bounds = array([x.upper_bound
                                  for x in reaction_list])
            objective_coefficients = array([x.objective_coefficient
                                            for x in reaction_list])
            self._lower_bounds = hstack((self._lower_bounds, lower_bounds)).view(LinkedArray)
            self._lower_bounds.__init__(self.reactions, "lower_bound")
            self._upper_bounds = hstack((self._upper_bounds, upper_bounds))
            self._upper_bounds.__init__(self.reactions, "upper_bound")
            self._objective_coefficients = hstack((self._objective_coefficients,
                                                  objective_coefficients))
            self._objective_coefficients.__init__(self.reactions,
                                                  "objective_coefficient")

        #Use dok format to speed up additions.
        coefficient_dictionary = {}
        #SPEED this up. This is the slow part.  Can probably use a dict to index.
        for the_reaction in reaction_list:
            reaction_index = self.reactions.index(the_reaction.id)
            for the_key, the_value in the_reaction._metabolites.items():
                coefficient_dictionary[(self.metabolites.index(the_key.id),
                                        reaction_index)] = the_value

        if not self.S.getformat() == 'dok':
            self.S = self.S.todok()
        self.S.update(coefficient_dictionary)
        if self.matrix_type == 'scipy.lil_matrix':
            try:
                self.S = self.S.tolil()
            except Exception, e:
                warn('Unable to convert S to lil_matrix maintaining in dok_matrix format: %s'%e)

    def update(self):
        """Regenerates the stoichiometric matrix and vectors
        
        """
        self._update_matrices()
        self._update_metabolite_vectors()


class LinkedArray(ndarray):
    """A ndarray which updates an attribute from a list"""
    def __new__(cls, list, attribute):
        # construct a new ndarray with the values from the list
        # For example, if the list if model.reactions and the attribute is
        # "lower_bound" create an array of [reaction.lower_bound for ... ]
        return array([getattr(i, attribute) for i in list]).view(cls)
    def __init__(self, list, attribute):
        self._list = list
        self._attr = attribute
        
    def __setitem__(self, index, value):
        numpy.ndarray.__setitem__(self, index, value)
        setattr(self._list[index], self._attr, value)
    def __setslice__(self, i, j, value):
        numpy.ndarray.__setslice__(self, i, j, value)
        if j == maxint:
            j = len(self)
        if hasattr(value, "__getitem__"):  # setting to a list
            for index in range(i, j):
                setattr(self._list[index], self._attr, value[index])
        else:
            for index in range(i, j):
                setattr(self._list[index], self._attr, value)


class SMatrix_dok(dok_matrix):
    """A 2D sparse dok matrix which maintains links to a cobra Model"""
    def __init__(self, *args, **kwargs):
        dok_matrix.__init__(self, *args)
        self.format = "dok"
        self._model = kwargs["model"] if "model" in kwargs else None

    def __setitem__(self, index, value):
        dok_matrix.__setitem__(self, index, value)
        if type(index[0]) is int and type(index[1]) is int:
            reaction = self._model.reactions[index[1]]
            if value != 0:
                reaction.add_metabolites({self._model.metabolites[index[0]]: value}, combine=False)
            else:  # setting 0 means metabolites should be removed
                metabolite = self._model.metabolites[index[0]]
                if metabolite in reaction._metabolites:
                    reaction.pop(metabolite)

    def tolil(self):
        new = SMatrix_lil(dok_matrix.tolil(self), model=self._model)
        return new

class SMatrix_lil(lil_matrix):
    """A 2D sparse lil matrix which maintains links to a cobra Model"""
    def __init__(self, *args, **kwargs):
        lil_matrix.__init__(self, *args)
        self.format = "lil"
        self._model = kwargs["model"] if "model" in kwargs else None
    
    def __setitem__(self, index, value):
        lil_matrix.__setitem__(self, index, value)
        if isinstance(index[0], int):
            metabolites = [self._model.metabolites[index[0]]]
        else:
            metabolites = self._model.metabolites[index[0]]
        
        if isinstance(index[1], int):
            reactions = [self._model.reactions[index[1]]]
        else:
            reactions = self._model.reactions[index[1]]
        
        if value == 0:  # remove_metabolites
            met_set = set(metabolites)
            for reaction in reactions:
                to_remove = met_set.intersection(reaction._metabolites)
                for i in to_remove:
                    reaction.pop(i)
        else:  # add metabolites
            met_dict = {met: value for met in metabolites}
            for reaction in reactions:
                reaction.add_metabolites(met_dict, combine=False)

    def todok(self):
        new = SMatrix_dok(lil_matrix.todok(self), model=self._model)
        return new


SMatrix_classes = {"scipy.dok_matrix": SMatrix_dok,
                   "scipy.lil_matrix": SMatrix_lil}
