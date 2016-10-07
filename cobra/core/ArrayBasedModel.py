from sys import maxsize
from warnings import warn
from six import iteritems

from numpy import array, ndarray
from scipy.sparse import lil_matrix, dok_matrix

from .Model import Model


class ArrayBasedModel(Model):
    """ArrayBasedModel is a class that adds arrays and vectors to
    a cobra.Model to make it easier to perform linear algebra operations.

    """

    def __init__(self, description=None, deepcopy_model=False,
                 matrix_type='scipy.lil_matrix'):
        """
        description: None | String | cobra.Model

        deepcopy_model: Boolean.  If True and description is
        a cobra.Model then make a deepcopy of the Model before
        creating the ArrayBasedModel.

        matrix_type: 'scipy.lil_matrix' or 'scipy.dok_matrix'
            Specifies which type of backend matrix to use for S.

        """
        if deepcopy_model and isinstance(description, Model):
            description = description.copy()
        Model.__init__(self, description)
        self._S = None
        self.matrix_type = matrix_type
        self.update()

    # no setter for S at the moment
    @property
    def S(self):
        """Stoichiometric matrix of the model

        This will be formatted as either :class:`~scipy.sparse.lil_matrix`
        or :class:`~scipy.sparse.dok_matrix`

        """
        return self._S

    @property
    def lower_bounds(self):
        return self._lower_bounds

    @lower_bounds.setter
    def lower_bounds(self, vector):
        self._update_from_vector("lower_bounds", vector)

    @property
    def upper_bounds(self):
        return self._upper_bounds

    @upper_bounds.setter
    def upper_bounds(self, vector):
        self._update_from_vector("upper_bounds", vector)

    @property
    def objective_coefficients(self):
        return self._objective_coefficients

    @objective_coefficients.setter
    def objective_coefficients(self, vector):
        self._update_from_vector("objective_coefficients", vector)

    @property
    def b(self):
        """bounds for metabolites as :class:`numpy.ndarray`"""
        return self._b

    @b.setter
    def b(self, vector):
        self._update_from_vector("b", vector)

    @property
    def constraint_sense(self):
        return self._constraint_sense

    @constraint_sense.setter
    def constraint_sense(self, vector):
        self._update_from_vector("_constraint_sense", vector)

    def copy(self):
        """Provides a partial 'deepcopy' of the Model.  All of the Metabolite,
        Gene, and Reaction objects are created anew but in a faster fashion
        than deepcopy

        """
        the_copy = Model.copy(self)
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
        Model.add_metabolites(self, metabolite_list)
        if self._S is not None and expand_stoichiometric_matrix:
            s_expansion = len(self.metabolites) - self._S.shape[0]
            if s_expansion > 0:
                self._S.resize((self._S.shape[0] + s_expansion,
                                self._S.shape[1]))
        self._update_metabolite_vectors()

    def _update_from_vector(self, attribute, vector):
        """convert from model.reactions = v to model.reactions[:] = v"""
        # this will fail if vector is the wrong length
        getattr(self, attribute)[:] = vector

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
                warn(the_reaction.id + ' is not in the model')
                continue

            # zero reaction stoichiometry column
            the_column = self._S[:, reaction_index]
            for nonzero_index in the_column.nonzero()[0]:
                the_column[nonzero_index, 0] = 0
            self._lower_bounds[reaction_index] = the_reaction.lower_bound
            self._upper_bounds[reaction_index] = the_reaction.upper_bound
            self.objective_coefficients[
                reaction_index] = the_reaction.objective_coefficient
            self.add_metabolites(the_reaction._metabolites)
            # Make sure that the metabolites are the ones contained in the
            # model
            the_reaction._metabolites = [self.metabolite.get_by_id(x.id)
                                         for x in the_reaction._metabolites]
            # Update the stoichiometric matrix
            metabolite_indices = map(
                self.metabolites.index,
                the_reaction._metabolites)
            for (index, metabolite_index) in enumerate(metabolite_indices):
                self._S[metabolite_index, reaction_index] = \
                    the_reaction.stoichiometric_coefficients[index]

    def add_reactions(self, reaction_list, update_matrices=True):
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
        Model.add_reactions(self, reaction_list)
        if update_matrices:
            self._update_matrices(reaction_list)

    def remove_reactions(self, reactions, update_matrices=True, **kwargs):
        """remove reactions from the model

        See :func:`cobra.core.Model.Model.remove_reactions`

        update_matrices:  Boolean
            If true populate / update matrices S, lower_bounds, upper_bounds.
            Note that this is slow to run for very large models, and using this
            option with repeated calls will degrade performance.

        """
        Model.remove_reactions(self, reactions, **kwargs)
        if update_matrices:
            self._update_matrices()

    def _construct_matrices(self):
        """Large sparse matrices take time to construct and to read / write.
        This function allows one to let the model exists without cobra_model.S
        and then generate it at needed.

        """
        self._update_matrices()  # This does basic construction as well.

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

    def _update_metabolite_vectors(self):
        """regenerates _b and _constraint_sense

        WARNING: This function is only used after the Model has been
        converted to matrices.  It is typically faster to access the objects
        in the Model directly.  This function will eventually moved to another
        module for advanced users due to the potential for mistakes.

        """
        self._b = LinkedArray(self.metabolites, "_bound")
        self._constraint_sense = LinkedArray(
            self.metabolites,
            "_constraint_sense")

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
        # no need to create matrix if there are no reactions or metabolites
        if len(self.reactions) == 0 and len(self.metabolites) == 0:
            return
        elif len(self.metabolites) == 0:
            self._update_reaction_vectors()
            return
        elif len(self.reactions) == 0:
            self._update_metabolite_vectors()
            return
        # Pretty much all of these things are unnecessary to use the objects
        # and interact with the optimization solvers.  It might be best to move
        # them to linear algebra modules.  If no reactions are present in the
        # Model, initialize the arrays
        if self._S is None or reaction_list is None:
            reaction_list = self.reactions
            SMatrix = SMatrix_classes[self.matrix_type]
            self._S = SMatrix((len(self.metabolites),
                               len(self.reactions)), model=self)
            self._update_reaction_vectors()
        else:  # Expand the arrays to accomodate the new reaction
            self._S.resize((len(self.metabolites),
                            len(self.reactions)))
            lower_bounds = array([x.lower_bound
                                  for x in reaction_list])
            upper_bounds = array([x.upper_bound
                                  for x in reaction_list])
            objective_coefficients = array([x.objective_coefficient
                                            for x in reaction_list])
            self._lower_bounds._extend(lower_bounds)
            self._upper_bounds._extend(upper_bounds)
            self._objective_coefficients._extend(objective_coefficients)

        coefficient_dictionary = {}
        for the_reaction in reaction_list:
            reaction_index = self.reactions.index(the_reaction.id)
            for the_key, the_value in the_reaction._metabolites.items():
                coefficient_dictionary[(self.metabolites.index(the_key.id),
                                        reaction_index)] = the_value

        self._S.update(coefficient_dictionary)

    def update(self):
        """Regenerates the stoichiometric matrix and vectors"""
        self._update_matrices()
        self._update_metabolite_vectors()


class LinkedArray(ndarray):
    """A :class:`numpy.ndarray` which updates an attribute from a list"""
    def __new__(cls, list, attribute):
        # construct a new ndarray with the values from the list
        # For example, if the list if model.reactions and the attribute is
        # "lower_bound" create an array of [reaction.lower_bound for ... ]
        x = array([getattr(i, attribute) for i in list]).view(cls)
        return x.copy()

    def __init__(self, list, attribute):
        self._list = list
        self._attr = attribute

    def __setitem__(self, index, value):
        ndarray.__setitem__(self, index, value)
        if isinstance(index, slice):
            # not sure why that is here
            if index.stop == maxsize:
                index = slice(index.start, len(self))
            if hasattr(value, "__getitem__"):
                for i, entry in enumerate(self._list[index]):
                    setattr(entry, self._attr, value[i])
            else:
                for i, entry in enumerate(self._list[index]):
                    setattr(entry, self._attr, value)
        else:
            setattr(self._list[index], self._attr, value)

    def __setslice__(self, i, j, value):
        self.__setitem__(slice(i, j), value)

    def _extend(self, other):
        old_size = len(self)
        new_size = old_size + len(other)
        self.resize(new_size, refcheck=False)
        ndarray.__setitem__(self, slice(old_size, new_size), other)


class SMatrix_dok(dok_matrix):
    """A 2D sparse dok matrix which maintains links to a cobra Model"""

    def __init__(self, *args, **kwargs):
        dok_matrix.__init__(self, *args)
        self.format = "dok"
        self._model = kwargs["model"] if "model" in kwargs else None

    def __setitem__(self, index, value):
        dok_matrix.__setitem__(self, index, value)
        if isinstance(index[0], int) and isinstance(index[1], int):
            reaction = self._model.reactions[index[1]]
            if value != 0:
                reaction.add_metabolites(
                    {self._model.metabolites[index[0]]: value}, combine=False)
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

    def update(self, value_dict):
        """update matrix without propagating to model"""
        if len(value_dict) < 100:  # TODO benchmark for heuristic
            for index, value in iteritems(value_dict):
                lil_matrix.__setitem__(self, index, value)
        else:
            matrix = lil_matrix.todok(self)
            matrix.update(value_dict)
            self = SMatrix_lil(matrix.tolil(), model=self._model)
            self._model._S = self

    def todok(self):
        new = SMatrix_dok(lil_matrix.todok(self), model=self._model)
        return new

    # TODO: check if implemented before using own function
    def resize(self, shape):
        matrix = lil_matrix.todok(self)
        matrix.resize(shape)
        self = SMatrix_lil(matrix.tolil(), model=self._model)
        self._model._S = self


SMatrix_classes = {"scipy.dok_matrix": SMatrix_dok,
                   "scipy.lil_matrix": SMatrix_lil}
