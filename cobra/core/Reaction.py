from __future__ import print_function

from ..external.six import string_types, iteritems


#Is it better to restrict a Reaction to a single model or
#should we allow a Reaction to be associated with multiple models?
#
from collections import defaultdict
import re
from copy import deepcopy
from .Object import Object
from .Metabolite import Metabolite
from .Gene import Gene

from warnings import warn

class Frozendict(dict):
    def __setitem__(self, key, value):
        raise NotImplementedError("read-only")

    def __delitem__(self, key):
        raise NotImplementedError("read-only")

    def pop(self, key, value):
        raise NotImplementedError("read-only")

    def popitem(self):
        raise NotImplementedError("read-only")


and_or_search = re.compile('\(| and| or|\+|\)', re.IGNORECASE)

class Reaction(Object):
    """Reaction is a class for holding information regarding
    a biochemical reaction in a cobra.Model object 

    """
    ## __slots__ = ['id', '_metabolites', '_gene_reaction_rule',
    ##              'subsystem', '_genes', '_model',
    ##              'name', 'lower_bound', 'upper_bound', 'objective_coefficient',
    ##              ]

    def __init__(self, name=None):
        """An object for housing reactions and associated information
        for cobra modeling..

        boundary: Boolean.  Indicates whether the reaction is at a boundary, i.e.,
        a source or a sink.
        
        """
        Object.__init__(self, name)
        self._gene_reaction_rule = ''
        self.subsystem = ''
        self._genes = set() #The cobra.Genes that are used to catalyze the reaction
        #reaction.  _ Indicates that it is not preferred to add a gene to a reaction
        #directly.
        #A dictionary of metabolites and their stoichiometric coefficients in
        #this reaction  _ Indicates that it is not preferred to add a metabolite to a reaction
        #directly.
        self._metabolites = {}
        self.name = name
        #self.model is None or refers to the cobra.Model that
        #contains self
        self._model = None

        self.objective_coefficient = self.lower_bound = 0.
        self.upper_bound = 1000.
        self.reflection = None #Either None or if this reaction is irreversible then
        #a reaction in the model that is essentially self * -1
        self.variable_kind = 'continuous' #Used during optimization.  Indicates whether the
        #variable is modeled as continuous, integer, binary, semicontinous, or semiinteger.

    # read-only
    @property
    def metabolites(self):
        return Frozendict(self._metabolites)

    @property
    def genes(self):
        return frozenset(self._genes)

    @property
    def gene_reaction_rule(self):
        return self._gene_reaction_rule

    @gene_reaction_rule.setter
    def gene_reaction_rule(self, new_rule):
        self._gene_reaction_rule = new_rule
        gene_names = set((re.compile(' {2,}').sub(' ', and_or_search.sub('', self._gene_reaction_rule))).split(' ' ))
        if '' in gene_names:
                gene_names.remove('')
        old_genes = self._genes
        if self._model is None:
            self._genes = {Gene(i) for i in gene_names}
        else:
            model_genes = self._model.genes
            self._genes = set()
            for id in gene_names:
                if model_genes.has_id(id):
                    self._genes.add(model_genes.get_by_id(id))
                else:
                    new_gene = Gene(id)
                    new_gene._model = self._model
                    self._genes.add(new_gene)
                    model_genes.append(new_gene)

        # Make the genes aware that it is involved in this reaction
        for g in self._genes:
            g._reaction.add(self)

        # make the old genes aware they are no longer involved in this reaction
        for g in old_genes:
            if g not in self._genes:  # if an old gene is not a new gene
                try:
                    g._reaction.remove(self)
                except:
                    warn("could not remove old gene %s from reaction %s" %
                         (g.id, self.id))



    @property
    def reversibility(self):
        """This property removes the independence of the reversibility attribute and the reaction's
        current upper and lower bounds.

        reversibility is defined in the context of the current instantiation.
        
        """
        return self.lower_bound < 0 and self.upper_bound > 0
    
    @property
    def boundary(self):
        # single metabolite implies it must be a boundary
        if len(self._metabolites) == 1:
            return "system_boundary"
        # if there is more than one metabolite, if it ONLY produces or ONLY
        # consumes, it is also a boundary.
        all_stoichiometry = self._metabolites.values()
        if not min(all_stoichiometry) < 0 < max(all_stoichiometry):
            return "system_boundary"
        return None

    def _update_awareness(self):
        """Make sure all metabolites and genes that are associated with
        this reaction are aware of it.
        
        """
        [x._reaction.add(self) for x in self._metabolites]
        [x._reaction.add(self) for x in self._genes]


    def get_model(self):
        """Returns the Model object that this Reaction is associated with.

        """
        return self._model

    def remove_from_model(self, model=None):
        """Removes the association

        model: deprecated argument, should be None
        
        """
        # why is model being taken in as a parameter? This plays
        #back to the question of allowing a Metabolite to be associated
        #with multiple Models
        if model != self._model and model is not None:
            raise Exception('%s not in %s ergo it cannot be removed. (%s)'%(self,
                                                                  model,
                                                                  self._model))
        if self._model is None:
            raise Exception("Reaction %s not in a model" % self.id)
        if model is not None:
            warn("model does not need to be passed into remove_from_model")
        new_metabolites = deepcopy(self._metabolites)
        new_genes = deepcopy(self._genes)
        self._model.reactions.remove(self)
        #Remove associations between the reaction and its container _model
        #and elements in the model
        self._model = None
        [x._reaction.remove(self)
         for x in self._metabolites.keys()]
        [x._reaction.remove(self)
         for x in self._genes]
        #Replace the model-linked metabolites with the new independent metabolites
        self._metabolites = {}
        self.add_metabolites(new_metabolites)
        #Replace the model-linked genes with new indepenent genes
        self._genes = set()
        for k in new_genes:
            self._associate_gene(k)

    def delete(self):
        """Removes all associations between a reaction and its container
        _model and metabolites and genes.
        
        TODO: Decide whether orphan metabolites should automatically be removed
        from the model.
        
        """
        self._model = None
        [x._reaction.remove(self)
         for x in self._metabolites.keys() if self in x._reaction]
        [x._reaction.remove(self)
         for x in self._genes if self in x._reaction]
        self._metabolites = {}
        self._genes = set()
        
        
    def __setstate__(self, state):
        """Probably not necessary to set _model as the cobra.Model that
        contains self sets the _model attribute for all metabolites and genes in the reaction.

        However, to increase performance speed we do want to let the metabolite and gene
        know that they are employed in this reaction

        """
        # These are necessary for old pickles which store attributes
        # which have since been superceded by properties.
        if "reaction" in state:
            state.pop("reaction")
        if "gene_reaction_rule" in state:
            state["_gene_reaction_rule"] = state.pop("gene_reaction_rule")

        self.__dict__.update(state)
        for x in state['_metabolites']:
            setattr(x, '_model', self._model)
            x._reaction.add(self)
        for x in state['_genes']:
            setattr(x, '_model', self._model)
            x._reaction.add(self)

    
    def copy(self):
        """When copying a reaction, it is necessary to deepcopy the
        components so the list references aren't carried over.

        Additionally, a copy of a reaction is no longer in a cobra.Model.

        This should be fixed with self.__deecopy__ if possible
        """
        ## the_model = self._model
        ## self._model = None
        new_reaction = deepcopy(self)
        ## self._model = the_model
        return new_reaction

    def guided_copy(self, the_model, metabolite_dict, gene_dict=None):
        """.. deprecated :: 0.3 use copy directly"""
        warn("deprecated")
        the_copy = Object.guided_copy(self)
        #Replace the complex items in a faster fashion
        the_copy._model = the_model
        if gene_dict:
            the_copy._genes = set([gene_dict[k.id]
                                for k in self._genes])
        the_copy._metabolites = {metabolite_dict[k.id]: v
                                 for k, v in iteritems(self._metabolites)}

        #make the metabolites and genes aware of the reaction
        [k._reaction.add(the_copy)
         for k in the_copy._genes]
        [k._reaction.add(the_copy)
         for k in the_copy._metabolites.keys()]

        return(the_copy)

    def pop(self, metabolite_id):
        """Remove a metabolite from the reaction and return the
        stoichiometric coefficient.

        metabolite_id: str or :class:`~cobra.core.Metabolite.Metabolite`

        """
        the_metabolite = metabolite_id
        if isinstance(the_metabolite, string_types):
            found_match = None
            for possible_match in self._metabolites:
                if possible_match.id == the_metabolite:
                    found_match = possible_match
                    break
            if found_match is None:
                raise KeyError("No metabolite named %s in the reaction" % the_metabolite)
            else:
                the_metabolite = found_match
        the_coefficient = self._metabolites.pop(the_metabolite)
        the_metabolite._reaction.remove(self)
        return the_coefficient
    
    def __add__(self, other_reaction):
        """Adds two reactions to each other.  Default behavior is
        to combine the metabolites but only use the remaining parameters
        from the first object.
        
        TODO: Either clean up metabolite associations or remove function

        TODO: Deal with gene association logic from adding reactions.

        TODO: Simplify and add in an __iadd__
        
        """
        new_reaction = deepcopy(self)
        new_reaction.id = self.id + '_' + other_reaction.id
        new_reaction.add_metabolites(deepcopy(other_reaction._metabolites))
        new_reaction._genes.update(deepcopy(other_reaction._genes))
        #Make all the genes aware of this reaction
        [x._reaction.add(new_reaction) for x in new_reaction._genes]
        gpr_1 = new_reaction.gene_reaction_rule
        gpr_2 = other_reaction.gene_reaction_rule
        if gpr_1 != '' and gpr_2 != '':
            new_reaction.gene_reaction_rule = '%s and %s'%(gpr_1, gpr_2)
        elif gpr_2 != '':
            new_reaction.gene_reaction_rule = gpr_2
        return new_reaction
    
    def __sub__(self, other_reaction):
        """Subtracts two reactions.  Default behavior is
        to combine the metabolites but only use the remaining parameters
        from the first object.

        Note: This is equivalent to adding reactions after changing the sign
        of the metabolites in other_reaction

        """
        new_reaction = deepcopy(self)
        if self is other_reaction:
            other_reaction = deepcopy(other_reaction)
        new_reaction.id = self.id + '_' + other_reaction.id
        new_reaction.subtract_metabolites(deepcopy(other_reaction._metabolites))
        return new_reaction

    def __imul__(self, the_coefficient):
        """Allows the reaction coefficients to be rapidly scaled.
        
        """
        [self._metabolites.update({k: the_coefficient * v})
         for k, v in self._metabolites.items()]
        return self

    def __mul__(self, the_coefficient):
        """Allows a reaction to be multiplied by a coefficient.
        
        TODO: this should return a new reaction.
        
        """
        [self._metabolites.update({k: the_coefficient * v})
         for k, v in self._metabolites.items()]
        return self
        

    def parse_gene_association(self, the_type='gene'):
        """.. deprecated :: 0.3 Set gene_reaction_rule directly"""
        warn("deprecated function")
        # trigger the update if that was the desired behavior for some reason
        self._gene_reaction_rule = self._gene_reaction_rule


    def add_gene_reaction_rule(self, the_rule):
        """.. deprecated :: 0.3 Set gene_reaction_rule directly"""
        self.gene_reaction_rule = the_rule
        warn("deprecated, assign to gene_reaction_rule directly")

    def get_reactants(self):
        """.. deprecated :: 0.3 use reactants property instead"""
        warn("deprecated, use the reactants property instead")
        return self.reactants


    @property
    def reactants(self):
        """Return a list of reactants for the reaction."""
        return [k for k, v in self._metabolites.items() if v < 0]

    def get_products(self):
        """.. deprecated :: 0.3 use products property instead"""
        warn("depreacated, use the products property instead")
        return self.products


    @property
    def products(self):
        """Return a list of products for the reaction"""
        return [k for k, v in self._metabolites.items() if v > 0]

    def get_gene(self):
        """.. deprecated :: 0.3 use genes property instead"""
        warn("deprecated, use the genes property instead")
        return list(self._genes)


    def get_coefficient(self, metabolite_id):
        """Return the stoichiometric coefficient for a metabolite in
        the reaction.

        metabolite_id: str or :class:`~cobra.core.Metabolite.Metabolite`
        
        """
        _id_to_metabolites = dict([(x.id, x)
                                        for x in self._metabolites])

        if hasattr(metabolite_id, 'id'):
            metabolite_id = metabolite_id.id
        return self._metabolites[_id_to_metabolites[metabolite_id]]
    
    def get_coefficients(self, metabolite_ids):
        """Return the stoichiometric coefficients for a list of
        metabolites in the reaction.

        metabolite_ids: iterable
            Containing str or :class:`~cobra.core.Metabolite.Metabolite`
        
        """
        return map(self.get_coefficient, metabolite_ids)
    
    def add_metabolites(self, metabolites, combine=True, add_to_container_model=True):
        """Add metabolites and stoichiometric coefficients to the reaction.
        If the final coefficient for a metabolite is 0 then it is removed
        from the reaction.

        metabolites: dict
            {:class:`~cobra.core.Metabolite.Metabolite`: coefficient}

        combine: Boolean.
            Describes behavior a metabolite already exists in the reaction.
            True causes the coefficients to be added.
            False causes the coefficient to be replaced.
            True and a metabolite already exists in the

        add_to_container_model: Boolean.
            Add the metabolite to the :class:`~cobra.core.Model.Model`
            the reaction is associated with (i.e. self.model)

        """
        the_metabolites = metabolites
        _id_to_metabolites = dict([(x.id, x)
                                        for x in self._metabolites])
        new_metabolites = []
        for the_metabolite, the_coefficient in the_metabolites.items():
            #If a metabolite already exists in the reaction then
            #just add them.
            if the_metabolite.id in _id_to_metabolites:
                reaction_metabolite = _id_to_metabolites[the_metabolite.id]
                if combine:
                    self._metabolites[reaction_metabolite] += the_coefficient
                else:
                    self._metabolites[reaction_metabolite] = the_coefficient
            else:
                self._metabolites[the_metabolite] = the_coefficient
                #make the metabolite aware that it is involved in this reaction
                the_metabolite._reaction.add(self)
                new_metabolites.append(the_metabolite)
        for the_metabolite, the_coefficient in self._metabolites.items():
            if the_coefficient == 0:
                #make the metabolite aware that it no longer participates
                #in this reaction
                the_metabolite._reaction.remove(self)
                self._metabolites.pop(the_metabolite)
        _id_to_metabolites = dict([(x.id, x)
                                        for x in self._metabolites])
        if add_to_container_model and hasattr(self._model, 'add_metabolites'):
            self._model.add_metabolites(new_metabolites)
            

    def subtract_metabolites(self, metabolites):
        """This function will 'subtract' cobra.metabolites from a reaction, which
        means add the metabolites with -1*coefficient.  If the final coefficient
        for a metabolite is 0 then the metabolite is removed from the reaction.

        metabolites: dict of {:class:`~cobra.core.Metabolite`: coefficient}
            These metabolites will be added to the reaction

        .. note:: A final coefficient < 0 implies a reactant.

        .. note:: This function uses deepcopy in case the reaction is being
                  subtracted from itself.
        
        """
        metabolites = deepcopy(metabolites)
        metabolites = dict([(k, -v)
                                for k, v in metabolites.items()])
        self.add_metabolites(metabolites)

    @property
    def reaction(self):
        """Human readable reaction string"""
        return self.build_reaction_string()


    def build_reaction_string(self, use_metabolite_names=False):
        """Generate a human readable reaction string"""
        def format(number):
            if number == 1:
                return ""
            if number == int(number):
                return str(int(number)) + " "
            return str(number) + " "
        reactant_dict = {}
        product_dict = {}
        id_type = 'id'
        if use_metabolite_names:
            id_type = 'name'
        reactant_bits = []
        product_bits = []
        for the_metabolite, coefficient in self._metabolites.items():
            name = str(getattr(the_metabolite, id_type))
            if coefficient > 0:
                product_bits.append(format(coefficient) + name)
            else:
                reactant_bits.append(format(abs(coefficient)) + name)

        reaction_string = ' + '.join(reactant_bits)
        if not self.reversibility:
            if self.lower_bound < 0 and self.upper_bound <=0:
                reaction_string += ' <-- '
            else:
                reaction_string += ' --> ' 
        else:
            reaction_string += ' <=> '
        reaction_string += ' + '.join(product_bits)
        return reaction_string


    def check_mass_balance(self):
        """Makes sure that the reaction is elementally-balanced."""
        reaction_element_dict = defaultdict(list)
        for the_metabolite, the_coefficient in self._metabolites.items():
            if the_metabolite.formula is not None:
                [reaction_element_dict[k].append(the_coefficient*v)
                 for k, v in the_metabolite.formula.elements.items()]
        reaction_element_dict = dict([(k, sum(v))
                                      for k, v in reaction_element_dict.items()])
        if sum(map(abs, reaction_element_dict.values())) != 0:
            return [self.id, reaction_element_dict]
        else:
            return []
        
    def print_values(self):
        """.. deprecated :: 0.3"""
        warn("deprecated")
        print("reaction:", self.id)
        print("subsystem", self.subsystem)
        print(self.reaction)
        print("bounds: (%.2f, %.2f)" % (self.lower_bound, self.upper_bound))
        print("objective_coefficient", self.objective_coefficient)
        print("gene reaction rule:", self.gene_reaction_rule)


    def get_compartments(self):
        """
        """
        return(list(set([x.compartment for x in self._metabolites])))


    def remove_gene(self, cobra_gene):
        """.. deprecated :: 0.3 update the gene_reaction_rule instead"""
        warn("deprecated: update the gene_reaction_rule instead")
        try:
            self._genes.remove(cobra_gene)
            cobra_gene._reaction.remove(self)
        except Exception as e:
            try:
                if hasattr(self._genes, 'keys'):
                    self._genes = set(self._genes.keys())
                    self.remove_gene(cobra_gene)
            except:
                raise Exception('Unable to remove gene %s from reaction %s: %s'%(cobra_gene.id, self.id, e))

    def add_gene(self, cobra_gene):
        """.. deprecated :: 0.3 update the gene_reaction_rule instead"""
        warn("deprecated: update the gene_reaction_rule instead")
        try:
            self._genes.add(cobra_gene)
            cobra_gene._reaction.add(self)
            cobra_gene._model = self._model
        except Exception as e:
            try:
                if hasattr(self._genes, 'keys'):
                    self._genes = set(self._genes.keys())
                    self.add_gene(cobra_gene)
            except:
                raise Exception('Unable to add gene %s to reaction %s: %s'%(cobra_gene.id, self.id, e))

    def _associate_gene(self, cobra_gene):
        """Associates a cobra.Gene object with a cobra.Reaction.

        cobra_gene : :class:`~cobra.core.Gene.Gene`

        """
        self._genes.add(cobra_gene)
        cobra_gene._reaction.add(self)
        cobra_gene._model = self._model

    def _dissociate_gene(self, cobra_gene):
        """Dissociates a cobra.Gene object with a cobra.Reaction.

        cobra_gene : :class:`~cobra.core.Gene.Gene`

        """
        self._genes.remove(cobra_gene)
        cobra_gene._reaction.remove(self)
                
    def knock_out(self):
        """Change the upper and lower bounds of the reaction to 0."""
        self.lower_bound = 0
        self.upper_bound = 0
