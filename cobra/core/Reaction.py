from __future__ import print_function

from collections import defaultdict
import re
from copy import copy, deepcopy
from warnings import warn

from six import string_types, iteritems

from .Object import Object
from .Gene import Gene, parse_gpr, ast2str
from .Metabolite import Metabolite


class Frozendict(dict):
    """Read-only dictionary view"""

    def __setitem__(self, key, value):
        raise NotImplementedError("read-only")

    def __delitem__(self, key):
        raise NotImplementedError("read-only")

    def pop(self, key, value):
        raise NotImplementedError("read-only")

    def popitem(self):
        raise NotImplementedError("read-only")


# precompiled regular expressions
# Matches and/or in a gene reaction rule
and_or_search = re.compile(r'\(| and| or|\+|\)', re.IGNORECASE)
gpr_clean = re.compile(' {2,}')
# This regular expression finds any single letter compartment enclosed in
# square brackets at the beginning of the string. For example [c] : foo --> bar
compartment_finder = re.compile("^\s*(\[[A-Za-z]\])\s*:*")
# Regular expressions to match the arrows
_reversible_arrow_finder = re.compile("<(-+|=+)>")
_forward_arrow_finder = re.compile("(-+|=+)>")
_reverse_arrow_finder = re.compile("<(-+|=+)")


class Reaction(Object):
    """Reaction is a class for holding information regarding
    a biochemical reaction in a cobra.Model object

    """

    def __init__(self, name=None):
        """An object for housing reactions and associated information
        for cobra modeling.

        """
        Object.__init__(self, name)
        self._gene_reaction_rule = ''
        self.subsystem = ''
        # The cobra.Genes that are used to catalyze the reaction
        self._genes = set()
        # A dictionary of metabolites and their stoichiometric coefficients in
        # this reaction.
        self._metabolites = {}
        self.name = name
        # self.model is None or refers to the cobra.Model that
        # contains self
        self._model = None

        self.objective_coefficient = self.lower_bound = 0.
        self.upper_bound = 1000.
        # Used during optimization.  Indicates whether the
        # variable is modeled as continuous, integer, binary, semicontinous, or
        # semiinteger.
        self.variable_kind = 'continuous'

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
        self._gene_reaction_rule = new_rule.strip()
        try:
            _, gene_names = parse_gpr(self._gene_reaction_rule)
        except (SyntaxError, TypeError) as e:
            warn("malformed gene_reaction_rule '%s' for %s" %
                 (new_rule, repr(self)))
            tmp_str = and_or_search.sub('', self._gene_reaction_rule)
            gene_names = set((gpr_clean.sub(' ', tmp_str).split(' ')))
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
    def gene_name_reaction_rule(self):
        """Display gene_reaction_rule with names intead.

        Do NOT use this string for computation. It is intended to give a
        representation of the rule using more familiar gene names instead of
        the often cryptic ids.

        """
        names = {i.id: i.name for i in self._genes}
        ast = parse_gpr(self._gene_reaction_rule)[0]
        return ast2str(ast, names=names)

    @property
    def x(self):
        """The flux through the reaction in the most recent solution

        Flux values are computed from the primal values of the variables in
        the solution.

        """
        try:
            return self._model.solution.x_dict[self.id]
        except Exception as e:
            if self._model is None:
                raise Exception("not part of a model")
            if not hasattr(self._model, "solution") or \
                    self._model.solution is None or \
                    self._model.solution.status == "NA":
                raise Exception("model has not been solved")
            if self._model.solution.status != "optimal":
                raise Exception("model solution was not optimal")
            raise e  # Not sure what the exact problem was

    @property
    def reversibility(self):
        """Whether the reaction can proceed in both directions (reversible)

        This is computed from the current upper and lower bounds.

        """
        return self.lower_bound < 0 and self.upper_bound > 0

    @reversibility.setter
    def reversibility(self, value):
        warn("Setting reaction reversibility is ignored")

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

    @property
    def model(self):
        """returns the model the reaction is a part of"""
        return self._model

    def _update_awareness(self):
        """Make sure all metabolites and genes that are associated with
        this reaction are aware of it.

        """
        for x in self._metabolites:
            x._reaction.add(self)
        for x in self._genes:
            x._reaction.add(self)

    def remove_from_model(self, model=None, remove_orphans=False):
        """Removes the reaction from the model while keeping it intact

        remove_orphans: Boolean
            Remove orphaned genes and metabolites from the model as well

        model: deprecated argument, must be None

        """
        if model is not None:
            warn("model does not need to be passed into remove_from_model")
            if model != self._model:
                raise Exception("Can not remove from a different model")
        if self._model is None:
            raise Exception("Reaction %s not in a model" % self.id)
        # preserve the original attributes (but as copies)
        model = self._model
        new_metabolites = {copy(met): value
                           for met, value in iteritems(self._metabolites)}
        new_genes = {copy(i) for i in self._genes}
        # Begin removing from the model
        self._model = None
        model.reactions.remove(self)
        for x in self._metabolites:
            x._reaction.remove(self)
            if remove_orphans and len(x._reaction) == 0:
                model.metabolites.remove(x)
        for x in self._genes:
            x._reaction.remove(self)
            if remove_orphans and len(x._reaction) == 0:
                model.genes.remove(x)
        # Rebuild the model with the new independent genes/metabolites
        self._metabolites = {}
        self.add_metabolites(new_metabolites)
        self._genes = set()
        for k in new_genes:
            self._associate_gene(k)

    def delete(self, remove_orphans=False):
        """Completely delete a reaction

        This removes all associations between a reaction the associated
        model, metabolites and genes (unlike remove_from_model which only
        dissociates the reaction from the model).

        remove_orphans: Boolean
            Remove orphaned genes and metabolites from the model as well

        """
        model = self._model
        if model is not None:
            self._model.reactions.remove(self)
        elif remove_orphans:
            # can't remove orphans if not part of a model
            remove_orphans = False
        self._model = None
        for x in self._metabolites:
            if self in x._reaction:
                x._reaction.remove(self)
                if remove_orphans and len(x._reaction) == 0:
                    model.metabolites.remove(x)
        for x in self._genes:
            if self in x._reaction:
                x._reaction.remove(self)
                if remove_orphans and len(x._reaction) == 0:
                    model.genes.remove(x)
        self._metabolites = {}
        self._genes = set()

    def __setstate__(self, state):

        """Probably not necessary to set _model as the cobra.Model that
        contains self sets the _model attribute for all metabolites and genes
        in the reaction.

        However, to increase performance speed we do want to let the metabolite
        and gene know that they are employed in this reaction

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
        # the_model = self._model
        # self._model = None
        new_reaction = deepcopy(self)
        # self._model = the_model
        return new_reaction

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
                raise KeyError(
                    "No metabolite named %s in the reaction" % the_metabolite)
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
        # Make all the genes aware of this reaction
        [x._reaction.add(new_reaction) for x in new_reaction._genes]
        gpr_1 = new_reaction.gene_reaction_rule
        gpr_2 = other_reaction.gene_reaction_rule
        if gpr_1 != '' and gpr_2 != '':
            new_reaction.gene_reaction_rule = '%s and %s' % (gpr_1, gpr_2)
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
        new_reaction.subtract_metabolites(
            deepcopy(other_reaction._metabolites))
        return new_reaction

    def __imul__(self, the_coefficient):
        """Allows the reaction coefficients to be rapidly scaled.

        """
        self._metabolites = {k: the_coefficient * v for k, v in
                             iteritems(self._metabolites)}
        return self

    def __mul__(self, the_coefficient):
        """Allows a reaction to be multiplied by a coefficient.

        TODO: this should return a new reaction.

        """
        self *= the_coefficient
        return self

    @property
    def reactants(self):
        """Return a list of reactants for the reaction."""
        return [k for k, v in self._metabolites.items() if v < 0]

    @property
    def products(self):
        """Return a list of products for the reaction"""
        return [k for k, v in self._metabolites.items() if v > 0]

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

    def add_metabolites(self, metabolites, combine=True,
                        add_to_container_model=True):
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
        _id_to_metabolites = {x.id: x for x in self._metabolites}
        new_metabolites = []
        for metabolite, coefficient in iteritems(metabolites):
            # If a metabolite already exists in the reaction then
            # just add them.
            if metabolite.id in _id_to_metabolites:
                reaction_metabolite = _id_to_metabolites[metabolite.id]
                if combine:
                    self._metabolites[reaction_metabolite] += coefficient
                else:
                    self._metabolites[reaction_metabolite] = coefficient
            else:
                # If the reaction is in a model, ensure we aren't using
                # a duplicate metabolite.
                try:
                    metabolite = \
                        self._model.metabolites.get_by_id(metabolite.id)
                except:
                    new_metabolites.append(metabolite)
                self._metabolites[metabolite] = coefficient
                # make the metabolite aware that it is involved in this
                # reaction
                metabolite._reaction.add(self)
        for metabolite, the_coefficient in list(self._metabolites.items()):
            if the_coefficient == 0:
                # make the metabolite aware that it no longer participates
                # in this reaction
                metabolite._reaction.remove(self)
                self._metabolites.pop(metabolite)
        if add_to_container_model and hasattr(self._model, 'add_metabolites'):
            self._model.add_metabolites(new_metabolites)

    def subtract_metabolites(self, metabolites):
        """This function will 'subtract' metabolites from a reaction, which
        means add the metabolites with -1*coefficient. If the final coefficient
        for a metabolite is 0 then the metabolite is removed from the reaction.

        metabolites: dict of {:class:`~cobra.core.Metabolite`: coefficient}
            These metabolites will be added to the reaction

        .. note:: A final coefficient < 0 implies a reactant.

        """
        self.add_metabolites({k: -v for k, v in iteritems(metabolites)})

    def clear_metabolites(self):
        """Remove all metabolites from the reaction"""
        for metabolite in list(self._metabolites.keys()):
            self.pop(metabolite)

    @property
    def reaction(self):
        """Human readable reaction string"""
        return self.build_reaction_string()

    @reaction.setter
    def reaction(self, value):
        return self.build_reaction_from_string(value)

    def build_reaction_string(self, use_metabolite_names=False):
        """Generate a human readable reaction string"""
        def format(number):
            return "" if number == 1 else str(number).rstrip(".") + " "
        reactant_dict = {}
        product_dict = {}
        id_type = 'id'
        if use_metabolite_names:
            id_type = 'name'
        reactant_bits = []
        product_bits = []
        for the_metabolite, coefficient in iteritems(self._metabolites):
            name = str(getattr(the_metabolite, id_type))
            if coefficient > 0:
                product_bits.append(format(coefficient) + name)
            else:
                reactant_bits.append(format(abs(coefficient)) + name)

        reaction_string = ' + '.join(reactant_bits)
        if not self.reversibility:
            if self.lower_bound < 0 and self.upper_bound <= 0:
                reaction_string += ' <-- '
            else:
                reaction_string += ' --> '
        else:
            reaction_string += ' <=> '
        reaction_string += ' + '.join(product_bits)
        return reaction_string

    def check_mass_balance(self):
        """Compute mass and charge balance for the reaction

        returns a dict of {element: amount} for unbalanced elements.
        "charge" is treated as an element in this dict
        This should be empty for balanced reactions.
        """
        reaction_element_dict = defaultdict(int)
        for metabolite, coefficient in self._metabolites.items():
            if metabolite.charge is not None:
                reaction_element_dict["charge"] += \
                    coefficient * metabolite.charge
            for element, amount in iteritems(metabolite.elements):
                reaction_element_dict[element] += coefficient * amount
        # filter out 0 values
        return {k: v for k, v in iteritems(reaction_element_dict) if v != 0}

    def get_compartments(self):
        """lists compartments the metabolites are in"""
        return list({x.compartment for x in self._metabolites})

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

    def build_reaction_from_string(self, reaction_str, verbose=True,
                                   fwd_arrow=None, rev_arrow=None,
                                   reversible_arrow=None, term_split="+"):
        # set the arrows
        forward_arrow_finder = _forward_arrow_finder if fwd_arrow is None \
            else re.compile(re.escape(fwd_arrow))
        reverse_arrow_finder = _reverse_arrow_finder if rev_arrow is None \
            else re.compile(re.escape(rev_arrow))
        reversible_arrow_finder = _reversible_arrow_finder \
            if reversible_arrow is None \
            else re.compile(re.escape(reversible_arrow))
        if self._model is None:
            warn("no model found")
            model = None
        else:
            model = self._model
        original_str = "" + reaction_str  # copy
        found_compartments = compartment_finder.findall(reaction_str)
        if len(found_compartments) == 1:
            compartment = found_compartments[0]
            reaction_str = compartment_finder.sub("", reaction_str)
        else:
            compartment = ""

        # reversible case
        arrow_match = reversible_arrow_finder.search(reaction_str)
        if arrow_match is not None:
            self.lower_bound = -1000
            self.upper_bound = 1000
        else:  # irreversible
            # try forward
            arrow_match = forward_arrow_finder.search(reaction_str)
            if arrow_match is not None:
                self.upper_bound = 1000
                self.lower_bound = 0
            else:
                # must be reverse
                arrow_match = reverse_arrow_finder.search(reaction_str)
                if arrow_match is None:
                    raise ValueError("no suitable arrow found in '%s'" %
                                     reaction_str)
                else:
                    self.upper_bound = 0
                    self.lower_bound = -1000
        reactant_str = reaction_str[:arrow_match.start()].strip()
        product_str = reaction_str[arrow_match.end():].strip()

        self.clear_metabolites()

        for substr, factor in ((reactant_str, -1), (product_str, 1)):
            if len(substr) == 0:
                continue
            for term in substr.split(term_split):
                term = term.strip()
                if term.lower() == "nothing":
                    continue
                if " " in term:
                    num_str, met_id = term.split()
                    num = float(num_str.lstrip("(").rstrip(")")) * factor
                else:
                    met_id = term
                    num = factor
                met_id += compartment
                try:
                    met = model.metabolites.get_by_id(met_id)
                except KeyError:
                    if verbose:
                        print("unknown metabolite '%s' created" % met_id)
                    met = Metabolite(met_id)
                self.add_metabolites({met: num})
