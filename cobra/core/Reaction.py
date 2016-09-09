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


def _is_positive(n):
    try:
        if n >= 0:
            return True
        else:
            return False
    except:
        return True

# precompiled regular expressions
# Matches and/or in a gene reaction rule
and_or_search = re.compile(r'\(| and| or|\+|\)', re.IGNORECASE)
uppercase_AND = re.compile(r'\bAND\b')
uppercase_OR = re.compile(r'\bOR\b')
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

    def __init__(self, id=None, name='', subsystem='', lower_bound=0.,
                 upper_bound=1000., objective_coefficient=0.):
        """An object for housing reactions and associated information
        for cobra modeling.

        """
        Object.__init__(self, id, name)
        self._gene_reaction_rule = ''
        self.subsystem = subsystem
        # The cobra.Genes that are used to catalyze the reaction
        self._genes = set()
        # A dictionary of metabolites and their stoichiometric coefficients in
        # this reaction.
        self._metabolites = {}
        # self.model is None or refers to the cobra.Model that
        # contains self
        self._model = None

        self.objective_coefficient = objective_coefficient
        self.upper_bound = upper_bound
        self.lower_bound = lower_bound
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
            if "AND" in new_rule or "OR" in new_rule:
                warn("uppercase AND/OR found in rule '%s' for '%s'" %
                     (new_rule, repr(self)))
                new_rule = uppercase_AND.sub("and", new_rule)
                new_rule = uppercase_OR.sub("or", new_rule)
                self.gene_reaction_rule = new_rule
                return
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
    def bounds(self):
        """ A more convienient bounds structure than seperate upper and lower
        bounds """

        return (self.lower_bound, self.upper_bound)

    @bounds.setter
    def bounds(self, value):
        """ Set the bounds directly from a tuple """

        self.lower_bound = value[0]
        self.upper_bound = value[1]

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
        """Copy a reaction

        The referenced metabolites and genes are also copied.

        """
        # no references to model when copying
        model = self._model
        self._model = None
        for i in self._metabolites:
            i._model = None
        for i in self._genes:
            i._model = None
        # now we can copy
        new_reaction = deepcopy(self)
        # restore the references
        self._model = model
        for i in self._metabolites:
            i._model = model
        for i in self._genes:
            i._model = model
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

    def __add__(self, other):
        """Add two reactions

        The stoichiometry will be the combined stoichiometry of the two
        reactions, and the gene reaction rule will be both rules combined by an
        and. All other attributes (i.e. reaction bounds) will match those of
        the first reaction

        """
        new_reaction = self.copy()
        new_reaction += other
        return new_reaction

    def __iadd__(self, other):
        self.add_metabolites(other._metabolites, combine=True)
        gpr1 = self.gene_reaction_rule.strip()
        gpr2 = other.gene_reaction_rule.strip()
        if gpr1 != '' and gpr2 != '':
            self.gene_reaction_rule = "(%s) and (%s)" % \
                (self.gene_reaction_rule, other.gene_reaction_rule)
        elif gpr1 != '' and gpr2 == '':
            self.gene_reaction_rule = gpr1
        elif gpr1 == '' and gpr2 != '':
            self.gene_reaction_rule = gpr2
        return self

    def __sub__(self, other):
        new = self.copy()
        new -= other
        return new

    def __isub__(self, other):
        self.subtract_metabolites(other._metabolites, combine=True)
        return self

    def __imul__(self, coefficient):
        """Scale coefficients in a reaction"""
        self._metabolites = {k: coefficient * v for k, v in
                             iteritems(self._metabolites)}
        return self

    def __mul__(self, coefficient):
        new = self.copy()
        new *= coefficient
        return new

    @property
    def reactants(self):
        """Return a list of reactants for the reaction."""
        return [k for k, v in self._metabolites.items() if not _is_positive(v)]

    @property
    def products(self):
        """Return a list of products for the reaction"""
        return [k for k, v in self._metabolites.items() if _is_positive(v)]

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
            {str or :class:`~cobra.core.Metabolite.Metabolite`: coefficient}

        combine: Boolean.
            Describes behavior a metabolite already exists in the reaction.
            True causes the coefficients to be added.
            False causes the coefficient to be replaced.
            True and a metabolite already exists in the

        add_to_container_model: Boolean.
            Add the metabolite to the :class:`~cobra.core.Model.Model`
            the reaction is associated with (i.e. self.model)

        """
        _id_to_metabolites = {str(x): x for x in self._metabolites}
        new_metabolites = []
        for metabolite, coefficient in iteritems(metabolites):
            met_id = str(metabolite)
            # If a metabolite already exists in the reaction then
            # just add them.
            if met_id in _id_to_metabolites:
                reaction_metabolite = _id_to_metabolites[met_id]
                if combine:
                    self._metabolites[reaction_metabolite] += coefficient
                else:
                    self._metabolites[reaction_metabolite] = coefficient
            else:
                # If the reaction is in a model, ensure we aren't using
                # a duplicate metabolite.
                if self._model:
                    try:
                        metabolite = \
                            self._model.metabolites.get_by_id(met_id)
                    except KeyError as e:
                        if isinstance(metabolite, Metabolite):
                            new_metabolites.append(metabolite)
                        else:
                            # do we want to handle creation here?
                            raise e
                elif isinstance(metabolite, string_types):
                    # if we want to handle creation, this should be changed
                    raise ValueError("reaction '%s' does not belong to a model"
                                     % self.id)
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

    def subtract_metabolites(self, metabolites, combine=True):
        """This function will 'subtract' metabolites from a reaction, which
        means add the metabolites with -1*coefficient. If the final coefficient
        for a metabolite is 0 then the metabolite is removed from the reaction.

        metabolites: dict of {:class:`~cobra.core.Metabolite`: coefficient}
            These metabolites will be added to the reaction

        .. note:: A final coefficient < 0 implies a reactant.

        """
        self.add_metabolites({k: -v for k, v in iteritems(metabolites)},
                             combine=combine)

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
        id_type = 'id'
        if use_metabolite_names:
            id_type = 'name'
        reactant_bits = []
        product_bits = []
        for the_metabolite, coefficient in sorted(
                iteritems(self._metabolites), key=lambda x: x[0].id):
            name = str(getattr(the_metabolite, id_type))
            if _is_positive(coefficient):
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
        self._genes.discard(cobra_gene)
        cobra_gene._reaction.discard(self)

    def knock_out(self):
        """Change the upper and lower bounds of the reaction to 0."""
        self.lower_bound = 0
        self.upper_bound = 0

    def build_reaction_from_string(self, reaction_str, verbose=True,
                                   fwd_arrow=None, rev_arrow=None,
                                   reversible_arrow=None, term_split="+"):
        """Builds reaction from reaction equation reaction_str using parser

        Takes a string and using the specifications supplied in the optional
        arguments infers a set of metabolites, metabolite compartments and
        stoichiometries for the reaction.  It also infers the reversibility
        of the reaction from the reaction arrow.

        Args:
            reaction_str: a string containing a reaction formula (equation)
            verbose: Boolean setting verbosity of function
                (optional, default=True)
            fwd_arrow: re.compile for forward irreversible reaction arrows
                (optional, default=_forward_arrow_finder)
            reverse_arrow: re.compile for backward irreversible reaction arrows
                (optional, default=_reverse_arrow_finder)
            fwd_arrow: re.compile for reversible reaction arrows
                (optional, default=_reversible_arrow_finder)
            term_split: String dividing individual metabolite entries
                (optional, default='+')
        """
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
