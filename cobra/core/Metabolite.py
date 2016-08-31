from warnings import warn
import re

from six import iteritems

from .Species import Species

# Numbers are not required because of the |(?=[A-Z])? block. See the
# discussion in https://github.com/opencobra/cobrapy/issues/128 for
# more details.
element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")


class Metabolite(Species):
    """Metabolite is a class for holding information regarding
    a metabolite in a cobra.Reaction object.

    """

    def __init__(self, id=None, formula=None, name="",
                 charge=None, compartment=None):
        """
        id: str

        formula: str
            Chemical formula (i.e. H2O)

        name: str
            A human readable name.

        compartment: str or None
            Compartment of metabolite.
        """
        Species.__init__(self, id, name)
        self.formula = formula
        # because in a Model a metabolite may participate in multiple Reactions
        self.compartment = compartment
        self.charge = charge

        self._constraint_sense = 'E'
        self._bound = 0.

    @property
    def elements(self):
        tmp_formula = self.formula
        if tmp_formula is None:
            return {}
        # necessary for some old pickles which use the deprecated
        # Formula class
        tmp_formula = str(self.formula)
        # commonly occuring characters in incorrectly constructed formulas
        if "*" in tmp_formula:
            warn("invalid character '*' found in formula '%s'" % self.formula)
            tmp_formula = tmp_formula.replace("*", "")
        if "(" in tmp_formula or ")" in tmp_formula:
            warn("invalid formula (has parenthesis) in '%s'" % self.formula)
            return None
        composition = {}
        parsed = element_re.findall(tmp_formula)
        for (element, count) in parsed:
            if count == '':
                count = 1
            else:
                try:
                    count = float(count)
                    int_count = int(count)
                    if count == int_count:
                        count = int_count
                    else:
                        warn("%s is not an integer (in formula %s)" %
                             (count, self.formula))
                except ValueError:
                    warn("failed to parse %s (in formula %s)" %
                         (count, self.formula))
                    return None
            if element in composition:
                composition[element] += count
            else:
                composition[element] = count
        return composition

    @elements.setter
    def elements(self, elements_dict):
        def stringify(element, number):
            return element if number == 1 else element + str(number)

        self.formula = ''.join(stringify(e, n) for e, n in
                               sorted(iteritems(elements_dict)))

    @property
    def formula_weight(self):
        """Calculate the formula weight"""
        try:
            return sum([count * elements_and_molecular_weights[element]
                        for element, count in self.elements.items()])
        except KeyError as e:
            warn("The element %s does not appear in the peridic table" % e)

    @property
    def y(self):
        """The shadow price for the metabolite in the most recent solution

        Shadow prices are computed from the dual values of the bounds in
        the solution.

        """
        try:
            return self._model.solution.y_dict[self.id]
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

    def remove_from_model(self, method='subtractive', **kwargs):
        """Removes the association from self.model

        method: 'subtractive' or 'destructive'.
            If 'subtractive' then the metabolite is removed from all
            associated reactions.  If 'destructive' then all associated
            reactions are removed from the Model.

        """
        # why is model being taken in as a parameter? This plays
        # back to the question of allowing a Metabolite to be associated
        # with multiple Models
        if "model" in kwargs:
            warn("model argument deprecated")

        self._model.metabolites.remove(self)
        self._model = None
        if method.lower() == 'subtractive':
            for the_reaction in list(self._reaction):
                the_coefficient = the_reaction._metabolites[self]
                the_reaction.subtract_metabolites({self: the_coefficient})
        elif method.lower() == 'destructive':
            for x in self._reaction:
                x.remove_from_model()
        else:
            raise Exception(method + " is not 'subtractive' or 'destructive'")

    def summary(self, **kwargs):
        """Print a summary of the reactions which produce and consume this
        metabolite. This method requires the model for which this metabolite is
        a part to be solved.

        threshold: float
            a value below which to ignore reaction fluxes

        fva: float (0->1), or None
            Whether or not to include flux variability analysis in the output.
            If given, fva should be a float between 0 and 1, representing the
            fraction of the optimum objective to be searched.

        floatfmt: string
            format method for floats, passed to tabulate. Default is '.3g'.

        """
        try:
            from ..flux_analysis.summary import metabolite_summary
            return metabolite_summary(self, **kwargs)
        except ImportError:
            warn('Summary methods require pandas/tabulate')

elements_and_molecular_weights = {
    'H':   1.007940,
    'He':  4.002602,
    'Li':  6.941000,
    'Be':  9.012182,
    'B':   10.811000,
    'C':   12.010700,
    'N':   14.006700,
    'O':   15.999400,
    'F':   18.998403,
    'Ne':  20.179700,
    'Na':  22.989770,
    'Mg':  24.305000,
    'Al':  26.981538,
    'Si':  28.085500,
    'P':   30.973761,
    'S':   32.065000,
    'Cl':  35.453000,
    'Ar':  39.948000,
    'K':   39.098300,
    'Ca':  40.078000,
    'Sc':  44.955910,
    'Ti':  47.867000,
    'V':   50.941500,
    'Cr':  51.996100,
    'Mn':  54.938049,
    'Fe':  55.845000,
    'Co':  58.933200,
    'Ni':  58.693400,
    'Cu':  63.546000,
    'Zn':  65.409000,
    'Ga':  69.723000,
    'Ge':  72.640000,
    'As':  74.921600,
    'Se':  78.960000,
    'Br':  79.904000,
    'Kr':  83.798000,
    'Rb':  85.467800,
    'Sr':  87.620000,
    'Y':   88.905850,
    'Zr':  91.224000,
    'Nb':  92.906380,
    'Mo':  95.940000,
    'Tc':  98.000000,
    'Ru':  101.070000,
    'Rh':  102.905500,
    'Pd':  106.420000,
    'Ag':  107.868200,
    'Cd':  112.411000,
    'In':  114.818000,
    'Sn':  118.710000,
    'Sb':  121.760000,
    'Te':  127.600000,
    'I':   126.904470,
    'Xe':  131.293000,
    'Cs':  132.905450,
    'Ba':  137.327000,
    'La':  138.905500,
    'Ce':  140.116000,
    'Pr':  140.907650,
    'Nd':  144.240000,
    'Pm':  145.000000,
    'Sm':  150.360000,
    'Eu':  151.964000,
    'Gd':  157.250000,
    'Tb':  158.925340,
    'Dy':  162.500000,
    'Ho':  164.930320,
    'Er':  167.259000,
    'Tm':  168.934210,
    'Yb':  173.040000,
    'Lu':  174.967000,
    'Hf':  178.490000,
    'Ta':  180.947900,
    'W':   183.840000,
    'Re':  186.207000,
    'Os':  190.230000,
    'Ir':  192.217000,
    'Pt':  195.078000,
    'Au':  196.966550,
    'Hg':  200.590000,
    'Tl':  204.383300,
    'Pb':  207.200000,
    'Bi':  208.980380,
    'Po':  209.000000,
    'At':  210.000000,
    'Rn':  222.000000,
    'Fr':  223.000000,
    'Ra':  226.000000,
    'Ac':  227.000000,
    'Th':  232.038100,
    'Pa':  231.035880,
    'U':   238.028910,
    'Np':  237.000000,
    'Pu':  244.000000,
    'Am':  243.000000,
    'Cm':  247.000000,
    'Bk':  247.000000,
    'Cf':  251.000000,
    'Es':  252.000000,
    'Fm':  257.000000,
    'Md':  258.000000,
    'No':  259.000000,
    'Lr':  262.000000,
    'Rf':  261.000000,
    'Db':  262.000000,
    'Sg':  266.000000,
    'Bh':  264.000000,
    'Hs':  277.000000,
    'Mt':  268.000000,
    'Ds':  281.000000,
    'Rg':  272.000000,
    'Cn':  285.000000,
    'Uuq': 289.000000,
    'Uuh': 292.000000
}
