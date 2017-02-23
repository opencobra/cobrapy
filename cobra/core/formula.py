# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
from warnings import warn

from cobra.core.object import Object

# Numbers are not required because of the |(?=[A-Z])? block. See the
# discussion in https://github.com/opencobra/cobrapy/issues/128 for
# more details.
element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")


class Formula(Object):
    """Describes a Chemical Formula

    Parameters
    ---------
    formula : string
        A legal formula string contains only letters and numbers.
    """
    def __init__(self, formula=None):
        Object.__init__(self, formula)
        self.formula = formula
        self.elements = {}
        if self.formula is not None:
            self.parse_composition()

    def __add__(self, other_formula):
        """Combine two molecular formulas.

        Parameters
        ----------
        other_formula : Formula, str
            string for a chemical formula

        Returns
        -------
        Formula
           The combined formula
        """
        return Formula(self.formula + other_formula.formula)

    def parse_composition(self):
        """Breaks the chemical formula down by element."""
        tmp_formula = self.formula
        # commonly occuring characters in incorrectly constructed formulas
        if "*" in tmp_formula:
            warn("invalid character '*' found in formula '%s'" % self.formula)
            tmp_formula = self.formula.replace("*", "")
        if "(" in tmp_formula or ")" in tmp_formula:
            warn("parenthesis found in formula '%s'" % self.formula)
            return
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
                    self.elements = {}
                    return
            if element in composition:
                composition[element] += count
            else:
                composition[element] = count
        self.elements = composition

    @property
    def weight(self):
        """Calculate the mol mass of the compound

        Returns
        -------
        float
            the mol mass
        """
        try:
            return sum([count * elements_and_molecular_weights[element]
                        for element, count in self.elements.items()])
        except KeyError as e:
            warn("The element %s does not appear in the periodic table" % e)


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
