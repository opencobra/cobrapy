# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pytest
from cobra.core import Metabolite, Model, Compartment


def test_compartment_setter():
    met = Metabolite('x', compartment=None)
    cytosol = Compartment('c')
    model = Model()
    model.add_compartments([cytosol])
    met2 = Metabolite("y")
    met5 = Metabolite("b")
    model.add_metabolites([met2, met5])
    met2.compartment = Compartment("c")
    met3 = Metabolite('z', compartment=cytosol)
    met4 = Metabolite('a', compartment='e')
    met5.compartment = 'c'

    assert met.compartment is None
    assert met2.compartment is cytosol
    assert met3.compartment is cytosol
    assert isinstance(met4.compartment, Compartment)
    assert met5.compartment is cytosol
    with pytest.raises(TypeError):
        Metabolite("c", compartment="Sane Compartment, Not!")
