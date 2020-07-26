# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pandas as pd
import pytest

import cobra.medium as medium
from cobra import Metabolite, Reaction


class TestModelMedium:
    def test_model_medium(self, model):
        # Add a dummy 'malformed' import reaction
        bad_import = Reaction('bad_import')
        bad_import.add_metabolites({model.metabolites.pyr_c: 1})
        bad_import.bounds = (0, 42)
        model.add_reaction(bad_import)

        # Test basic setting and getting methods
        medium = model.medium
        model.medium = medium
        assert model.medium == medium

        # Test context management
        with model:
            # Ensure the bounds are correct beforehand
            assert model.reactions.EX_glc__D_e.lower_bound == -10
            assert model.reactions.bad_import.upper_bound == 42
            assert model.reactions.EX_co2_e.lower_bound == -1000

            # Make changes to the media
            new_medium = model.medium
            new_medium['EX_glc__D_e'] = 20
            new_medium['bad_import'] = 24
            del new_medium['EX_co2_e']

            # Change the medium, make sure changes work
            model.medium = new_medium
            assert model.reactions.EX_glc__D_e.lower_bound == -20
            assert model.reactions.bad_import.upper_bound == 24
            assert model.reactions.EX_co2_e.lower_bound == 0

        # Make sure changes revert after the contex
        assert model.reactions.EX_glc__D_e.lower_bound == -10
        assert model.reactions.bad_import.upper_bound == 42
        assert model.reactions.EX_co2_e.lower_bound == -1000

        new_medium['bogus_rxn'] = 0
        with pytest.raises(KeyError):
            model.medium = new_medium


class TestTypeDetection:

    def test_external_compartment(self, model):
        # by name
        assert medium.find_external_compartment(model) == "e"
        # from boundary counts
        for m in model.metabolites:
            if m.compartment == "e":
                m.compartment = "outside"
        for r in model.reactions:
            r._compartments = None
        assert medium.find_external_compartment(model) == "outside"
        # names are always right
        model.exchanges[0].reactants[0].compartment = "extracellular"
        assert medium.find_external_compartment(model) == "extracellular"

    def test_multi_external(self, model):
        for r in model.reactions:
            r._compartments = None
        model.exchanges[0].reactants[0].compartment = "extracellular"
        # still works due to different boundary numbers
        assert medium.find_external_compartment(model) == "e"
        model.exchanges[1].reactants[0].compartment = "extra cellular"
        model.remove_reactions(model.exchanges)
        # Now fails because same boundary count
        with pytest.raises(RuntimeError):
            medium.find_external_compartment(model)

    def test_exchange(self, model):
        ex = model.exchanges
        assert all(r.id.startswith("EX_") for r in ex)
        ex = medium.find_boundary_types(model, "exchange", "e")
        assert all(r.id.startswith("EX_") for r in ex)

    def test_demand(self, model):
        dm = Reaction("demand")
        model.add_reaction(dm)
        dm.build_reaction_from_string("atp_c ->")
        dm = model.demands
        assert len(dm) == 1
        assert "demand" in [r.id for r in dm]

    def test_sink(self, model):
        sn = Reaction("sink")
        model.add_reaction(sn)
        sn.build_reaction_from_string("atp_c <->")
        sn.bounds = -1000, 1000
        sn = model.sinks
        assert len(sn) == 1
        assert "sink" in [r.id for r in sn]

    def test_sbo_terms(self, model):
        assert not medium.is_boundary_type(
            model.reactions.ATPM, "exchange", "e")
        model.reactions.ATPM.annotation["sbo"] = "SBO:0000627"
        assert medium.is_boundary_type(model.reactions.ATPM, "exchange", "bla")
        model.reactions.ATPM.annotation["sbo"] = "SBO:0000632"
        assert not medium.is_boundary_type(
            model.reactions.ATPM, "exchange", "e")


class TestMinimalMedia:

    def test_medium_linear(self, model):
        med = medium.minimal_medium(model, 0.8)
        assert len(med) <= 4
        assert all(med > 1e-6)

    def test_medium_mip(self, model):
        med = medium.minimal_medium(model, 0.8, minimize_components=True)
        assert len(med) <= 4
        assert all(med > 1e-6)

        # Anaerobic growth
        med = medium.minimal_medium(model, 0.1, minimize_components=True)
        assert len(med) <= 3
        assert all(med > 1e-6)

    def test_medium_alternative_mip(self, model):
        med = medium.minimal_medium(model, 0.8, minimize_components=5,
                                    open_exchanges=True)
        assert isinstance(med, pd.DataFrame)
        assert med.shape[0] >= 5
        assert med.shape[1] == 5
        assert all((med > 0).sum() == 3)
        assert all(med.sum(axis=1) > 1e-6)

    def test_benchmark_medium_linear(self, model, benchmark):
        benchmark(medium.minimal_medium, model, 0.8)

    def test_benchmark_medium_mip(self, model, benchmark):
        benchmark(medium.minimal_medium, model, 0.8, True)

    def test_medium_exports(self, model):
        med = medium.minimal_medium(model, 0.8, exports=True,
                                    minimize_components=True)
        assert len(med) > 4
        assert any(med < -1e-6)

    def test_open_exchanges(self, model):
        model.reactions.EX_glc__D_e.bounds = 0, 0
        med = medium.minimal_medium(model, 0.8)
        assert med is None
        med = medium.minimal_medium(model, 0.8, minimize_components=True)
        assert med is None

        med = medium.minimal_medium(model, 0.8, open_exchanges=True)
        assert len(med) >= 3
        med = medium.minimal_medium(model, 0.8, open_exchanges=100)
        assert len(med) >= 3


class TestErrorsAndExceptions:

    def test_no_boundary_reactions(self, empty_model):
        assert medium.find_boundary_types(empty_model, 'e', None) == []

    def test_no_names_or_boundary_reactions(self, empty_model):
        with pytest.raises(RuntimeError):
            medium.find_external_compartment(empty_model)

    def test_bad_exchange(self, model):
        with pytest.raises(ValueError):
            m = Metabolite("baddy", compartment="nonsense")
            model.add_boundary(m, type="exchange")
        m = Metabolite("goody", compartment="e")
        rxn = model.add_boundary(m, type="exchange")
        assert isinstance(rxn, Reaction)
