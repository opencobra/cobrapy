from copy import deepcopy
import warnings
import pytest
from cobra.core import Model, Metabolite, Reaction
from cobra.solvers import solver_dict
from .conftest import model

try:
    import numpy
    from cobra.util import create_stoichiometric_array
except ImportError:
    numpy = None

try:
    import scipy
except ImportError:
    scipy = None

class TestReactions:
    def test_gpr(self):
        model = Model()
        reaction = Reaction("test")
        # set a gpr to  reaction not in a model
        reaction.gene_reaction_rule = "(g1 or g2) and g3"
        assert reaction.gene_reaction_rule == "(g1 or g2) and g3"
        assert len(reaction.genes) == 3
        # adding reaction with a GPR propagates to the model
        model.add_reaction(reaction)
        assert len(model.genes) == 3
        # ensure the gene objects are the same in the model and reaction
        reaction_gene = list(reaction.genes)[0]
        model_gene = model.genes.get_by_id(reaction_gene.id)
        assert reaction_gene is model_gene
        # test ability to handle uppercase AND/OR
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            reaction.gene_reaction_rule = "(b1 AND b2) OR (b3 and b4)"
        assert reaction.gene_reaction_rule == "(b1 and b2) or (b3 and b4)"
        assert len(reaction.genes) == 4
        # ensure regular expressions correctly extract genes from malformed
        # GPR string
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            reaction.gene_reaction_rule = "(a1 or a2"
            assert len(reaction.genes) == 2
            reaction.gene_reaction_rule = "(forT or "
            assert len(reaction.genes) == 1

    def test_gpr_modification(self, model):
        reaction = model.reactions.get_by_id("PGI")
        old_gene = list(reaction.genes)[0]
        new_gene = model.genes.get_by_id("s0001")
        # add an existing 'gene' to the gpr
        reaction.gene_reaction_rule = 's0001'
        assert new_gene in reaction.genes
        assert reaction in new_gene.reactions
        # removed old gene correctly
        assert old_gene not in reaction.genes
        assert reaction not in old_gene.reactions
        # add a new 'gene' to the gpr
        reaction.gene_reaction_rule = 'fake_gene'
        assert model.genes.has_id("fake_gene")
        fake_gene = model.genes.get_by_id("fake_gene")
        assert fake_gene in reaction.genes
        assert reaction in fake_gene.reactions
        fake_gene.name = "foo_gene"
        assert reaction.gene_name_reaction_rule == fake_gene.name

    @pytest.mark.parametrize("solver", list(solver_dict))
    def test_add_metabolite_benchmark(self, model, benchmark, solver):
        reaction = model.reactions.get_by_id("PGI")
        many_metabolites = dict((m, 1) for m in model.metabolites[0:50])

        def add_remove_metabolite():
            reaction.add_metabolites(many_metabolites)
            if not getattr(model, 'solver', None):
                solver_dict[solver].create_problem(model)
            for m, c in many_metabolites.items():
                try:
                    reaction.pop(m.id)
                except KeyError:
                    pass
        benchmark(add_remove_metabolite)

    def test_add_metabolite(self, model):
        reaction = model.reactions.get_by_id("PGI")
        reaction.add_metabolites({model.metabolites[0]: 1})
        assert model.metabolites[0] in reaction._metabolites
        fake_metabolite = Metabolite("fake")
        reaction.add_metabolites({fake_metabolite: 1})
        assert fake_metabolite in reaction._metabolites
        assert model.metabolites.has_id("fake")
        assert model.metabolites.get_by_id("fake") is fake_metabolite

        # test adding by string
        reaction.add_metabolites({"g6p_c": -1})  # already in reaction
        assert reaction._metabolites[
                   model.metabolites.get_by_id("g6p_c")] == -2
        reaction.add_metabolites({"h_c": 1})
        assert reaction._metabolites[model.metabolites.get_by_id("h_c")] == 1
        with pytest.raises(KeyError):
            reaction.add_metabolites({"missing": 1})

        # test adding to a new Reaction
        reaction = Reaction("test")
        assert len(reaction._metabolites) == 0
        reaction.add_metabolites({Metabolite("test_met"): -1})
        assert len(reaction._metabolites) == 1

    @pytest.mark.parametrize("solver", list(solver_dict))
    def test_subtract_metabolite_benchmark(self, model, benchmark, solver):
        benchmark(self.test_subtract_metabolite, model, solver)

    @pytest.mark.parametrize("solver", list(solver_dict))
    def test_subtract_metabolite(self, model, solver):
        reaction = model.reactions.get_by_id("PGI")
        reaction.subtract_metabolites(reaction.metabolites)
        if not getattr(model, 'solver', None):
            solver_dict[solver].create_problem(model)
        assert len(reaction.metabolites) == 0

    def test_mass_balance(self, model):
        reaction = model.reactions.get_by_id("PGI")
        # should be balanced now
        assert len(reaction.check_mass_balance()) == 0
        # should not be balanced after adding a hydrogen
        reaction.add_metabolites({model.metabolites.get_by_id("h_c"): 1})
        imbalance = reaction.check_mass_balance()
        assert "H" in imbalance
        assert imbalance["H"] == 1

    def test_build_from_string(self, model):
        m = len(model.metabolites)
        pgi = model.reactions.get_by_id("PGI")
        pgi.reaction = "g6p_c --> f6p_c"
        assert pgi.lower_bound == 0
        pgi.bounds = (0, 1000)
        assert pgi.bounds == (0, 1000)
        assert not pgi.reversibility
        pgi.reaction = "g6p_c <== f6p_c"
        assert pgi.upper_bound == 0
        assert pgi.reaction.strip() == "g6p_c <-- f6p_c"
        pgi.reaction = "g6p_c --> f6p_c + h2o_c"
        assert model.metabolites.h2o_c, pgi._metabolites
        pgi.build_reaction_from_string("g6p_c --> f6p_c + foo", verbose=False)
        assert model.metabolites.h2o_c not in pgi._metabolites
        assert "foo" in model.metabolites
        assert model.metabolites.foo in pgi._metabolites
        assert len(model.metabolites) == m + 1

    def test_bounds_setter(self, model):
        rxn = model.reactions.get_by_id("PGI")
        with pytest.raises(AssertionError):
            rxn.bounds = (1, 0)

    def test_copy(self, model):
        PGI = model.reactions.PGI
        copied = PGI.copy()
        assert PGI is not copied
        assert PGI._model is model
        assert copied._model is not model
        # the copy should refer to different metabolites and genes
        for met in copied.metabolites:
            assert met is not model.metabolites.get_by_id(met.id)
            assert met.model is not model
        for gene in copied.genes:
            assert gene is not model.genes.get_by_id(gene.id)
            assert gene.model is not model

    def test_iadd(self, model):
        PGI = model.reactions.PGI
        EX_h2o = model.reactions.EX_h2o_e
        original_PGI_gpr = PGI.gene_reaction_rule
        PGI += EX_h2o
        assert PGI.gene_reaction_rule == original_PGI_gpr
        assert PGI.metabolites[model.metabolites.h2o_e] == -1.0
        # original should not have changed
        assert EX_h2o.gene_reaction_rule == ''
        assert EX_h2o.metabolites[model.metabolites.h2o_e] == -1.0
        # what about adding a reaction not in the model
        new_reaction = Reaction("test")
        new_reaction.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
        PGI += new_reaction
        assert PGI.gene_reaction_rule == original_PGI_gpr
        assert len(PGI.gene_reaction_rule) == 5
        # and vice versa
        new_reaction += PGI
        assert len(new_reaction.metabolites) == 5  # not
        assert len(new_reaction.genes) == 1
        assert new_reaction.gene_reaction_rule == original_PGI_gpr
        # what about combining 2 gpr's
        model.reactions.ACKr += model.reactions.ACONTa
        expected_rule = '(b2296 or b3115 or b1849) and (b0118 or b1276)'
        assert model.reactions.ACKr.gene_reaction_rule == expected_rule
        assert len(model.reactions.ACKr.genes) == 5

    def test_add(self, model):
        # not in place addition should work on a copy
        new = model.reactions.PGI + model.reactions.EX_h2o_e
        assert new._model is not model
        assert len(new.metabolites) == 3
        # the copy should refer to different metabolites and genes
        # This currently fails because add_metabolites does not copy.
        # Should that be changed?
        # for met in new.metabolites:
        #    assert met is not model.metabolites.get_by_id(met.id)
        #    assert met.model is not model
        for gene in new.genes:
            assert gene is not model.genes.get_by_id(gene.id)
            assert gene.model is not model

    def test_mul(self, model):
        new = model.reactions.PGI * 2
        assert set(new.metabolites.values()) == {-2, 2}

    def test_sub(self, model):
        new = model.reactions.PGI - model.reactions.EX_h2o_e
        assert new._model is not model
        assert len(new.metabolites) == 3


class TestCobraMetabolites:
    def test_metabolite_formula(self):
        met = Metabolite("water")
        met.formula = "H2O"
        assert met.elements == {"H": 2, "O": 1}
        assert met.formula_weight == 18.01528

    def test_formula_element_setting(self, model):
        met = model.metabolites[1]
        orig_formula = str(met.formula)
        orig_elements = dict(met.elements)
        met.formula = ''
        assert met.elements == {}
        met.elements = orig_elements
        assert met.formula == orig_formula


class TestCobraModel:
    """test core cobra functions"""

    @pytest.mark.parametrize("solver", list(solver_dict))
    def test_add_remove_reaction_benchmark(self, model, benchmark, solver):
        metabolite_foo = Metabolite("test_foo")
        metabolite_bar = Metabolite("test_bar")
        metabolite_baz = Metabolite("test_baz")
        actual_metabolite = model.metabolites[0]
        dummy_reaction = Reaction("test_foo_reaction")
        dummy_reaction.add_metabolites({metabolite_foo: -1,
                                        metabolite_bar: 1,
                                        metabolite_baz: -2,
                                        actual_metabolite: 1})

        def benchmark_add_reaction():
            model.add_reaction(dummy_reaction)
            if not getattr(model, 'solver', None):
                solver_dict[solver].create_problem(model)
            model.remove_reactions([dummy_reaction], delete=False)
        benchmark(benchmark_add_reaction)

    def test_add_reaction(self, model):
        old_reaction_count = len(model.reactions)
        old_metabolite_count = len(model.metabolites)
        dummy_metabolite_1 = Metabolite("test_foo_1")
        dummy_metabolite_2 = Metabolite("test_foo_2")
        actual_metabolite = model.metabolites[0]
        copy_metabolite = model.metabolites[1].copy()
        dummy_reaction = Reaction("test_foo_reaction")
        dummy_reaction.add_metabolites({dummy_metabolite_1: -1,
                                        dummy_metabolite_2: 1,
                                        copy_metabolite: -2,
                                        actual_metabolite: 1})

        model.add_reaction(dummy_reaction)
        assert model.reactions.get_by_id(dummy_reaction.id) == dummy_reaction
        for x in [dummy_metabolite_1, dummy_metabolite_2]:
            assert model.metabolites.get_by_id(x.id) == x
        # should have added 1 reaction and 2 metabolites
        assert len(model.reactions) == old_reaction_count + 1
        assert len(model.metabolites) == old_metabolite_count + 2
        # tests on the added reaction
        reaction_in_model = model.reactions.get_by_id(dummy_reaction.id)
        assert type(reaction_in_model) is Reaction
        assert reaction_in_model is dummy_reaction
        assert len(reaction_in_model._metabolites) == 4
        for i in reaction_in_model._metabolites:
            assert type(i) == Metabolite
        # tests on the added metabolites
        met1_in_model = model.metabolites.get_by_id(dummy_metabolite_1.id)
        assert met1_in_model is dummy_metabolite_1
        copy_in_model = model.metabolites.get_by_id(copy_metabolite.id)
        assert copy_metabolite is not copy_in_model
        assert type(copy_in_model) is Metabolite
        assert dummy_reaction in actual_metabolite._reaction
        # test adding a different metabolite with the same name as an
        # existing one uses the metabolite in the model
        r2 = Reaction("test_foo_reaction2")
        model.add_reaction(r2)
        r2.add_metabolites({Metabolite(model.metabolites[0].id): 1})
        assert model.metabolites[0] is list(r2._metabolites)[0]

    def test_add_reaction_from_other_model(self, model):
        other = model.copy()
        for i in other.reactions:
            i.id += "_other"
        other.repair()
        model.add_reactions(other.reactions)
        # what if the other reaction has an error in its GPR
        m1 = model.copy()
        m2 = model.copy()
        m1.reactions.PGI.remove_from_model()
        m2.genes.b4025._reaction.clear()
        m1.add_reaction(m2.reactions.PGI)

    def test_model_remove_reaction(self, model):
        old_reaction_count = len(model.reactions)
        model.remove_reactions(["PGI"])
        assert len(model.reactions) == old_reaction_count - 1
        with pytest.raises(KeyError):
            model.reactions.get_by_id("PGI")
        model.remove_reactions(model.reactions[:1])
        assert len(model.reactions) == old_reaction_count - 2
        tmp_metabolite = Metabolite("testing")
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        assert tmp_metabolite in model.metabolites
        model.remove_reactions(model.reactions[:1],
                               remove_orphans=True)
        assert tmp_metabolite not in model.metabolites

    def test_reaction_remove(self, model):
        old_reaction_count = len(model.reactions)
        tmp_metabolite = Metabolite("testing")
        # Delete without removing orphan
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        assert len(tmp_metabolite.reactions) == 1
        # esnsure the stoichiometry is still the same using different objects
        removed_reaction = model.reactions[0]
        original_stoich = {i.id: value for i, value
                           in removed_reaction._metabolites.items()}
        model.reactions[0].remove_from_model(remove_orphans=False)
        assert len(original_stoich) == len(removed_reaction._metabolites)
        for met in removed_reaction._metabolites:
            assert original_stoich[met.id] == removed_reaction._metabolites[
                met]
            assert met is not model.metabolites
        # make sure it's still in the model
        assert tmp_metabolite in model.metabolites
        assert len(tmp_metabolite.reactions) == 0
        assert len(model.reactions) == old_reaction_count - 1

        # Now try it with removing orphans
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        assert len(tmp_metabolite.reactions) == 1
        model.reactions[0].remove_from_model(remove_orphans=True)
        assert tmp_metabolite not in model.metabolites
        assert len(tmp_metabolite.reactions) == 0
        assert len(model.reactions) == old_reaction_count - 2

        # It shouldn't remove orphans if it's in 2 reactions however
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        model.reactions[1].add_metabolites({tmp_metabolite: 1})
        assert len(tmp_metabolite.reactions) == 2
        model.reactions[0].remove_from_model(remove_orphans=False)
        assert tmp_metabolite in model.metabolites
        assert len(tmp_metabolite.reactions) == 1
        assert len(model.reactions) == old_reaction_count - 3

    def test_reaction_delete(self, model):
        old_reaction_count = len(model.reactions)
        tmp_metabolite = Metabolite("testing")
        # Delete without removing orphan
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        assert len(tmp_metabolite.reactions) == 1
        model.reactions[0].delete(remove_orphans=False)
        # make sure it's still in the model
        assert tmp_metabolite in model.metabolites
        assert len(tmp_metabolite.reactions) == 0
        assert len(model.reactions) == old_reaction_count - 1

        # Now try it with removing orphans
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        assert len(tmp_metabolite.reactions) == 1
        model.reactions[0].delete(remove_orphans=True)
        assert tmp_metabolite not in model.metabolites
        assert len(tmp_metabolite.reactions) == 0
        assert len(model.reactions) == old_reaction_count - 2

        # It shouldn't remove orphans if it's in 2 reactions however
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        model.reactions[1].add_metabolites({tmp_metabolite: 1})
        assert len(tmp_metabolite.reactions) == 2
        model.reactions[0].delete(remove_orphans=False)
        assert tmp_metabolite in model.metabolites
        assert len(tmp_metabolite.reactions) == 1
        assert len(model.reactions) == old_reaction_count - 3

    def test_remove_gene(self, model):
        target_gene = model.genes[0]
        gene_reactions = list(target_gene.reactions)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            target_gene.remove_from_model()
        assert target_gene.model is None
        # make sure the reaction was removed from the model
        assert target_gene not in model.genes
        # ensure the old reactions no longer have a record of the gene
        for reaction in gene_reactions:
            assert target_gene not in reaction.genes

    @pytest.mark.parametrize("solver", list(solver_dict))
    def test_copy_benchmark(self, model, solver, benchmark):
        def _():
            model.copy()
            if not getattr(model, 'solver', None):
                solver_dict[solver].create_problem(model)
        benchmark(_)

    @pytest.mark.parametrize("solver", list(solver_dict))
    def test_copy_benchmark_large_model(self, large_model, solver, benchmark):
        def _():
            large_model.copy()
            if not getattr(large_model, 'solver', None):
                solver_dict[solver].create_problem(large_model)
        benchmark(_)

    def test_copy(self, model):
        """modifying copy should not modify the original"""
        # test that deleting reactions in the copy does not change the
        # number of reactions in the original model
        model_copy = model.copy()
        old_reaction_count = len(model.reactions)
        assert model_copy.notes is not model.notes
        assert model_copy.annotation is not model.annotation
        assert len(model.reactions) == len(model_copy.reactions)
        assert len(model.metabolites) == len(model_copy.metabolites)
        model_copy.remove_reactions(model_copy.reactions[0:5])
        assert old_reaction_count == len(model.reactions)
        assert len(model.reactions) != len(model_copy.reactions)

    def test_deepcopy_benchmark(self, model, benchmark):
        benchmark(deepcopy, model)

    def test_deepcopy(self, model):
        """Reference structures are maintained when deepcopying"""
        model_copy = deepcopy(model)
        for gene, gene_copy in zip(model.genes, model_copy.genes):
            assert gene.id == gene_copy.id
            reactions = sorted(i.id for i in gene.reactions)
            reactions_copy = sorted(i.id for i in gene_copy.reactions)
            assert reactions == reactions_copy
        for reaction, reaction_copy in zip(model.reactions,
                                           model_copy.reactions):
            assert reaction.id == reaction_copy.id
            metabolites = sorted(i.id for i in reaction._metabolites)
            metabolites_copy = sorted(i.id for i in reaction_copy._metabolites)
            assert metabolites == metabolites_copy

    def test_add_reaction_orphans(self, model):
        """test reaction addition

        Need to verify that no orphan genes or metabolites are
        contained in reactions after adding them to the model.
        """
        model = model.__class__('test')
        model.add_reactions((x.copy() for x in model.reactions))
        genes = []
        metabolites = []
        for x in model.reactions:
            genes.extend(x.genes)
            metabolites.extend(x._metabolites)
        orphan_genes = [x for x in genes if x.model is not model]
        orphan_metabolites = [x for x in metabolites if x.model is not model]
        # check not dangling genes when running Model.add_reactions
        assert len(orphan_genes) == 0
        # 'check not dangling metabolites when running Model.add_reactions
        assert len(orphan_metabolites) == 0

    @pytest.mark.parametrize("solver", list(solver_dict))
    def test_change_objective_benchmark(self, model, benchmark, solver):
        atpm = model.reactions.get_by_id("ATPM")

        def benchmark_change_objective():
            model.objective = atpm.id
            if not getattr(model, 'solver', None):
                solver_dict[solver].create_problem(model)
        benchmark(benchmark_change_objective)

    def test_change_objective(self, model):
        biomass = model.reactions.get_by_id("Biomass_Ecoli_core")
        atpm = model.reactions.get_by_id("ATPM")
        model.objective = atpm.id
        assert atpm.objective_coefficient == 1.
        assert biomass.objective_coefficient == 0.
        assert model.objective == {atpm: 1.}
        # change it back using object itself
        model.objective = biomass
        assert atpm.objective_coefficient == 0.
        assert biomass.objective_coefficient == 1.
        # set both to 1 with a list
        model.objective = [atpm, biomass]
        assert atpm.objective_coefficient == 1.
        assert biomass.objective_coefficient == 1.
        # set both using a dict
        model.objective = {atpm: 0.2, biomass: 0.3}
        assert abs(atpm.objective_coefficient - 0.2) < 10 ** -9
        assert abs(biomass.objective_coefficient - 0.3) < 10 ** -9
        # test setting by index
        model.objective = model.reactions.index(atpm)
        assert model.objective == {atpm: 1.}
        # test by setting list of indexes
        model.objective = [model.reactions.index(reaction) for
                           reaction in [atpm, biomass]]
        assert model.objective == {atpm: 1., biomass: 1.}


    def test_context_manager(self, model):
        bounds0 = model.reactions[0].bounds
        bounds1 = (1, 2)
        bounds2 = (3, 4)

        # Trigger a nested model context, ensuring that bounds are
        # preserved at each level
        with model:
            model.reactions[0].bounds = bounds1
            with model:
                model.reactions[0].bounds = bounds2

                assert model.reactions[0].bounds == bounds2
            assert model.reactions[0].bounds == bounds1
        assert model.reactions[0].bounds == bounds0


class TestStoichiometricMatrix:
    """Test the simple replacement for ArrayBasedModel"""

    @pytest.mark.skipif(not numpy, reason='Array methods require numpy')
    def test_dense_matrix(self, model):
        S = create_stoichiometric_array(model, array_type='dense', dtype=int)
        assert S.dtype == int
        assert numpy.allclose(S.max(), 59)

        S = create_stoichiometric_array(model, array_type='dense', dtype=float)
        model.optimize()
        # Is this really the best way to get a vector of fluxes?
        mass_balance = S.dot(numpy.array(list(model.solution.fluxes.values())))
        assert numpy.allclose(mass_balance, 0)

        # Test model property
        assert numpy.allclose(model.S, S)

    @pytest.mark.skipif(not scipy, reason='Sparse array methods require scipy')
    def test_sparse_matrix(self, model):
        sparse_types = ['dok', 'lil']

        model.optimize()
        fluxes = numpy.array(list(model.solution.fluxes.values()))
        for sparse_type in sparse_types:
            S = create_stoichiometric_array(model, array_type=sparse_type)
            mass_balance = S.dot(fluxes)
            assert numpy.allclose(mass_balance, 0)

        # Is this really the best way to get a vector of fluxes?
