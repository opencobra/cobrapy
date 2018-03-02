# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
import sys
import warnings
import math
import pytest
import numpy
from contextlib import contextmanager
from optlang.interface import OPTIMAL, INFEASIBLE
from six import StringIO, iteritems

import cobra.util.solver as sutil
from cobra.core import Metabolite, Model, Reaction
from cobra.flux_analysis import *
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.flux_analysis.sampling import ACHRSampler, OptGPSampler
from cobra.flux_analysis.reaction import assess
from cobra.exceptions import Infeasible
from cobra.flux_analysis.moma import add_moma

# The scipy interface is currently unstable and may yield errors or infeasible
# solutions.
stable_optlang = ["glpk", "cplex", "gurobi"]
optlang_solvers = ["optlang-" + s for s in stable_optlang if s in
                   sutil.solvers]
all_solvers = ["optlang-" + s for s in stable_optlang if
               s in sutil.solvers]
qp_solvers = ["cplex"]


def construct_ll_test_model():
    test_model = Model()
    test_model.add_metabolites(Metabolite("A"))
    test_model.add_metabolites(Metabolite("B"))
    test_model.add_metabolites(Metabolite("C"))
    EX_A = Reaction("EX_A")
    EX_A.add_metabolites({test_model.metabolites.A: 1})
    DM_C = Reaction("DM_C")
    DM_C.add_metabolites({test_model.metabolites.C: -1})
    v1 = Reaction("v1")
    v1.add_metabolites({test_model.metabolites.A: -1,
                        test_model.metabolites.B: 1})
    v2 = Reaction("v2")
    v2.add_metabolites({test_model.metabolites.B: -1,
                        test_model.metabolites.C: 1})
    v3 = Reaction("v3")
    v3.add_metabolites({test_model.metabolites.C: -1,
                        test_model.metabolites.A: 1})
    test_model.add_reactions([EX_A, DM_C, v1, v2, v3])
    DM_C.objective_coefficient = 1
    return test_model


@pytest.fixture(scope="function", params=optlang_solvers)
def ll_test_model(request):
    test_model = construct_ll_test_model()
    test_model.solver = request.param
    return test_model


@contextmanager
def captured_output():
    """A context manager to test the IO summary methods."""
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class TestCobraFluxAnalysis:
    """Test the simulation functions in cobra.flux_analysis."""

    @pytest.mark.parametrize("solver", all_solvers)
    def test_pfba_benchmark(self, large_model, benchmark, solver):
        large_model.solver = solver
        benchmark(pfba, large_model)

    @pytest.mark.parametrize("solver", all_solvers)
    def test_pfba(self, model, solver):
        model.solver = solver
        with model:
            add_pfba(model)
            with pytest.raises(ValueError):
                add_pfba(model)

        expression = model.objective.expression
        n_constraints = len(model.constraints)
        solution = pfba(model)
        assert solution.status == "optimal"
        assert numpy.isclose(solution.x_dict["Biomass_Ecoli_core"],
                             0.8739, atol=1e-4, rtol=0.0)
        abs_x = [abs(i) for i in solution.x]
        assert numpy.isclose(sum(abs_x), 518.4221, atol=1e-4, rtol=0.0)
        # test changes to model reverted
        assert expression == model.objective.expression
        assert len(model.constraints) == n_constraints

        # needed?
        # Test desired_objective_value
        # desired_objective = 0.8
        # pfba(model, solver=solver,
        #                       desired_objective_value=desired_objective)
        # abs_x = [abs(i) for i in model.solution.x]
        # assert model.solution.status == "optimal"
        # assert abs(model.solution.f - desired_objective) < 0.001
        # assert abs(sum(abs_x) - 476.1594) < 0.001

        # TODO: parametrize fraction (DRY it up)
        # Test fraction_of_optimum
        solution = pfba(model, fraction_of_optimum=0.95)
        assert solution.status == "optimal"
        assert numpy.isclose(solution.x_dict["Biomass_Ecoli_core"],
                             0.95 * 0.8739, atol=1e-4, rtol=0.0)
        abs_x = [abs(i) for i in solution.x]
        assert numpy.isclose(sum(abs_x), 493.4400, atol=1e-4, rtol=0.0)

        # Infeasible solution
        model.reactions.ATPM.lower_bound = 500
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            with pytest.raises((UserWarning, Infeasible, ValueError)):
                pfba(model)

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_single_gene_deletion_fba_benchmark(self, model, benchmark,
                                                solver):
        model.solver = solver
        benchmark(single_gene_deletion, model)

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_single_gene_deletion_fba(self, model, solver):
        # expected knockouts for textbook model
        model.solver = solver
        growth_dict = {"b0008": 0.87, "b0114": 0.80, "b0116": 0.78,
                       "b2276": 0.21, "b1779": 0.00}
        rates = single_gene_deletion(model,
                                     gene_list=growth_dict.keys(),
                                     method="fba")["growth"]
        for gene, expected_value in iteritems(growth_dict):
            assert abs(rates[frozenset([gene])] - expected_value) < 0.01

    def test_single_gene_deletion_moma_benchmark(self, model, benchmark):
        try:
            solver = sutil.get_solver_name(qp=True)
            model.solver = solver
        except sutil.SolverNotFound:
            pytest.skip("no qp support")

        genes = ['b0008', 'b0114', 'b2276', 'b1779']
        benchmark(single_gene_deletion, model, gene_list=genes,
                  method="moma")

    def test_single_deletion_linear_moma_benchmark(self, model, benchmark):
        try:
            solver = sutil.get_solver_name(qp=True)
            model.solver = solver
        except sutil.SolverNotFound:
            pytest.skip("no qp support")

        genes = ['b0008', 'b0114', 'b2276', 'b1779']
        benchmark(single_gene_deletion, model, gene_list=genes,
                  method="linear moma")

    def test_moma_sanity(self, model):
        """Test optimization criterion and optimality."""
        try:
            solver = sutil.get_solver_name(qp=True)
            model.solver = solver
        except sutil.SolverNotFound:
            pytest.skip("no qp support")

        sol = model.optimize()
        with model:
            model.reactions.PFK.knock_out()
            knock_sol = model.optimize()
            ssq = (knock_sol.fluxes - sol.fluxes).pow(2).sum()

        with model:
            add_moma(model)
            model.reactions.PFK.knock_out()
            moma_sol = model.optimize()
            moma_ssq = (moma_sol.fluxes - sol.fluxes).pow(2).sum()

        assert numpy.allclose(moma_sol.objective_value, moma_ssq)
        assert moma_ssq < ssq

    def test_single_gene_deletion_moma(self, model):
        try:
            solver = sutil.get_solver_name(qp=True)
            model.solver = solver
        except sutil.SolverNotFound:
            pytest.skip("no qp support")

        # expected knockouts for textbook model
        growth_dict = {"b0008": 0.87, "b0114": 0.71, "b0116": 0.56,
                       "b2276": 0.11, "b1779": 0.00}

        result = single_gene_deletion(
            model, gene_list=growth_dict.keys(), method="moma")["growth"]
        for gene, expected_value in iteritems(growth_dict):
            assert abs(result[frozenset([gene])] - expected_value) < 0.01

    def test_linear_moma_sanity(self, model):
        """Test optimization criterion and optimality."""
        try:
            solver = sutil.get_solver_name(qp=True)
            model.solver = solver
        except sutil.SolverNotFound:
            pytest.skip("no qp support")

        sol = model.optimize()
        with model:
            model.reactions.PFK.knock_out()
            knock_sol = model.optimize()
            sabs = (knock_sol.fluxes - sol.fluxes).abs().sum()

        with model:
            add_moma(model, linear=True)
            model.reactions.PFK.knock_out()
            moma_sol = model.optimize()
            moma_sabs = (moma_sol.fluxes - sol.fluxes).abs().sum()

        assert numpy.allclose(moma_sol.objective_value, moma_sabs)
        assert moma_sabs < sabs

    def test_single_gene_deletion_linear_moma(self, model):
        try:
            solver = sutil.get_solver_name(qp=True)
            model.solver = solver
        except sutil.SolverNotFound:
            pytest.skip("no qp support")

        # expected knockouts for textbook model
        growth_dict = {"b0008": 0.87, "b0114": 0.76, "b0116": 0.65,
                       "b2276": 0.08, "b1779": 0.00}

        result = single_gene_deletion(
            model, gene_list=growth_dict.keys(),
            method="linear moma")['growth']
        assert all(abs(result[frozenset([gene])] - expected) < 0.01
                   for gene, expected in iteritems(growth_dict))
        with model:
            add_moma(model, linear=True)
            with pytest.raises(ValueError):
                add_moma(model)

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_single_gene_deletion_benchmark(self, model, benchmark,
                                            solver):
        model.solver = solver
        benchmark(single_reaction_deletion, model)

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_single_reaction_deletion(self, model, solver):
        expected_results = {'FBA': 0.70404, 'FBP': 0.87392, 'CS': 0,
                            'FUM': 0.81430, 'GAPD': 0, 'GLUDy': 0.85139}

        model.solver = solver
        results = single_reaction_deletion(
            model, reaction_list=expected_results.keys()).to_dict()['growth']
        for reaction, value in iteritems(expected_results):
            assert abs(results[frozenset([reaction])] - value) < 1E-05

    @classmethod
    def compare_matrices(cls, matrix1, matrix2, places=3):
        nrows = len(matrix1)
        ncols = len(matrix1[0])
        assert nrows == len(matrix2)
        assert ncols == len(matrix2[0])
        for i in range(nrows):
            for j in range(ncols):
                assert abs(matrix1[i][j] - matrix2[i][j]) < 10 ** -places

    def test_double_gene_deletion_benchmark(self, large_model, benchmark):
        genes = ["b0726", "b4025", "b0724", "b0720", "b2935", "b2935",
                 "b1276",
                 "b1241"]
        benchmark(double_gene_deletion, large_model, gene_list1=genes)

    def test_double_gene_deletion(self, model):
        genes = ["b0726", "b4025", "b0724", "b0720", "b2935", "b2935",
                 "b1276",
                 "b1241"]
        growth_dict = {'b0720': {'b0720': 0.0,
                                 'b0724': 0.0,
                                 'b0726': 0.0,
                                 'b1241': 0.0,
                                 'b1276': 0.0,
                                 'b2935': 0.0,
                                 'b4025': 0.0},
                       'b0724': {'b0720': 0.0,
                                 'b0724': 0.814,
                                 'b0726': 0.814,
                                 'b1241': 0.814,
                                 'b1276': 0.814,
                                 'b2935': 0.814,
                                 'b4025': 0.739},
                       'b0726': {'b0720': 0.0,
                                 'b0724': 0.814,
                                 'b0726': 0.858,
                                 'b1241': 0.858,
                                 'b1276': 0.858,
                                 'b2935': 0.858,
                                 'b4025': 0.857},
                       'b1241': {'b0720': 0.0,
                                 'b0724': 0.814,
                                 'b0726': 0.858,
                                 'b1241': 0.874,
                                 'b1276': 0.874,
                                 'b2935': 0.874,
                                 'b4025': 0.863},
                       'b1276': {'b0720': 0.0,
                                 'b0724': 0.814,
                                 'b0726': 0.858,
                                 'b1241': 0.874,
                                 'b1276': 0.874,
                                 'b2935': 0.874,
                                 'b4025': 0.863},
                       'b2935': {'b0720': 0.0,
                                 'b0724': 0.814,
                                 'b0726': 0.858,
                                 'b1241': 0.874,
                                 'b1276': 0.874,
                                 'b2935': 0.874,
                                 'b4025': 0.863},
                       'b4025': {'b0720': 0.0,
                                 'b0724': 0.739,
                                 'b0726': 0.857,
                                 'b1241': 0.863,
                                 'b1276': 0.863,
                                 'b2935': 0.863,
                                 'b4025': 0.863}}
        solution = double_gene_deletion(
            model, gene_list1=genes, processes=4)['growth']
        solution_one_process = double_gene_deletion(
            model, gene_list1=genes, processes=1)['growth']
        for (rxn_a, sub) in iteritems(growth_dict):
            for rxn_b, growth in iteritems(sub):
                sol = solution[frozenset([rxn_a, rxn_b])]
                sol_one = solution_one_process[frozenset([rxn_a, rxn_b])]
                assert round(sol, 3) == growth
                assert round(sol_one, 3) == growth

    def test_double_reaction_deletion(self, model):
        reactions = ['FBA', 'ATPS4r', 'ENO', 'FRUpts2']
        growth_dict = {
            "FBA": {
                "ATPS4r": 0.135,
                "ENO": float('nan'),
                "FRUpts2": 0.704
            },
            "ATPS4r": {
                "ENO": float('nan'),
                "FRUpts2": 0.374
            },
            "ENO": {
                "FRUpts2": 0.0
            },
        }

        solution = double_reaction_deletion(
            model, reaction_list1=reactions, processes=3)['growth']
        solution_one_process = double_reaction_deletion(
            model, reaction_list1=reactions, processes=1)['growth']
        for (rxn_a, sub) in iteritems(growth_dict):
            for rxn_b, growth in iteritems(sub):
                sol = solution[frozenset([rxn_a, rxn_b])]
                sol_one = solution_one_process[frozenset([rxn_a, rxn_b])]
                if math.isnan(growth):
                    assert math.isnan(sol)
                    assert math.isnan(sol_one)
                else:
                    assert round(sol, 3) == growth
                    assert round(sol_one, 3) == growth

    def test_double_reaction_deletion_benchmark(self, large_model,
                                                benchmark):
        reactions = large_model.reactions[1::100]
        benchmark(double_reaction_deletion, large_model,
                  reaction_list1=reactions)

    @pytest.mark.parametrize("solver", all_solvers)
    def test_flux_variability_benchmark(self, large_model, benchmark,
                                        solver):
        large_model.solver = solver
        benchmark(flux_variability_analysis, large_model,
                  reaction_list=large_model.reactions[1::3])

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_flux_variability_loopless_benchmark(self, model, benchmark,
                                                 solver):
        model.solver = solver
        benchmark(flux_variability_analysis, model, loopless=True,
                  reaction_list=model.reactions[1::3])

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_pfba_flux_variability(self, model, pfba_fva_results,
                                   fva_results, solver):
        model.solver = solver
        with pytest.warns(UserWarning):
            flux_variability_analysis(
                model, pfba_factor=0.1, reaction_list=model.reactions[1::3])
        fva_out = flux_variability_analysis(
            model, pfba_factor=1.1, reaction_list=model.reactions)
        for name, result in iteritems(fva_out.T):
            for k, v in iteritems(result):
                assert abs(pfba_fva_results[k][name] - v) < 0.00001
                assert abs(pfba_fva_results[k][name]) <= abs(
                    fva_results[k][name])
        loop_reactions = [model.reactions.get_by_id(rid)
                          for rid in ("FRD7", "SUCDi")]
        fva_loopless = flux_variability_analysis(
            model, pfba_factor=1.1, reaction_list=loop_reactions,
            loopless=True)
        assert numpy.allclose(fva_loopless["maximum"],
                              fva_loopless["minimum"])

    @pytest.mark.parametrize("solver", all_solvers)
    def test_flux_variability(self, model, fva_results, solver):
        model.solver = solver
        fva_out = flux_variability_analysis(
            model, reaction_list=model.reactions)
        for name, result in iteritems(fva_out.T):
            for k, v in iteritems(result):
                assert abs(fva_results[k][name] - v) < 0.00001

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_flux_variability_loopless(self, model, solver):
        model.solver = solver
        loop_reactions = [model.reactions.get_by_id(rid)
                          for rid in ("FRD7", "SUCDi")]
        fva_normal = flux_variability_analysis(
            model, reaction_list=loop_reactions)
        fva_loopless = flux_variability_analysis(
            model, reaction_list=loop_reactions, loopless=True)

        assert not numpy.allclose(fva_normal["maximum"],
                                  fva_normal["minimum"])
        assert numpy.allclose(fva_loopless["maximum"],
                              fva_loopless["minimum"])

    def test_fva_data_frame(self, model):
        df = flux_variability_analysis(model)
        assert numpy.all([df.columns.values == ['maximum', 'minimum']])

    def test_fva_infeasible(self, model):
        infeasible_model = model.copy()
        infeasible_model.reactions.get_by_id("EX_glc__D_e").lower_bound = 0
        # ensure that an infeasible model does not run FVA
        with pytest.raises(Infeasible):
            flux_variability_analysis(infeasible_model)

    def test_fva_minimization(self, model):
        model.objective = model.reactions.EX_glc__D_e
        model.objective_direction = 'min'
        solution = flux_variability_analysis(model, fraction_of_optimum=.95)
        assert solution.at['EX_glc__D_e', 'minimum'] == -10.0
        assert solution.at['EX_glc__D_e', 'maximum'] == -9.5

    def test_find_blocked_reactions_solver_none(self, model):
        result = find_blocked_reactions(model, model.reactions[40:46])
        assert result == ['FRUpts2']

    def test_essential_genes(self, model):
        essential_genes = {'b2779', 'b1779', 'b0720', 'b2416',
                           'b2926', 'b1136', 'b2415'}
        observed_essential_genes = {g.id for g in
                                    find_essential_genes(model)}
        assert observed_essential_genes == essential_genes

    def test_essential_reactions(self, model):
        essential_reactions = {'GLNS', 'Biomass_Ecoli_core', 'PIt2r',
                               'GAPD',
                               'ACONTb', 'EX_nh4_e', 'ENO', 'EX_h_e',
                               'EX_glc__D_e', 'ICDHyr', 'CS', 'NH4t',
                               'GLCpts',
                               'PGM', 'EX_pi_e', 'PGK', 'RPI', 'ACONTa'}
        observed_essential_reactions = {r.id for r in
                                        find_essential_reactions(model)}
        assert observed_essential_reactions == essential_reactions

    @pytest.mark.parametrize("solver", all_solvers)
    def test_find_blocked_reactions(self, model, solver):
        model.solver = solver
        result = find_blocked_reactions(model, model.reactions[40:46])
        assert result == ['FRUpts2']

        result = find_blocked_reactions(model, model.reactions[42:48])
        assert set(result) == {'FUMt2_2', 'FRUpts2'}

        result = find_blocked_reactions(model, model.reactions[30:50],
                                        open_exchanges=True)
        assert result == []

    def test_loopless_benchmark_before(self, benchmark):
        test_model = construct_ll_test_model()

        def _():
            with test_model:
                add_loopless(test_model)
                test_model.optimize()

        benchmark(_)

    def test_loopless_benchmark_after(self, benchmark):
        test_model = construct_ll_test_model()
        benchmark(loopless_solution, test_model)

    def test_loopless_solution(self, ll_test_model):
        solution_feasible = loopless_solution(ll_test_model)
        ll_test_model.reactions.v3.lower_bound = 1
        ll_test_model.optimize()
        solution_infeasible = loopless_solution(ll_test_model)
        assert solution_feasible.fluxes["v3"] == 0.0
        assert solution_infeasible.fluxes["v3"] == 1.0

    def test_loopless_solution_fluxes(self, model):
        fluxes = model.optimize().fluxes
        ll_solution = loopless_solution(model, fluxes=fluxes)
        assert len(ll_solution.fluxes) == len(model.reactions)

    def test_add_loopless(self, ll_test_model):
        add_loopless(ll_test_model)
        feasible_status = ll_test_model.optimize().status
        ll_test_model.reactions.v3.lower_bound = 1
        ll_test_model.slim_optimize()
        infeasible_status = ll_test_model.solver.status
        assert feasible_status == OPTIMAL
        assert infeasible_status == INFEASIBLE

    def test_gapfilling(self, salmonella):
        m = Model()
        m.add_metabolites([Metabolite(m_id) for m_id in ["a", "b", "c"]])
        exa = Reaction("EX_a")
        exa.add_metabolites({m.metabolites.a: 1})
        b2c = Reaction("b2c")
        b2c.add_metabolites({m.metabolites.b: -1, m.metabolites.c: 1})
        dmc = Reaction("DM_c")
        dmc.add_metabolites({m.metabolites.c: -1})
        m.add_reactions([exa, b2c, dmc])
        m.objective = 'DM_c'

        universal = Model()
        a2b = Reaction("a2b")
        a2d = Reaction("a2d")
        universal.add_reactions([a2b, a2d])
        a2b.build_reaction_from_string("a --> b", verbose=False)
        a2d.build_reaction_from_string("a --> d", verbose=False)

        # # GrowMatch
        # result = gapfilling.growMatch(m, universal)[0]
        result = gapfilling.gapfill(m, universal)[0]
        assert len(result) == 1
        assert result[0].id == "a2b"

        # # SMILEY
        # result = gapfilling.SMILEY(m, "b", universal)[0]
        with m:
            m.objective = m.add_boundary(m.metabolites.b, type='demand')
            result = gapfilling.gapfill(m, universal)[0]
            assert len(result) == 1
            assert result[0].id == "a2b"

        # # 2 rounds of GrowMatch with exchange reactions
        # result = gapfilling.growMatch(m, None, ex_rxns=True, iterations=2)
        result = gapfilling.gapfill(m, None, exchange_reactions=True,
                                    iterations=2)
        assert len(result) == 2
        assert len(result[0]) == 1
        assert len(result[1]) == 1
        assert {i[0].id for i in result} == {"EX_b", "EX_c"}

        # somewhat bigger model
        universal = Model("universal_reactions")
        with salmonella as model:
            for i in [i.id for i in model.metabolites.f6p_c.reactions]:
                reaction = model.reactions.get_by_id(i)
                universal.add_reactions([reaction.copy()])
                model.remove_reactions([reaction])
            gf = gapfilling.GapFiller(model, universal,
                                      penalties={'TKT2': 1e3},
                                      demand_reactions=False)
            solution = gf.fill()
            assert 'TKT2' not in {r.id for r in solution[0]}
            assert gf.validate(solution[0])

    def check_line(self, output, expected_entries,
                   pattern=re.compile(r"\s")):
        """Ensure each expected entry is in the output."""
        output_set = set(
            pattern.sub("", line) for line in output.splitlines())
        for elem in expected_entries:
            assert pattern.sub("", elem) in output_set

    def check_in_line(self, output, expected_entries,
                      pattern=re.compile(r"\s")):
        """Ensure each expected entry is contained in the output."""
        output_strip = [pattern.sub("", line) for line in
                        output.splitlines()]
        for elem in expected_entries:
            assert any(
                pattern.sub("", elem) in line for line in output_strip), \
                "Not found: {} in:\n{}".format(pattern.sub("", elem),
                                               "\n".join(output_strip))

    def test_model_summary_previous_solution(self, model, opt_solver):
        model.solver = opt_solver
        solution = model.optimize()
        rxn_test = model.exchanges[0]
        met_test = list(rxn_test.metabolites.keys())[0].id

        solution.fluxes[rxn_test.id] = 321

        with captured_output() as (out, err):
            model.summary(solution)
        self.check_in_line(out.getvalue(), [met_test + '321'])

    def test_model_summary(self, model, opt_solver):
        model.solver = opt_solver
        # test non-fva version (these should be fixed for textbook model
        expected_entries = [
            'o2_e      21.8',
            'glc__D_e  10',
            'nh4_e      4.77',
            'pi_e       3.21',
            'h2o_e  29.2',
            'co2_e  22.8',
            'h_e    17.5',
            'Biomass_Ecol...  0.874',
        ]
        # Need to use a different method here because
        # there are multiple entries per line.
        model.optimize()
        with captured_output() as (out, err):
            model.summary()
        self.check_in_line(out.getvalue(), expected_entries)

        with model:
            model.objective = model.exchanges[0]
            model.summary()

    @pytest.mark.parametrize("fraction", [0.95])
    def test_model_summary_with_fva(self, model, opt_solver, fraction):
        if opt_solver == "optlang-gurobi":
            pytest.xfail("FVA currently buggy")
        # test non-fva version (these should be fixed for textbook model
        expected_entries = [
            'idFluxRangeidFluxRangeBiomass_Ecol...0.874',
            'o2_e       21.8   [19.9, 23.7]'
            'h2o_e       29.2  [25, 30.7]',
            'glc__D_e   10     [9.52, 10]'
            'co2_e       22.8  [18.9, 24.7]',
            'nh4_e       4.77  [4.53, 5.16]'
            'h_e         17.5  [16.7, 22.4]',
            'pi_e        3.21  [3.05, 3.21]'
            'for_e        0    [0, 5.72]',
            'ac_e         0    [0, 1.91]',
            'pyr_e        0    [0, 1.27]',
            'lac__D_e     0    [0, 1.07]',
            'succ_e       0    [0, 0.837]',
            'glu__L_e     0    [0, 0.636]',
            'akg_e        0    [0, 0.715]',
            'etoh_e       0    [0, 1.11]',
            'acald_e      0    [0, 1.27]',
        ]
        # Need to use a different method here because
        # there are multiple entries per line.
        model.solver = opt_solver
        solution = model.optimize()
        with captured_output() as (out, err):
            model.summary(solution, fva=fraction)
        self.check_in_line(out.getvalue(), expected_entries)

    @pytest.mark.parametrize("met", ["q8_c"])
    def test_metabolite_summary_previous_solution(
            self, model, opt_solver, met):
        model.solver = opt_solver
        solution = pfba(model)
        model.metabolites.get_by_id(met).summary(solution)

    @pytest.mark.parametrize("met", ["q8_c"])
    def test_metabolite_summary(self, model, opt_solver, met):
        model.solver = opt_solver
        model.optimize()
        with captured_output() as (out, err):
            model.metabolites.get_by_id(met).summary()

        expected_entries = [
            'PRODUCING REACTIONS -- Ubiquinone-8 (q8_c)',
            '%       FLUX  RXN ID    REACTION',
            '100%   43.6   CYTBD     '
            '2.0 h_c + 0.5 o2_c + q8h2_c --> h2o_c + 2.0 h_e...',
            'CONSUMING REACTIONS -- Ubiquinone-8 (q8_c)',
            '%       FLUX  RXN ID    REACTION',
            '88%    38.5   NADH16    '
            '4.0 h_c + nadh_c + q8_c --> 3.0 h_e + nad_c + q...',
            '12%     5.06  SUCDi     q8_c + succ_c --> fum_c + q8h2_c',
        ]

        self.check_in_line(out.getvalue(), expected_entries)

    @pytest.mark.parametrize("fraction, met", [(0.99, "fdp_c")])
    def test_metabolite_summary_with_fva(self, model, opt_solver, fraction,
                                         met):
        if opt_solver in (
                "optlang-glpk", "optlang-cplex", "optlang-gurobi"):
            pytest.xfail("FVA currently buggy")

        model.solver = opt_solver
        model.optimize()
        with captured_output() as (out, err):
            model.metabolites.get_by_id(met).summary(fva=fraction)

        expected_entries = [
            'PRODUCING REACTIONS -- D-Fructose 1,6-bisphosphate (fdp_c)',
            '%       FLUX  RANGE         RXN ID    REACTION',
            '100%    7.48  [6.17, 9.26]  PFK       '
            'atp_c + f6p_c --> adp_c + fdp_c + h_c',
            'CONSUMING REACTIONS -- D-Fructose 1,6-bisphosphate (fdp_c)',
            '%       FLUX  RANGE         RXN ID    REACTION',
            '100%    7.48  [6.17, 8.92]  FBA       fdp_c <=> dhap_c + g3p_c',
            '0%      0     [0, 1.72]     FBP       '
            'fdp_c + h2o_c --> f6p_c + pi_c',
        ]

        self.check_line(out.getvalue(), expected_entries)


class TestCobraFluxSampling:
    """Tests and benchmark flux sampling."""

    def test_single_achr(self, model):
        s = sample(model, 10, method="achr")
        assert s.shape == (10, len(model.reactions))

    def test_single_optgp(self, model):
        s = sample(model, 10, processes=1)
        assert s.shape == (10, len(model.reactions))

    def test_multi_optgp(self, model):
        s = sample(model, 10, processes=2)
        assert s.shape == (10, len(model.reactions))

    def test_wrong_method(self, model):
        with pytest.raises(ValueError):
            sample(model, 1, method="schwupdiwupp")

    def test_validate_wrong_sample(self, model):
        s = self.achr.sample(10)
        s["hello"] = 1
        with pytest.raises(ValueError):
            self.achr.validate(s)

    def test_fixed_seed(self, model):
        s = sample(model, 1, seed=42)
        assert numpy.allclose(s.TPI[0], 9.12037487)

    def test_equality_constraint(self, model):
        model.reactions.ACALD.bounds = (-1.5, -1.5)
        s = sample(model, 10)
        assert numpy.allclose(s.ACALD, -1.5, atol=1e-6, rtol=0)
        s = sample(model, 10, method="achr")
        assert numpy.allclose(s.ACALD, -1.5, atol=1e-6, rtol=0)

    def test_inequality_constraint(self, model):
        co = model.problem.Constraint(
            model.reactions.ACALD.flux_expression, lb=-0.5)
        model.add_cons_vars(co)
        s = sample(model, 10)
        assert all(s.ACALD > -0.5 - 1e-6)
        s = sample(model, 10, method="achr")
        assert all(s.ACALD > -0.5 - 1e-6)

    def setup_class(self):
        from . import create_test_model
        model = create_test_model("textbook")
        achr = ACHRSampler(model, thinning=1)
        assert ((achr.n_warmup > 0) and
                (achr.n_warmup <= 2 * len(model.variables)))
        assert all(achr.validate(achr.warmup) == "v")
        self.achr = achr

        optgp = OptGPSampler(model, processes=1, thinning=1)
        assert ((optgp.n_warmup > 0) and
                (optgp.n_warmup <= 2 * len(model.variables)))
        assert all(optgp.validate(optgp.warmup) == "v")
        self.optgp = optgp

    def test_achr_init_benchmark(self, model, benchmark):
        benchmark(lambda: ACHRSampler(model))

    def test_optgp_init_benchmark(self, model, benchmark):
        benchmark(lambda: OptGPSampler(model, processes=2))

    def test_sampling(self):
        s = self.achr.sample(10)
        assert all(self.achr.validate(s) == "v")

        s = self.optgp.sample(10)
        assert all(self.optgp.validate(s) == "v")

    def test_achr_sample_benchmark(self, benchmark):
        benchmark(self.achr.sample, 1)

    def test_optgp_sample_benchmark(self, benchmark):
        benchmark(self.optgp.sample, 1)

    def test_batch_sampling(self):
        for b in self.achr.batch(5, 4):
            assert all(self.achr.validate(b) == "v")

        for b in self.optgp.batch(5, 4):
            assert all(self.optgp.validate(b) == "v")

    def test_variables_samples(self):
        vnames = numpy.array([v.name for v in self.achr.model.variables])
        s = self.achr.sample(10, fluxes=False)
        assert s.shape == (10, self.achr.warmup.shape[1])
        assert (s.columns == vnames).all()
        assert (self.achr.validate(s) == "v").all()
        s = self.optgp.sample(10, fluxes=False)
        assert s.shape == (10, self.optgp.warmup.shape[1])
        assert (s.columns == vnames).all()
        assert (self.optgp.validate(s) == "v").all()

    def test_inhomogeneous_sanity(self, model):
        """Test whether inhomogeneous sampling gives approximately the same
           standard deviation as a homogeneous version."""
        model.reactions.ACALD.bounds = (-1.5, -1.5)
        s_inhom = sample(model, 64)
        model.reactions.ACALD.bounds = (-1.5 - 1e-3, -1.5 + 1e-3)
        s_hom = sample(model, 64)
        relative_diff = (s_inhom.std() + 1e-12) / (s_hom.std() + 1e-12)
        assert 0.5 < relative_diff.abs().mean() < 2

        model.reactions.ACALD.bounds = (-1.5, -1.5)
        s_inhom = sample(model, 64, method="achr")
        model.reactions.ACALD.bounds = (-1.5 - 1e-3, -1.5 + 1e-3)
        s_hom = sample(model, 64, method="achr")
        relative_diff = (s_inhom.std() + 1e-12) / (s_hom.std() + 1e-12)
        assert 0.5 < relative_diff.abs().mean() < 2

    def test_reproject(self):
        s = self.optgp.sample(10, fluxes=False).as_matrix()
        proj = numpy.apply_along_axis(self.optgp._reproject, 1, s)
        assert all(self.optgp.validate(proj) == "v")
        s = numpy.random.rand(10, self.optgp.warmup.shape[1])
        proj = numpy.apply_along_axis(self.optgp._reproject, 1, s)
        assert all(self.optgp.validate(proj) == "v")

    def test_complicated_model(self):
        """Difficult model since the online mean calculation is numerically
        unstable so many samples weakly violate the equality constraints."""
        model = Model('flux_split')
        reaction1 = Reaction('V1')
        reaction2 = Reaction('V2')
        reaction3 = Reaction('V3')
        reaction1.lower_bound = 0
        reaction2.lower_bound = 0
        reaction3.lower_bound = 0
        reaction1.upper_bound = 6
        reaction2.upper_bound = 8
        reaction3.upper_bound = 10
        A = Metabolite('A')
        reaction1.add_metabolites({A: -1})
        reaction2.add_metabolites({A: -1})
        reaction3.add_metabolites({A: 1})
        model.add_reactions([reaction1])
        model.add_reactions([reaction2])
        model.add_reactions([reaction3])

        optgp = OptGPSampler(model, 1, seed=42)
        achr = ACHRSampler(model, seed=42)
        optgp_samples = optgp.sample(100)
        achr_samples = achr.sample(100)
        assert any(optgp_samples.corr().abs() < 1.0)
        assert any(achr_samples.corr().abs() < 1.0)
        # > 95% are valid
        assert(sum(optgp.validate(optgp_samples) == "v") > 95)
        assert(sum(achr.validate(achr_samples) == "v") > 95)

    def test_single_point_space(self, model):
        """Model where constraints reduce the sampling space to one point."""
        pfba_sol = pfba(model)
        pfba_const = model.problem.Constraint(
            sum(model.variables), ub=pfba_sol.objective_value)
        model.add_cons_vars(pfba_const)
        model.reactions.Biomass_Ecoli_core.lower_bound = \
            pfba_sol.fluxes.Biomass_Ecoli_core
        with pytest.raises(ValueError):
            s = sample(model, 1)


class TestProductionEnvelope:
    """Test the production envelope."""

    def test_envelope_one(self, model):
        df = production_envelope(model, ["EX_o2_e"])
        assert numpy.isclose(df["flux_maximum"].sum(), 9.342, atol=1e-3)

    def test_envelope_multi_reaction_objective(self, model):
        obj = {model.reactions.EX_ac_e: 1,
               model.reactions.EX_co2_e: 1}
        with pytest.raises(ValueError):
            production_envelope(model, "EX_o2_e", obj)

    @pytest.mark.parametrize("variables, num", [
        (["EX_glc__D_e"], 30),
        (["EX_glc__D_e", "EX_o2_e"], 20),
        (["EX_glc__D_e", "EX_o2_e", "EX_ac_e"], 10)
    ])
    def test_multi_variable_envelope(self, model, variables, num):
        df = production_envelope(model, variables, points=num)
        assert len(df) == num ** len(variables)

    def test_envelope_two(self, model):
        df = production_envelope(model, ["EX_glc__D_e", "EX_o2_e"],
                                 objective="EX_ac_e")
        assert numpy.isclose(df["flux_maximum"].sum(), 1737.466, atol=1e-3)
        assert numpy.isclose(df["carbon_yield_maximum"].sum(), 83.579,
                             atol=1e-3)
        assert numpy.isclose(df["mass_yield_maximum"].sum(), 82.176,
                             atol=1e-3)


class TestReactionUtils:
    """Test the assess_ functions in reactions.py."""

    @pytest.mark.parametrize("solver", all_solvers)
    def test_assess(self, model, solver):
        with model:
            assert assess(model, model.reactions.GLCpts,
                          solver=solver) is True
            pyr = model.metabolites.pyr_c
            a = Metabolite('a')
            b = Metabolite('b')
            model.add_metabolites([a, b])
            pyr_a2b = Reaction('pyr_a2b')
            pyr_a2b.add_metabolites({pyr: -1, a: -1, b: 1})
            model.add_reactions([pyr_a2b])
            res = assess(model, pyr_a2b, 0.01, solver=solver)
            expected = {
                'precursors': {a: {'required': 0.01, 'produced': 0.0}},
                'products': {b: {'required': 0.01, 'capacity': 0.0}}}
            assert res == expected
