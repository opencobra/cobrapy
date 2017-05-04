# -*- coding: utf-8 -*-
from __future__ import absolute_import

import re
import sys
import warnings
from contextlib import contextmanager
from os import name

import pytest
import numpy
from optlang.interface import OPTIMAL, INFEASIBLE
from six import StringIO, iteritems

import cobra.util.solver as sutil
from cobra.core import Metabolite, Model, Reaction
from cobra.flux_analysis import *
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.flux_analysis.sampling import ACHRSampler, OptGPSampler
from cobra.manipulation import convert_to_irreversible
from cobra.solvers import SolverNotFound, get_solver_name, solver_dict
from cobra.exceptions import Infeasible

try:
    import scipy
    from cobra.flux_analysis.moma import add_moma
except ImportError:
    scipy = None
    add_moma = None
try:
    import matplotlib
except (ImportError, RuntimeError):
    matplotlib = None
try:
    from matplotlib import pyplot
    from mpl_toolkits.mplot3d import axes3d
except (ImportError, RuntimeError):
    pyplot = None
    axes3d = None

# The scipt interface is currently unstable and may yield errors or infeasible
# solutions
stable_optlang = ["glpk", "cplex", "gurobi"]
optlang_solvers = ["optlang-" + s for s in stable_optlang if s in
                   sutil.solvers]
all_solvers = optlang_solvers + list(solver_dict)


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
    """ A context manager to test the IO summary methods """
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class TestCobraFluxAnalysis:
    """Test the simulation functions in cobra.flux_analysis"""

    @pytest.mark.parametrize("solver", all_solvers)
    def test_pfba_benchmark(self, large_model, benchmark, solver):
        convert_to_irreversible(large_model)

        def do_pfba(solver):
            pfba(large_model, solver=solver,
                 already_irreversible=True)

        benchmark(do_pfba, solver)

    @pytest.mark.parametrize("solver", all_solvers)
    def test_pfba(self, model, solver):
        with model:
            add_pfba(model)
            with pytest.raises(ValueError):
                add_pfba(model)

        if solver in optlang_solvers:
            model.solver = solver
        expression = model.objective.expression
        n_constraints = len(model.constraints)
        solution = pfba(model, solver=solver)
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
        solution = pfba(model, solver=solver,
                        fraction_of_optimum=0.95)
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
                pfba(model, solver=solver)

    @pytest.mark.parametrize("solver", all_solvers)
    def test_single_gene_deletion_fba_benchmark(self, model, benchmark,
                                                solver):
        benchmark(single_gene_deletion, model, solver=solver)

    @pytest.mark.parametrize("solver", all_solvers)
    def test_single_gene_deletion_fba(self, model, solver):
        # expected knockouts for textbook model
        growth_dict = {"b0008": 0.87, "b0114": 0.80, "b0116": 0.78,
                       "b2276": 0.21, "b1779": 0.00}
        df = single_gene_deletion(model, gene_list=growth_dict.keys(),
                                  method="fba", solver=solver)
        assert numpy.all([df.status == 'optimal'])
        assert all(abs(df.flux[gene] - expected) < 0.01 for
                   gene, expected in iteritems(growth_dict))

    @pytest.mark.skipif(scipy is None,
                        reason="moma gene deletion requires scipy")
    def test_single_gene_deletion_moma_benchmark(self, model, benchmark):
        try:
            sutil.get_solver_name(qp=True)
        except sutil.SolverNotFound:
            pytest.skip("no qp support")
        genes = ['b0008', 'b0114', 'b2276', 'b1779']
        benchmark(single_gene_deletion, model, gene_list=genes, method="moma")

    @pytest.mark.skipif(scipy is None,
                        reason="moma gene deletion requires scipy")
    def test_single_gene_deletion_moma(self, model):
        try:
            sutil.get_solver_name(qp=True)
        except sutil.SolverNotFound:
            pytest.skip("no qp support")

        # expected knockouts for textbook model
        growth_dict = {"b0008": 0.87, "b0114": 0.71, "b0116": 0.56,
                       "b2276": 0.11, "b1779": 0.00}

        df = single_gene_deletion(model, gene_list=growth_dict.keys(),
                                  method="moma")
        assert numpy.all([df.status == 'optimal'])
        assert all(abs(df.flux[gene] - expected) < 0.01
                   for gene, expected in iteritems(growth_dict))
        with model:
            add_moma(model)
            with pytest.raises(ValueError):
                add_moma(model)

    @pytest.mark.parametrize("solver", all_solvers)
    def test_single_gene_deletion_benchmark(self, model, benchmark,
                                            solver):
        benchmark(single_reaction_deletion, model, solver=solver)

    @pytest.mark.parametrize("solver", all_solvers)
    def test_single_reaction_deletion(self, model, solver):
        expected_results = {'FBA': 0.70404, 'FBP': 0.87392, 'CS': 0,
                            'FUM': 0.81430, 'GAPD': 0, 'GLUDy': 0.85139}

        df = single_reaction_deletion(
            model, reaction_list=expected_results.keys(), solver=solver)
        assert len(df) == 6
        assert numpy.all([df.status == 'optimal'])
        assert all(abs(df.flux[gene] - expected) < 0.00001 for
                   gene, expected in iteritems(expected_results))

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
        genes = ["b0726", "b4025", "b0724", "b0720", "b2935", "b2935", "b1276",
                 "b1241"]
        benchmark(double_gene_deletion, large_model, gene_list1=genes)

    def test_double_gene_deletion(self, model):
        genes = ["b0726", "b4025", "b0724", "b0720", "b2935", "b2935", "b1276",
                 "b1241"]
        growth_list = [
            [0.858, 0.857, 0.814, 0.000, 0.858, 0.858, 0.858, 0.858],
            [0.857, 0.863, 0.739, 0.000, 0.863, 0.863, 0.863, 0.863],
            [0.814, 0.739, 0.814, 0.000, 0.814, 0.814, 0.814, 0.814],
            [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
            [0.858, 0.863, 0.814, 0.000, 0.874, 0.874, 0.874, 0.874],
            [0.858, 0.863, 0.814, 0.000, 0.874, 0.874, 0.874, 0.874],
            [0.858, 0.863, 0.814, 0.000, 0.874, 0.874, 0.874, 0.874],
            [0.858, 0.863, 0.814, 0.000, 0.874, 0.874, 0.874, 0.874]]
        opts = {"number_of_processes": 1} if name == "nt" else {}
        solution = double_gene_deletion(model, gene_list1=genes, **opts)
        assert solution["x"] == genes
        assert solution["y"] == genes
        self.compare_matrices(growth_list, solution["data"])
        # test when lists differ slightly
        solution = double_gene_deletion(model, gene_list1=genes[:-1],
                                        gene_list2=genes,
                                        number_of_processes=1)
        assert solution["x"] == genes[:-1]
        assert solution["y"] == genes
        self.compare_matrices(growth_list[:-1], solution["data"])

    def test_double_reaction_deletion(self, model):
        reactions = ['FBA', 'ATPS4r', 'ENO', 'FRUpts2']
        growth_list = [[0.704, 0.135, 0.000, 0.704],
                       [0.135, 0.374, 0.000, 0.374],
                       [0.000, 0.000, 0.000, 0.000],
                       [0.704, 0.374, 0.000, 0.874]]

        solution = double_reaction_deletion(model,
                                            reaction_list1=reactions,
                                            number_of_processes=1)
        assert solution["x"] == reactions
        assert solution["y"] == reactions
        self.compare_matrices(growth_list, solution["data"])

    @pytest.mark.parametrize("solver", all_solvers)
    def test_flux_variability_benchmark(self, large_model, benchmark, solver):
        benchmark(flux_variability_analysis, large_model, solver=solver,
                  reaction_list=large_model.reactions[1::3])

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_flux_variability_loopless_benchmark(self, model, benchmark,
                                                 solver):
        benchmark(flux_variability_analysis, model, loopless=True,
                  solver=solver, reaction_list=model.reactions[1::3])

    @pytest.mark.parametrize("solver", all_solvers)
    def test_flux_variability(self, model, fva_results, solver):
        if solver == "esolver":
            pytest.skip("esolver too slow...")
        fva_out = flux_variability_analysis(
            model, solver=solver, reaction_list=model.reactions)
        for name, result in iteritems(fva_out.T):
            for k, v in iteritems(result):
                assert abs(fva_results[k][name] - v) < 0.00001

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_flux_variability_loopless(self, model, fva_results, solver):
        loop_reactions = [model.reactions.get_by_id(rid)
                          for rid in ("FRD7", "SUCDi")]
        fva_normal = flux_variability_analysis(
            model, solver=solver, reaction_list=loop_reactions)
        fva_loopless = flux_variability_analysis(
            model, solver=solver, reaction_list=loop_reactions, loopless=True)

        assert not numpy.allclose(fva_normal["maximum"], fva_normal["minimum"])
        assert numpy.allclose(fva_loopless["maximum"], fva_loopless["minimum"])

    def test_fva_data_frame(self, model):
        df = flux_variability_analysis(model, return_frame=True)
        assert numpy.all([df.columns.values == ['maximum', 'minimum']])

    def test_fva_infeasible(self, model):
        infeasible_model = model.copy()
        infeasible_model.reactions.get_by_id("EX_glc__D_e").lower_bound = 0
        # ensure that an infeasible model does not run FVA
        with pytest.raises(Infeasible):
            flux_variability_analysis(infeasible_model)

    def test_find_blocked_reactions_solver_none(self, model):
        result = find_blocked_reactions(model, model.reactions[40:46])
        assert result == ['FRUpts2']

    @pytest.mark.parametrize("solver", all_solvers)
    def test_find_blocked_reactions(self, model, solver):
        result = find_blocked_reactions(model, model.reactions[40:46],
                                        solver=solver)
        assert result == ['FRUpts2']

        result = find_blocked_reactions(model, model.reactions[42:48],
                                        solver=solver)
        assert set(result) == {'FUMt2_2', 'FRUpts2'}

        result = find_blocked_reactions(model, model.reactions[30:50],
                                        solver=solver,
                                        open_exchanges=True)
        assert result == []

    def test_legacy_loopless_benchmark(self, benchmark):
        test_model = construct_ll_test_model()
        benchmark(lambda: construct_loopless_model(test_model).optimize(
            solver="cglpk"))

    def test_loopless_benchmark_before(self, benchmark):
        test_model = construct_ll_test_model()

        def _():
            with test_model:
                add_loopless(test_model)
                test_model.optimize(solver="optlang-glpk")

        benchmark(_)

    def test_loopless_benchmark_after(self, benchmark):
        test_model = construct_ll_test_model()
        benchmark(loopless_solution, test_model)

    def test_legacy_loopless(self):
        try:
            get_solver_name(mip=True)
        except SolverNotFound:
            pytest.skip("no MILP solver found")
        test_model = construct_ll_test_model()
        feasible_sol = construct_loopless_model(test_model).optimize(
            solver="cglpk")
        assert feasible_sol.status == "optimal"
        test_model.reactions.v3.lower_bound = 1
        infeasible_mod = construct_loopless_model(test_model)

        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            with pytest.raises(UserWarning):
                infeasible_mod.optimize(solver="cglpk")

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
        fluxes["Biomass_Ecoli_core"] = 1
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            with pytest.raises(UserWarning):
                loopless_solution(model, fluxes=fluxes)

    def test_add_loopless(self, ll_test_model):
        add_loopless(ll_test_model)
        feasible_status = ll_test_model.optimize().status
        ll_test_model.reactions.v3.lower_bound = 1
        ll_test_model.solver.optimize()
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

    def test_phenotype_phase_plane_benchmark(self, model, benchmark):
        benchmark(calculate_phenotype_phase_plane,
                  model, "EX_glc__D_e", "EX_o2_e",
                  reaction1_npoints=20, reaction2_npoints=20)

    def test_phenotype_phase_plane(self, model):
        data = calculate_phenotype_phase_plane(
            model, "EX_glc__D_e", "EX_o2_e",
            reaction1_npoints=20, reaction2_npoints=20)
        assert data.growth_rates.shape == (20, 20)
        assert abs(data.growth_rates.max() - 1.20898) < 0.0001
        assert abs(data.growth_rates[0, :].max()) < 0.0001
        if pyplot is None or axes3d is None:
            pytest.skip("can't test plots without 3D plotting")
        data.plot()

    def check_line(self, output, expected_entries, pattern=re.compile(r"\s")):
        """Ensure each expected entry is in the output."""
        output_set = set(pattern.sub("", line) for line in output.splitlines())
        for elem in expected_entries:
            assert pattern.sub("", elem) in output_set

    def check_in_line(self, output, expected_entries,
                      pattern=re.compile(r"\s")):
        """Ensure each expected entry is contained in the output."""
        output_strip = [pattern.sub("", line) for line in output.splitlines()]
        for elem in expected_entries:
            assert any(
                pattern.sub("", elem) in line for line in output_strip), \
                "Not found: {} in:\n{}".format(pattern.sub("", elem),
                                               "\n".join(output_strip))

    def test_model_summary_unoptimized(self, model, opt_solver):
        model.solver = opt_solver
        with pytest.raises(RuntimeError):
            model.summary()

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
        model.optimize()
        with captured_output() as (out, err):
            model.summary(fva=fraction)
        self.check_in_line(out.getvalue(), expected_entries)

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
        if opt_solver in ("optlang-glpk", "optlang-cplex", "optlang-gurobi"):
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
    """Test and benchmark flux sampling"""

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
        assert numpy.allclose(s.TPI[0], 8.94197599534)

    def test_equality_constraint(self, model):
        model.reactions.ACALD.bounds = (-1.5, -1.5)
        s = sample(model, 10)
        assert numpy.allclose(s.ACALD, -1.5, atol=1e-6, rtol=0)

    def test_inequality_constraint(self, model):
        co = model.problem.Constraint(
            model.reactions.ACALD.flux_expression, lb=-0.5)
        model.add_cons_vars(co)
        s = sample(model, 10)
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


class TestProductionEnvelope:
    """Test the production envelope"""

    def test_envelope_one(self, model):
        df = production_envelope(model, ["EX_o2_e"])
        assert abs(sum(df.flux) - 9.342) < 0.001

    def test_envelope_two(self, model):
        df = production_envelope(model, ["EX_glc__D_e", "EX_o2_e"],
                                 objective="EX_ac_e",
                                 c_source="EX_glc__D_e")
        assert abs(numpy.sum(df.carbon_yield) - 83.579) < 0.001
        assert abs(numpy.sum(df.flux) - 1737.466) < 0.001
        assert abs(numpy.sum(df.mass_yield) - 82.176) < 0.001
