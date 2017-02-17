import pytest
import sys
from os import name
from contextlib import contextmanager
import re
from six import iteritems, StringIO
from cobra.core import Model, Reaction, Metabolite, Solution
from cobra.core.solution import LegacySolution
from cobra.solvers import solver_dict, get_solver_name
import cobra.util.solver as sutil
from cobra.flux_analysis import *
from cobra.solvers import SolverNotFound
from .conftest import model, large_model, solved_model, fva_results
from cobra.manipulation import convert_to_irreversible
from cobra.exceptions import SolveError

try:
    import numpy
    from cobra.flux_analysis.sampling import ARCHSampler, OptGPSampler
except ImportError:
    numpy = None
try:
    import scipy
except ImportError:
    scipy = None
try:
    import matplotlib
except ImportError:
    matplotlib = None
try:
    import pandas
except ImportError:
    pandas = None
try:
    import tabulate
except ImportError:
    tabulate = None

# The scipt interface is currently unstable and may yield errors or infeasible
# solutions
stable_optlang = ["glpk", "cplex", "gurobi"]
all_solvers = ["optlang-" + s for s in stable_optlang if s in
               sutil.solvers] + list(solver_dict)


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
            optimize_minimal_flux(large_model, solver=solver,
                                  already_irreversible=True)

        benchmark(do_pfba, solver)

    @pytest.mark.parametrize("solver", all_solvers)
    def test_pfba(self, model, solver):
        expression = model.objective.expression
        n_constraints = len(model.solver.constraints)
        optimize_minimal_flux(model, solver=solver)
        abs_x = [abs(i) for i in model.solution.x]
        assert model.solution.status == "optimal"
        if isinstance(model.solution, LegacySolution):
            assert abs(model.solution.f - 0.8739) < 0.001
        else:
            assert abs(
                model.solution.fluxes['Biomass_Ecoli_core'] - 0.8739) < 0.001
        assert abs(sum(abs_x) - 518.4221) < 0.001

        # test changes to model reverted
        assert expression == model.objective.expression
        assert len(model.solver.constraints) == n_constraints

        # needed?
        # Test desired_objective_value
        # desired_objective = 0.8
        # optimize_minimal_flux(model, solver=solver,
        #                       desired_objective_value=desired_objective)
        # abs_x = [abs(i) for i in model.solution.x]
        # assert model.solution.status == "optimal"
        # assert abs(model.solution.f - desired_objective) < 0.001
        # assert abs(sum(abs_x) - 476.1594) < 0.001

        # Test fraction_of_optimum
        optimize_minimal_flux(model, solver=solver,
                              fraction_of_optimum=0.95)
        abs_x = [abs(i) for i in model.solution.x]
        assert model.solution.status == "optimal"
        if isinstance(model.solution, LegacySolution):
            assert abs(model.solution.f - 0.95 * 0.8739) < 0.001
        else:
            assert abs(
                model.solution.fluxes[
                    'Biomass_Ecoli_core'] - 0.95 * 0.8739) < 0.001
        assert abs(sum(abs_x) - 493.4400) < 0.001

        # Infeasible solution
        model.reactions.ATPM.lower_bound = 500
        with pytest.raises((SolveError, ValueError)):
            optimize_minimal_flux(model, solver=solver)

    @pytest.mark.parametrize("solver", all_solvers)
    def test_single_gene_deletion_fba_benchmark(self, model, benchmark,
                                                solver):
        benchmark(single_gene_deletion, model, solver=solver)

    @pytest.mark.parametrize("solver", all_solvers)
    def test_single_gene_deletion_fba(self, model, solver):
        # expected knockouts for textbook model
        growth_dict = {"b0008": 0.87, "b0114": 0.80, "b0116": 0.78,
                       "b2276": 0.21, "b1779": 0.00}
        rates, statuses = single_gene_deletion(model,
                                               gene_list=growth_dict.keys(),
                                               method="fba",
                                               solver=solver)
        for gene, expected_value in iteritems(growth_dict):
            assert statuses[gene] == 'optimal'
            assert abs(rates[gene] - expected_value) < 0.01

    def test_single_gene_deletion_moma_benchmark(self, model, benchmark):
        try:
            sutil.get_solver_name(qp=True)
        except sutil.SolverNotFound:
            pytest.skip("no qp support")
        genes = ['b0008', 'b0114', 'b2276', 'b1779']
        benchmark(single_gene_deletion, model, gene_list=genes,
                  method="moma")

    def test_single_gene_deletion_moma(self, model):
        try:
            sutil.get_solver_name(qp=True)
        except sutil.SolverNotFound:
            pytest.skip("no qp support")

        # expected knockouts for textbook model
        growth_dict = {"b0008": 0.87, "b0114": 0.71, "b0116": 0.56,
                       "b2276": 0.11, "b1779": 0.00}

        rates, statuses = single_gene_deletion(model,
                                               gene_list=growth_dict.keys(),
                                               method="moma")
        for gene, expected_value in iteritems(growth_dict):
            assert statuses[gene] == 'optimal'
            assert abs(rates[gene] - expected_value) < 0.01

    @pytest.mark.parametrize("solver", all_solvers)
    def test_single_gene_deletion_benchmark(self, model, benchmark,
                                            solver):
        benchmark(single_reaction_deletion, model, solver=solver)

    @pytest.mark.parametrize("solver", all_solvers)
    def test_single_reaction_deletion(self, model, solver):
        expected_results = {'FBA': 0.70404, 'FBP': 0.87392, 'CS': 0,
                            'FUM': 0.81430, 'GAPD': 0, 'GLUDy': 0.85139}

        results, status = single_reaction_deletion(
            model, reaction_list=expected_results.keys(), solver=solver)
        assert len(results) == 6
        assert len(status) == 6
        for status_value in status.values():
            assert status_value == "optimal"
        for reaction, value in results.items():
            assert abs(value - expected_results[reaction]) < 0.00001

    @classmethod
    def compare_matrices(cls, matrix1, matrix2, places=3):
        nrows = len(matrix1)
        ncols = len(matrix1[0])
        assert nrows == len(matrix2)
        assert ncols == len(matrix2[0])
        for i in range(nrows):
            for j in range(ncols):
                assert abs(matrix1[i][j] - matrix2[i][j]) < 10 ** -places

    @pytest.mark.skipif(numpy is None, reason="double deletions require numpy")
    def test_double_gene_deletion_benchmark(self, large_model, benchmark):
        genes = ["b0726", "b4025", "b0724", "b0720", "b2935", "b2935", "b1276",
                 "b1241"]
        benchmark(double_gene_deletion, large_model, gene_list1=genes)

    @pytest.mark.skipif(numpy is None, reason="double deletions require numpy")
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

    @pytest.mark.skipif(numpy is None, reason="double deletions require numpy")
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

    @pytest.mark.parametrize("solver", all_solvers)
    def test_flux_variability(self, model, fva_results, solver):
        if solver == "esolver":
            pytest.skip("esolver too slow...")
        fva_out = flux_variability_analysis(
            model, solver=solver, reaction_list=model.reactions[1::3])
        for name, result in iteritems(fva_out):
            for k, v in iteritems(result):
                assert abs(fva_results[name][k] - v) < 0.00001

    def test_fva_infeasible(self, model):
        infeasible_model = model.copy()
        infeasible_model.reactions.get_by_id("EX_glc__D_e").lower_bound = 0
        # ensure that an infeasible model does not run FVA
        with pytest.raises(ValueError):
            flux_variability_analysis(infeasible_model)

    def test_find_blocked_reactions(self, model):
        result = find_blocked_reactions(model, model.reactions[40:46])
        assert result == ['FRUpts2']

        result = find_blocked_reactions(model, model.reactions[42:48])
        assert set(result) == {'FUMt2_2', 'FRUpts2'}

        result = find_blocked_reactions(model, model.reactions[30:50],
                                        open_exchanges=True)
        assert result == []

    @classmethod
    def construct_ll_test_model(cls):
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

    def test_loopless_benchmark(self, benchmark):
        test_model = self.construct_ll_test_model()
        benchmark(lambda: construct_loopless_model(test_model).optimize())

    def test_loopless(self):
        try:
            get_solver_name(mip=True)
        except SolverNotFound:
            pytest.skip("no MILP solver found")
        test_model = self.construct_ll_test_model()
        feasible_sol = construct_loopless_model(test_model).optimize()
        test_model.reactions.get_by_id('v3').lower_bound = 1
        infeasible_mod = construct_loopless_model(test_model)
        infeasible_mod.solver.optimize()
        assert feasible_sol.status == "optimal"
        assert infeasible_mod.solver.status == "infeasible"

    @pytest.mark.xfail(reason='broken / needs re-implementation')
    def test_gapfilling(self):
        try:
            get_solver_name(mip=True)
        except SolverNotFound:
            pytest.skip("no MILP solver found")
        m = Model()
        m.add_metabolites(map(Metabolite, ["a", "b", "c"]))
        r = Reaction("EX_A")
        m.add_reaction(r)
        r.add_metabolites({m.metabolites.a: 1})
        r = Reaction("r1")
        m.add_reaction(r)
        r.add_metabolites({m.metabolites.b: -1, m.metabolites.c: 1})
        r = Reaction("DM_C")
        m.add_reaction(r)
        r.add_metabolites({m.metabolites.c: -1})
        r.objective_coefficient = 1
        U = Model()
        r = Reaction("a2b")
        U.add_reaction(r)

        r.build_reaction_from_string("a --> b", verbose=False)
        r = Reaction("a2d")
        U.add_reaction(r)
        r.build_reaction_from_string("a --> d", verbose=False)

        # GrowMatch
        result = gapfilling.growMatch(m, U)[0]
        assert len(result) == 1
        assert result[0].id == "a2b"
        # SMILEY
        result = gapfilling.SMILEY(m, "b", U)[0]
        assert len(result) == 1
        assert result[0].id == "a2b"

        # 2 rounds of GrowMatch with exchange reactions
        result = gapfilling.growMatch(m, None, ex_rxns=True, iterations=2)
        assert len(result) == 2
        assert len(result[0]) == 1
        assert len(result[1]) == 1
        assert {i[0].id for i in result} == {"SMILEY_EX_b", "SMILEY_EX_c"}

    @pytest.mark.skipif(numpy is None, reason="phase plane require numpy")
    def test_phenotype_phase_plane_benchmark(self, model, benchmark):
        benchmark(calculate_phenotype_phase_plane,
                  model, "EX_glc__D_e", "EX_o2_e",
                  reaction1_npoints=20, reaction2_npoints=20)

    @pytest.mark.skipif(numpy is None, reason="phase plane require numpy")
    def test_phenotype_phase_plane(self, model):
        data = calculate_phenotype_phase_plane(
            model, "EX_glc__D_e", "EX_o2_e",
            reaction1_npoints=20, reaction2_npoints=20)
        assert data.growth_rates.shape == (20, 20)
        assert abs(data.growth_rates.max() - 1.20898) < 0.0001
        assert abs(data.growth_rates[0, :].max()) < 0.0001
        if matplotlib is None:
            pytest.skip("can't test plots without matplotlib")
        data.plot()

    def check_entries(self, out, desired_entries):
        """ensure each entry in desired_entries appears in output"""
        output = out.getvalue().strip()
        output_set = set((re.sub('\s', '', l) for l in output.splitlines()))
        for item in desired_entries:
            assert re.sub('\s', '', item) in output_set

    @pytest.mark.skipif((pandas is None) or (tabulate is None),
                        reason="summary methods require pandas and tabulate")
    def test_summary_methods(self, model, solved_model):
        # Test model summary methods
        with pytest.raises(Exception):
            model.summary()

        desired_entries = [
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
        for solver in solver_dict:
            with captured_output() as (out, err):
                solved_model.summary(fva=0.95, solver=solver)
            self.check_entries(out, desired_entries)

        # test non-fva version (these should be fixed for textbook model
        desired_entries = [
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
        for solver in solver_dict:
            with captured_output() as (out, err):
                solved_model.summary()

            s = out.getvalue()
            for i in desired_entries:
                assert i in s

        # Test metabolite summary methods
        desired_entries = [
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

        for solver in solver_dict:
            with captured_output() as (out, err):
                solved_model.metabolites.q8_c.summary()
            self.check_entries(out, desired_entries)

        desired_entries = [
            'PRODUCING REACTIONS -- D-Fructose 1,6-bisphosphate (fdp_c)',
            '----------------------------------------------------------',
            '%       FLUX  RANGE         RXN ID    REACTION',
            '100%    7.48  [6.17, 9.26]  PFK       '
            'atp_c + f6p_c --> adp_c + fdp_c + h_c',
            'CONSUMING REACTIONS -- D-Fructose 1,6-bisphosphate (fdp_c)',
            '----------------------------------------------------------',
            '%       FLUX  RANGE         RXN ID    REACTION',
            '100%    7.48  [6.17, 8.92]  FBA       fdp_c <=> dhap_c + g3p_c',
            '0%      0     [0, 1.72]     FBP       '
            'fdp_c + h2o_c --> f6p_c + pi_c',
        ]

        for solver in solver_dict:
            with captured_output() as (out, err):
                solved_model.metabolites.fdp_c.summary(fva=0.99, solver=solver)
            self.check_entries(out, desired_entries)


@pytest.mark.skipif(numpy is None, reason="flux sampling requires numpy")
class TestCobraFluxSampling:
    """Test and benchmark flux sampling"""

    def test_single_arch(self, model):
        s = sample(model, 10, method="arch")
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

    def test_fixed_seed(self, model):
        s = sample(model, 1, seed=42)
        assert numpy.allclose(s[0, 94], [8.38570846])

    def setup_class(self):
        from . import create_test_model
        model = create_test_model("textbook")
        arch = ARCHSampler(model, thinning=1)
        assert ((arch.n_warmup > 0) and
                (arch.n_warmup <= 2 * len(model.reactions)))
        assert all(arch.validate(arch.warmup) == "v")
        self.arch = arch

        optgp = OptGPSampler(model, processes=1, thinning=1)
        assert ((optgp.n_warmup > 0) and
                (optgp.n_warmup <= 2 * len(model.reactions)))
        assert all(optgp.validate(optgp.warmup) == "v")
        self.optgp = optgp

    def test_arch_init_benchmark(self, model, benchmark):
        benchmark(lambda: ARCHSampler(model))

    def test_optgp_init_benchmark(self, model, benchmark):
        benchmark(lambda: OptGPSampler(model, processes=2))

    def test_sampling(self):
        s = self.arch.sample(10)
        assert all(self.arch.validate(s) == "v")

        s = self.optgp.sample(10)
        assert all(self.optgp.validate(s) == "v")

    def test_arch_sample_benchmark(self, benchmark):
        benchmark(self.arch.sample, 1)

    def test_optgp_sample_benchmark(self, benchmark):
        benchmark(self.optgp.sample, 1)

    def test_batch_sampling(self):
        for b in self.arch.batch(5, 4):
            assert all(self.arch.validate(b) == "v")

        for b in self.optgp.batch(5, 4):
            assert all(self.optgp.validate(b) == "v")
