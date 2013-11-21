from unittest import TestCase, TestLoader, TextTestRunner
import sys
# deal with absolute imports by adding the appropriate directory to the path
if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.manipulation import initialize_growth_medium
    from cobra.test import create_test_model
    from cobra import Model, Reaction, Metabolite
    from cobra import solvers
    from cobra.solvers import __legacy_solver
    sys.path.pop(0)  # remove the added directory to the path
else:
    from ..manipulation import initialize_growth_medium
    from . import create_test_model
    from .. import Model, Reaction, Metabolite
    from .. import solvers
    __legacy_solver = solvers.__legacy_solver

class TestCobraSolver(TestCase):
    def setUp(self):
        self.model = create_test_model()
        initialize_growth_medium(self.model, 'MgM')
        self.old_solution = 0.320064
        self.infeasible_model = Model()
        metabolite_1 = Metabolite("met1")
        metabolite_2 = Metabolite("met2")
        reaction_1 = Reaction("rxn1")
        reaction_2 = Reaction("rxn2")
        reaction_1.add_metabolites({metabolite_1: 1})
        reaction_2.add_metabolites({metabolite_2: 1})
        reaction_1.lower_bound = 1
        reaction_2.upper_bound = 2
        self.infeasible_model.add_reactions([reaction_1, reaction_2])
        #self.infeasible_model.update()


def add_new_test(TestCobraSolver, solver_name, solver):
    """Creates a test set for each of the solvers that are installed
    using the modular interface.

    """
    def attributes(self):
        self.assertTrue(hasattr(solver, "create_problem"))
        self.assertTrue(hasattr(solver, "solve_problem"))
        # self.assertTrue(hasattr(solver, "update_problem"))

    def creation(self):
        solver.create_problem(self.model)

    def solve_feasible(self):
        solver.solve(self.model)
        solution = self.model.solution        
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(self.old_solution, \
            solution.f, places=4)

    def solve_minimize(self):
        solver.solve(self.model, objective_sense='minimize')
        solution = self.model.solution        
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(0, \
            solution.f, places=4)

    def set_objective_sense(self):
        maximize = solver.create_problem(self.model, objective_sense="maximize")
        minimize = solver.create_problem(self.model, objective_sense="minimize")
        solver.solve_problem(maximize)
        solver.solve_problem(minimize)
        max_solution = solver.format_solution(maximize, self.model)
        min_solution = solver.format_solution(minimize, self.model)
        self.assertAlmostEqual(0, min_solution.f, places=4)
        self.assertAlmostEqual(self.old_solution, max_solution.f, places=4)
        # if we set minimize at creation, can we override it at solve
        solver.solve_problem(minimize, objective_sense="maximize")
        override_minimize = solver.format_solution(minimize, self.model)
        self.assertAlmostEqual(max_solution.f, override_minimize.f, places=4)

    def solve_mip(self):
        cone_selling_price = 7.
        cone_production_cost = 3.
        popsicle_selling_price = 2.
        popsicle_production_cost = 1.
        starting_budget = 100.
        cobra_model = Model('MILP_implementation_test')
        cone_out = Metabolite(id='cone_out', compartment='c')
        cone_in = Metabolite(id='cone_in', compartment='c')
        cone_consumed = Metabolite(id='cone_consumed', compartment='c')

        popsicle_out = Metabolite(id='popsicle_out', compartment='c')
        popsicle_in = Metabolite(id='popsicle_in', compartment='c')
        popsicle_consumed = Metabolite(id='popsicle_consumed', compartment='c')

        the_reactions = []

        # SOURCE
        Cone_source = Reaction(name='Cone_source')
        temp_metabolite_dict = {cone_out: 1}
        Cone_source.add_metabolites(temp_metabolite_dict)
        the_reactions.append(Cone_source)

        Popsicle_source = Reaction(name='Popsicle_source')
        temp_metabolite_dict = {popsicle_out: 1}
        Popsicle_source.add_metabolites(temp_metabolite_dict)
        the_reactions.append(Popsicle_source)


        ## PRODUCTION
        Cone_production = Reaction(name='Cone_production')
        temp_metabolite_dict = {cone_out: -1,
                                cone_in: 1}
        Cone_production.add_metabolites(temp_metabolite_dict)
        the_reactions.append(Cone_production)


        Popsicle_production = Reaction(name='Popsicle_production')
        temp_metabolite_dict = {popsicle_out: -1,
                                popsicle_in: 1}
        Popsicle_production.add_metabolites(temp_metabolite_dict)
        the_reactions.append(Popsicle_production)

        ## CONSUMPTION
        Cone_consumption = Reaction(name='Cone_consumption')
        temp_metabolite_dict = {cone_in: -1,
                                cone_consumed: 1}
        Cone_consumption.add_metabolites(temp_metabolite_dict)
        the_reactions.append(Cone_consumption)

        Popsicle_consumption = Reaction(name='Popsicle_consumption')
        temp_metabolite_dict = {popsicle_in: -1,
                                popsicle_consumed: 1}
        Popsicle_consumption.add_metabolites(temp_metabolite_dict)
        the_reactions.append(Popsicle_consumption)

        # SINK
        Cone_consumed_sink = Reaction(name='Cone_consumed_sink')
        temp_metabolite_dict = {cone_consumed: -1}
        Cone_consumed_sink.add_metabolites(temp_metabolite_dict)
        the_reactions.append(Cone_consumed_sink)

        Popsicle_consumed_sink = Reaction(name='Popsicle_consumed_sink')
        temp_metabolite_dict = {popsicle_consumed: -1}
        Popsicle_consumed_sink.add_metabolites(temp_metabolite_dict)
        the_reactions.append(Popsicle_consumed_sink)

        ## add all reactions
        cobra_model.add_reactions(the_reactions)

        # set objective coefficients
        Cone_consumption.objective_coefficient = cone_selling_price
        Popsicle_consumption.objective_coefficient = popsicle_selling_price

        Cone_production.objective_coefficient = -1*cone_production_cost
        Popsicle_production.objective_coefficient = -1*popsicle_production_cost
        

        #Make sure we produce whole cones
        Cone_production.variable_kind = 'integer'
        Popsicle_production.variable_kind = 'integer'


        production_capacity_constraint = Metabolite(id='production_capacity_constraint')
        production_capacity_constraint._constraint_sense = 'L'
        production_capacity_constraint._bound = starting_budget;

        Cone_production.add_metabolites({production_capacity_constraint: cone_production_cost })

        Popsicle_production.add_metabolites({production_capacity_constraint: popsicle_production_cost })
        cobra_model.optimize()
        self.assertEqual(133, cobra_model.solution.f)

    def solve_infeasible(self):
        solver.solve(self.infeasible_model)
        solution = self.infeasible_model.solution
        self.assertEqual(solution.status, "infeasible")

    def independent_creation(self):
        feasible_lp = solver.create_problem(self.model)
        infeasible_lp = solver.create_problem(self.infeasible_model)
        solver.solve_problem(feasible_lp)
        solver.solve_problem(infeasible_lp)
        feasible_solution = solver.format_solution(feasible_lp, self.model)
        infeasible_solution = solver.format_solution(infeasible_lp, self.infeasible_model)
        self.assertEqual(feasible_solution.status, "optimal")
        self.assertAlmostEqual(self.old_solution, \
            feasible_solution.f, places=4)
        self.assertEqual(infeasible_solution.status, "infeasible")

    for tester in [attributes, creation, solve_feasible, solve_minimize,
                    set_objective_sense, solve_mip, solve_infeasible,
                    independent_creation]:
        setattr(TestCobraSolver, "test_%s_%s" %
                (solver_name, tester.func_name), tester)


def add_legacy_test(TestCobraSolver, solver_name, solver_function):
    """Creates a test set for each of the installed solvers using the
    legacy interface.

    """
    def test_solve_feasible(self):

        the_solution = solver_function(self.model)['the_solution']
        self.assertEqual(the_solution.status, 'optimal')
        self.assertAlmostEqual(self.old_solution, the_solution.f, places=4)
    setattr(TestCobraSolver, "test_%s_feasible_solve" % solver_name, \
        test_solve_feasible)
        


if not __legacy_solver:
    add_test = add_new_test
else:
    add_test = add_legacy_test

    
for solver_name, solver in solvers.solver_dict.iteritems():
    add_test(TestCobraSolver, solver_name, solver)
# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
