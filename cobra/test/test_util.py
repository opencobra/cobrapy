from __future__ import absolute_import, print_function

from copy import deepcopy, copy
import re
from cobra import DictList, Object
from cobra.core import LazySolution
from pickle import loads, dumps, HIGHEST_PROTOCOL
from functools import partial
from itertools import chain
from six.moves import range
from cobra.util import TimeMachine, generate_colors, Singleton, partition, \
    frozendict, ProblemCache
from . import create_test_model
import pytest
from .conftest import model


@pytest.fixture(scope="session")
def tm():
    return TimeMachine()


@pytest.fixture(scope="session")
def seed():
    return 1234


@pytest.fixture(scope="session")
def frozen_dict():
    return frozendict({"A": 1, "B": 2, "C": 3, "D": 4, "E": [2, 3, 4, 5]})


@pytest.fixture(scope="function")
def problem_cache_trial():
    model = create_test_model("textbook")
    reference = model.optimize(solution_type=LazySolution).fluxes
    n_constraints = len(model.solver.constraints)
    n_variables = len(model.solver.variables)
    return model, reference, n_constraints, n_variables


@pytest.fixture(scope="function")
def dict_list():
    obj = Object("test1")
    test_list = DictList()
    test_list.append(obj)
    return obj, test_list


class TestDictList:
    def test_contains(self, dict_list):
        obj, test_list = dict_list
        assert obj in test_list
        assert obj.id in test_list
        assert Object("not_in") not in test_list
        assert 'not_in' not in test_list

    def test_index(self, dict_list):
        obj, test_list = dict_list
        assert test_list.index("test1") == 0
        assert test_list.index(obj) == 0
        with pytest.raises(ValueError):
            test_list.index("f")
        with pytest.raises(ValueError):
            test_list.index(Object("f"))
        # ensure trying to index with an object that is a different object
        # also raises an error
        with pytest.raises(ValueError):
            test_list.index(Object("test1"))

    def test_independent(self):
        a = DictList([Object("o1"), Object("o2")])
        b = DictList()
        assert "o1" in a
        assert "o1" not in b
        b.append(Object("o3"))
        assert "o3" not in a
        assert "o3" in b

    def test_append(self, dict_list):
        obj, test_list = dict_list
        obj2 = Object("test2")
        test_list.append(obj2)
        with pytest.raises(ValueError):
            test_list.append(Object("test1"))
        assert test_list.index(obj2) == 1
        assert test_list[1] == obj2
        assert test_list.get_by_id("test2") is obj2
        assert len(test_list) == 2

    def test_insert(self, dict_list):
        obj, test_list = dict_list
        obj2 = Object("a")
        test_list.insert(0, obj2)
        assert test_list.index(obj2) == 0
        assert test_list.index("test1") == 1
        assert test_list.get_by_id("a") is obj2
        assert len(test_list) == 2
        with pytest.raises(ValueError):
            test_list.append(obj2)

    def test_extend(self, dict_list):
        obj, test_list = dict_list
        obj_list = [Object("test%d" % (i)) for i in range(2, 10)]
        test_list.extend(obj_list)
        assert test_list[1].id == "test2"
        assert test_list.get_by_id("test2") == obj_list[0]
        assert test_list[8].id == "test9"
        assert len(test_list) == 9
        with pytest.raises(ValueError):
            test_list.extend([Object("test1")])
        # Even if the object is unique, if it is present twice in the new
        # list, it should still raise an exception
        with pytest.raises(ValueError):
            test_list.extend([Object("testd"), Object("testd")])

    def test_iadd(self, dict_list):
        obj, test_list = dict_list
        obj_list = [Object("test%d" % (i)) for i in range(2, 10)]
        test_list += obj_list
        assert test_list[1].id == "test2"
        assert test_list.get_by_id("test2") == obj_list[0]
        assert test_list[8].id == "test9"
        assert len(test_list) == 9

    def test_add(self, dict_list):
        obj, test_list = dict_list
        obj_list = [Object("test%d" % (i)) for i in range(2, 10)]
        sum = test_list + obj_list
        assert sum is not test_list
        assert sum is not obj_list
        assert test_list[0].id == "test1"
        assert sum[1].id == "test2"
        assert sum.get_by_id("test2") == obj_list[0]
        assert sum[8].id == "test9"
        assert len(test_list) == 1
        assert len(sum) == 9

    def test_init_copy(self, dict_list):
        obj, test_list = dict_list
        test_list.append(Object("test2"))
        copied = DictList(test_list)
        assert test_list is not copied
        assert isinstance(copied, test_list.__class__)
        assert len(test_list) == len(copied)
        for i, v in enumerate(test_list):
            assert test_list[i].id == copied[i].id
            assert i == copied.index(v.id)
            assert test_list[i] is copied[i]
            assert v is copied.get_by_id(v.id)

    def test_slice(self, dict_list):
        obj, test_list = dict_list
        test_list.append(Object("test2"))
        test_list.append(Object("test3"))
        sliced = test_list[:-1]
        assert test_list is not sliced
        assert isinstance(sliced, test_list.__class__)
        assert len(test_list) == len(sliced) + 1
        for i, v in enumerate(sliced):
            assert test_list[i].id == sliced[i].id
            assert i == sliced.index(v.id)
            assert test_list[i] is sliced[i]
            assert test_list[i] is sliced.get_by_id(v.id)

    def test_copy(self, dict_list):
        obj, test_list = dict_list
        test_list.append(Object("test2"))
        copied = copy(test_list)
        assert test_list is not copied
        assert isinstance(copied, test_list.__class__)
        assert len(test_list) == len(copied)
        for i, v in enumerate(test_list):
            assert test_list[i].id == copied[i].id
            assert i == copied.index(v.id)
            assert test_list[i] is copied[i]
            assert v is copied.get_by_id(v.id)

    def test_deepcopy(self, dict_list):
        obj, test_list = dict_list
        test_list.append(Object("test2"))
        copied = deepcopy(test_list)
        assert test_list is not copied
        assert isinstance(copied, test_list.__class__)
        assert len(test_list) == len(copied)
        for i, v in enumerate(test_list):
            assert test_list[i].id == copied[i].id
            assert i == copied.index(v.id)
            assert test_list[i] is not copied[i]
            assert v is not copied.get_by_id(v.id)

    def test_pickle(self, dict_list):
        obj, test_list = dict_list
        test_list.append(Object("test2"))
        for protocol in range(HIGHEST_PROTOCOL):
            pickle_str = dumps(test_list, protocol=protocol)
            copied = loads(pickle_str)
            assert test_list is not copied
            assert isinstance(copied, test_list.__class__)
            assert len(test_list) == len(copied)
            for i, v in enumerate(test_list):
                assert test_list[i].id == copied[i].id
                assert i == copied.index(v.id)
                assert test_list[i] is not copied[i]
                assert v is not copied.get_by_id(v.id)

    def test_query(self, dict_list):
        obj, test_list = dict_list
        obj2 = Object("test2")
        obj2.name = "foobar1"
        test_list.append(obj2)
        result = test_list.query("test1")  # matches only test1
        assert len(result) == 1
        result = test_list.query(u"test1")  # matches with unicode
        assert len(result) == 1
        assert result[0] == obj
        result = test_list.query("foo", "name")  # matches only test2
        assert len(result) == 1
        assert result[0] == obj2
        result = test_list.query("test")  # matches test1 and test2
        assert len(result) == 2
        # test with a regular expression
        result = test_list.query(re.compile("test[0-9]"))
        assert len(result) == 2
        result = test_list.query(re.compile("test[29]"))
        assert len(result) == 1
        # test query of name
        result = test_list.query(re.compile("foobar."), "name")
        assert len(result) == 1

    def test_removal(self):
        obj_list = DictList(Object("test%d" % (i)) for i in range(2, 10))
        del obj_list[3]
        assert "test5" not in obj_list
        assert obj_list.index(obj_list[-1]) == len(obj_list) - 1
        assert len(obj_list) == 7
        del obj_list[3:5]
        assert "test6" not in obj_list
        assert "test7" not in obj_list
        assert obj_list.index(obj_list[-1]) == len(obj_list) - 1
        assert len(obj_list) == 5
        removed = obj_list.pop(1)
        assert obj_list.index(obj_list[-1]) == len(obj_list) - 1
        assert removed.id == "test3"
        assert "test3" not in obj_list
        assert len(obj_list) == 4
        removed = obj_list.pop()
        assert removed.id == "test9"
        assert removed.id not in obj_list
        assert len(obj_list) == 3

    def test_set(self):
        obj_list = DictList(Object("test%d" % (i)) for i in range(10))
        obj_list[4] = Object("testa")
        assert obj_list.index("testa") == 4
        assert obj_list[4].id == "testa"
        obj_list[5:7] = [Object("testb"), Object("testc")]
        assert obj_list.index("testb") == 5
        assert obj_list[5].id == "testb"
        assert obj_list.index("testc") == 6
        assert obj_list[6].id == "testc"
        # Even if the object is unique, if it is present twice in the new
        # list, it should still raise an exception
        with pytest.raises(ValueError):
            obj_list.__setitem__(slice(5, 7),
                                 [Object("testd"), Object("testd")])

    def test_sort_and_reverse(self):
        dl = DictList(Object("test%d" % (i)) for i in reversed(range(10)))
        assert dl[0].id == "test9"
        dl.sort()
        assert len(dl) == 10
        assert dl[0].id == "test0"
        assert dl.index("test0") == 0
        dl.reverse()
        assert dl[0].id == "test9"
        assert dl.index("test0") == 9

    def test_dir(self, dict_list):
        obj, test_list = dict_list
        """makes sure tab complete will work"""
        attrs = dir(test_list)
        assert "test1" in attrs
        assert "_dict" in attrs  # attribute of DictList

    def test_union(self, dict_list):
        obj, test_list = dict_list
        test_list.union([Object("test1"), Object("test2")])
        # should only add 1 element
        assert len(test_list) == 2
        assert test_list.index("test2") == 1


class TimeMachineTestCase:
    def test_one_change_list(self, tm):
        l = [1, 2, 3, 4]
        tm(do=partial(l.append, 5), undo=l.pop)
        assert l == [1, 2, 3, 4, 5]
        tm.reset()
        assert l == [1, 2, 3, 4]

    def test_str_handles_different_types_of_stored_operations(self, tm):
        def normal_function():
            pass

        partial_function = partial(str, 1)
        tm(do=normal_function, undo=partial_function)
        assert tm.__str__().split('\n')[2:-1] == [
            "undo: " + str(str) + " (1,) {}", 'redo: normal_function']

    def test_with_statement(self):
        l = [1, 2, 3, 4]
        with TimeMachine() as tm:
            tm(do=partial(l.append, 33), undo=partial(l.pop))
            tm(do=partial(l.append, 66), undo=partial(l.pop))
            tm(do=partial(l.append, 99), undo=partial(l.pop))
        assert l == [1, 2, 3, 4]


def some_method_that_adds_stuff(model, cache):
    assert isinstance(cache, ProblemCache)

    def create_variable(model, var_id, lb, ub):
        return model.solver.interface.Variable(var_id, ub=ub, lb=lb)

    def update_variable(model, var, lb, ub):
        var.lb = lb
        var.ub = ub

    def create_constraint(model, cid, vars, lb, ub):
        return model.solver.interface.Constraint(vars[0] + vars[1], ub=ub,
                                                 lb=lb)

    def update_constraint(model, constraint, vars, lb, ub):
        constraint.lb = lb
        constraint.ub = ub

    for i in range(10):
        cache.add_variable("var_%i" % (i + 1), create_variable,
                           update_variable, 10, 15)

    for i in range(9):
        v1 = cache.variables["var_%i" % (i + 1)]
        v2 = cache.variables["var_%i" % (i + 2)]
        cache.add_constraint("c_%i" % (i + 1), create_constraint,
                             update_constraint, [v1, v2], -20, 100)


class TestProblemCache:
    def test_add_variable(self, model):
        cache = ProblemCache(model)

        def add_var(model, var_id):
            return model.solver.interface.Variable(var_id, ub=0)

        def update_var(model, var):
            return setattr(var, "ub", 1000)

        for i in range(10):
            cache.add_variable("%i" % i, add_var, update_var)

        for i in range(10):
            assert cache.variables["%i" % i] in model.solver.variables
            assert cache.variables["%i" % i].ub == 0
            assert model.solver.variables["%i" % i].ub == 0

        for i in range(10):
            cache.add_variable("%i" % i, add_var, update_var)
            assert cache.variables["%i" % i].ub == 1000
            assert model.solver.variables["%i" % i].ub == 1000

        cache.reset()

        for i in range(10):
            with pytest.raises(KeyError):
                model.solver.variables.__getitem__("%i" % i)

    def test_add_constraint(self, model):
        cache = ProblemCache(model)

        def add_var(model, var_id):
            return model.solver.interface.Variable(var_id, ub=0)

        def add_constraint(m, const_id, var):
            return m.solver.interface.Constraint(var, lb=-10, ub=10,
                                                 name=const_id)

        def update_constraint(model, const, var):
            return setattr(const, "ub", 1000)

        for i in range(10):
            cache.add_variable("%i" % i, add_var, None)
            cache.add_constraint("c%i" % i, add_constraint,
                                 update_constraint,
                                 cache.variables["%i" % i])

        for i in range(10):
            assert cache.constraints[
                       "c%i" % i] in model.solver.constraints
            assert cache.constraints["c%i" % i].ub == 10
            assert cache.constraints["c%i" % i].lb == -10
            assert model.solver.constraints["c%i" % i].ub == 10
            assert model.solver.constraints["c%i" % i].lb == -10

        for i in range(10):
            cache.add_constraint("c%i" % i, add_constraint,
                                 update_constraint,
                                 cache.variables["%i" % i])
            assert model.solver.constraints["c%i" % i].ub == 1000

        cache.reset()

        for i in range(10):
            with pytest.raises(KeyError):
                model.solver.variables.__getitem__("%i" % i)
            with pytest.raises(KeyError):
                model.solver.constraints.__getitem__("c%i" % i)

    def test_cache_problem(self, problem_cache_trial):
        model, reference, n_constraints, n_variables = problem_cache_trial
        # After the number of variables and constraints remains the same if
        # nothing happens
        assert n_constraints == len(model.solver.constraints)
        assert n_variables == len(model.solver.variables)

        cache = ProblemCache(model)
        some_method_that_adds_stuff(model, cache)
        # After running some_method_that_adds_stuff with cache, problem has
        # 10 more variables
        assert n_variables + 10 == len(model.solver.variables)
        # And has 9 more more constraints
        assert n_constraints + 9 == len(model.solver.constraints)

        cache.reset()
        # After reset cache, the problem should return to its original size
        assert n_constraints == len(model.solver.constraints)
        assert n_variables == len(model.solver.variables)

    def test_with(self, problem_cache_trial):
        model, reference, n_constraints, n_variables = problem_cache_trial
        with ProblemCache(model) as cache:
            some_method_that_adds_stuff(model, cache)
            # After running some_method_that_adds_stuff with cache, problem
            # has 10 more variables
            assert n_variables + 10 == len(model.solver.variables)
            # And has 9 more more constraints
            assert n_constraints + 9 == len(model.solver.constraints)

            # If the method runs again, it does not add repeated variables
            some_method_that_adds_stuff(model, cache)
            # After running some_method_that_adds_stuff with cache, problem
            # has 10 more variables
            assert n_variables + 10 == len(model.solver.variables)
            # And has 9 more more constraints
            assert n_constraints + 9 == len(model.solver.constraints)

        # After reset cache, the problem should return to its original size
        assert n_constraints == len(model.solver.constraints)
        assert n_variables == len(model.solver.variables)


class TestUtils:
    def test_color_generation(self):
        for i in range(1, 100):
            color_map = generate_colors(i)
            assert len(color_map) == i
            assert len(color_map) == len(set(color_map.values()))

    def test_partition(self):
        chunks = 3
        iterables = [
            [1, 2, 3, 4, 5, 6, 7, 8, 9],
            {5, 3, 8, 3, 8, 5, 8, 0, 10, 11, 15},
            range(29)
        ]
        for fixture in iterables:
            test_output = partition(fixture, chunks)
            assert len(fixture) == sum(map(len, test_output))
            assert len(test_output) == chunks
            assert list(fixture) == list(chain(*test_output))
            for out_chunk in test_output:
                assert set(out_chunk).issubset(set(fixture))

        bad_input = 5
        with pytest.raises(TypeError):
            partition(bad_input, chunks)


class FrozendictTestCase:
    def test_frozen_attributes(self, frozen_dict):
        with pytest.raises(AttributeError):
            frozen_dict.popitem()
        with pytest.raises(AttributeError):
            frozen_dict.pop("A")
        with pytest.raises(AttributeError):
            frozen_dict.__setitem__("C", 1)
        with pytest.raises(AttributeError):
            frozen_dict.setdefault("K")
        with pytest.raises(AttributeError):
            frozen_dict.__delitem__("A")
        with pytest.raises(AttributeError):
            frozen_dict.update()

        assert hasattr(frozen_dict, "__hash__")


class TestSingleton:
    def test_singleton(self):
        s1 = Singleton()
        s2 = Singleton()
        assert s1 is s2
