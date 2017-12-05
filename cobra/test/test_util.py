# -*- coding: utf-8 -*-
from __future__ import absolute_import

import re
from copy import copy, deepcopy
from pickle import HIGHEST_PROTOCOL, dumps, loads

import pytest
from six.moves import range

from cobra import DictList, Object
from cobra.util.version_info import (
    get_sys_info, get_pkg_info, show_versions, SYS_ORDER, PKG_ORDER)


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

    def test_get_by_any(self, dict_list):
        obj, test_list = dict_list
        assert test_list.get_by_any(0) == [obj]
        assert test_list.get_by_any('test1') == [obj]
        with pytest.raises(KeyError):
            test_list.get_by_any('not-in-list')
        with pytest.raises(TypeError):
            test_list.get_by_any(1.1)
        assert test_list.get_by_any(obj) == [obj]

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

    def test_sub(self, dict_list):
        obj, test_list = dict_list
        obj_list = [Object("test%d" % i) for i in range(2, 10)]
        sum = test_list + obj_list
        sub = sum - test_list
        assert test_list[0].id == "test1"
        assert sub[0].id == "test2"
        assert len(sub) == 8
        assert sum - obj_list == test_list

    def test_isub(self, dict_list):
        obj, test_list = dict_list
        obj_list = [Object("test%d" % i) for i in range(2, 10)]
        sum = test_list + obj_list
        sum -= obj_list[2:4]
        assert len(sum) == 7
        with pytest.raises(ValueError):
            sum -= [Object('bogus')]

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
        result = test_list.query(u"test1", "id")  # matches with unicode
        assert len(result) == 1
        assert result[0] == obj
        result = test_list.query("foo", "name")  # matches only test2
        assert len(result) == 1
        assert result[0] == obj2
        result = test_list.query("test", "id")  # matches test1 and test2
        assert len(result) == 2
        # test with a regular expression
        result = test_list.query(re.compile("test[0-9]"), "id")
        assert len(result) == 2
        result = test_list.query(re.compile("test[29]"), "id")
        assert len(result) == 1
        # test query of name
        result = test_list.query(re.compile("foobar."), "name")
        assert len(result) == 1
        # test query with lambda function
        result = test_list.query(lambda x: x.id == 'test1')
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


class TestVersionInfo:
    SKIP_OPTIONAL = frozenset([
        "cobra", "python-libsbml", "lxml", "matplotlib", "palettable", "scipy",
        "pymatbridge"])

    @pytest.fixture(scope="module")
    def sys_info(self):
        return get_sys_info()

    @pytest.fixture(scope="module")
    def pkg_info(self):
        return get_pkg_info()

    @pytest.mark.parametrize("key", SYS_ORDER)
    def test_sys_info_key(self, key, sys_info):
        assert key in sys_info

    @pytest.mark.parametrize("key", SYS_ORDER)
    def test_sys_info_value(self, key, sys_info):
        assert len(sys_info[key]) > 0

    @pytest.mark.parametrize("key", PKG_ORDER)
    def test_pkg_info_key(self, key, pkg_info):
        if key in self.SKIP_OPTIONAL:
            pytest.skip()
        assert key in pkg_info

    @pytest.mark.parametrize("key", PKG_ORDER)
    def test_pkg_info_value(self, key, pkg_info):
        if key in self.SKIP_OPTIONAL:
            pytest.skip()
        assert len(pkg_info[key]) > 0

    def test_show_versions(self, sys_info, pkg_info, capsys):
        show_versions()
        out, err = capsys.readouterr()
        lines = out.split("\n")
        i = 3
        for key in SYS_ORDER:
            line = lines[i]
            assert line.startswith(key)
            assert line.endswith(sys_info[key])
            i += 1
        i += 3
        for key in PKG_ORDER:
            line = lines[i]
            if key in self.SKIP_OPTIONAL:
                if line.startswith(key):
                    i += 1
                continue
            assert line.startswith(key)
            assert line.endswith(pkg_info[key])
            i += 1
