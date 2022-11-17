"""Test functions of dictlist.py."""

import re
from copy import copy, deepcopy
from pickle import HIGHEST_PROTOCOL, dumps, loads
from typing import Tuple

import pandas as pd
import pytest

from cobra.core import DictList, Object


@pytest.fixture(scope="function")
def dict_list() -> Tuple[Object, DictList]:
    """Provide function-level fixture for a filled dictlist.

    Returns
    -------
    tuple of Object and DictList
        The tuple with an Object and a filled DictList.

    """
    obj = Object("test1")
    test_list = DictList()
    test_list.append(obj)
    return obj, test_list


def test_contains(dict_list: Tuple[Object, DictList]) -> None:
    """Test containment check for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    assert obj in test_list
    assert obj.id in test_list
    assert Object("not_in") not in test_list
    assert "not_in" not in test_list


def test_index(dict_list: Tuple[Object, DictList]) -> None:
    """Test indexing for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    assert test_list.index("test1") == 0
    assert test_list.index(obj) == 0
    with pytest.raises(ValueError):
        test_list.index("f")
    with pytest.raises(ValueError):
        test_list.index(Object("f"))
    # Ensure indexing with an object that is a different object
    # also raises an error
    with pytest.raises(ValueError):
        test_list.index(Object("test1"))


def test_independent() -> None:
    """Test proper instance creation for dictlist."""
    a = DictList([Object("o1"), Object("o2")])
    b = DictList()
    assert "o1" in a
    assert "o1" not in b
    b.append(Object("o3"))
    assert "o3" not in a
    assert "o3" in b


def test_get_by_any(dict_list: Tuple[Object, DictList]) -> None:
    """Test get_by_any() for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    assert test_list.get_by_any(0) == [obj]
    assert test_list.get_by_any("test1") == [obj]
    with pytest.raises(KeyError):
        test_list.get_by_any("not-in-list")
    with pytest.raises(TypeError):
        test_list.get_by_any(1.1)
    assert test_list.get_by_any(obj) == [obj]


def test_append(dict_list: Tuple[Object, DictList]) -> None:
    """Test append() for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    obj2 = Object("test2")
    test_list.append(obj2)
    with pytest.raises(ValueError):
        test_list.append(Object("test1"))
    assert test_list.index(obj2) == 1
    assert test_list[1] == obj2
    assert test_list.get_by_id("test2") is obj2
    assert len(test_list) == 2


def test_insert(dict_list: Tuple[Object, DictList]) -> None:
    """Test insert() for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    obj2 = Object("a")
    test_list.insert(0, obj2)
    assert test_list.index(obj2) == 0
    assert test_list.index("test1") == 1
    assert test_list.get_by_id("a") is obj2
    assert len(test_list) == 2
    with pytest.raises(ValueError):
        test_list.append(obj2)


def test_extend(dict_list: Tuple[Object, DictList]) -> None:
    """Test extend() for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    obj_list = [Object(f"test{i:d}") for i in range(2, 10)]
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


def test_iadd(dict_list: Tuple[Object, DictList]) -> None:
    """Test in-place addition for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    obj_list = [Object(f"test{i:d}") for i in range(2, 10)]
    test_list += obj_list
    assert test_list[1].id == "test2"
    assert test_list.get_by_id("test2") == obj_list[0]
    assert test_list[8].id == "test9"
    assert len(test_list) == 9


def test_add(dict_list: Tuple[Object, DictList]) -> None:
    """Test addition for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    obj_list = [Object(f"test{i:d}") for i in range(2, 10)]
    sum_ = test_list + obj_list
    assert sum_ is not test_list
    assert sum_ is not obj_list
    assert test_list[0].id == "test1"
    assert sum_[1].id == "test2"
    # noinspection PyUnresolvedReferences
    assert sum_.get_by_id("test2") == obj_list[0]
    assert sum_[8].id == "test9"
    assert len(test_list) == 1
    assert len(sum_) == 9


def test_sub(dict_list: Tuple[Object, DictList]) -> None:
    """Test subtraction for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    obj_list = [Object("test%d" % i) for i in range(2, 10)]
    sum_ = test_list + obj_list
    sub = sum_ - test_list
    assert test_list[0].id == "test1"
    assert sub[0].id == "test2"
    assert len(sub) == 8
    assert sum_ - obj_list == test_list


def test_isub(dict_list: Tuple[Object, DictList]) -> None:
    """Test in-place subtraction for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    obj_list = [Object("test%d" % i) for i in range(2, 10)]
    sum_ = test_list + obj_list
    sum_ -= obj_list[2:4]
    assert len(sum_) == 7
    with pytest.raises(ValueError):
        sum_ -= [Object("bogus")]


def test_init_copy(dict_list: Tuple[Object, DictList]) -> None:
    """Test instance comparison for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
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


def test_slice(dict_list: Tuple[Object, DictList]) -> None:
    """Test slicing for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
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


def test_copy(dict_list: Tuple[Object, DictList]) -> None:
    """Test soft-copy for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
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


def test_deepcopy(dict_list: Tuple[Object, DictList]) -> None:
    """Test deep-copy for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
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


def test_pickle(dict_list: Tuple[Object, DictList]) -> None:
    """Test pickling for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
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


def test_query(dict_list: Tuple[Object, DictList]) -> None:
    """Test query() for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    obj2 = Object("test2")
    obj2.name = "foobar1"
    test_list.append(obj2)
    result = test_list.query("test1")  # matches only test1
    assert len(result) == 1
    result = test_list.query("test1", "id")
    assert len(result) == 1
    assert result[0] == obj
    result = test_list.query("foo", "name")  # matches only test2
    assert len(result) == 1
    assert result[0] == obj2
    result = test_list.query("test", "id")  # matches test1 and test2
    assert len(result) == 2
    # Test with a regular expression
    result = test_list.query(re.compile("test[0-9]"), "id")
    assert len(result) == 2
    result = test_list.query(re.compile("test[29]"), "id")
    assert len(result) == 1
    # Test query of name
    result = test_list.query(re.compile("foobar."), "name")
    assert len(result) == 1
    # Test query with lambda function
    result = test_list.query(lambda x: x.id == "test1")
    assert len(result) == 1


def test_removal() -> None:
    """Test pop() for dictlist."""
    obj_list = DictList(Object(f"test{i:d}") for i in range(2, 10))
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


def test_set() -> None:
    """Test set item for dictlist."""
    obj_list = DictList(Object(f"test{i:d}") for i in range(10))
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
        obj_list.__setitem__(slice(5, 7), [Object("testd"), Object("testd")])


def test_sort_and_reverse() -> None:
    """Test sort() and reverse() for dictlist."""
    dl = DictList(Object(f"test{i:d}") for i in reversed(range(10)))
    assert dl[0].id == "test9"
    dl.sort()
    assert len(dl) == 10
    assert dl[0].id == "test0"
    assert dl.index("test0") == 0
    dl.reverse()
    assert dl[0].id == "test9"
    assert dl.index("test0") == 9


def test_dir(dict_list: Tuple[Object, DictList]) -> None:
    """Test local scope item listings for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    # Make sure tab completion works
    attrs = dir(test_list)
    assert "test1" in attrs
    assert "_dict" in attrs  # attribute of DictList


def test_union(dict_list: Tuple[Object, DictList]) -> None:
    """Test union() for dictlist.

    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    obj, test_list = dict_list
    test_list.union([Object("test1"), Object("test2")])
    # Add only 1 element
    assert len(test_list) == 2
    assert test_list.index("test2") == 1

def test_to_df(dict_list: Tuple[Object, DictList]) -> None:
    """Test to_df for dictlist.
    
    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    _, test_list = dict_list
    test_list.name = "foo"
    df = test_list.to_df()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert "id" in df.columns
    assert "test1" in df["id"].values
    assert "name" in df.columns
    assert "foo" in df["name"].values

def test__repr_html(dict_list: Tuple[Object, DictList]) -> None:
    """Test _repr_html_ for dictlist.
    
    Parameters
    ----------
    dict_list : tuple
        The fixture for filled dictlist.

    """
    _, test_list = dict_list
    test_list.name = "foo"
    html = test_list._repr_html_()
    assert isinstance(html, str)
    assert "test1" in html
    assert "foo" in html