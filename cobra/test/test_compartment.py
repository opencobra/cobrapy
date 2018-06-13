# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pytest
from cobra.core import Compartment


def test__eq__():
    a = b = Compartment("x")
    c = "x"
    d = Compartment("x")
    e = Compartment("z")
    assert a.__eq__(b) is True
    assert a.__eq__(c) is True
    assert a.__eq__(d) is True
    assert a.__eq__(e) is False


def test__ne__():
    a = b = Compartment("x")
    c = "x"
    d = Compartment("x")
    e = Compartment("z")
    assert a.__ne__(b) is False
    assert a.__ne__(c) is False
    assert a.__ne__(d) is False
    assert a.__ne__(e) is True


def test_error_with_non_sane_id():
    with pytest.raises(TypeError):
        Compartment("Trust Me This ID Is Sane")
