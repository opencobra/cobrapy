# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from operator import attrgetter
from os.path import join
from random import sample

import pytest

import cobra.io as cio
from cobra import DictList, Model

LOGGER = logging.getLogger(__name__)


@pytest.fixture(scope="module")
def tmp_path(tmpdir_factory):
    return str(tmpdir_factory.mktemp("model_order"))


@pytest.fixture(scope="module")
def minimized_shuffle(small_model):
    model = small_model.copy()
    chosen = sample(list(set(model.reactions) - set(model.exchanges)), 10)
    new = Model("minimized_shuffle")
    new.add_reactions(chosen)
    LOGGER.debug("'%s' has %d metabolites, %d reactions, and %d genes.",
                 new.id, new.metabolites, new.reactions, new.genes)
    return new


@pytest.fixture(scope="module")
def minimized_sorted(minimized_shuffle):
    model = minimized_shuffle.copy()
    model.id = "minimized_sorted"
    model.metabolites = DictList(
        sorted(model.metabolites, key=attrgetter("id")))
    model.genes = DictList(sorted(model.genes, key=attrgetter("id")))
    model.reactions = DictList(sorted(model.reactions, key=attrgetter("id")))
    return model


@pytest.fixture(scope="module")
def minimized_reverse(minimized_shuffle):
    model = minimized_shuffle.copy()
    model.id = "minimized_reverse"
    model.metabolites = DictList(
        sorted(model.metabolites, key=attrgetter("id"), reverse=True))
    model.genes = DictList(
        sorted(model.genes, key=attrgetter("id"), reverse=True))
    model.reactions = DictList(
        sorted(model.reactions, key=attrgetter("id"), reverse=True))
    return model


@pytest.fixture(scope="module", params=[
    "minimized_shuffle", "minimized_reverse", "minimized_sorted"])
def template(request, minimized_shuffle, minimized_reverse, minimized_sorted):
    return locals()[request.param]


@pytest.fixture(scope="module", params=["metabolites", "reactions", "genes"])
def attribute(request):
    return request.param


def get_ids(iterable):
    return [x.id for x in iterable]


@pytest.mark.parametrize("read, write, ext", [
    ("read_sbml_model", "write_sbml_model", ".xml"),
    pytest.mark.skip(("read_legacy_sbml", "write_legacy_sbml", ".xml"),
                     reason="Order for legacy SBML I/O is uninteresting."),
    pytest.mark.skip(("load_matlab_model", "save_matlab_model", ".mat"),
                     reason="Order for Matlab model I/O is uninteresting."),
    ("load_json_model", "save_json_model", ".json"),
    ("load_yaml_model", "save_yaml_model", ".yml"),
])
def test_io_order(attribute, read, write, ext, template, tmp_path):
    read = getattr(cio, read)
    write = getattr(cio, write)
    filename = join(tmp_path, "template" + ext)
    write(template, filename)
    model = read(filename)
    model_elements = get_ids(getattr(model, attribute))
    template_elements = get_ids(getattr(template, attribute))
    assert len(model_elements) == len(template_elements)
    assert set(model_elements) == set(template_elements)
    assert model_elements == template_elements
