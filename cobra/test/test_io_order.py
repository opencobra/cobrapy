# -*- coding: utf-8 -*-

from __future__ import absolute_import

from operator import attrgetter
from os.path import join, dirname
from random import sample

import pytest

import cobra.io as cio
from cobra import DictList


@pytest.fixture(scope="module")
def tmp_path(tmpdir_factory):
    return str(tmpdir_factory.mktemp("model_order"))


@pytest.fixture(scope="module")
def minimized_core():
    model = cio.read_sbml_model(join(dirname(__file__), "data",
                                "textbook.xml.gz"))
    model.id = "minimized_core"
    chosen = sample(model.reactions, 10)
    model.remove_reactions(set(model.reactions).difference(chosen),
                           remove_orphans=True)
    return model


@pytest.fixture(scope="module")
def minimized_sorted(minimized_core):
    model = minimized_core.copy()
    model.id = "minimized_sorted"
    model.metabolites = DictList(
        sorted(model.metabolites, key=attrgetter("id")))
    model.genes = DictList(sorted(model.genes, key=attrgetter("id")))
    model.reactions = DictList(sorted(model.reactions, key=attrgetter("id")))
    return model


@pytest.fixture(scope="module")
def minimized_reverse(minimized_core):
    model = minimized_core.copy()
    model.id = "minimized_reverse"
    model.metabolites = DictList(
        sorted(model.metabolites, key=attrgetter("id"), reverse=True))
    model.genes = DictList(
        sorted(model.genes, key=attrgetter("id"), reverse=True))
    model.reactions = DictList(
        sorted(model.reactions, key=attrgetter("id"), reverse=True))
    return model


@pytest.fixture(scope="module", params=[
    "minimized_core", "minimized_reverse", "minimized_sorted"])
def template(request, minimized_core, minimized_reverse, minimized_sorted):
    return locals()[request.param]


@pytest.fixture(scope="module", params=["metabolites", "reactions", "genes"])
def attribute(request):
    return request.param


def get_ids(iterable):
    return [x.id for x in iterable]


@pytest.mark.parametrize("read, write, ext", [
    ("read_sbml_model", "write_sbml_model", ".xml"),
    ("read_legacy_sbml", "write_legacy_sbml", ".xml"),
    ("load_matlab_model", "save_matlab_model", ".mat"),
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

