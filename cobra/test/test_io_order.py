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


def get_attr(iterable, attr="id"):
    get = attrgetter(attr)
    return [get(x) for x in iterable]


@pytest.mark.parametrize("read, write, ext", [
    (cio.read_sbml_model, cio.write_sbml_model, ".xml"),
    (cio.read_legacy_sbml, cio.write_legacy_sbml, ".xml"),
    (cio.load_matlab_model, cio.save_matlab_model, ".mat"),
    (cio.load_json_model, cio.save_json_model, ".json"),
    (cio.load_yaml_model, cio.save_yaml_model, ".yml"),
])
def test_io_order(read, write, ext, template, tmp_path):
    filename = join(tmp_path, "template" + ext)
    write(template, filename)
    model = read(filename)
    # test metabolite order
    assert len(model.metabolites) == len(template.metabolites)
    assert set(get_attr(model.metabolites)) == set(
        get_attr(template.metabolites))
    assert get_attr(model.metabolites) == get_attr(template.metabolites)
    # test reaction order
    assert len(model.reactions) == len(template.reactions)
    assert set(get_attr(model.reactions)) == set(get_attr(template.reactions))
    assert get_attr(model.reactions) == get_attr(template.reactions)
    # test gene order
    assert len(model.genes) == len(template.genes)
    assert set(get_attr(model.genes)) == set(get_attr(template.genes))
    assert get_attr(model.genes) == get_attr(template.genes)

