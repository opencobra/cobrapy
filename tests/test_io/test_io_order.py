"""Test functionalities of I/O in an usage order."""

import logging
from operator import attrgetter
from pathlib import Path
from random import sample
from typing import Iterable, List

import pytest

from cobra import DictList, Model
from cobra import io as cio


logger = logging.getLogger(__name__)


@pytest.fixture(scope="module")
def tmp_path_order(tmp_path_factory: Path) -> Path:
    """Temporary path for I/O order tests."""
    return tmp_path_factory.mktemp("model_order")


@pytest.fixture(scope="module")
def minimized_shuffle(small_model: Model) -> Model:
    """Generate a minimal shuffled model for I/O order tests."""
    model = small_model.copy()
    chosen = sample(list(set(model.reactions) - set(model.exchanges)), 10)
    new = Model("minimized_shuffle")
    new.add_reactions(chosen)
    logger.debug(
        f"'{new.id}' has {new.metabolites} metabolites, {new.reactions} reactions, and "
        f"{new.genes} genes."
    )
    return new


@pytest.fixture(scope="module")
def minimized_sorted(minimized_shuffle: Model) -> Model:
    """Generate a minimal sorted model for I/O order tests."""
    model = minimized_shuffle.copy()
    model.id = "minimized_sorted"
    model.metabolites = DictList(sorted(model.metabolites, key=attrgetter("id")))
    model.genes = DictList(sorted(model.genes, key=attrgetter("id")))
    model.reactions = DictList(sorted(model.reactions, key=attrgetter("id")))
    return model


@pytest.fixture(scope="module")
def minimized_reverse(minimized_shuffle: Model) -> Model:
    """Generate a minimal reversed model for I/O order tests."""
    model = minimized_shuffle.copy()
    model.id = "minimized_reverse"
    model.metabolites = DictList(
        sorted(model.metabolites, key=attrgetter("id"), reverse=True)
    )
    model.genes = DictList(sorted(model.genes, key=attrgetter("id"), reverse=True))
    model.reactions = DictList(
        sorted(model.reactions, key=attrgetter("id"), reverse=True)
    )
    return model


@pytest.fixture(
    scope="module",
    params=["minimized_shuffle", "minimized_reverse", "minimized_sorted"],
)
def template(
    request: pytest.FixtureRequest,
    minimized_shuffle: Model,
    minimized_reverse: Model,
    minimized_sorted: Model,
) -> Model:
    """Return the cobra Model instances found in the current local symbol table."""
    return locals()[request.param]


@pytest.fixture(scope="module", params=["metabolites", "reactions", "genes"])
def attribute(request: pytest.FixtureRequest) -> str:
    """Return the parameter passed."""
    return request.param


def _get_ids(iterable: Iterable) -> List[str]:
    """Get IDs for elements in `iterable`."""
    return [x.id for x in iterable]


@pytest.mark.parametrize(
    "read, write, ext",
    [
        ("read_sbml_model", "write_sbml_model", ".xml"),
        ("load_json_model", "save_json_model", ".json"),
        ("load_yaml_model", "save_yaml_model", ".yml"),
    ],
)
def test_io_order(
    attribute: str,
    read: str,
    write: str,
    ext: str,
    template: Model,
    tmp_path_order: Path,
) -> None:
    """Test loading and saving of models in order.

    Parameters
    ----------
    attribute : str
        The attribute of cobra Model to access.
    read : str
        The function name for loading model, defined as string.
    write : str
        The function name for saving model, defined as string.
    ext : str
        The extension of the file format for loading and saving model.
    template : cobra.Model
        The cobra Model instance to load and save.
    tmp_path_order : pathlib.Path
        The folder path for storing I/O order test files.

    """
    read = getattr(cio, read)
    write = getattr(cio, write)
    file_path = tmp_path_order / f"template{ext}"
    write(template, str(file_path.resolve()))
    model = read(str(file_path.resolve()))
    model_elements = _get_ids(getattr(model, attribute))
    template_elements = _get_ids(getattr(template, attribute))
    assert len(model_elements) == len(template_elements)
    assert set(model_elements) == set(template_elements)
    assert model_elements == template_elements
