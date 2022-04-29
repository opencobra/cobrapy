"""Test data storage and recovery using pickle."""

from pathlib import Path
from pickle import dump, load
from typing import Callable

import pytest

from cobra import Model


@pytest.mark.parametrize("load_function", [load])
def test_read_pickle(
    compare_models: Callable,
    data_directory: Path,
    mini_model: Model,
    load_function: Callable,
):
    """Test the reading of model from pickle."""
    if load_function is None:
        pytest.skip()

    with open(data_directory.joinpath("mini.pickle"), "rb") as infile:
        pickle_model = load_function(infile)

    assert compare_models(mini_model, pickle_model) is None


@pytest.mark.parametrize("dump_function", [dump])
def test_write_pickle(tmp_path: Path, mini_model: Model, dump_function: Callable):
    """Test the writing of model to pickle."""
    if dump_function is None:
        pytest.skip()

    output_file = tmp_path.joinpath("mini.pickle")
    with open(str(output_file), "wb") as outfile:
        dump_function(mini_model, outfile)

    assert output_file.exists()
