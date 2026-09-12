"""Test model loading from local and remote model repositores."""

import gzip
import os
from pathlib import Path
from typing import TYPE_CHECKING, Callable
from unittest.mock import Mock

import pytest

from cobra import Configuration
from cobra.io import BiGGModels, BioModels, load_model
from cobra.io.web import biomodels_repository


if TYPE_CHECKING:
    from pytest_mock import MockerFixture

    from cobra import Model


@pytest.fixture(scope="module")
def mini_sbml(cobra_data_directory: Path) -> bytes:
    """Provide a gzip-compressed SBML document."""
    with (cobra_data_directory / "mini_cobra.xml").open(mode="rb") as handle:
        return gzip.compress(handle.read())


@pytest.fixture
def bigg_models(mini_sbml: bytes, mocker: "MockerFixture") -> Mock:
    """Provide a mocked BiGG Models repository interface."""
    result = mocker.Mock(spec_set=BiGGModels)
    result.get_sbml.return_value = mini_sbml
    return result


@pytest.fixture
def biomodels(mini_sbml: bytes, mocker: "MockerFixture") -> Mock:
    """Provide a mocked BioModels repository interface."""
    result = mocker.Mock(spec_set=BioModels)
    result.get_sbml.return_value = mini_sbml
    return result


def test_bigg_access(bigg_models: Mock) -> None:
    """Test that SBML would be retrieved from the BiGG Models repository.

    Parameters
    ----------
    bigg_models : unittest.mock.Mock
        The mocked object for BiGG model respository.

    """
    load_model("e_coli_core", cache=False, repositories=[bigg_models])
    bigg_models.get_sbml.assert_called_once_with(model_id="e_coli_core")


@pytest.mark.xfail(
    condition=os.environ.get("GITHUB_ACTIONS"),
    reason="BioModels often blocks Github Actions",
)
def test_biomodels_access(biomodels: Mock) -> None:
    """Test that SBML would be retrieved from the BioModels repository.

    Parameters
    ----------
    biomodels : unittest.mock.Mock
        The mocked object for BioModels model respository.

    """
    load_model("BIOMD0000000633", cache=False, repositories=[biomodels])
    biomodels.get_sbml.assert_called_once_with(model_id="BIOMD0000000633")


def test_unknown_model() -> None:
    """Expect that a not found error is raised (e2e)."""
    with pytest.raises(RuntimeError):
        load_model("MODELWHO?", cache=False)


@pytest.mark.parametrize(
    "main_files",
    [
        pytest.param(
            [{"name": "253_2019_9630_MOESM3_ESM.omex", "fileSize": "106359"}],
            id="omex-only",
        ),
        pytest.param([], id="empty-listing"),
    ],
)
def test_biomodels_without_sbml_document(
    mocker: "MockerFixture", main_files: list
) -> None:
    """Expect a clear error when BioModels lists no SBML document.

    Some models are deposited only as COMBINE (OMEX) archives, so the file
    listing contains no ``.xml`` entry. The download must not be attempted,
    because the archive bytes are not an SBML document and would surface much
    later as a ``UnicodeDecodeError``.

    Parameters
    ----------
    mocker : pytest_mock.MockerFixture
        The mocking fixture.
    main_files : list
        The ``main`` section of the mocked BioModels files response.

    """
    response = mocker.Mock()
    response.json.return_value = {"additional": [], "main": main_files}
    mocker.patch.object(biomodels_repository.httpx, "get", return_value=response)
    stream = mocker.patch.object(biomodels_repository.httpx, "stream")

    with pytest.raises(RuntimeError, match="Could not find an SBML document"):
        BioModels().get_sbml(model_id="BIOMD0000001089")

    stream.assert_not_called()


@pytest.mark.parametrize(
    "model_id, num_metabolites, num_reactions",
    [("e_coli_core", 72, 95), ("BIOMD0000000633", 50, 35)],
)
def test_remote_load(model_id: str, num_metabolites: int, num_reactions: int) -> None:
    """Test that sample models can be loaded from remote repositories (e2e).

    Parameters
    ----------
    model_id : str
        The ID of the model.
    num_metabolites : int
        The total number of metabolites in the model having ID `model_id`.
    num_reactions : int
        The total number of reactions in the model having ID `model_id`.

    """
    if model_id.startswith("BIOMD"):
        pytest.xfail(reason="BioModels often blocks Github Actions.")

    model = load_model(model_id, cache=False)
    assert len(model.metabolites) == num_metabolites
    assert len(model.reactions) == num_reactions


@pytest.mark.skip("caching currently removed because of a CVE in diskcache")
def test_cache(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, bigg_models: Mock, biomodels: Mock
) -> None:
    """Test that remote models are properly cached.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        The monkeypatch-ing object.
    tmp_path : pathlib.Path
        The path to the temporary test assets store.
    bigg_models : Mock
        The mocked object for BiGG model respository.
    biomodels : unittest.mock.Mock
        The mocked object for BioModels model respository.

    """
    config = Configuration()
    monkeypatch.setattr(config, "cache_directory", tmp_path)
    remote_model = load_model("e_coli_core")
    cached_model = load_model("e_coli_core", repositories=[bigg_models, biomodels])
    bigg_models.get_sbml.assert_not_called()
    biomodels.get_sbml.assert_not_called()
    assert len(cached_model.metabolites) == len(remote_model.metabolites)
    assert len(cached_model.reactions) == len(remote_model.reactions)


def test_local_load(model: "Model", compare_models: Callable) -> None:
    """Test model loading from local repository.

    Parameters
    ----------
    model : cobra.Model
        The model to compare local loading against.
    compare_models : Callable
        A callable object to compare local loading.

    """
    model_local = load_model("textbook")
    compare_models(model, model_local)
