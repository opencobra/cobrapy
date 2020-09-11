import gzip
import pathlib

import pytest

from cobra import Configuration
from cobra.io import BiGGModels, BioModels, load_model


@pytest.fixture(scope="module")
def mini_sbml(data_directory):
    """Provide a gzip-compressed SBML document."""
    with (pathlib.Path(data_directory) / "mini_cobra.xml").open(mode="rb") as handle:
        return gzip.compress(handle.read())


@pytest.fixture
def bigg_models(mini_sbml, mocker):
    """Provide a mocked BiGG Models repository interface."""
    result = mocker.Mock(spec_set=BiGGModels)
    result.get_sbml.return_value = mini_sbml
    return result


@pytest.fixture
def biomodels(mini_sbml, mocker):
    """Provide a mocked BioModels repository interface."""
    result = mocker.Mock(spec_set=BioModels)
    result.get_sbml.return_value = mini_sbml
    return result


def test_bigg_access(bigg_models):
    """Test that SBML would be retrieved from the BiGG Models repository."""
    load_model("e_coli_core", cache=False, repositories=[bigg_models])
    bigg_models.get_sbml.assert_called_once_with(model_id="e_coli_core")


def test_biomodels_access(biomodels):
    """Test that SBML would be retrieved from the BioModels repository."""
    load_model("BIOMD0000000633", cache=False, repositories=[biomodels])
    biomodels.get_sbml.assert_called_once_with(model_id="BIOMD0000000633")


@pytest.mark.raises(exception=RuntimeError, message="could not be found")
def test_unknown_model():
    """Expect that a not found error is raised (e2e)."""
    load_model("MODELWHO?", cache=False)


@pytest.mark.parametrize(
    "model_id, num_metabolites, num_reactions",
    [("e_coli_core", 72, 95), ("BIOMD0000000633", 50, 35)],
)
def test_remote_load(model_id: str, num_metabolites: int, num_reactions: int):
    """Test that sample models can be loaded from remote repositories (e2e)."""
    model = load_model(model_id, cache=False)
    assert len(model.metabolites) == num_metabolites
    assert len(model.reactions) == num_reactions


def test_cache(monkeypatch, tmp_path, bigg_models, biomodels):
    """Test that remote models are properly cached."""
    config = Configuration()
    monkeypatch.setattr(config, "cache_directory", tmp_path)
    remote_model = load_model("e_coli_core")
    cached_model = load_model("e_coli_core", repositories=[bigg_models, biomodels])
    bigg_models.get_sbml.assert_not_called()
    biomodels.get_sbml.assert_not_called()
    assert len(cached_model.metabolites) == len(remote_model.metabolites)
    assert len(cached_model.reactions) == len(remote_model.reactions)
