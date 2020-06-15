from cobra.core.metadata import *
from cobra.core.species import Species
import pytest


def test_annotation():
    s = Species()
    print(s.annotation)
    s.annotation["chebi"] = ["1234", "23423432"]
    s.annotation["sbo"] = ["SBO123"]
    print(s.annotation)

    assert "chebi" in s.annotation
    assert "sbo" in s.annotation
    assert len(s.annotation) == 2
    for key in ["keys", "items", "values"]:
        assert hasattr(s.annotation, key)

    # assert 0 == 1


def test_metadata():
    s = Species()

    meta = MetaData()

    # history

    # keyValuePair

    # cvterms
    meta["cvterms"] = {
        "is": [{"resources": ["http://identifiers.org/chebi/CHEBI:17847"]}]}
    print(meta["annotation"])
    print(meta["annotation"]["chebi"])
    meta["cvterms"]["isDescribedBy"] = [{"resources": [
        "https://identifiers.org/pubmed/1111111",
        "https://identifiers.org/pubmed/111321"]}]
    meta["cvterms"]["is"].append({"resources": [
        "https://identifiers.org/pubmed/1111111",
        "https://identifiers.org/pubmed/111321"]})
    print(meta["annotation"])
    meta["cvterms"]["is"].append(
        {"resources": ["https://someotherurl.org/pubmed/1111111"]})

    s.annotation = meta

# if __name__ == "__main__":
    # example_metadata()
    # example_annotation()
#    pass
