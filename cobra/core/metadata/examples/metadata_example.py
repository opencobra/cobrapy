from pathlib import Path
from cobra.core.metadata import *
from cobra.core.species import Species
import pytest
import json
from pprint import pprint
from cobra.core.metadata.cvterm import CVTerms


def example_metadata():
    s = Species()

    for json_example in ["cvterms_flat", "cvterms_nested"]:

        with open(Path(__file__).parent / f"{json_example}.json", "r") as f_cvterms:
            cvterms_data = json.load(f_cvterms)
            print("-" * 80)
            pprint(cvterms_data)
            print("-" * 80)

            # FIXME:
            cvterms = CVTerms(cvterms_data)
            print(cvterms)
            meta = MetaData(cvterms=cvterms)
    return


    # history

    # keyValuePair

    # cvterms
    # FIXME: update example
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


if __name__ == "__main__":
    example_metadata()
    # example_annotation()
