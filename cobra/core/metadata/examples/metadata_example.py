from pathlib import Path
from cobra.core.metadata import *
from cobra.core.species import Species
import pytest
import json
from pprint import pprint
from cobra.core.metadata.cvterm import CVTerms


def example_metadata():
    s = Species()

    for json_example in ["cvterms_alternative", "cvterms_nested"]:

        with open(Path(__file__).parent / f"{json_example}.json", "r") as f_cvterms:
            cvterms_data = json.load(f_cvterms)
            print("-" * 80)
            print("{} Annotation: ".format(json_example))
            pprint(cvterms_data)


            # s.annotation = cvterms_data
            # print("Using species: ")
            # print("Direct annotation:")
            # print(s.annotation)
            # print("CVTerms")
            # pprint(s.annotation.cvterms)

            # print("Using CVTerm:")
            # cvterms = CVTerms(cvterms_data)
            # pprint(dict(cvterms), width=1)

            print("Reading using Metadata Class:\n")
            meta = MetaData(cvterms_data)
            print("Printing the direct annotation\n")
            print(meta,"\n")
            print("Printing the collection of cvterms\n")
            print(meta.cvterms,"\n")
            print("adding new valid cvterm: {}\n".format("http://identifiers.org/chebi/CHEBI:12newchebiterm"))
            cvt = CVTerm(Qualifier["bqb_hasPart"], "http://identifiers.org/chebi/CHEBI:12newchebiterm")
            meta.add_cvterm(cvt, 0)
            print("Printing the direct annotation\n")
            print(meta,"\n")
            print("Printing the collection of cvterms\n")
            print(meta.cvterms,"\n")

            cvt2 = CVTerm(Qualifier["bqb_hasPart"], "http://notidentifiers.org/uniprot/newuniport")
            meta.add_cvterm(cvt2, 0)
            print("Printing the direct annotation\n")
            print(meta,"\n")
            print("Printing the collection of cvterms\n")
            print(meta.cvterms,"\n")
            cvt2 = CVTerm(Qualifier["bqb_hasPart"], "http://notidentifiers.org/uniprot/newuniport")
            print("-" * 80)

    return


    # history

    # keyValuePair

    # cvterms
    # FIXME: update example

    # meta["cvterms"] = {
    #     "is": [{"resources": ["http://identifiers.org/chebi/CHEBI:17847"]}]}
    # print(meta["annotation"])
    # print(meta["annotation"]["chebi"])
    # meta["cvterms"]["isDescribedBy"] = [{"resources": [
    #     "https://identifiers.org/pubmed/1111111",
    #     "https://identifiers.org/pubmed/111321"]}]
    # meta["cvterms"]["is"].append({"resources": [
    #     "https://identifiers.org/pubmed/1111111",
    #     "https://identifiers.org/pubmed/111321"]})
    # print(meta["annotation"])
    # meta["cvterms"]["is"].append(
    #     {"resources": ["https://someotherurl.org/pubmed/1111111"]})
    #
    # s.annotation = meta


if __name__ == "__main__":
    example_metadata()
    # example_annotation()
