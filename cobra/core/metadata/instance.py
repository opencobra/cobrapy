from cobra.core.metadata import *
meta = MetaData()
meta["cvterms"] = {"is":[{"resources":["http://identifiers.org/chebi/CHEBI:17847"]}]}
print(meta["annotation"])
print(meta["annotation"]["chebi"])
meta["cvterms"]["isDescribedBy"] = [{"resources":["https://identifiers.org/pubmed/1111111", "https://identifiers.org/pubmed/111321"]}]
print(meta["annotation"])
meta["cvterms"]["is"].append({"resources":["https://identifiers.org/pubmed/1111111", "https://identifiers.org/pubmed/111321"]})
print(meta["annotation"])
meta["cvterms"]["is"].append({"resources":["https://someotherurl.org/pubmed/1111111"]})
