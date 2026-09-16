"""SBML-named exchanges must not be occluded by the boundary-species fill."""

from tempfile import NamedTemporaryFile

import libsbml
from cobra.io import read_sbml_model


def _document_with_named_and_unnamed_boundary() -> str:
    xmlns = libsbml.SBMLNamespaces(3, 1, "fbc", 2)
    doc = libsbml.SBMLDocument(xmlns)
    doc.setPackageRequired("fbc", False)
    model = doc.createModel()
    model.setId("named_boundary")
    fbc = model.getPlugin("fbc")
    if fbc is not None:
        fbc.setStrict(False)
    compartment = model.createCompartment()
    compartment.setId("e")
    compartment.setConstant(True)
    compartment.setSpatialDimensions(3)
    compartment.setSize(1.0)

    def add_species(sid: str, boundary: bool) -> None:
        species = model.createSpecies()
        species.setId(sid)
        species.setCompartment("e")
        species.setHasOnlySubstanceUnits(False)
        species.setBoundaryCondition(boundary)
        species.setConstant(False)

    add_species("A_e", True)
    add_species("B_e", True)

    lower = model.createParameter()
    lower.setId("lb")
    lower.setValue(-10.0)
    lower.setConstant(True)
    upper = model.createParameter()
    upper.setId("ub")
    upper.setValue(10.0)
    upper.setConstant(True)

    reaction = model.createReaction()
    reaction.setId("EX_A_e")
    reaction.setReversible(True)
    reaction.setFast(False)
    reaction_fbc = reaction.getPlugin("fbc")
    if reaction_fbc is not None:
        reaction_fbc.setLowerFluxBound("lb")
        reaction_fbc.setUpperFluxBound("ub")
    reactant = reaction.createReactant()
    reactant.setSpecies("A_e")
    reactant.setStoichiometry(1.0)
    reactant.setConstant(True)
    return libsbml.writeSBMLToString(doc)


def test_named_sbml_exchange_keeps_document_bounds() -> None:
    """A ListOfReactions exchange wins over a later synthetic EX_ id.

    A second boundary metabolite with no SBML exchange still receives
    the extra exchange suggested for unsupported boundary species.
    """
    xml = _document_with_named_and_unnamed_boundary()
    with NamedTemporaryFile("w", suffix=".xml", delete=False) as handle:
        handle.write(xml)
        path = handle.name
    model = read_sbml_model(path)
    named = model.reactions.get_by_id("EX_A_e")
    unnamed = model.reactions.get_by_id("EX_B_e")
    assert named.bounds == (-10.0, 10.0)
    assert unnamed.bounds == (-1000.0, 1000.0)
