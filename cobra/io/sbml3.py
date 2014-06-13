from cobra import Metabolite, Reaction, DictList, Formula, Model

# import xml parsing libraries
from xml.dom import minidom  # only used for prettyprinting
try:
    from xml.etree import cElementTree as ET
except ImportError:
    from xml.etree import ElementTree as ET
try:
    from lxml.etree import parse
except ImportError:
    from xml.etree.ElementTree import parse
    # cElementTree does not take namespaces argument to .find()

ns = {"fbc": "http://www.sbml.org/sbml/level3/version1/fbc/version1",
      "sbml": "http://www.sbml.org/sbml/level3/version1/core",
      "cobra": "http://opencobra.github.io/cobra_sbml"}


def strnum(number):
    """Utility function to convert a number to a string"""
    return str(int(number)) if number == int(number) else str(number)


def metabolite_from_species(species):
    """create a metabolite from an sbml species"""
    met = Metabolite(get_attribute(species, "id", str).lstrip("_"))
    met.name = get_attribute(species, "name", str)
    met.compartment = get_attribute(species, "compartment", str)
    met.charge = get_attribute(species, "fbc:charge", int)
    met.formula = get_attribute(species, "fbc:chemicalFormula", Formula)
    return met


def parse_xml_into_model(xml):
    """create cobra model from xml etree"""
    xml_model = get_model(xml)
    metabolites = DictList(metabolite_from_species(i)
                           for i in list_metabolites(xml_model))

    fbc_operation_mapping = {"lessEqual": "upper_bound", "greaterEqual": "lower_bound"}
    fbc_bounds = {}
    for fbc_bound in list_flux_bounds(xml_model):
        reaction_id = get_attribute(fbc_bound, "fbc:reaction").lstrip("_")
        if reaction_id not in fbc_bounds:
            fbc_bounds[reaction_id] = {}
        bound_type = fbc_operation_mapping[get_attribute(fbc_bound, "fbc:operation")]
        fbc_bounds[reaction_id][bound_type] = get_attribute(fbc_bound, "fbc:value", float)

    reactions = []
    for sbml_reaction in list_reactions(xml_model):
        reaction = Reaction(get_attribute(sbml_reaction, "id", str).lstrip("_"))
        reaction.name = get_attribute(sbml_reaction, "name", str)
        if reaction.id in fbc_bounds:
            # todo: handle infinity
            bounds = fbc_bounds[reaction.id]
            if "lower_bound" in bounds:
                reaction.lower_bound = bounds["lower_bound"]
            if "upper_bound" in bounds:
                reaction.upper_bound = bounds["upper_bound"]
        # get annotation information
        extra = "sbml:annotation/cobra:extension/"
        subsystem = sbml_reaction.find(extra + "cobra:subsystem", namespaces=ns)
        if subsystem is not None:
            reaction.subsystem = subsystem.text
        gpr = sbml_reaction.find(extra + "cobra:geneAssociation", namespaces=ns)
        if gpr is not None:
            reaction.gene_reaction_rule = gpr.text
        reactions.append(reaction)

        stoichiometry = {}
        for species_reference in list_reactants(sbml_reaction):
            stoichiometry[get_attribute(species_reference, "species").lstrip("_")] = get_attribute(species_reference, "stoichiometry", float) * -1
        for species_reference in list_products(sbml_reaction):
            stoichiometry[get_attribute(species_reference, "species").lstrip("_")] = get_attribute(species_reference, "stoichiometry", float)

        object_stoichiometry = {}  # needs to have keys of metabolite objects, not ids
        for met_id in stoichiometry:
            object_stoichiometry[metabolites.get_by_id(met_id)] = stoichiometry[met_id]

        reaction.add_metabolites(object_stoichiometry)

    model = Model()
    model.id = get_attribute(xml_model, "id")
    # objective coefficients
    model.add_reactions(reactions)
    for sbml_objective in list_objectives(xml_model):
        rxn_id = get_attribute(sbml_objective, "fbc:reaction").lstrip("_")
        model.reactions.get_by_id(rxn_id).objective_coefficient = \
            get_attribute(sbml_objective, "fbc:coefficient")

    return model


def get_attribute(tag, attribute, type=lambda x: x):
    if ":" in attribute:
        split = attribute.split(":")
        attribute = "{%s}%s" % (ns[split[0]], split[1])
    return type(tag.attrib[attribute]) if attribute in tag.attrib else None


def get_model(xml):
    return xml.find("sbml:model", namespaces=ns)


def list_metabolites(xml_model):
    return xml_model.findall("sbml:listOfSpecies/sbml:species", namespaces=ns)


def list_objectives(xml_model):
    obj_list = xml_model.find("fbc:listOfObjectives", namespaces=ns)
    target_objective = get_attribute(obj_list, "fbc:activeObjective")
    all_objectives = [i for i in obj_list if get_attribute(i, "fbc:id") == target_objective]
    return all_objectives[0].findall("fbc:listOfFluxObjectives/fbc:fluxObjective", namespaces=ns)


def list_flux_bounds(xml_model):
    return xml_model.findall("fbc:listOfFluxBounds/fbc:fluxBound", namespaces=ns)


def list_reactions(xml_model):
    return xml_model.findall("sbml:listOfReactions/sbml:reaction", namespaces=ns)


def list_reactants(sbml_reaction):
    return sbml_reaction.findall("sbml:listOfReactants/sbml:speciesReference", namespaces=ns)


def list_products(sbml_reaction):
    return sbml_reaction.findall("sbml:listOfProducts/sbml:speciesReference", namespaces=ns)


def read_sbml_model(filename):
    xmlfile = parse(filename)
    xml = xmlfile.getroot()
    return parse_xml_into_model(xml)


def write_sbml_model(cobra_model, filename):
    xml = ET.Element("sbml", xmlns=ns["sbml"], level="3", version="1")
    xml.attrib["xmlns:fbc"] = ns["fbc"]
    xml.attrib["xmlns:cobra"] = ns["cobra"]

    xml.attrib["fbc:required"] = "false"
    xml_model = ET.SubElement(xml, "model")
    if cobra_model.id is not None:
        xml_model.attrib["id"] = cobra_model.id

    flux_bound_list = ET.SubElement(xml_model, "fbc:listOfFluxBounds")
    objectives_list_tmp = ET.SubElement(xml_model, "fbc:listOfObjectives")
    objectives_list_tmp.attrib["fbc:activeObjective"] = "obj"
    objectives_list_tmp = ET.SubElement(objectives_list_tmp, "fbc:objective")
    objectives_list_tmp.attrib.update({"fbc:id": "obj", "fbc:type": "maximize"})
    flux_objectives_list = ET.SubElement(objectives_list_tmp, "fbc:listOfFluxObjectives")

    # add in compartments
    compartmenst_list = ET.SubElement(xml_model, "listOfCompartments")
    compartments = cobra_model.compartments
    for compartment, name in compartments.items():
        ET.SubElement(compartmenst_list, "compartment", id=compartment, name=name, constant="true")

    # add in metabolites as species
    species_list = ET.SubElement(xml_model, "listOfSpecies")
    # Useless required SBML params
    extra_attributes = {"constant": "false", "boundaryCondition": "false",
                        "hasOnlySubstanceUnits": "false", "initialAmount": "NaN"}
    for met in cobra_model.metabolites:
        species = ET.SubElement(species_list, "species")
        attributes = species.attrib
        attributes["id"] = met.id if not met.id[0].isdigit() else "_" + met.id
        attributes["fbc:charge"] = str(met.charge)
        attributes["fbc:chemicalFormula"] = str(met.formula)
        attributes["name"] = met.name
        attributes["compartment"] = met.compartment
        attributes.update(extra_attributes)

    # add in reactions
    reactions_list = ET.SubElement(xml_model, "listOfReactions")
    for reaction in cobra_model.reactions:
        sbml_reaction = ET.SubElement(reactions_list, "reaction")
        attributes = sbml_reaction.attrib
        id = reaction.id if not reaction.id[0].isdigit() else "_" + reaction.id
        attributes["id"] = id
        attributes["name"] = reaction.name
        # Useless required SBML params
        attributes["fast"] = "false"
        attributes["reversible"] = str(reaction.lower_bound <= 0 and
                                       reaction.upper_bound >= 0).lower()
        # Add in annotation information (gpr and subsystem)
        # prevent empty strings
        strcheck = lambda x: None if x is None or len(x) == 0 else x
        gpr = strcheck(reaction.gene_reaction_rule)
        subsystem = strcheck(reaction.subsystem)
        if gpr is not None or subsystem is not None:
            annotation = ET.SubElement(sbml_reaction, "annotation")
            annotation = ET.SubElement(annotation, "cobra:extension")
        if gpr is not None:
            ET.SubElement(annotation, "cobra:geneAssociation").text = gpr
        if subsystem is not None:
            ET.SubElement(annotation, "cobra:subsystem").text = subsystem
        # add in bounds
        lb = ET.SubElement(flux_bound_list, "fbc:fluxBound")
        lb.attrib["fbc:reaction"] = id
        lb.attrib["fbc:value"] = strnum(reaction.lower_bound)
        lb.attrib["fbc:operation"] = "greaterEqual"
        ub = ET.SubElement(flux_bound_list, "fbc:fluxBound")
        ub.attrib["fbc:reaction"] = id
        ub.attrib["fbc:value"] = strnum(reaction.upper_bound)
        ub.attrib["fbc:operation"] = "lessEqual"
        # objective coefficient
        if reaction.objective_coefficient != 0:
            attr = {"fbc:reaction": id,
                    "fbc:coefficient": strnum(reaction.objective_coefficient)}
            ET.SubElement(flux_objectives_list, "fbc:fluxObjective").attrib.update(attr)

        # stoichiometry
        reactants = {}
        products = {}
        for metabolite, stoichiomety in reaction._metabolites.items():
            met_id = metabolite.id if not metabolite.id[0].isdigit() else "_" + metabolite.id
            if stoichiomety > 0:
                products[met_id] = strnum(stoichiomety)
            else:
                reactants[met_id] = strnum(-stoichiomety)
        if len(reactants) > 0:
            reactant_list = ET.SubElement(sbml_reaction, "listOfReactants")
            for met_id, stoichiomety in reactants.items():
                ET.SubElement(reactant_list, "speciesReference", species=met_id, stoichiometry=stoichiomety, constant="true")
        if len(products) > 0:
            product_list = ET.SubElement(sbml_reaction, "listOfProducts")
            for met_id, stoichiomety in products.items():
                ET.SubElement(product_list, "speciesReference", species=met_id, stoichiometry=stoichiomety, constant="true")
    # write out file
    minidom_xml = minidom.parseString(ET.tostring(xml))
    with open(filename, "w") as outfile:
        outfile.write(minidom_xml.toprettyxml(indent="  ", encoding="utf-8"))

if __name__ == "__main__":
    from cobra.test import create_test_model
    m = create_test_model()
    write_sbml_model(m, "test.xml")
    m2 = read_sbml_model("test.xml")
    assert len(m2.reactions) == len(m.reactions)
    assert len(m2.metabolites) == len(m.metabolites)
    from IPython import embed; embed()
