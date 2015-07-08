from collections import defaultdict
from warnings import warn
from decimal import Decimal
from ast import parse as ast_parse, Name, Or, And, BoolOp
from gzip import GzipFile
from bz2 import BZ2File

from six import iteritems

from .. import Metabolite, Reaction, Gene, Model
from ..core.Gene import parse_gpr


try:
    from xml.etree.cElementTree import parse, Element, SubElement, \
        ElementTree, tostring, register_namespace
except ImportError:
    from xml.etree.ElementTree import parse, Element, SubElement, \
        ElementTree, tostring, register_namespace

# use sbml level 2 from sbml.py (which uses libsbml). Eventually, it would
# be nice to use the libSBML converters directly instead.
try:
    import libsbml
except ImportError:
    libsbml = None
else:
    from .sbml import create_cobra_model_from_sbml_file as read_sbml2
    from .sbml import write_cobra_model_to_sbml_file as write_sbml2

try:
    from sympy import Basic
except:
    class Basic:
        pass


# deal with namespaces
namespaces = {"fbc": "http://www.sbml.org/sbml/level3/version1/fbc/version2",
              "sbml": "http://www.sbml.org/sbml/level3/version1/core"}
for key in namespaces:
    register_namespace(key, namespaces[key])

# XPATH query wrappers
fbc_prefix = "{" + namespaces["fbc"] + "}"
sbml_prefix = "{" + namespaces["sbml"] + "}"

SBML_DOT = "__SBML_DOT__"
# FBC TAGS
OR_TAG = "fbc:or"
AND_TAG = "fbc:and"
GENEREF_TAG = "fbc:geneProductRef"
GPR_TAG = "fbc:geneProductAssociation"
GENELIST_TAG = "fbc:listOfGeneProducts"
GENE_TAG = "fbc:geneProduct"
# XPATHS
BOUND_XPATH = "sbml:listOfParameters/sbml:parameter"  # TODO only get constant
COMPARTMENT_XPATH = "sbml:listOfCompartments/sbml:compartment"
GENES_XPATH = GENELIST_TAG + "/" + GENE_TAG
SPECIES_XPATH = "sbml:listOfSpecies/sbml:species[@boundaryCondition='%s']"
OBJECTIVES_XPATH = ("fbc:objective[@fbc:id='%s']/"
                    "fbc:listOfFluxObjectives/"
                    "fbc:fluxObjective")


def ns(query):
    """replace prefixes with namespace"""
    return query.replace("fbc:", fbc_prefix).replace("sbml:", sbml_prefix)


def get_attrib(tag, attribute, type=lambda x: x, require=False):
    value = tag.get(ns(attribute))
    if require and value is None:
        raise Exception("required attribute '%s' not found in tag '%s'" %
                        (attribute, tag.tag))
    return type(value) if value is not None else None


def set_attrib(xml, attribute_name, value):
    if value is None or value == "":
        return
    xml.set(ns(attribute_name), str(value))


# string utility functions
def clip(string, prefix):
    """clips a prefix from the beginning of a string if it exists

    >>> clip("R_pgi", "R_")
    "pgi"

    """
    return string[len(prefix):] if string.startswith(prefix) else string


def strnum(number):
    """Utility function to convert a number to a string"""
    if isinstance(number, (Decimal, Basic, str)):
        return str(number)
    s = "%.15g" % number
    return s.rstrip(".")


def construct_gpr_xml(parent, expression):
    """create gpr xml under parent node"""
    if isinstance(expression, BoolOp):
        op = expression.op
        if isinstance(op, And):
            new_parent = SubElement(parent, ns(AND_TAG))
        elif isinstance(op, Or):
            new_parent = SubElement(parent, ns(OR_TAG))
        else:
            raise Exception("unsupported operation " + op.__class__)
        for arg in expression.values:
            construct_gpr_xml(new_parent, arg)
    elif isinstance(expression, Name):
        gene_elem = SubElement(parent, ns(GENEREF_TAG))
        set_attrib(gene_elem, "fbc:geneProduct", "G_" + expression.id)
    else:
        raise Exception("unsupported operation  " + repr(expression))


def parse_xml_into_model(xml, number=float):
    xml_model = xml.find(ns("sbml:model"))
    if get_attrib(xml_model, "fbc:strict") != "true":
        warn('loading SBML model without fbc:strict="true"')

    model_id = get_attrib(xml_model, "id")
    model = Model(model_id)

    model.compartments = {c.get("id"): c.get("name") for c in
                          xml.findall(ns(COMPARTMENT_XPATH))}

    # add metabolites
    for species in xml_model.findall(ns(SPECIES_XPATH) % 'false'):
        met = Metabolite(clip(species.get("id"), "M_"))
        met.name = species.get("name")
        met.compartment = species.get("compartment")
        met.charge = get_attrib(species, "fbc:charge", int)
        met.formula = get_attrib(species, "fbc:chemicalFormula")
        model.add_metabolites([met])
    # Detect boundary metabolites - In case they have been mistakenly
    # added. They should not actually appear in a model
    boundary_metabolites = {clip(i.get("id"), "M_")
                            for i in xml_model.findall(ns(SPECIES_XPATH)
                                                       % 'true')}

    # add genes
    for sbml_gene in xml_model.findall(ns(GENES_XPATH)):
        gene_id = get_attrib(sbml_gene, "fbc:id").replace(SBML_DOT, ".")
        gene = Gene(clip(gene_id, "G_"))
        gene.name = get_attrib(sbml_gene, "fbc:label")
        model.genes.append(gene)

    or_tag = ns(OR_TAG)
    and_tag = ns(AND_TAG)
    generef_tag = ns(GENEREF_TAG)

    def process_gpr(sub_xml):
        """recursively convert gpr xml to a gpr string"""
        if sub_xml.tag == or_tag:
            return "( " + ' or '.join(process_gpr(i) for i in sub_xml) + " )"
        elif sub_xml.tag == and_tag:
            return "( " + ' and '.join(process_gpr(i) for i in sub_xml) + " )"
        elif sub_xml.tag == generef_tag:
            gene_id = get_attrib(sub_xml, "fbc:geneProduct", require=True)
            return clip(gene_id, "G_")
        else:
            raise Exception("unsupported tag " + sub_xml.tag)

    # options for get_attrib for numbers
    opts = {"type": number, "require": True}
    # TODO handle infinity
    bounds = {bound.get("id"): get_attrib(bound, "value", **opts)
              for bound in xml_model.findall(ns(BOUND_XPATH))}
    # add reactions
    reactions = []
    for sbml_reaction in xml_model.findall(
            ns("sbml:listOfReactions/sbml:reaction")):
        reaction = Reaction(clip(sbml_reaction.get("id"), "R_"))
        reaction.name = sbml_reaction.get("name")
        lb_id = get_attrib(sbml_reaction, "fbc:lowerFluxBound", require=True)
        reaction.lower_bound = bounds[lb_id]
        ub_id = get_attrib(sbml_reaction, "fbc:upperFluxBound", require=True)
        reaction.upper_bound = bounds[ub_id]
        reactions.append(reaction)

        stoichiometry = defaultdict(lambda: 0)
        for species_reference in sbml_reaction.findall(
                ns("sbml:listOfReactants/sbml:speciesReference")):
            met_name = clip(species_reference.get("species"), "M_")
            stoichiometry[met_name] -= \
                number(species_reference.get("stoichiometry"))
        for species_reference in sbml_reaction.findall(
                ns("sbml:listOfProducts/sbml:speciesReference")):
            met_name = clip(species_reference.get("species"), "M_")
            stoichiometry[met_name] += \
                get_attrib(species_reference, "stoichiometry", **opts)
        # needs to have keys of metabolite objects, not ids
        object_stoichiometry = {}
        for met_id in stoichiometry:
            if met_id in boundary_metabolites:
                continue
            try:
                metabolite = model.metabolites.get_by_id(met_id)
            except KeyError:
                warn("ignoring unknown metabolite '%s' in %s" %
                     (met_id, repr(reaction)))
                continue
            object_stoichiometry[metabolite] = stoichiometry[met_id]
        reaction.add_metabolites(object_stoichiometry)
        # set gene reaction rule
        gpr_xml = sbml_reaction.find(ns(GPR_TAG))
        if gpr_xml is not None and len(gpr_xml) != 1:
            warn("ignoring invalid geneAssocation for " + repr(reaction))
            gpr_xml = None
        gpr = process_gpr(gpr_xml[0]) if gpr_xml is not None else ''
        # remove outside parenthesis, if any
        if gpr.startswith("(") and gpr.endswith(")"):
            gpr = gpr[1:-1].strip()
        gpr = gpr.replace(SBML_DOT, ".")
        reaction.gene_reaction_rule = gpr
    model.add_reactions(reactions)

    # objective coefficients are handled after all reactions are added
    obj_list = xml_model.find(ns("fbc:listOfObjectives"))
    target_objective = get_attrib(obj_list, "fbc:activeObjective")
    obj_query = ns(OBJECTIVES_XPATH) % target_objective
    for sbml_objective in obj_list.findall(obj_query):
        rxn_id = clip(get_attrib(sbml_objective, "fbc:reaction"), "R_")
        model.reactions.get_by_id(rxn_id).objective_coefficient = \
            get_attrib(sbml_objective, "fbc:coefficient", type=number)

    return model


def model_to_xml(cobra_model, units=True):
    xml = Element("sbml", xmlns=namespaces["sbml"], level="3", version="1",
                  sboTerm="SBO:0000624")
    set_attrib(xml, "fbc:required", "false")
    xml_model = SubElement(xml, "model")
    set_attrib(xml_model, "fbc:strict", "true")
    if cobra_model.id is not None:
        xml_model.set("id", cobra_model.id)

    # if using units, add in mmol/gdw/hr
    if units:
        unit_def = SubElement(
            SubElement(xml_model, "listOfUnitDefinitions"),
            "unitDefinition", id="mmol_per_gDW_per_hr")
        list_of_units = SubElement(unit_def, "listOfUnits")
        SubElement(list_of_units, "unit", kind="mole", scale="-3",
                   multiplier="1", exponent="1")
        SubElement(list_of_units, "unit", kind="gram", scale="0",
                   multiplier="1", exponent="-1")
        SubElement(list_of_units, "unit", kind="second", scale="0",
                   multiplier="3600", exponent="-1")

    # create the element for the flux objective
    obj_list_tmp = SubElement(xml_model, ns("fbc:listOfObjectives"))
    set_attrib(obj_list_tmp, "fbc:activeObjective", "obj")
    obj_list_tmp = SubElement(obj_list_tmp, ns("fbc:objective"))
    set_attrib(obj_list_tmp, "fbc:id", "obj")
    set_attrib(obj_list_tmp, "fbc:type", "maximize")
    flux_objectives_list = SubElement(obj_list_tmp,
                                      ns("fbc:listOfFluxObjectives"))

    # create the element for the flux bound parameters
    parameter_list = SubElement(xml_model, "listOfParameters")
    param_attr = {"constant": "true"}
    if units:
        param_attr["units"] = "mmol_per_gDW_per_hr"
    # the most common bounds are the minimum, maxmium, and 0
    min_value = min(cobra_model.reactions.list_attr("lower_bound"))
    max_value = max(cobra_model.reactions.list_attr("upper_bound"))
    SubElement(parameter_list, "parameter", value=strnum(min_value),
               id="cobra_default_lb", sboTerm="SBO:0000626", **param_attr)
    SubElement(parameter_list, "parameter", value=strnum(max_value),
               id="cobra_default_ub", sboTerm="SBO:0000626", **param_attr)
    SubElement(parameter_list, "parameter", value="0",
               id="cobra_0_bound", sboTerm="SBO:0000626", **param_attr)

    def create_bound(reaction, bound_type):
        """returns the str id of the appropriate bound for the reaction

        The bound will also be created if necessary"""
        value = getattr(reaction, bound_type)
        if value == min_value:
            return "cobra_default_lb"
        elif value == 0:
            return "cobra_0_bound"
        elif value == max_value:
            return "cobra_default_ub"
        else:
            param_id = "R_" + reaction.id + "_" + bound_type
            SubElement(parameter_list, "parameter", id=param_id,
                       value=strnum(value), sboTerm="SBO:0000625",
                       **param_attr)
            return param_id

    # add in compartments
    compartmenst_list = SubElement(xml_model, "listOfCompartments")
    compartments = cobra_model.compartments
    for compartment, name in iteritems(compartments):
        SubElement(compartmenst_list, "compartment", id=compartment, name=name,
                   constant="true")

    # add in metabolites
    species_list = SubElement(xml_model, "listOfSpecies")
    for met in cobra_model.metabolites:
        species = SubElement(species_list, "species",
                             id="M_" + met.id,
                             compartment=str(met.compartment),
                             # Useless required SBML parameters
                             constant="false",
                             boundaryCondition="false",
                             hasOnlySubstanceUnits="false")
        set_attrib(species, "name", met.name)
        set_attrib(species, "fbc:charge", met.charge)
        set_attrib(species, "fbc:chemicalFormula", met.formula)

    # add in genes
    if len(cobra_model.genes) > 0:
        genes_list = SubElement(xml_model, ns(GENELIST_TAG))
    for gene in cobra_model.genes:
        gene_id = gene.id.replace(".", SBML_DOT)
        sbml_gene = SubElement(genes_list, ns(GENE_TAG))
        set_attrib(sbml_gene, "fbc:id", "G_" + gene_id)
        name = gene.name
        if name is None or len(name) == 0:
            name = gene.id
        set_attrib(sbml_gene, "fbc:label", gene_id)
        set_attrib(sbml_gene, "fbc:name", gene.name)

    # add in reactions
    reactions_list = SubElement(xml_model, "listOfReactions")
    for reaction in cobra_model.reactions:
        id = "R_" + reaction.id
        sbml_reaction = SubElement(
            reactions_list, "reaction",
            id=id,
            # Useless required SBML parameters
            fast="false",
            reversible=str(reaction.lower_bound < 0).lower())
        set_attrib(sbml_reaction, "name", reaction.name)
        # add in bounds
        set_attrib(sbml_reaction, "fbc:upperFluxBound",
                   create_bound(reaction, "upper_bound"))
        set_attrib(sbml_reaction, "fbc:lowerFluxBound",
                   create_bound(reaction, "lower_bound"))

        # objective coefficient
        if reaction.objective_coefficient != 0:
            objective = SubElement(flux_objectives_list,
                                   ns("fbc:fluxObjective"))
            set_attrib(objective, "fbc:reaction", id)
            set_attrib(objective, "fbc:coefficient",
                       strnum(reaction.objective_coefficient))

        # stoichiometry
        reactants = {}
        products = {}
        for metabolite, stoichiomety in iteritems(reaction._metabolites):
            met_id = "M_" + metabolite.id
            if stoichiomety > 0:
                products[met_id] = strnum(stoichiomety)
            else:
                reactants[met_id] = strnum(-stoichiomety)
        if len(reactants) > 0:
            reactant_list = SubElement(sbml_reaction, "listOfReactants")
            for met_id, stoichiomety in iteritems(reactants):
                SubElement(reactant_list, "speciesReference", species=met_id,
                           stoichiometry=stoichiomety, constant="true")
        if len(products) > 0:
            product_list = SubElement(sbml_reaction, "listOfProducts")
            for met_id, stoichiomety in iteritems(products):
                SubElement(product_list, "speciesReference", species=met_id,
                           stoichiometry=stoichiomety, constant="true")

        # gene reaction rule
        gpr = reaction.gene_reaction_rule
        if gpr is not None and len(gpr) > 0:
            gpr = gpr.replace(".", SBML_DOT)
            gpr_xml = SubElement(sbml_reaction, ns(GPR_TAG))
            try:
                parsed = parse_gpr(gpr)[0]
                construct_gpr_xml(gpr_xml, parsed.body)
            except Exception as e:
                print("failed on '%s' in %s" %
                      (reaction.gene_reaction_rule, repr(reaction)))
                raise e

    return xml


def read_sbml_model(filename, number=float, **kwargs):
    if filename.endswith(".gz"):
        with GzipFile(filename) as infile:
            xmlfile = parse(infile)
    elif filename.endswith(".bz2"):
        with BZ2File(filename) as infile:
            xmlfile = parse(infile)
    else:
        xmlfile = parse(filename)
    xml = xmlfile.getroot()
    if xml.get("level") != "3" or xml.get("version") != "1" or \
            get_attrib(xml, "fbc:required") is None:
        if libsbml is None:
            raise Exception("libSBML required for fbc < 2")
        return read_sbml2(filename, **kwargs)
    return parse_xml_into_model(xml, number=number)


def write_sbml_model(cobra_model, filename, use_fbc_package=True, **kwargs):
    if not use_fbc_package:
        write_sbml2(cobra_model, filename, use_fbc_package=False, **kwargs)
        return
    xml = model_to_xml(cobra_model, **kwargs)
    indent_xml(xml)
    if filename.endswith(".gz"):
        xmlfile = GzipFile(filename, "wb")
    elif filename.endswith(".bz2"):
        xmlfile = BZ2File(filename, "wb")
    else:
        xmlfile = open(filename, "wb")
    ElementTree(xml).write(xmlfile, encoding="UTF-8")
    xmlfile.close()


# inspired by http://effbot.org/zone/element-lib.htm#prettyprint
def indent_xml(elem, level=0):
    """indent xml for pretty printing"""
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent_xml(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
