# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
from ast import And, BoolOp, Name, Or
from bz2 import BZ2File
from collections import defaultdict
from decimal import Decimal
from gzip import GzipFile
from tempfile import NamedTemporaryFile
from warnings import catch_warnings, simplefilter, warn

from six import iteritems, string_types

from cobra.core import Configuration, Gene, Metabolite, Model, Reaction
from cobra.core.gene import parse_gpr
from cobra.manipulation.modify import _renames
from cobra.manipulation.validate import check_metabolite_compartment_formula
from cobra.util.solver import set_objective


try:
    from lxml.etree import (
        parse, Element, SubElement, ElementTree, register_namespace,
        ParseError, XPath)

    _with_lxml = True
except ImportError:
    _with_lxml = False
    try:
        from xml.etree.cElementTree import (
            parse, Element, SubElement, ElementTree, register_namespace,
            ParseError)
    except ImportError:
        XPath = None
        from xml.etree.ElementTree import (
            parse, Element, SubElement, ElementTree, register_namespace,
            ParseError)

# use sbml level 2 from sbml.py (which uses libsbml). Eventually, it would
# be nice to use the libSBML converters directly instead.
try:
    import libsbml
except ImportError:
    libsbml = None
else:
    from cobra.io.sbml import create_cobra_model_from_sbml_file as read_sbml2
    from cobra.io.sbml import write_cobra_model_to_sbml_file as write_sbml2

try:
    from optlang.symbolics import Basic
except ImportError:
    class Basic:
        pass


CONFIGURATION = Configuration()

# deal with namespaces
namespaces = {"fbc": "http://www.sbml.org/sbml/level3/version1/fbc/version2",
              "sbml": "http://www.sbml.org/sbml/level3/version1/core",
              "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
              "bqbiol": "http://biomodels.net/biology-qualifiers/"}
for key in namespaces:
    register_namespace(key, namespaces[key])


def ns(query):
    """replace prefixes with namespace"""
    for prefix, uri in iteritems(namespaces):
        query = query.replace(prefix + ":", "{" + uri + "}")
    return query


# XPATH query wrappers
fbc_prefix = "{" + namespaces["fbc"] + "}"
sbml_prefix = "{" + namespaces["sbml"] + "}"

SBML_DOT = "__SBML_DOT__"
# FBC TAGS
OR_TAG = ns("fbc:or")
AND_TAG = ns("fbc:and")
GENEREF_TAG = ns("fbc:geneProductRef")
GPR_TAG = ns("fbc:geneProductAssociation")
GENELIST_TAG = ns("fbc:listOfGeneProducts")
GENE_TAG = ns("fbc:geneProduct")
# XPATHS
BOUND_XPATH = ns("sbml:listOfParameters/sbml:parameter[@value]")
COMPARTMENT_XPATH = ns("sbml:listOfCompartments/sbml:compartment")
GENES_XPATH = GENELIST_TAG + "/" + GENE_TAG
SPECIES_XPATH = ns("sbml:listOfSpecies/sbml:species[@boundaryCondition='%s']")
OBJECTIVES_XPATH = ns("fbc:objective[@fbc:id='%s']/"
                      "fbc:listOfFluxObjectives/"
                      "fbc:fluxObjective")
LONG_SHORT_DIRECTION = {'maximize': 'max', 'minimize': 'min'}
SHORT_LONG_DIRECTION = {'min': 'minimize', 'max': 'maximize'}

if _with_lxml:
    RDF_ANNOTATION_XPATH = ("sbml:annotation/rdf:RDF/"
                            "rdf:Description[@rdf:about=$metaid]/"
                            "*[self::bqbiol:isEncodedBy or self::bqbiol:is]/"
                            "rdf:Bag/rdf:li/@rdf:resource")
    extract_rdf_annotation = XPath(RDF_ANNOTATION_XPATH, namespaces=namespaces,
                                   smart_strings=False)
else:
    RDF_ANNOTATION_XPATH = ns("sbml:annotation/rdf:RDF/"
                              "rdf:Description[@rdf:about='%s']/"
                              "bqbiol:isEncodedBy/"
                              "rdf:Bag/rdf:li[@rdf:resource]")

    def extract_rdf_annotation(sbml_element, metaid):
        search_xpath = RDF_ANNOTATION_XPATH % metaid
        for i in sbml_element.iterfind(search_xpath):
            yield get_attrib(i, "rdf:resource")
        for i in sbml_element.iterfind(search_xpath.replace(
                "isEncodedBy", "is")):
            yield get_attrib(i, "rdf:resource")


class CobraSBMLError(Exception):
    pass


def get_attrib(tag, attribute, type=lambda x: x, require=False):
    value = tag.get(ns(attribute))
    if require and value is None:
        msg = "required attribute '%s' not found in tag '%s'" % \
              (attribute, tag.tag)
        if tag.get("id") is not None:
            msg += " with id '%s'" % tag.get("id")
        elif tag.get("name") is not None:
            msg += " with name '%s'" % tag.get("name")
        raise CobraSBMLError(msg)
    return type(value) if value is not None else None


def set_attrib(xml, attribute_name, value):
    if value is None or value == "":
        return
    xml.set(ns(attribute_name), str(value))


def parse_stream(filename):
    """parses filename or compressed stream to xml"""
    try:
        if hasattr(filename, "read"):
            return parse(filename)
        elif filename.endswith(".gz"):
            with GzipFile(filename) as infile:
                return parse(infile)
        elif filename.endswith(".bz2"):
            with BZ2File(filename) as infile:
                return parse(infile)
        else:
            return parse(filename)
    except ParseError as e:
        raise CobraSBMLError("Malformed XML file: " + str(e))


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
            new_parent = SubElement(parent, AND_TAG)
        elif isinstance(op, Or):
            new_parent = SubElement(parent, OR_TAG)
        else:
            raise Exception("unsupported operation " + op.__class__)
        for arg in expression.values:
            construct_gpr_xml(new_parent, arg)
    elif isinstance(expression, Name):
        gene_elem = SubElement(parent, GENEREF_TAG)
        set_attrib(gene_elem, "fbc:geneProduct", "G_" + expression.id)
    else:
        raise Exception("unsupported operation  " + repr(expression))


def annotate_cobra_from_sbml(cobra_element, sbml_element):
    sbo_term = sbml_element.get("sboTerm")
    if sbo_term is not None:
        cobra_element.annotation["sbo"] = sbo_term
    meta_id = get_attrib(sbml_element, "metaid")
    if meta_id is None:
        return
    annotation = cobra_element.annotation
    for uri in extract_rdf_annotation(sbml_element, metaid="#" + meta_id):
        if not uri.startswith("http://identifiers.org/"):
            warn("%s does not start with http://identifiers.org/" % uri)
            continue
        try:
            provider, identifier = uri[23:].split("/", 1)
        except ValueError:
            warn("%s does not conform to http://identifiers.org/provider/id"
                 % uri)
            continue
        # handle multiple id's in the same database
        if provider in annotation:
            # make into a list if necessary
            if isinstance(annotation[provider], string_types):
                annotation[provider] = [annotation[provider]]
            annotation[provider].append(identifier)
        else:
            cobra_element.annotation[provider] = identifier


def annotate_sbml_from_cobra(sbml_element, cobra_element):
    if len(cobra_element.annotation) == 0:
        return
    # get the id so we can set the metaid
    tag = sbml_element.tag
    if tag.startswith(sbml_prefix) or tag[0] != "{":
        prefix = ""
    elif tag.startswith(fbc_prefix):
        prefix = fbc_prefix
    else:
        raise ValueError("Can not annotate " + repr(sbml_element))
    id = sbml_element.get(prefix + "id")
    if len(id) == 0:
        raise ValueError("%s does not have id set" % repr(sbml_element))
    set_attrib(sbml_element, "metaid", id)
    annotation = SubElement(sbml_element, ns("sbml:annotation"))
    rdf_desc = SubElement(SubElement(annotation, ns("rdf:RDF")),
                          ns("rdf:Description"))
    set_attrib(rdf_desc, "rdf:about", "#" + id)
    bag = SubElement(SubElement(rdf_desc, ns("bqbiol:is")),
                     ns("rdf:Bag"))
    for provider, identifiers in sorted(iteritems(cobra_element.annotation)):
        if provider == "sbo":
            set_attrib(sbml_element, "sboTerm", identifiers)
            continue
        if isinstance(identifiers, string_types):
            identifiers = (identifiers,)
        for identifier in identifiers:
            li = SubElement(bag, ns("rdf:li"))
            set_attrib(li, "rdf:resource", "http://identifiers.org/%s/%s" %
                       (provider, identifier))


def parse_xml_into_model(xml, number=float):
    xml_model = xml.find(ns("sbml:model"))
    if get_attrib(xml_model, "fbc:strict") != "true":
        warn('loading SBML model without fbc:strict="true"')

    model_id = get_attrib(xml_model, "id")
    model = Model(model_id)
    model.name = xml_model.get("name")

    model.compartments = {c.get("id"): c.get("name") for c in
                          xml_model.findall(COMPARTMENT_XPATH)}
    # add metabolites
    for species in xml_model.findall(SPECIES_XPATH % 'false'):
        met = get_attrib(species, "id", require=True)
        met = Metabolite(clip(met, "M_"))
        met.name = species.get("name")
        annotate_cobra_from_sbml(met, species)
        met.compartment = species.get("compartment")
        met.charge = get_attrib(species, "fbc:charge", int)
        met.formula = get_attrib(species, "fbc:chemicalFormula")
        model.add_metabolites([met])
    # Detect boundary metabolites - In case they have been mistakenly
    # added. They should not actually appear in a model
    boundary_metabolites = {clip(i.get("id"), "M_")
                            for i in xml_model.findall(SPECIES_XPATH % 'true')}

    # add genes
    for sbml_gene in xml_model.iterfind(GENES_XPATH):
        gene_id = get_attrib(sbml_gene, "fbc:id").replace(SBML_DOT, ".")
        gene = Gene(clip(gene_id, "G_"))
        gene.name = get_attrib(sbml_gene, "fbc:name")
        if gene.name is None:
            gene.name = get_attrib(sbml_gene, "fbc:label")
        annotate_cobra_from_sbml(gene, sbml_gene)
        model.genes.append(gene)

    def process_gpr(sub_xml):
        """recursively convert gpr xml to a gpr string"""
        if sub_xml.tag == OR_TAG:
            return "( " + ' or '.join(process_gpr(i) for i in sub_xml) + " )"
        elif sub_xml.tag == AND_TAG:
            return "( " + ' and '.join(process_gpr(i) for i in sub_xml) + " )"
        elif sub_xml.tag == GENEREF_TAG:
            gene_id = get_attrib(sub_xml, "fbc:geneProduct", require=True)
            return clip(gene_id, "G_")
        else:
            raise Exception("unsupported tag " + sub_xml.tag)

    bounds = {bound.get("id"): get_attrib(bound, "value", type=number)
              for bound in xml_model.iterfind(BOUND_XPATH)}
    # add reactions
    reactions = []
    for sbml_reaction in xml_model.iterfind(
            ns("sbml:listOfReactions/sbml:reaction")):
        reaction = get_attrib(sbml_reaction, "id", require=True)
        reaction = Reaction(clip(reaction, "R_"))
        reaction.name = sbml_reaction.get("name")
        annotate_cobra_from_sbml(reaction, sbml_reaction)
        lb_id = get_attrib(sbml_reaction, "fbc:lowerFluxBound", require=True)
        ub_id = get_attrib(sbml_reaction, "fbc:upperFluxBound", require=True)
        try:
            reaction.upper_bound = bounds[ub_id]
            reaction.lower_bound = bounds[lb_id]
        except KeyError as e:
            raise CobraSBMLError("No constant bound with id '%s'" % str(e))
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
                get_attrib(species_reference, "stoichiometry",
                           type=number, require=True)
        # needs to have keys of metabolite objects, not ids
        object_stoichiometry = {}
        for met_id in stoichiometry:
            if met_id in boundary_metabolites:
                warn("Boundary metabolite '%s' used in reaction '%s'" %
                     (met_id, reaction.id))
                continue
            try:
                metabolite = model.metabolites.get_by_id(met_id)
            except KeyError:
                warn("ignoring unknown metabolite '%s' in reaction %s" %
                     (met_id, reaction.id))
                continue
            object_stoichiometry[metabolite] = stoichiometry[met_id]
        reaction.add_metabolites(object_stoichiometry)
        # set gene reaction rule
        gpr_xml = sbml_reaction.find(GPR_TAG)
        if gpr_xml is not None and len(gpr_xml) != 1:
            warn("ignoring invalid geneAssociation for " + repr(reaction))
            gpr_xml = None
        gpr = process_gpr(gpr_xml[0]) if gpr_xml is not None else ''
        # remove outside parenthesis, if any
        if gpr.startswith("(") and gpr.endswith(")"):
            gpr = gpr[1:-1].strip()
        gpr = gpr.replace(SBML_DOT, ".")
        reaction.gene_reaction_rule = gpr
    try:
        model.add_reactions(reactions)
    except ValueError as e:
        warn(str(e))

    # objective coefficients are handled after all reactions are added
    obj_list = xml_model.find(ns("fbc:listOfObjectives"))
    if obj_list is None:
        warn("listOfObjectives element not found")
        return model
    target_objective_id = get_attrib(obj_list, "fbc:activeObjective")
    target_objective = obj_list.find(
        ns("fbc:objective[@fbc:id='{}']".format(target_objective_id)))
    obj_direction_long = get_attrib(target_objective, "fbc:type")
    obj_direction = LONG_SHORT_DIRECTION[obj_direction_long]

    obj_query = OBJECTIVES_XPATH % target_objective_id
    coefficients = {}
    for sbml_objective in obj_list.findall(obj_query):
        rxn_id = clip(get_attrib(sbml_objective, "fbc:reaction"), "R_")
        try:
            objective_reaction = model.reactions.get_by_id(rxn_id)
        except KeyError:
            raise CobraSBMLError("Objective reaction '%s' not found" % rxn_id)
        try:
            coefficients[objective_reaction] = get_attrib(
                sbml_objective, "fbc:coefficient", type=number)
        except ValueError as e:
            warn(str(e))
    set_objective(model, coefficients)
    model.solver.objective.direction = obj_direction
    return model


def model_to_xml(cobra_model, units=True):
    xml = Element("sbml", xmlns=namespaces["sbml"], level="3", version="1",
                  sboTerm="SBO:0000624")
    set_attrib(xml, "fbc:required", "false")
    xml_model = SubElement(xml, "model")
    set_attrib(xml_model, "fbc:strict", "true")
    if cobra_model.id is not None:
        xml_model.set("id", cobra_model.id)
    if cobra_model.name is not None:
        xml_model.set("name", cobra_model.name)

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
    set_attrib(obj_list_tmp, "fbc:type",
               SHORT_LONG_DIRECTION[cobra_model.objective.direction])
    flux_objectives_list = SubElement(obj_list_tmp,
                                      ns("fbc:listOfFluxObjectives"))

    # create the element for the flux bound parameters
    parameter_list = SubElement(xml_model, "listOfParameters")
    param_attr = {"constant": "true"}
    if units:
        param_attr["units"] = "mmol_per_gDW_per_hr"
    # the most common bounds are the minimum, maximum, and 0
    if len(cobra_model.reactions) > 0:
        min_value = min(cobra_model.reactions.list_attr("lower_bound"))
        max_value = max(cobra_model.reactions.list_attr("upper_bound"))
    else:
        min_value = CONFIGURATION.lower_bound
        max_value = CONFIGURATION.upper_bound

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
    compartments_list = SubElement(xml_model, "listOfCompartments")
    compartments = cobra_model.compartments
    for compartment, name in iteritems(compartments):
        SubElement(compartments_list, "compartment", id=compartment, name=name,
                   constant="true")

    # add in metabolites
    species_list = SubElement(xml_model, "listOfSpecies")
    for met in cobra_model.metabolites:
        species = SubElement(species_list, "species",
                             id="M_" + met.id,
                             # Useless required SBML parameters
                             constant="false",
                             boundaryCondition="false",
                             hasOnlySubstanceUnits="false")
        set_attrib(species, "name", met.name)
        annotate_sbml_from_cobra(species, met)
        set_attrib(species, "compartment", met.compartment)
        set_attrib(species, "fbc:charge", met.charge)
        set_attrib(species, "fbc:chemicalFormula", met.formula)

    # add in genes
    if len(cobra_model.genes) > 0:
        genes_list = SubElement(xml_model, GENELIST_TAG)
        for gene in cobra_model.genes:
            gene_id = gene.id.replace(".", SBML_DOT)
            sbml_gene = SubElement(genes_list, GENE_TAG)
            set_attrib(sbml_gene, "fbc:id", "G_" + gene_id)
            name = gene.name
            if name is None or len(name) == 0:
                name = gene.id
            set_attrib(sbml_gene, "fbc:label", gene_id)
            set_attrib(sbml_gene, "fbc:name", name)
            annotate_sbml_from_cobra(sbml_gene, gene)

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
        annotate_sbml_from_cobra(sbml_reaction, reaction)
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
            for met_id, stoichiomety in sorted(iteritems(reactants)):
                SubElement(reactant_list, "speciesReference", species=met_id,
                           stoichiometry=stoichiomety, constant="true")
        if len(products) > 0:
            product_list = SubElement(sbml_reaction, "listOfProducts")
            for met_id, stoichiomety in sorted(iteritems(products)):
                SubElement(product_list, "speciesReference", species=met_id,
                           stoichiometry=stoichiomety, constant="true")

        # gene reaction rule
        gpr = reaction.gene_reaction_rule
        if gpr is not None and len(gpr) > 0:
            gpr = gpr.replace(".", SBML_DOT)
            gpr_xml = SubElement(sbml_reaction, GPR_TAG)
            try:
                parsed, _ = parse_gpr(gpr)
                construct_gpr_xml(gpr_xml, parsed.body)
            except Exception as e:
                print("failed on '%s' in %s" %
                      (reaction.gene_reaction_rule, repr(reaction)))
                raise e

    return xml


def read_sbml_model(filename, number=float, **kwargs):
    if not _with_lxml:
        warn("Install lxml for faster SBML I/O", ImportWarning)
    xmlfile = parse_stream(filename)
    xml = xmlfile.getroot()
    # use libsbml if not l3v1 with fbc v2
    use_libsbml = (xml.get("level") != "3" or xml.get("version") != "1" or
                   get_attrib(xml, "fbc:required") is None)
    if use_libsbml:
        if libsbml is None:
            raise ImportError("libSBML required for fbc < 2")
        # libsbml needs a file string, so write to temp file if a file handle
        if hasattr(filename, "read"):
            with NamedTemporaryFile(suffix=".xml", delete=False) as outfile:
                xmlfile.write(outfile, encoding="UTF-8", xml_declaration=True)
            filename = outfile.name
        return read_sbml2(filename, **kwargs)
    try:
        return parse_xml_into_model(xml, number=number)
    except Exception:
        raise CobraSBMLError(
            "Something went wrong reading the model. You can get a detailed "
            "report using the `cobra.io.sbml3.validate_sbml_model` function "
            "or using the online validator at http://sbml.org/validator")


id_required = {ns(i) for i in ("sbml:model", "sbml:reaction:", "sbml:species",
                               "fbc:geneProduct", "sbml:compartment",
                               "sbml:parameter", "sbml:UnitDefinition",
                               "fbc:objective")}
invalid_id_detector = re.compile("|".join(re.escape(i[0]) for i in _renames))


def validate_sbml_model(filename, check_model=True):
    """Returns the model along with a list of errors.

    Parameters
    ----------
    filename : str
        The filename of the SBML model to be validated.
    check_model: bool, optional
        Whether to also check some basic model properties such as reaction
        boundaries and compartment formulas.

    Returns
    -------
    model : :class:`~cobra.core.Model.Model` object
        The cobra model if the file could be read succesfully or None
        otherwise.
    errors : dict
        Warnings and errors grouped by their respective types.

    Raises
    ------
    CobraSBMLError
        If the file is not a valid SBML Level 3 file with FBC.
    """
    xmlfile = parse_stream(filename)
    xml = xmlfile.getroot()
    # use libsbml if not l3v1 with fbc v2
    use_libsbml = (xml.get("level") != "3" or xml.get("version") != "1" or
                   get_attrib(xml, "fbc:required") is None)
    if use_libsbml:
        raise CobraSBMLError("XML is not SBML level 3 v1 with fbc v2")

    errors = {k: [] for k in ("validator", "warnings", "SBML errors", "other")}

    def err(err_msg, group="validator"):
        errors[group].append(err_msg)

    # make sure there is exactly one model
    xml_models = xml.findall(ns("sbml:model"))
    if len(xml_models) == 0:
        raise CobraSBMLError("No SBML model detected in file")
    if len(xml_models) > 1:
        err("More than 1 SBML model detected in file")
    xml_model = xml_models[0]

    # make sure all sbml id's are valid
    all_ids = set()
    for element in xmlfile.iter():
        if element.tag.startswith(fbc_prefix):
            prefix = fbc_prefix
        elif element.tag.startswith(sbml_prefix):
            prefix = ""
        else:
            continue
        str_id = element.get(prefix + "id", None)
        element_name = element.tag.split("}")[-1]
        id_repr = "%s id '%s' " % (element_name, str_id)
        if str_id is None or len(str_id) == 0:
            if element.tag in id_required:
                err(element_name + " missing id")
        else:
            if str_id in all_ids:
                err("duplicate id for " + id_repr)
            all_ids.add(str_id)
            try:
                str_id.encode("ascii")
            except UnicodeEncodeError as e:
                err("invalid character '%s' found in %s" %
                    (str_id[e.start:e.end], id_repr))
            if invalid_id_detector.search(str_id):
                bad_chars = "".join(invalid_id_detector.findall(str_id))
                err("invalid character%s '%s' found in %s" %
                    ("s" if len(bad_chars) > 1 else "", bad_chars, id_repr))
            if not str_id[0].isalpha():
                err("%s does not start with alphabet character" % id_repr)

    # check SBO terms
    for element in xml.findall(".//*[@sboTerm]"):
        sbo_term = element.get("sboTerm")
        if not sbo_term.startswith("SBO:"):
            err("sboTerm '%s' does not begin with 'SBO:'" % sbo_term)

    # ensure can be made into model
    # all warnings generated while loading will be logged as errors
    with catch_warnings(record=True) as warning_list:
        simplefilter("always")
        try:
            model = parse_xml_into_model(xml)
        except CobraSBMLError as e:
            err(str(e), "SBML errors")
            return None, errors
        except Exception as e:
            err(str(e), "other")
            return None, errors
    errors["warnings"].extend(str(i.message) for i in warning_list)

    # check genes
    xml_genes = {
        get_attrib(i, "fbc:id").replace(SBML_DOT, ".")
        for i in xml_model.iterfind(GENES_XPATH)}
    for gene in model.genes:
        if "G_" + gene.id not in xml_genes and gene.id not in xml_genes:
            err("No gene specfied with id 'G_%s'" % gene.id)

    if check_model:
        errors["validator"].extend(check_metabolite_compartment_formula(model))

    return model, errors


def write_sbml_model(cobra_model, filename, use_fbc_package=True, **kwargs):
    if not _with_lxml:
        warn("Install lxml for faster SBML I/O", ImportWarning)
    if not use_fbc_package:
        if libsbml is None:
            raise ImportError("libSBML required to write non-fbc models")
        write_sbml2(cobra_model, filename, use_fbc_package=False, **kwargs)
        return
    # create xml
    xml = model_to_xml(cobra_model, **kwargs)
    write_args = {"encoding": "UTF-8", "xml_declaration": True}
    if _with_lxml:
        write_args["pretty_print"] = True
        write_args["pretty_print"] = True
    else:
        indent_xml(xml)
    # write xml to file
    should_close = True
    if hasattr(filename, "write"):
        xmlfile = filename
        should_close = False
    elif filename.endswith(".gz"):
        xmlfile = GzipFile(filename, "wb")
    elif filename.endswith(".bz2"):
        xmlfile = BZ2File(filename, "wb")
    else:
        xmlfile = open(filename, "wb")
    ElementTree(xml).write(xmlfile, **write_args)
    if should_close:
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
