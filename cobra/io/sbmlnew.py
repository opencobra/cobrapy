"""
SBML import and export using libsbml.
"""
from __future__ import absolute_import

import tempfile
from warnings import warn
from six import string_types

import libsbml
from cobra.core import Gene, Metabolite, Model, Reaction


'''
try:
    from optlang.symbolics import Basic
except ImportError:
    class Basic:
        pass
'''


class CobraSBMLError(Exception):
    pass


def read_sbml_model(filename):
    """ Reads model from given filename.

    :param filename: path to SBML file or SBML string
    :return:
    """

    # FIXME: support file handles
    # libsbml needs a file string, so write to temp file if a file handle
    '''
    if hasattr(filename, "read"):
        with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as outfile:
            xmlfile.write(outfile, encoding="UTF-8")
        filename = outfile.name
    '''

    try:

        doc = libsbml.readSBML(filename)  # type: libsbml.SBMLDocument

        return _parse_sbml_into_model(doc)

    except Exception:
        raise CobraSBMLError(
            "Something went wrong reading the model. You can get a detailed "
            "report using the `cobra.io.sbml3.validate_sbml_model` function "
            "or using the online validator at http://sbml.org/validator")

    return None


def _get_required_attrib(sbase, attribute):
    value = getattr(sbase, attribute, None)
    if value is None:
        msg = "required attribute '%s' not found in '%s'" % \
              (attribute, sbase)
        if sbase.id is not None:
            msg += " with id '%s'" % sbase.id
        elif sbase.name is not None:
            msg += " with name '%s'" % tag.get("name")
        raise CobraSBMLError(msg)
    return value

def annotate_cobra_from_sbase(cobj, sbase):

    # TODO: implement
    pass



# clip ids
# TODO: clip(met, "M_")
# clip_prefixes = {'compartment': None, 'specie': 'M_', 'gene': 'G_'}
# replace ids
# get_attrib(sbml_gene, "fbc:id").replace(SBML_DOT, ".")


def _parse_sbml_into_model(doc):
    """

    :param doc: libsbml.SBMLDocument
    :return: cobrapy model
    """
    # SBML model
    doc_fbc = doc.getPlugin("fbc")  # type: libsbml.FbcSBMLDocumentPlugin
    model = doc.getModel()  # type: libsbml.Model
    model_fbc = model.getPlugin("fbc")  # type: libsbml.FbcModelPlugin

    if not model_fbc:
        warn("Model does not contain FBC information.")
    else:
        if not model_fbc.isSetStrict():
            warn('Loading SBML model without fbc:strict="true"')

    # Model
    cmodel = Model(model.id)
    cmodel.name = model.name

    # Compartments
    cmodel.compartments = {c.id: c.name for c in model.compartments}

    # Species
    boundary_ids = set()
    for s in model.species:  # type: libsbml.Species
        sid = _get_required_attrib(s, "id")
        met = Metabolite(sid)
        met.name = s.name
        met.compartment = s.compartment
        s_fbc = s.getPlugin("fbc")
        if s_fbc:
            met.charge = s_fbc.getCharge()
            met.formula = s_fbc.getChemicalFormula()

        # Detect boundary metabolites - In case they have been mistakenly
        # added. They should not actually appear in a model
        if s.getBoundaryCondition() is True:
            boundary_ids.add(s.id)

        annotate_cobra_from_sbase(met, s)

        cmodel.add_metabolites([met])

    # Genes
    for gp in model_fbc.getListOfGeneProducts():  # type: libsbml.GeneProduct
        gene = Gene(gp.id)
        gene.name = gp.name
        if gene.name is None:
            gene.name = gp.get
        annotate_cobra_from_sbase(gene, gp)
        cmodel.genes.append(gene)


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


    # Reactions
    reactions = []
    for r in model.reactions:  # type: libsbml.Reaction
        rid = _get_required_attrib(r, "id")
        reaction = Reaction(rid)
        reaction.name = r.name
        annotate_cobra_from_sbml(reaction, sbml_reaction)

        # set bounds
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


def write_sbml_model(path, legacy=True):
    """ Reads model from given path.

    :param path:
    :param legacy: write legacy format
    :return:
    """
    return None


def validate_sbml_model(path):
    """ Validate given SBML model.

    :param path:
    :return:
    """
    assert 0 == 1