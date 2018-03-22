"""
SBML import and export using libsbml.


TODO: converters
- COBRA to FBCV2
- FBCV2 to COBRA
- FBCV1 to FBCV2

- SBMLIdConverter

"""
from __future__ import absolute_import

import tempfile
from warnings import warn
from six import string_types
from collections import defaultdict

import libsbml
from cobra.core import Gene, Metabolite, Model, Reaction
from cobra.util.solver import set_objective

LONG_SHORT_DIRECTION = {'maximize': 'max', 'minimize': 'min'}
SHORT_LONG_DIRECTION = {'min': 'minimize', 'max': 'maximize'}

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


# clip ids
# TODO: clip(met, "M_")
# clip_prefixes = {'compartment': None, 'specie': 'M_', 'gene': 'G_'}
# replace ids
# get_attrib(sbml_gene, "fbc:id").replace(SBML_DOT, ".")


def _parse_sbml_into_model(doc, number=float):
    """

    :param doc: libsbml.SBMLDocument
    'param number: data type of stoichiometry
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
        sid = _check_required(s, s.id, "id")
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

    # Reactions
    reactions = []
    for r in model.reactions:  # type: libsbml.Reaction
        rid = _check_required(r, r.id, "id")
        reaction = Reaction(rid)
        reaction.name = r.name
        annotate_cobra_from_sbase(reaction, r)

        # set bounds
        r_fbc = r.getPlugin("fbc")  # type: libsbml.FbcReactionPlugin
        if r_fbc is None:
            raise CobraSBMLError("No flux bounds on reaction '%s'" % r)
        else:
            # FIXME: remove code duplication in this section
            lb_id = _check_required(r_fbc, r_fbc.getLowerFluxBound(), "lowerFluxBound")
            ub_id = _check_required(r_fbc, r_fbc.getUpperFluxBound(), "upperFluxBound")
            p_lb = model.getParameter(lb_id)
            p_ub = model.getParameter(ub_id)

            if p_lb.constant and (p_lb.value is not None):
                reaction.lower_bound = p_lb.value
            else:
                raise CobraSBMLError("No constant bound '%s' for reaction '%s" % (p_lb, r))

            if p_ub.constant and (p_ub.value is not None):
                reaction.upper_bound = p_ub.value
            else:
                raise CobraSBMLError("No constant bound '%s' for reaction '%s" % (p_ub, r))
                bounds.append(p.value)

        reactions.append(reaction)

        # parse equation
        stoichiometry = defaultdict(lambda: 0)
        for sref in r.getListOfReactants():  # type: libsbml.SpeciesReference
            sid = sref.getSpecies()
            # FIXME: clip
            stoichiometry[sid] -= number(_check_required(sref, sref.stoichiometry, "stoichiometry"))

        for sref in r.getListOfProducts():  # type: libsbml.SpeciesReference
            sid = sref.getSpecies()
            stoichiometry[sid] += number(_check_required(sref, sref.stoichiometry, "stoichiometry"))

        # needs to have keys of metabolite objects, not ids
        object_stoichiometry = {}
        for met_id in stoichiometry:
            if met_id in boundary_ids:
                warn("Boundary metabolite '%s' used in reaction '%s'" %
                     (met_id, reaction.id))
                continue
            try:
                metabolite = cmodel.metabolites.get_by_id(met_id)
            except KeyError:
                warn("ignoring unknown metabolite '%s' in reaction %s" %
                     (met_id, reaction.id))
                continue
            object_stoichiometry[metabolite] = stoichiometry[met_id]
        reaction.add_metabolites(object_stoichiometry)


        # GPR rules
        # TODO
        '''
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
        

        def process_association(association):
            """ Recursively convert gpr xml to a gpr string. """
            type_code = association.getTypeCode()
            if association.isFbcOr():
                association.get

                return "( " + ' or '.join(process_gpa(i) for i in gpa.getCh) + " )"
            elif sub_xml.tag == AND_TAG:
                return "( " + ' and '.join(process_gpr(i) for i in sub_xml) + " )"
            elif sub_xml.tag == GENEREF_TAG:
                gene_id = get_attrib(sub_xml, "fbc:geneProduct", require=True)
                return clip(gene_id, "G_")
            else:
                raise Exception("unsupported tag " + sub_xml.tag)
        '''
        gpa = r_fbc.getGeneProductAssociation()  # type: libsbml.GeneProductAssociation
        print(gpa)

        association = None
        if gpa is not None:
            association = gpa.getAssociation()  # type: libsbml.FbcAssociation
            print(association)
            print(association.getListOfAllElements())
            print(gpa.)
            print(association.getListOfFbcAssociations())



        # gpr = process_association(association) if association is not None else ''

        # remove outside parenthesis, if any
        if gpr.startswith("(") and gpr.endswith(")"):
            gpr = gpr[1:-1].strip()

        # gpr = gpr.replace(SBML_DOT, ".")
        reaction.gene_reaction_rule = gpr

        try:
            cmodel.add_reactions(reactions)
        except ValueError as e:
            warn(str(e))

    # Objective
    obj_list = model_fbc.getListOfObjectives()  # type: libsbml.ListOfObjectives
    if obj_list is None:
        warn("listOfObjectives element not found")
    else:
        obj_id = obj_list.getActiveObjective()
        obj = model_fbc.getObjective(obj_id)  # type: libsbml.Objective
        obj_direction = LONG_SHORT_DIRECTION[obj.getType()]

        coefficients = {}

        for flux_obj in obj.getListOfFluxObjectives():  # type: libsbml.FluxObjective
            # FIXME: clip id
            rid = flux_obj.getReaction()
            try:
                objective_reaction = cmodel.reactions.get_by_id(rid)
            except KeyError:
                raise CobraSBMLError("Objective reaction '%s' not found" % rid)
            try:
                coefficients[objective_reaction] = number(flux_obj.getCoefficient())
            except ValueError as e:
                warn(str(e))
        set_objective(cmodel, coefficients)
        cmodel.solver.objective.direction = obj_direction

    return cmodel


def _check_required(sbase, value, attribute):
    """ Get required attribute from the SBase.

    :param sbase:
    :param attribute:
    :return:
    """
    if value is None:
        msg = "required attribute '%s' not found in '%s'" % \
              (attribute, sbase)
        if sbase.id is not None:
            msg += " with id '%s'" % sbase.id
        elif sbase.name is not None:
            msg += " with name '%s'" % sbase.get("name")
        raise CobraSBMLError(msg)
    return value


def annotate_cobra_from_sbase(cobj, sbase):
    """ Read annotations from SBase into dictionary.

    :param cobj:
    :param sbase:
    :return:
    """
    annotation = cobj.annotation

    # SBO term
    if sbase.isSetSBOTerm():
        annotation["SBO"] = sbase.getSBOTerm()

    # RDF annotation

    cvterms = sbase.getCVTerms()
    if cvterms is None:
        return

    for cvterm in cvterms:  # type: libsbml.CVTerm
        # FIXME: currently only the terms, but not the qualifier
        # are stored (only subset of identifiers.org parsed)
        for k in range(cvterm.getNumResources()):
            uri = cvterm.getResourceURI(k)
            if not uri.startswith("http://identifiers.org/"):
                warn("%s does not start with http://identifiers.org/" % uri)
                continue
            try:
                provider, identifier = uri[23:].split("/", 1)
            except ValueError:
                warn("%s does not conform to http://identifiers.org/provider/id"
                     % uri)
                continue

            # handle multiple by same provider (create list)
            if provider in annotation:
                if isinstance(annotation[provider], string_types):
                    annotation[provider] = [annotation[provider]]
                annotation[provider].append(identifier)
            else:
                annotation[provider] = identifier


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