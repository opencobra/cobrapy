"""
SBML import and export using libsbml.


TODO: converters
- COBRA to FBCV2
- FBCV2 to COBRA
- FBCV1 to FBCV2

- SBMLIdConverter

"""
# -------------------------------
# TODO
# ------------------------------
# [1] Replacing/Changing of identifiers between SBML and cobra formats
# clip ids
# clip(met, "M_")
# clip_prefixes = {'compartment': None, 'specie': 'M_', 'gene': 'G_'}
# replace ids
# get_attrib(sbml_gene, "fbc:id").replace(SBML_DOT, ".")


# [2] Legacy format and COBRA format support
# [3] Conversion of FBCv1 to FBCv2


from __future__ import absolute_import

import os
import re
from warnings import catch_warnings, simplefilter, warn
from six import string_types, iteritems
from collections import defaultdict, namedtuple

import libsbml
from cobra.core import Gene, Metabolite, Model, Reaction
from cobra.util.solver import set_objective
from cobra.manipulation.validate import check_metabolite_compartment_formula

from .sbml import write_cobra_model_to_sbml_file as write_sbml2


LONG_SHORT_DIRECTION = {'maximize': 'max', 'minimize': 'min'}
SHORT_LONG_DIRECTION = {'min': 'minimize', 'max': 'maximize'}

# ----------------------------------------------------------
# Defaults for writing SBML
# ----------------------------------------------------------
LOWER_BOUND = -1000
UPPER_BOUND = 1000
LOWER_BOUND_ID = "cobra_default_lb"
UPPER_BOUND_ID = "cobra_default_ub"
ZERO_BOUND_ID = "cobra_0_bound"

SBO_FBA_FRAMEWORK = "SBO:0000624"
SBO_DEFAULT_FLUX_BOUND = "SBO:0000626"
SBO_FLUX_BOUND = "SBO:0000625"

Unit = namedtuple('Unit', ['kind', 'scale', 'multiplier', 'exponent'])
UNITS_FLUX = ("mmol_per_gDW_per_hr",
              [Unit(kind=libsbml.UNIT_KIND_MOLE, scale=-3, multiplier=1, exponent=1),
               Unit(kind=libsbml.UNIT_KIND_GRAM, scale=0, multiplier=1, exponent=-1),
               Unit(kind=libsbml.UNIT_KIND_SECOND, scale=0, multiplier=3600, exponent=-1)]
              )
# ----------------------------------------------------------


class CobraSBMLError(Exception):
    """ SBML error class. """
    pass


def read_sbml_model(filename):
    """ Reads model from given filename.

    If the given filename ends with the suffix ''.gz'' (for example,
    ''myfile.xml.gz'),' the file is assumed to be compressed in gzip
    format and will be automatically decompressed upon reading. Similarly,
    if the given filename ends with ''.zip'' or ''.bz2',' the file is
    assumed to be compressed in zip or bzip2 format (respectively).  Files
    whose names lack these suffixes will be read uncompressed.  Note that
    if the file is in zip format but the archive contains more than one
    file, only the first file in the archive will be read and the rest
    ignored.

    To read a gzip/zip file, libSBML needs to be configured and linked
    with the zlib library at compile time.  It also needs to be linked
    with the bzip2 library to read files in bzip2 format.  (Both of these
    are the default configurations for libSBML.)

    :param filename: path to SBML file or SBML string
    :param validate: validate the file on reading (additional overhead)
    :return:
    """
    try:
        doc = _get_doc_from_filename(filename)
        return _sbml_to_model(doc)
    except Exception:
        raise CobraSBMLError(
            "Something went wrong reading the model. You can get a detailed "
            "report using the `cobra.io.sbml3.validate_sbml_model` function "
            "or using the online validator at http://sbml.org/validator")


def _get_doc_from_filename(filename):
    """ SBMLDocument from given filename.

    :param filename:
    :return:
    """
    if os.path.exists(filename):
        doc = libsbml.readSBMLFromFile(filename)  # type: libsbml.SBMLDocument
    elif isinstance(filename, string_types):
        # SBML as string representation
        doc = libsbml.readSBMLFromString(filename)
    elif hasattr(filename, "read"):
        # File handle
        doc = libsbml.readSBMLFromString(filename.read())
    else:
        raise CobraSBMLError("Input format is not supported.")
    return doc


def write_sbml_model(cobra_model, filename, use_fbc_package, **kwargs):
    """ Writes cobra model to filename.

    The created model is SBML level 3 version 1 (L1V3) with
    fbc package v2 (fbc-v2).

    If the given filename ends with the suffix ".gz" (for example,
    "myfile.xml.gz"), libSBML assumes the caller wants the file to be
    written compressed in gzip format. Similarly, if the given filename
    ends with ".zip" or ".bz2", libSBML assumes the caller wants the
    file to be compressed in zip or bzip2 format (respectively). Files
    whose names lack these suffixes will be written uncompressed. Special
    considerations for the zip format: If the given filename ends with
    ".zip", the file placed in the zip archive will have the suffix
    ".xml" or ".sbml".  For example, the file in the zip archive will
    be named "test.xml" if the given filename is "test.xml.zip" or
    "test.zip". Similarly, the filename in the archive will be
    "test.sbml" if the given filename is "test.sbml.zip".

    :param cobra_model:
    :param filename:
    :param use_fbc_package:
    :param kwargs:
    :return:
    """
    if not use_fbc_package:
        # legacy cobra without fbc
        write_sbml2(cobra_model, filename, use_fbc_package=False, **kwargs)

    # create xml
    doc = _model_to_sbml(cobra_model, **kwargs)
    libsbml.writeSBMLToFile(doc, filename)


def _sbml_to_model(doc, number=float):
    """ Creates cobra model from SBMLDocument.

    :param doc: libsbml.SBMLDocument
    'param number: data type of stoichiometry
    :return: cobrapy model
    """
    # SBML model
    model = doc.getModel()  # type: libsbml.Model
    if model is None:
        raise CobraSBMLError("No SBML model detected in file.")
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
        sid = _check_required(s, s.id, "id")  # FIXME: S_ prefix (clip)
        met = Metabolite(sid)
        met.name = s.name
        met.compartment = s.compartment

        # parse notes in notes dictionary
        notes = s.getNotesString()
        if notes and len(notes) > 0:
            met.notes = _parse_notes(notes)

        s_fbc = s.getPlugin("fbc")
        if s_fbc:
            met.charge = s_fbc.getCharge()
            met.formula = s_fbc.getChemicalFormula()
        else:
            if 'CHARGE' in met.notes:
                met.charge = met.notes['CHARGE']
            if 'FORMULA' in met.notes:
                met.formula = met.notes['FORMULA']

        # Detect boundary metabolites - In case they have been mistakenly
        # added. They should not actually appear in a model
        if s.getBoundaryCondition() is True:
            boundary_ids.add(s.id)

        annotate_cobra_from_sbase(met, s)

        cmodel.add_metabolites([met])

    # Genes
    if model_fbc:
        for gp in model_fbc.getListOfGeneProducts():  # type: libsbml.GeneProduct
            gid = gp.id   # FIXME: G_ prefix (clip), DOT replacements
            gene = Gene(gid)
            gene.name = gp.name
            if gene.name is None:
                gene.name = gp.get
            annotate_cobra_from_sbase(gene, gp)
            cmodel.genes.append(gene)

    # Reactions
    reactions = []
    for r in model.reactions:  # type: libsbml.Reaction
        rid = _check_required(r, r.id, "id")  # FIXME: R_ prefix (clip)
        reaction = Reaction(rid)
        reaction.name = r.name
        annotate_cobra_from_sbase(reaction, r)

        # parse notes in notes dictionary
        notes = r.getNotesString()
        if notes and len(notes) > 0:
            reaction.notes = _parse_notes(notes)

        # set bounds
        r_fbc = r.getPlugin("fbc")  # type: libsbml.FbcReactionPlugin
        if r_fbc:
            # information in fbc
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

        elif r.isSetKineticLaw():
            # sometime information encoded in kinetic laws
            klaw = r.getKineticLaw()  # type: libsbml.KineticLaw
            p_lb = klaw.getParameter("LOWER_BOUND")
            if p_lb:
                reaction.lower_bound = p_lb.value
            else:
                raise CobraSBMLError("Missing flux bounds on reaction '%s'" % r)
            p_ub = klaw.getParameter("UPPER_BOUND")
            if p_ub:
                reaction.upper_bound = p_ub.value
            else:
                raise CobraSBMLError("Missing flux bounds on reaction '%s'" % r)
            warn("Bounds encoded in KineticLaw for '%s" % r)
        else:
            raise CobraSBMLError("No flux bounds on reaction '%s'" % r)


        reactions.append(reaction)

        # parse equation
        stoichiometry = defaultdict(lambda: 0)
        for sref in r.getListOfReactants():  # type: libsbml.SpeciesReference
            sid = sref.getSpecies()    # FIXME: M_ prefix (clip)
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
        if r_fbc:
            gpa = r_fbc.getGeneProductAssociation()  # type: libsbml.GeneProductAssociation

            # print(gpa)

            association = None
            if gpa is not None:
                association = gpa.getAssociation()  # type: libsbml.FbcAssociation
                # print(association)
                # print(association.getListOfAllElements())

            # gpr = process_association(association) if association is not None else ''
            gpr = ''
        else:
            # fallback to notes information
            gpr = reaction.notes.get('GENE_ASSOCIATION', '')

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
    if model_fbc:
        obj_list = model_fbc.getListOfObjectives()  # type: libsbml.ListOfObjectives
        if obj_list is None:
            warn("listOfObjectives element not found")
        else:
            obj_id = obj_list.getActiveObjective()
            obj = model_fbc.getObjective(obj_id)  # type: libsbml.Objective
            obj_direction = LONG_SHORT_DIRECTION[obj.getType()]

            coefficients = {}

            for flux_obj in obj.getListOfFluxObjectives():  # type: libsbml.FluxObjective
                rid = flux_obj.getReaction()  # FIXME: R_ prefix (clip)
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


def _model_to_sbml(cobra_model, units=True):
    """

    :param cobra_model:
    :param units: boolean, if true the FLUX_UNITS are written
    :return:
    """
    sbmlns = libsbml.SBMLNamespaces(3, 1)
    sbmlns.addPackageNamespace("fbc", 2)

    doc = libsbml.SBMLDocument(sbmlns)  # type: libsbml.SBMLDocument
    doc.setPackageRequired("fbc", False)
    doc.setSBOTerm(SBO_FBA_FRAMEWORK)
    model = doc.createModel()  # type: libsbml.Model
    model_fbc = model.getPlugin("fbc")  # type: libsbml.FbcModelPlugin
    model_fbc.setStrict(True)

    if cobra_model.id is not None:
        model.setId(cobra_model.id)
    if cobra_model.name is not None:
        model.setName(cobra_model.name)

    # Units
    if units:
        flux_udef = model.createUnitDefinition()  # type: libsbml.UnitDefinition
        flux_udef.setId(UNITS_FLUX[0])
        for u in UNITS_FLUX[1]:
            unit = flux_udef.createUnit()  # type: libsbml.Unit
            unit.setKind(u.kind)
            unit.setExponent(u.exponent)
            unit.setScale(u.scale)
            unit.setMultiplier(u.multiplier)

    # Flux bounds
    def _create_bound(model, reaction, bound_type):
        """ Creates bound parameter.

        :param model: SBML model
        :param reaction: cobra reaction
        :param bound_type:
        :return: parameter id of bound
        """
        value = getattr(reaction, bound_type)
        if value == LOWER_BOUND:
            return LOWER_BOUND_ID
        elif value == 0:
            return ZERO_BOUND_ID
        elif value == UPPER_BOUND:
            return UPPER_BOUND_ID
        else:
            # new parameter
            pid = reaction.id + "_" + bound_type  # FIXME: R_ prefix
            _create_parameter(model, pid=pid, value=value, sbo=SBO_FLUX_BOUND)
            return pid

    def _create_parameter(model, pid, value, sbo=None, constant=True):
        """ Create parameter in SBML model. """
        p = model.createParameter()  # type: libsbml.Parameter
        p.setId(pid)
        p.setValue(value)
        p.setConstant(constant)
        if sbo:
            p.setSBOTerm(sbo)
        if units:
            p.setUnits(flux_udef.getId())

    # minimum and maximum from model
    if len(cobra_model.reactions) > 0:
        min_value = min(cobra_model.reactions.list_attr("lower_bound"))
        max_value = max(cobra_model.reactions.list_attr("upper_bound"))
    else:
        min_value = LOWER_BOUND
        max_value = UPPER_BOUND

    _create_parameter(model, pid=LOWER_BOUND_ID, value=min_value, sbo=SBO_DEFAULT_FLUX_BOUND)
    _create_parameter(model, pid=UPPER_BOUND_ID, value=max_value, sbo=SBO_DEFAULT_FLUX_BOUND)
    _create_parameter(model, pid=ZERO_BOUND_ID, value=0, sbo=SBO_DEFAULT_FLUX_BOUND)

    # Compartments
    for cid, name in iteritems(cobra_model.compartments):
        c = model.createCompartment()  # type: libsbml.Compartment
        c.setId(cid)
        c.setName(name)
        c.setConstant(True)

    # Species
    for met in cobra_model.metabolites:
        s = model.createSpecies()  # type: libsbml.Species
        mid = met.id
        s.setId(mid)  # FIXME: id replacement (R_prefix)
        s.setConstant(True)
        s.setBoundaryCondition(True)
        s.setHasOnlySubstanceUnits(False)
        s.setName(met.name)
        s.setCompartment(met.compartment)
        s_fbc = s.getPlugin("fbc")  # type: libsbml.FbcSpeciesPlugin
        if met.charge is not None:
            s_fbc.setCharge(met.charge)
        if met.formula is not None:
            s_fbc.setChemicalFormula(met.formula)

        annotate_sbase_from_cobra(s, met)

    # Genes
    for gene in cobra_model.genes:
        gp = model_fbc.createGeneProduct()  # type: libsbml.GeneProduct
        gid = gene.id  # FIXME: id replacement (SBML_DOT, G_ prefix)
        gp.setId(gid)
        gname = gene.name
        if gname is None or len(gname) == 0:
            gname = gid
        gp.setName(gname)
        gp.setLabel(gid)

        annotate_sbase_from_cobra(gp, gene)

    # Objective
    objective = model_fbc.createObjective()  # type: libsbml.Objective
    objective.setId("obj")
    objective.setType(SHORT_LONG_DIRECTION[cobra_model.objective.direction])
    model_fbc.setActiveObjectiveId("obj")

    # Reactions
    for reaction in cobra_model.reactions:
        rid = reaction.id  # FIXME: id replacement (R_prefix)
        r = model.createReaction()  # type: libsbml.Reaction
        r.setId(rid)
        r.setName(reaction.name)
        r.setFast(False)
        r.setReversible((reaction.lower_bound < 0))

        annotate_sbase_from_cobra(r, reaction)

        # stoichiometry
        for metabolite, stoichiomety in iteritems(reaction._metabolites):
            sid = metabolite.id  # FIXME: id replacement (M_ prefix)
            if stoichiomety < 0:
                sref = r.createReactant()  # type: libsbml.SpeciesReference
                sref.setSpecies(sid)
                sref.setStoichiometry(-stoichiomety)
                sref.setConstant(True)
            else:
                sref = r.createProduct()  # type: libsbml.SpeciesReference
                sref.setSpecies(sid)
                sref.setStoichiometry(stoichiomety)
                sref.setConstant(True)

        # bounds
        r_fbc = r.getPlugin("fbc")  # type: libsbml.FbcReactionPlugin

        r_fbc.setLowerFluxBound(_create_bound(model, reaction, "lower_bound"))
        r_fbc.setUpperFluxBound(_create_bound(model, reaction, "upper_bound"))

        # GPR
        gpr = reaction.gene_reaction_rule
        if gpr is not None and len(gpr) > 0:
            gpa = r_fbc.createGeneProductAssociation()  # type: libsbml.GeneProductAssociation
            # This is a helper method that allows a user to set the
            # GeneProductAssociation via a string such as "a1 AND b1 OR C2" and
            # have the method work out the correct XML structure.
            gpa.setAssociation(gpr)

        # objective coefficient
        if reaction.objective_coefficient != 0:
            flux_obj = objective.createFluxObjective()  # type: libsbml.FluxObjective
            flux_obj.setReaction(rid)
            flux_obj.setCoefficient(reaction.objective_coefficient)

    return doc


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


def _parse_notes(notes):
    """ Creates dictionary of notes.

    :param notes:
    :return:
    """
    pattern = r"<p>\s*(\w+)\s*:\s*([\w|\s]+)<"
    matches = re.findall(pattern, notes)
    return {k.strip(): v.strip() for (k, v) in matches}

# ----------------------
# Annotations
# ----------------------
# FIXME: currently only the terms, but not the qualifier are parsed
# FIXME: migration to https, both should be supported
# (better parsing of collection & id via regular expression)

URL_IDENTIFIERS = "http://identifiers.org/"

QualifierType = {
    0: "MODEL_QUALIFIER",
    1: "BIOLOGICAL_QUALIFIER",
    2: "UNKNOWN_QUALIFIER"
}

ModelQualifierType = {
    0: "BQM_IS",
    1: "BQM_IS_DESCRIBED_BY",
    2: "BQM_IS_DERIVED_FROM",
    3: "BQM_IS_INSTANCE_OF",
    4: "BQM_HAS_INSTANCE",
    5: "BQM_UNKNOWN",
}

BiologicalQualifierType = {
    0: "BQB_IS",
    1: "BQB_HAS_PART",
    2: "BQB_IS_PART_OF",
    3: "BQB_IS_VERSION_OF",
    4: "BQB_HAS_VERSION",
    5: "BQB_IS_HOMOLOG_TO",
    6: "BQB_IS_DESCRIBED_BY",
    7: "BQB_IS_ENCODED_BY",
    8: "BQB_ENCODES",
    9: "BQB_OCCURS_IN",
    10: "BQB_HAS_PROPERTY",
    11: "BQB_IS_PROPERTY_OF",
    12: "BQB_HAS_TAXON",
    13: "BQB_UNKNOWN",
}


def annotate_cobra_from_sbase(cobj, sbase):
    """ Read annotations from SBase into dictionary.

    :param cobj:
    :param sbase:
    :return:
    """
    annotation = cobj.annotation

    # SBO term
    if sbase.isSetSBOTerm():
        annotation["sbo"] = sbase.getSBOTerm()

    # RDF annotation
    cvterms = sbase.getCVTerms()
    if cvterms is None:
        return

    for cvterm in cvterms:  # type: libsbml.CVTerm
        for k in range(cvterm.getNumResources()):
            uri = cvterm.getResourceURI(k)
            if not uri.startswith(URL_IDENTIFIERS):
                warn("%s does not start with %s" % (uri, URL_IDENTIFIERS))
                continue
            try:
                provider, identifier = uri[23:].split("/", 1)
            except ValueError:
                warn("%s does not conform to http://identifiers.org/collection/id"
                     % uri)
                continue

            # handle multiple by same provider (create list)
            if provider in annotation:
                if isinstance(annotation[provider], string_types):
                    annotation[provider] = [annotation[provider]]
                annotation[provider].append(identifier)
            else:
                annotation[provider] = identifier


def annotate_sbase_from_cobra(sbase, cobj):
    """ Set cobra annotations on SBase into.

    :param sbase:
    :param cobj: cobra object

    :return:
    """
    if len(cobj.annotation) == 0:
        return

    sbase.setMetaId("meta_{}".format(sbase.id))
    for provider, identifiers in sorted(iteritems(cobj.annotation)):
        if provider == "SBO":
            sbase.setSBOTerm(identifiers)
        else:
            if isinstance(identifiers, string_types):
                identifiers = (identifiers,)

            for identifier in identifiers:
                _add_cv_to_sbase(sbase, qualifier="BQB_IS",
                                 resource="%s/%s/%s" % (URL_IDENTIFIERS, provider, identifier))


def _add_cv_to_sbase(sbase, qualifier, resource):
    """ Adds RDF information to given element.

    :param sbase:
    :param qualifier:
    :param resource:
    :return:
    """
    cv = libsbml.CVTerm()  # type: libsbml.CVTerm

    if qualifier.startswith('BQB'):
        cv.setQualifierType(libsbml.BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(libsbml.__dict__.get(qualifier))
    elif qualifier.startswith('BQM'):
        cv.setQualifierType(libsbml.MODEL_QUALIFIER)
        cv.setModelQualifierType(libsbml.__dict__.get(qualifier))
    else:
        raise CobraSBMLError('Unsupported qualifier: {}'.format(qualifier))

    cv.addResource(resource)
    sbase.addCVTerm(cv)


# -----------------------------------
# Validation
# -----------------------------------
def validate_sbml_model(filename, use_libsbml=False, check_model=True, ucheck=False, internalConsistency=True,
                        check_units_consistency=False,
                        check_modeling_practice=False):
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
    # store errors
    errors = {key: [] for key in ("validator", "warnings", "other", "SBML errors")}
    if use_libsbml:
        for key in ["SBML_FATAL", "SBML ERROR", "SBML_SCHEMA_ERROR", "SBML_WARNING"]:
            errors[key] = []

    def err(err_msg, group="validator"):
        errors[group].append(err_msg)

    # make sure there is exactly one model
    doc = _get_doc_from_filename(filename)
    model = doc.getModel()  # type: libsbml.Model
    if model is None:
        raise CobraSBMLError("No SBML model detected in file.")

    if use_libsbml:
        # set the unit checking, similar for the other settings
        doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, check_units_consistency)
        doc.setConsistencyChecks(libsbml.LIBSBML_CAT_MODELING_PRACTICE, check_modeling_practice)

        # validate the document
        if internalConsistency:
            doc.checkInternalConsistency()
        doc.checkConsistency()

        for k in range(doc.getNumErrors()):
            e = doc.getError(k)
            sev = e.getSeverity()
            if sev == libsbml.LIBSBML_SEV_FATAL:
                err(error_string(e), "SBML_FATAL")
            elif sev == libsbml.LIBSBML_SEV_ERROR:
                err(error_string(e), "SBML_ERROR")
            elif sev == libsbml.LIBSBML_SEV_SCHEMA_ERROR:
                err(error_string(e), "SBML_SCHEMA_ERROR")
            elif sev == libsbml.LIBSBML_SEV_WARNING:
                err(error_string(e), "SBML_WARNING")

    # ensure can be made into model
    # all warnings generated while loading will be logged as errors
    with catch_warnings(record=True) as warning_list:
        simplefilter("always")
        try:
            model = _sbml_to_model(doc)
        except CobraSBMLError as e:
            err(str(e), "SBML errors")
            return None, errors
        except Exception as e:
            err(str(e), "other")
            return None, errors
    errors["warnings"].extend(str(i.message) for i in warning_list)

    if check_model:
        errors["validator"].extend(check_metabolite_compartment_formula(model))

    return model, errors


def error_string(error, k=None):
    """ String representation of SBMLError.

    :param error:
    :return:
    """
    package = error.getPackage()
    if package == '':
        package = 'core'

    error_str = 'E{}: {} ({}, L{}, {})  \n' \
                '{}\n' \
                '[{}] {}\n' \
                '{}\n'.format(
                    k, error.getCategoryAsString(), package, error.getLine(), 'code',
                    '-' * 60,
                    error.getSeverityAsString(), error.getShortMessage(),
                    error.getMessage())
    return error_str