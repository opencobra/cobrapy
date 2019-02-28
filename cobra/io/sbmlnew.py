"""
SBML import and export using python-libsbml(-experimental).

- The SBML importer supports all versions of SBML and the fbc package.
- The SBML exporter writes SBML L3 models.
- Annotation information is stored on the cobrapy objects
- Information from the group package is read

"""

from __future__ import absolute_import

import os
import re
import traceback
import logging
from warnings import catch_warnings, simplefilter
from six import string_types, iteritems
from collections import defaultdict, namedtuple

import libsbml
from cobra.core import Gene, Metabolite, Model, Reaction
from cobra.util.solver import set_objective, linear_reaction_coefficients
from cobra.manipulation.validate import check_metabolite_compartment_formula

from .sbml import write_cobra_model_to_sbml_file as write_sbml2


class CobraSBMLError(Exception):
    """ SBML error class. """
    pass


LOGGER = logging.getLogger(__name__)
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
              [
                  Unit(kind=libsbml.UNIT_KIND_MOLE, scale=-3, multiplier=1,
                       exponent=1),
                  Unit(kind=libsbml.UNIT_KIND_GRAM, scale=0, multiplier=1,
                       exponent=-1),
                  Unit(kind=libsbml.UNIT_KIND_SECOND, scale=0, multiplier=3600,
                       exponent=-1)
              ])
# ----------------------------------------------------------
# Functions for id replacements (import/export)
# ----------------------------------------------------------
F_GENE = "F_GENE"
F_GENE_REV = "F_GENE_REV"
F_SPECIE = "F_SPECIE"
F_SPECIE_REV = "F_SPECIE_REV"
F_REACTION = "F_REACTION"
F_REACTION_REV = "F_REACTION_REV"

SBML_DOT = "__SBML_DOT__"


def _clip(s, prefix):
    """Clips a prefix from the beginning of a string if it exists."""
    return s[len(prefix):] if s.startswith(prefix) else s


def _f_gene(s, prefix="G_"):
    s = s.replace(SBML_DOT, ".")
    return _clip(s, prefix)


def _f_gene_rev(s, prefix="G_"):
    return prefix + s.replace(".", SBML_DOT)


def _f_specie(s, prefix="M_"):
    return _clip(s, prefix)


def _f_specie_rev(s, prefix="M_"):
    return prefix + s


def _f_reaction(s, prefix="R_"):
    return _clip(s, prefix)


def _f_reaction_rev(s, prefix="R_"):
    return prefix + s


F_REPLACE = {
    F_GENE: _f_gene,
    F_GENE_REV: _f_gene_rev,
    F_SPECIE: _f_specie,
    F_SPECIE_REV: _f_specie_rev,
    F_REACTION: _f_reaction,
    F_REACTION_REV: _f_reaction_rev,
}


# ----------------------------------------------------------


def read_sbml_model(filename, number=float, f_replace=F_REPLACE, **kwargs):
    """Reads SBML model from given filename.

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

    This function supports SBML with FBC-v1 and FBC-v2. FBC-v1 models
    are converted to FBC-v2 models before reading.

    The parser tries to fall back to information in notes dictionaries
    if information is not available in the FBC packages, e.g.,
    CHARGE, FORMULA on species, or GENE_ASSOCIATION, SUBSYSTEM on reactions.

    Parameters
    ----------
    filename : path to SBML file, or SBML string, or SBML file handle
        SBML which is read into cobra model
    number: data type of stoichiometry: {float, int}
        In which data type should the stoichiometry be parsed.
    f_replace : dict of replacement functions for id replacement
        Dictionary of replacement functions for gene, specie, and reaction.
        By default the following id changes are performed on import:
        clip G_ from genes, clip M_ from species, clip R_ from reactions
        If no replacements should be performed, set f_replace={}, None

    Returns
    -------
    cobra.core.Model

    Notes
    -----
    Provided file handles cannot be opened in binary mode, i.e., use
        with open(path, "r" as f):
            read_sbml_model(f)
    File handles to compressed files are not supported yet.
    """
    try:
        doc = _get_doc_from_filename(filename)
        return _sbml_to_model(doc, number=number,
                              f_replace=f_replace, **kwargs)
    except Exception:
        print(traceback.print_exc())
        raise CobraSBMLError(
            "Something went wrong reading the SBML model. You can get a "
            "detailed report using the `cobra.io.sbml3.validate_sbml_model` "
            "function or using the online validator at "
            "http://sbml.org/validator")


def _get_doc_from_filename(filename):
    """Get SBMLDocument from given filename.

    Parameters
    ----------
    filename : path to SBML, or SBML string, or filehandle

    Returns
    -------
    libsbml.SBMLDocument
    """
    if isinstance(filename, string_types):
        if os.path.exists(filename):
            # SBML as file
            doc = libsbml.readSBMLFromFile(filename)  # type: libsbml.SBMLDocument
        else:
            # SBML as string representation
            doc = libsbml.readSBMLFromString(filename)
    elif hasattr(filename, "read"):
        # File handle
        doc = libsbml.readSBMLFromString(filename.read())
    else:
        raise CobraSBMLError("Input format for 'filename' is not supported.")

    return doc


def write_sbml_model(cobra_model, filename, use_fbc_package=True,
                     f_replace=F_REPLACE, **kwargs):
    """Writes cobra model to filename.

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

    Parameters
    ----------
    cobra_model : cobra.core.Model
        Model instance which is written to SBML
    filename : string
        path to which the model is written
    use_fbc_package : boolean {True, False}
        should the fbc package be used
    f_replace: dict of replacement functions for id replacement
    """
    if not use_fbc_package:
        # FIXME: get completely rid of the legacy non-sense
        # legacy cobra without fbc
        LOGGER.warning("Constrained based models should be stored with fbc-v2,"
                       "By setting 'use_fbc_package' using legacy writer. No support from here on.")
        write_sbml2(cobra_model, filename, use_fbc_package=False, **kwargs)
    else:
        doc = _model_to_sbml(cobra_model, f_replace=f_replace, **kwargs)
        libsbml.writeSBMLToFile(doc, filename)


def _sbml_to_model(doc, number=float, f_replace=None, **kwargs):
    """Creates cobra model from SBMLDocument.

    Parameters
    ----------
    doc: libsbml.SBMLDocument
    number: data type of stoichiometry: {float, int}
        In which data type should the stoichiometry be parsed.
    f_replace : dict of replacement functions for id replacement

    Returns
    -------
    cobra.core.Model
    """
    if f_replace is None:
        f_replace = {}

    # SBML model
    model = doc.getModel()  # type: libsbml.Model
    if model is None:
        raise CobraSBMLError("No SBML model detected in file.")
    model_fbc = model.getPlugin("fbc")  # type: libsbml.FbcModelPlugin

    if not model_fbc:
        LOGGER.warning("Model does not contain FBC information.")
    else:
        if not model_fbc.isSetStrict():
            LOGGER.warning('Loading SBML model without fbc:strict="true"')

        # fbc-v1 (legacy)
        doc_fbc = doc.getPlugin("fbc")  # type: libsbml.FbcSBMLDocumentPlugin
        fbc_version = doc_fbc.getPackageVersion()

        if fbc_version == 1:
            LOGGER.warning("Loading SBML with fbc-v1 (models should be encoded"
                           " using fbc-v2)")
            conversion_properties = libsbml.ConversionProperties()
            conversion_properties.addOption("convert fbc v1 to fbc v2", True,
                                            "Convert FBC-v1 model to FBC-v2")
            result = doc.convert(conversion_properties)
            if result != libsbml.LIBSBML_OPERATION_SUCCESS:
                raise Exception("Conversion of SBML fbc v1 to fbc v2 failed")

    # Model
    cmodel = Model(model.id)
    cmodel.name = model.name

    # meta information
    meta = {

        "level": model.getLevel(),
        "version": model.getVersion(),
        "packages": []
    }
    # History
    creators = []
    created = None
    if model.isSetModelHistory():
        history = model.getModelHistory()  # type: libsbml.ModelHistory

        if history.isSetCreatedDate():
            created = history.getCreatedDate()

        for c in history.getListCreators():  # type: libsbml.ModelCreator
            creators.append({
             "familyName": c.getFamilyName() if c.isSetFamilyName() else None,
             "givenName": c.getGivenName() if c.isSetGivenName() else None,
             "organisation": c.getOrganisation() if c.isSetOrganisation() else None,
             "email": c.getEmail() if c.isSetEmail() else None,
            })

    meta["creators"] = creators
    meta["created"] = created
    meta["notes"] = _parse_notes(doc)

    info = "SBML L{}V{}".format(model.getLevel(), model.getVersion())
    packages = {}
    for k in range(doc.getNumPlugins()):
        plugin = doc.getPlugin(k)  # type:libsbml.SBasePlugin
        key, value = plugin.getPackageName(), plugin.getPackageVersion()
        packages[key] = value
        info += ", {}-v{}".format(key, value)

    meta["info"] = info
    meta["packages"] = packages
    meta["notes"] = _parse_notes(doc)
    meta["annotation"] = _parse_annotations(doc)
    cmodel._sbml = meta

    # notes and annotations
    cmodel.notes = _parse_notes(model)
    cmodel.annotation = _parse_annotations(model)

    # Compartments
    cmodel.compartments = {c.id: c.name for c in model.compartments}

    # Species
    boundary_ids = set()
    metabolites = []
    for s in model.getListOfSpecies():  # type: libsbml.Species
        sid = _check_required(s, s.id, "id")
        if f_replace and F_SPECIE in f_replace:
            sid = f_replace[F_SPECIE](sid)

        met = Metabolite(sid)
        met.name = s.name
        met.notes = _parse_notes(s)
        met.annotation = _parse_annotations(s)
        met.compartment = s.compartment

        s_fbc = s.getPlugin("fbc")  # type: libsbml.FbcSpeciesPlugin
        if s_fbc:
            met.charge = s_fbc.getCharge()
            met.formula = s_fbc.getChemicalFormula()
        else:
            if s.isSetCharge():
                LOGGER.warning("Use of charge attribute is highly "
                               "discouraged '%s" % s)
                met.charge = s.getCharge()
            else:
                if met.notes and 'CHARGE' in met.notes:
                    LOGGER.warning("Use of CHARGE note is discouraged '%s" % s)
                    try:
                        met.charge = int(met.notes['CHARGE'])
                    except ValueError:
                        # handle nan, na, NA, ...
                        pass

            if met.notes and 'FORMULA' in met.notes:
                LOGGER.warning("Use of FORMULA note is discouraged '%s" % s)
                met.formula = met.notes['FORMULA']

        # Detect boundary metabolites - In case they have been mistakenly
        # added. They should not actually appear in a model
        if s.getBoundaryCondition() is True:
            boundary_ids.add(s.id)

        metabolites.append(met)

    cmodel.add_metabolites(metabolites)

    # Genes
    if model_fbc:
        for gp in model_fbc.getListOfGeneProducts():  # noqa: E501 type: libsbml.GeneProduct
            gid = gp.id
            if f_replace and F_GENE in f_replace:
                gid = f_replace[F_GENE](gid)
            gene = Gene(gid)
            gene.name = gp.name
            if gene.name is None:
                gene.name = gp.getId()
            gene.annotation = _parse_annotations(gp)
            gene.notes = _parse_notes(gp)

            cmodel.genes.append(gene)

    # GPR rules
    def process_association(ass):
        """ Recursively convert gpr association to a gpr string. """
        if ass.isFbcOr():
            return " ".join(
                ["(", ' or '.join(process_association(c)
                                  for c in ass.getListOfAssociations()), ")"]
            )
        elif ass.isFbcAnd():
            return " ".join(
                ["(", ' and '.join(process_association(c)
                                   for c in ass.getListOfAssociations()), ")"])
        elif ass.isGeneProductRef():
            gid = ass.getGeneProduct()
            if f_replace and F_GENE in f_replace:
                return f_replace[F_GENE](gid)
            else:
                return gid

    # Reactions
    reactions = []
    for r in model.getListOfReactions():  # type: libsbml.Reaction
        rid = _check_required(r, r.id, "id")
        if f_replace and F_REACTION in f_replace:
            rid = f_replace[F_REACTION](rid)
        reaction = Reaction(rid)
        reaction.name = r.name
        reaction.annotation = _parse_annotations(r)
        reaction.notes = _parse_notes(r)

        # set bounds
        r_fbc = r.getPlugin("fbc")  # type: libsbml.FbcReactionPlugin
        if r_fbc:
            # bounds in fbc
            lb_id = _check_required(r_fbc, r_fbc.getLowerFluxBound(),
                                    "lowerFluxBound")
            ub_id = _check_required(r_fbc, r_fbc.getUpperFluxBound(),
                                    "upperFluxBound")
            p_lb = model.getParameter(lb_id)
            p_ub = model.getParameter(ub_id)

            if p_lb.constant and (p_lb.value is not None):
                reaction.lower_bound = p_lb.value
            else:
                raise CobraSBMLError("No constant bound '%s' for "
                                     "reaction '%s" % (p_lb, r))

            if p_ub.constant and (p_ub.value is not None):
                reaction.upper_bound = p_ub.value
            else:
                raise CobraSBMLError("No constant bound '%s' for "
                                     "reaction '%s" % (p_ub, r))

        elif r.isSetKineticLaw():
            # some legacy models encode bounds in kinetic laws
            klaw = r.getKineticLaw()  # type: libsbml.KineticLaw
            p_lb = klaw.getParameter("LOWER_BOUND")
            if p_lb:
                reaction.lower_bound = p_lb.value
            else:
                raise CobraSBMLError("Missing flux bounds on reaction %s" % r)
            p_ub = klaw.getParameter("UPPER_BOUND")
            if p_ub:
                reaction.upper_bound = p_ub.value
            else:
                raise CobraSBMLError("Missing flux bounds on reaction %s" % r)
            LOGGER.warning("Bounds encoded in KineticLaw for '%s" % r)
        else:
            raise CobraSBMLError("No flux bounds on reaction '%s'" % r)

        reactions.append(reaction)

        # parse equation
        stoichiometry = defaultdict(lambda: 0)
        for sref in r.getListOfReactants():  # type: libsbml.SpeciesReference
            sid = sref.getSpecies()
            if f_replace and F_SPECIE in f_replace:
                sid = f_replace[F_SPECIE](sid)
            stoichiometry[sid] -= number(_check_required(sref,
                                                         sref.stoichiometry,
                                                         "stoichiometry"))

        for sref in r.getListOfProducts():  # type: libsbml.SpeciesReference
            sid = sref.getSpecies()
            if f_replace and F_SPECIE in f_replace:
                sid = f_replace[F_SPECIE](sid)
            stoichiometry[sid] += number(_check_required(sref,
                                                         sref.stoichiometry,
                                                         "stoichiometry"))

        # needs to have keys of metabolite objects, not ids
        object_stoichiometry = {}
        for met_id in stoichiometry:

            # FIXME: THIS IS INCORRECT BEHAVIOUR
            # boundary species must be created but additional exchange
            # reactions must be added to the model
            if met_id in boundary_ids:
                LOGGER.warning("Boundary metabolite '%s' used in "
                               "reaction '%s'" % (met_id, reaction.id))
                continue
            try:
                metabolite = cmodel.metabolites.get_by_id(met_id)
            except KeyError:
                LOGGER.warning("Ignoring unknown metabolite '%s' in "
                               "reaction %s" % (met_id, reaction.id))
                continue
            object_stoichiometry[metabolite] = stoichiometry[met_id]
        reaction.add_metabolites(object_stoichiometry)

        # GPR
        if r_fbc:
            gpr = ''
            gpa = r_fbc.getGeneProductAssociation()  # noqa: E501 type: libsbml.GeneProductAssociation
            if gpa is not None:
                # type: libsbml.FbcAssociation
                association = gpa.getAssociation()
                gpr = process_association(association)
        else:
            # fallback to notes information
            gpr = reaction.notes.get('GENE_ASSOCIATION', '')
            if f_replace and F_GENE in f_replace:
                gpr = " ".join(f_replace[F_GENE](t) for t in gpr.split(' '))

        # remove outside parenthesis, if any
        if gpr.startswith("(") and gpr.endswith(")"):
            gpr = gpr[1:-1].strip()

        reaction.gene_reaction_rule = gpr

    cmodel.add_reactions(reactions)

    # Objective
    if model_fbc:
        obj_list = model_fbc.getListOfObjectives()  # noqa: E501 type: libsbml.ListOfObjectives
        if obj_list is None:
            LOGGER.warning("listOfObjectives element not found")
        elif obj_list.size() == 0:
            LOGGER.warning("No objective in listOfObjectives")
        elif not obj_list.getActiveObjective():
            LOGGER.warning("No active objective in listOfObjectives")
        else:
            obj_id = obj_list.getActiveObjective()
            obj = model_fbc.getObjective(obj_id)  # type: libsbml.Objective
            obj_direction = LONG_SHORT_DIRECTION[obj.getType()]

            coefficients = {}

            for flux_obj in obj.getListOfFluxObjectives():  # noqa: E501 type: libsbml.FluxObjective
                rid = flux_obj.getReaction()
                if f_replace and F_REACTION in f_replace:
                    rid = f_replace[F_REACTION](rid)
                try:
                    objective_reaction = cmodel.reactions.get_by_id(rid)
                except KeyError:
                    raise CobraSBMLError("Objective reaction '%s' "
                                         "not found" % rid)
                try:
                    coefficients[objective_reaction] = number(
                        flux_obj.getCoefficient()
                    )
                except ValueError as e:
                    LOGGER.warning(str(e))
            set_objective(cmodel, coefficients)
            cmodel.solver.objective.direction = obj_direction

    return cmodel


def _model_to_sbml(cobra_model, f_replace=None, units=True):
    """Convert Cobra model to SBMLDocument.

    Parameters
    ----------
    cobra_model : cobra.core.Model
        Cobra model instance
    f_replace : dict of replacement functions
        Replacement to apply on identifiers.
    units : boolean
        Should the FLUX_UNITS be written in the SBMLDocument.

    Returns
    -------
    libsbml.SBMLDocument
    """
    if f_replace is None:
        f_replace = {}

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
        # type: libsbml.UnitDefinition
        flux_udef = model.createUnitDefinition()
        flux_udef.setId(UNITS_FLUX[0])
        for u in UNITS_FLUX[1]:
            unit = flux_udef.createUnit()  # type: libsbml.Unit
            unit.setKind(u.kind)
            unit.setExponent(u.exponent)
            unit.setScale(u.scale)
            unit.setMultiplier(u.multiplier)

    # Flux bounds
    def _create_bound(model, reaction, bound_type):
        """Creates bound in model for given reaction.

        Adds the parameters for the bounds to the SBML model.

        Parameters
        ----------
        model : libsbml.Model
            SBML model instance
        reaction : cobra.core.Reaction
            Cobra reaction instance from which the bounds are read.
        bound_type : {LOWER_BOUND, UPPER_BOUND}
            Type of bound

        Returns
        -------
        Id of bound parameter.
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
            rid = reaction.id
            if f_replace and F_REACTION_REV in f_replace:
                rid = f_replace[F_REACTION_REV](rid)
            pid = rid + "_" + bound_type
            _create_parameter(model, pid=pid, value=value, sbo=SBO_FLUX_BOUND)
            return pid

    def _create_parameter(model, pid, value, sbo=None, constant=True):
        """Create parameter in SBML model."""
        p = model.createParameter()  # type: libsbml.Parameter
        p.setId(pid)
        p.setValue(value)
        p.setConstant(constant)
        if sbo:
            p.setSBOTerm(sbo)
        if units:
            p.setUnits(flux_udef.getId())

    # minimum and maximum value from model
    if len(cobra_model.reactions) > 0:
        min_value = min(cobra_model.reactions.list_attr("lower_bound"))
        max_value = max(cobra_model.reactions.list_attr("upper_bound"))
    else:
        min_value = LOWER_BOUND
        max_value = UPPER_BOUND

    _create_parameter(model, pid=LOWER_BOUND_ID,
                      value=min_value, sbo=SBO_DEFAULT_FLUX_BOUND)
    _create_parameter(model, pid=UPPER_BOUND_ID,
                      value=max_value, sbo=SBO_DEFAULT_FLUX_BOUND)
    _create_parameter(model, pid=ZERO_BOUND_ID,
                      value=0, sbo=SBO_DEFAULT_FLUX_BOUND)

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
        if f_replace and F_SPECIE_REV in f_replace:
            mid = f_replace[F_SPECIE_REV](mid)
        s.setId(mid)
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
        gid = gene.id
        if f_replace and F_GENE_REV in f_replace:
            gid = f_replace[F_GENE_REV](gid)
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
    reaction_coefficients = linear_reaction_coefficients(cobra_model)
    for reaction in cobra_model.reactions:
        rid = reaction.id
        if f_replace and F_REACTION_REV in f_replace:
            rid = f_replace[F_REACTION_REV](rid)
        r = model.createReaction()  # type: libsbml.Reaction
        r.setId(rid)
        r.setName(reaction.name)
        r.setFast(False)
        r.setReversible((reaction.lower_bound < 0))

        annotate_sbase_from_cobra(r, reaction)

        # stoichiometry
        for metabolite, stoichiomety in iteritems(reaction._metabolites):
            sid = metabolite.id
            if f_replace and F_SPECIE_REV in f_replace:
                sid = f_replace[F_SPECIE_REV](sid)
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
            # type: libsbml.GeneProductAssociation
            gpa = r_fbc.createGeneProductAssociation()
            # replace ids
            if f_replace and F_GENE_REV in f_replace:
                tokens = gpr.split(' ')
                for k in range(len(tokens)):
                    if tokens[k] not in ['and', 'or', '(', ')']:
                        tokens[k] = f_replace[F_GENE_REV](tokens[k])
                gpr = " ".join(tokens)


            res = gpa.setAssociation(gpr)
            print("write", gpr, res)


        # objective coefficients
        if reaction_coefficients.get(reaction, 0) != 0:
            # type: libsbml.FluxObjective
            flux_obj = objective.createFluxObjective()
            flux_obj.setReaction(rid)
            flux_obj.setCoefficient(reaction.objective_coefficient)

    return doc


def _check_required(sbase, value, attribute):
    """Get required attribute from the SBase.

    Parameters
    ----------
    sbase : libsbml.SBase
    value : existing value
    attribute: name of attribute

    Returns
    -------
    attribute value (or value if already set)
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


# ----------------------
# Notes & Annotations
# ----------------------
def _parse_notes(sbase):
    """Creates dictionary of notes.

    Parameters
    ----------
    sbase : libsbml.SBase

    Returns
    -------
    dict of notes
    """
    notes = sbase.getNotesString()
    if notes and len(notes) > 0:
        pattern = r"<p>\s*(\w+)\s*:\s*([\w|\s]+)<"
        matches = re.findall(pattern, notes)
        d = {k.strip(): v.strip() for (k, v) in matches}
        return {k: v for k, v in d.items() if len(v) > 0}
    else:
        return {}


# FIXME: currently only the terms, but not the qualifier are parsed
URL_IDENTIFIERS_PATTERN = r"^http[s]{0,1}://identifiers.org/(.+)/(.+)"
URL_IDENTIFIERS_PREFIX = r"http://identifiers.org"
BIOLOGICAL_QUALIFIER_TYPES = set(["BQB_IS", "BQB_HAS_PART", "BQB_IS_PART_OF",
                                 "BQB_IS_VERSION_OF", "BQB_HAS_VERSION",
                                 "BQB_IS_HOMOLOG_TO", "BQB_IS_DESCRIBED_BY",
                                 "BQB_IS_ENCODED_BY", "BQB_ENCODES",
                                 "BQB_OCCURS_IN", "BQB_HAS_PROPERTY",
                                 "BQB_IS_PROPERTY_OF", "BQB_HAS_TAXON",
                                 "BQB_UNKNOWN"])


def _parse_annotations(sbase):
    """Parses cobra annotations from a given SBase object.

    Annotations are dictionaries with the providers as keys.

    Parameters
    ----------
    sbase : libsbml.SBase
        SBase from which the SBML annotations are read

    Returns
    -------
    dict (annotation dictionary)
    """
    annotation = {}

    # SBO term
    if sbase.isSetSBOTerm():
        annotation["sbo"] = sbase.getSBOTermID()

    # RDF annotation
    cvterms = sbase.getCVTerms()
    if cvterms is None:
        return

    for cvterm in cvterms:  # type: libsbml.CVTerm
        for k in range(cvterm.getNumResources()):
            uri = cvterm.getResourceURI(k)

            # FIXME: read and store the qualifier

            tokens = uri.split('/')
            if len(tokens) != 5 or not tokens[2] == "identifiers.org":
                LOGGER.warning("%s does not conform to "
                               "http(s)://identifiers.org/collection/id" % uri)
                continue

            provider, identifier = tokens[3], tokens[4]
            if provider in annotation:
                if isinstance(annotation[provider], string_types):
                    annotation[provider] = [annotation[provider]]
                annotation[provider].append(identifier)
            else:
                # FIXME: always in list
                annotation[provider] = identifier

    return annotation


def annotate_sbase_from_cobra(sbase, cobj):
    """Annotate SBase based on cobra object annotations.

    Parameters
    ----------
    sbase : libsbml.SBase
        SBML object to annotate
    cobj : cobra object
        cobra object with annotation information
    """

    if len(cobj.annotation) == 0:
        return

    # FIXME: currently no support for qualifiers
    qualifier_type = libsbml.BIOLOGICAL_QUALIFIER
    qualifier = libsbml.BQB_IS

    meta_id = "meta_{}".format(sbase.id)
    sbase.setMetaId(meta_id)

    # rdf_items = []
    for provider, identifiers in iteritems(cobj.annotation):
        if provider in ["SBO", "sbo"]:
            if provider == "SBO":
                logging.warning("'SBO' provider is deprecated, "
                                "use 'sbo' provider instead")
            sbase.setSBOTerm(identifiers)
        else:
            if isinstance(identifiers, string_types):
                identifiers = (identifiers,)

            for identifier in identifiers:
                cv = libsbml.CVTerm()  # type: libsbml.CVTerm
                cv.setQualifierType(qualifier_type)
                if qualifier_type == libsbml.BIOLOGICAL_QUALIFIER:
                    cv.setBiologicalQualifierType(qualifier)
                elif qualifier_type == libsbml.MODEL_QUALIFIER:
                    cv.setModelQualifierType(qualifier)
                else:
                    raise CobraSBMLError('Unsupported qualifier: '
                                         '%s' % qualifier)
                cv.addResource("%s/%s/%s" % (URL_IDENTIFIERS_PREFIX, provider,
                                             identifier))
                sbase.addCVTerm(cv)


# -----------------------------------
# Validation
# -----------------------------------
def validate_sbml_model(filename, use_libsbml=False, check_model=True,
                        internal_consistency=True,
                        check_units_consistency=False,
                        check_modeling_practice=False):
    """Validate SBML model and returns the model along with a list of errors.

    Parameters
    ----------
    filename : str
        The filename (or SBML string) of the SBML model to be validated.
    use_libsbml : boolean {True, False}
        Perform SBML validation via libsbml. This can take some time.
    internal_consistency: boolean {True, False}
        Check internal consistency.
    check_units_consistency: boolean {True, False}
        Check consistency of units.
    check_modeling_practice: boolean {True, False}
        Check modeling practise.
    check_model: boolean {True, False}
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
    errors = {key: [] for key in ("validator", "warnings", "other",
                                  "SBML errors")}
    if use_libsbml:
        for key in ["SBML_FATAL", "SBML ERROR", "SBML_SCHEMA_ERROR",
                    "SBML_WARNING"]:
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
        doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY,
                                 check_units_consistency)
        doc.setConsistencyChecks(libsbml.LIBSBML_CAT_MODELING_PRACTICE,
                                 check_modeling_practice)

        # validate the document
        if internal_consistency:
            doc.checkInternalConsistency()
        doc.checkConsistency()

        for k in range(doc.getNumErrors()):
            e = doc.getError(k)
            sev = e.getSeverity()
            if sev == libsbml.LIBSBML_SEV_FATAL:
                err(_error_string(e), "SBML_FATAL")
            elif sev == libsbml.LIBSBML_SEV_ERROR:
                err(_error_string(e), "SBML_ERROR")
            elif sev == libsbml.LIBSBML_SEV_SCHEMA_ERROR:
                err(_error_string(e), "SBML_SCHEMA_ERROR")
            elif sev == libsbml.LIBSBML_SEV_WARNING:
                err(_error_string(e), "SBML_WARNING")

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


def _error_string(error, k=None):
    """String representation of SBMLError.

    Parameters
    ----------
    error : libsbml.SBMLError
    k : index of error

    Returns
    -------
    string representation of error
    """
    package = error.getPackage()
    if package == '':
        package = 'core'

    error_str = 'E{}: {} ({}, L{}, {})  \n' \
                '{}\n' \
                '[{}] {}\n' \
                '{}\n'.format(
                    k, error.getCategoryAsString(), package, error.getLine(),
                    'code',
                    '-' * 60,
                    error.getSeverityAsString(), error.getShortMessage(),
                    error.getMessage())
    return error_str
