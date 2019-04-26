"""
SBML import and export using python-libsbml(-experimental).

- The SBML importer supports all versions of SBML and the fbc package.
- The SBML exporter writes SBML L3 models.
- Annotation information is stored on the cobrapy objects
- Information from the group package is read

Parsing of fbc models was implemented as efficient as possible, whereas
(discouraged) fallback solutions are not optimized for efficiency.

Notes are only supported in a minimal way relevant for constraint-based
models. I.e., structured information from notes in the form
   <p>key: value</p>
is read into the Object.notes dictionary when reading SBML files.
On writing the Object.notes dictionary is serialized to the SBML
notes information.

Annotations are read in the Object.annotation fields.

Some SBML related issues are still open, please refer to the respective issue:
- update annotation format and support qualifiers (depends on decision
    for new annotation format; https://github.com/opencobra/cobrapy/issues/684)
- write compartment annotations and notes (depends on updated first-class
    compartments; see https://github.com/opencobra/cobrapy/issues/760)
- support compression on file handles (depends on solution for
    https://github.com/opencobra/cobrapy/issues/812)
"""

from __future__ import absolute_import

import datetime
import logging
import os
import re
import traceback
from collections import defaultdict, namedtuple
from copy import deepcopy
from sys import platform

import libsbml
from six import iteritems, string_types

import cobra
from cobra.core import Gene, Group, Metabolite, Model, Reaction
from cobra.manipulation.validate import check_metabolite_compartment_formula
from cobra.util.solver import linear_reaction_coefficients, set_objective


try:
    from cStringIO import StringIO  # Python 2
except ImportError:
    from io import StringIO


class CobraSBMLError(Exception):
    """ SBML error class. """
    pass


LOGGER = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Defaults and constants for writing SBML
# -----------------------------------------------------------------------------
config = cobra.Configuration()  # for default bounds
LOWER_BOUND_ID = "cobra_default_lb"
UPPER_BOUND_ID = "cobra_default_ub"
ZERO_BOUND_ID = "cobra_0_bound"

BOUND_MINUS_INF = "minus_inf"
BOUND_PLUS_INF = "plus_inf"


SBO_FBA_FRAMEWORK = "SBO:0000624"
SBO_DEFAULT_FLUX_BOUND = "SBO:0000626"
SBO_FLUX_BOUND = "SBO:0000625"
SBO_EXCHANGE_REACTION = "SBO:0000627"

LONG_SHORT_DIRECTION = {'maximize': 'max', 'minimize': 'min'}
SHORT_LONG_DIRECTION = {'min': 'minimize', 'max': 'maximize'}

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

# -----------------------------------------------------------------------------
# Functions for id replacements (import/export)
# -----------------------------------------------------------------------------
SBML_DOT = "__SBML_DOT__"


def _clip(sid, prefix):
    """Clips a prefix from the beginning of a string if it exists."""
    return sid[len(prefix):] if sid.startswith(prefix) else sid


def _f_gene(sid, prefix="G_"):
    """Clips gene prefix from id."""
    sid = sid.replace(SBML_DOT, ".")
    return _clip(sid, prefix)


def _f_gene_rev(sid, prefix="G_"):
    """Adds gene prefix to id."""
    return prefix + sid.replace(".", SBML_DOT)


def _f_specie(sid, prefix="M_"):
    """Clips specie/metabolite prefix from id."""
    return _clip(sid, prefix)


def _f_specie_rev(sid, prefix="M_"):
    """Adds specie/metabolite prefix to id."""
    return prefix + sid


def _f_reaction(sid, prefix="R_"):
    """Clips reaction prefix from id."""
    return _clip(sid, prefix)


def _f_reaction_rev(sid, prefix="R_"):
    """Adds reaction prefix to id."""
    return prefix + sid


F_GENE = "F_GENE"
F_GENE_REV = "F_GENE_REV"
F_SPECIE = "F_SPECIE"
F_SPECIE_REV = "F_SPECIE_REV"
F_REACTION = "F_REACTION"
F_REACTION_REV = "F_REACTION_REV"

F_REPLACE = {
    F_GENE: _f_gene,
    F_GENE_REV: _f_gene_rev,
    F_SPECIE: _f_specie,
    F_SPECIE_REV: _f_specie_rev,
    F_REACTION: _f_reaction,
    F_REACTION_REV: _f_reaction_rev,
}


# -----------------------------------------------------------------------------
# Read SBML
# -----------------------------------------------------------------------------
def read_sbml_model(filename, number=float, f_replace=F_REPLACE,
                    set_missing_bounds=False, **kwargs):
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
    set_missing_bounds : boolean flag to set missing bounds
        Missing bounds are set to default bounds in configuration.

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
        return _sbml_to_model(doc,
                              number=number,
                              f_replace=f_replace,
                              set_missing_bounds=set_missing_bounds,
                              **kwargs)
    except IOError as e:
        raise e

    except Exception:
        LOGGER.error(traceback.print_exc())
        raise CobraSBMLError(
            "Something went wrong reading the SBML model. Most likely the SBML"
            " model is not valid. Please check that your model is valid using "
            "the `cobra.io.sbml.validate_sbml_model` function or via the "
            "online validator at http://sbml.org/validator .\n"
            "\t`(model, errors) = validate_sbml_model(filename)`"
            "\nIf the model is valid and cannot be read please open an issue "
            "at https://github.com/opencobra/cobrapy/issues .")


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
        if ("win" in platform) and (len(filename) < 260) \
                and os.path.exists(filename):
            # path (win)
            doc = libsbml.readSBMLFromFile(filename)  # noqa: E501 type: libsbml.SBMLDocument
        elif ("win" not in platform) and os.path.exists(filename):
            # path other
            doc = libsbml.readSBMLFromFile(filename)  # noqa: E501 type: libsbml.SBMLDocument
        else:
            # string representation
            if "<sbml" not in filename:
                raise IOError("The file with 'filename' does not exist, "
                              "or is not an SBML string. Provide the path to "
                              "an existing SBML file or a valid SBML string "
                              "representation: \n%s", filename)

            doc = libsbml.readSBMLFromString(filename)  # noqa: E501 type: libsbml.SBMLDocument

    elif hasattr(filename, "read"):
        # file handle
        doc = libsbml.readSBMLFromString(filename.read())  # noqa: E501 type: libsbml.SBMLDocument
    else:
        raise CobraSBMLError("Input type '%s' for 'filename' is not supported."
                             " Provide a path, SBML str, "
                             "or file handle.", type(filename))

    return doc


def _sbml_to_model(doc, number=float, f_replace=F_REPLACE,
                   set_missing_bounds=False, **kwargs):
    """Creates cobra model from SBMLDocument.

    Parameters
    ----------
    doc: libsbml.SBMLDocument
    number: data type of stoichiometry: {float, int}
        In which data type should the stoichiometry be parsed.
    f_replace : dict of replacement functions for id replacement
    set_missing_bounds : flag to set missing bounds

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
        LOGGER.warning("Model does not contain SBML fbc package information.")
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
    cobra_model = Model(model.getId())
    cobra_model.name = model.getName()

    # meta information
    meta = {
        "model.id": model.getId(),
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
                "familyName": c.getFamilyName() if c.isSetFamilyName() else None,  # noqa: E501
                "givenName": c.getGivenName() if c.isSetGivenName() else None,
                "organisation": c.getOrganisation() if c.isSetOrganisation() else None,  # noqa: E501
                "email": c.getEmail() if c.isSetEmail() else None,
            })

    meta["creators"] = creators
    meta["created"] = created
    meta["notes"] = _parse_notes_dict(doc)
    meta["annotation"] = _parse_annotations(doc)

    info = "<{}> SBML L{}V{}".format(model.getId(),
                                     model.getLevel(), model.getVersion())
    packages = {}
    for k in range(doc.getNumPlugins()):
        plugin = doc.getPlugin(k)  # type:libsbml.SBasePlugin
        key, value = plugin.getPackageName(), plugin.getPackageVersion()
        packages[key] = value
        info += ", {}-v{}".format(key, value)
        if key not in ["fbc", "groups"]:
            LOGGER.warning("SBML package '%s' not supported by cobrapy,"
                           "information is not parsed", key)
    meta["info"] = info
    meta["packages"] = packages
    cobra_model._sbml = meta

    # notes and annotations
    cobra_model.notes = _parse_notes_dict(model)
    cobra_model.annotation = _parse_annotations(model)

    # Compartments
    # FIXME: update with new compartments
    cobra_model.compartments = {c.getId(): c.getName()
                                for c in model.getListOfCompartments()}

    # Species
    metabolites = []
    boundary_metabolites = []
    if model.getNumSpecies() == 0:
        LOGGER.warning("No metabolites in model")

    for specie in model.getListOfSpecies():  # type: libsbml.Species
        sid = _check_required(specie, specie.getId(), "id")
        if f_replace and F_SPECIE in f_replace:
            sid = f_replace[F_SPECIE](sid)

        met = Metabolite(sid)
        met.name = specie.getName()
        met.notes = _parse_notes_dict(specie)
        met.annotation = _parse_annotations(specie)
        met.compartment = specie.getCompartment()

        specie_fbc = specie.getPlugin("fbc")  # type: libsbml.FbcSpeciesPlugin
        if specie_fbc:
            met.charge = specie_fbc.getCharge()
            met.formula = specie_fbc.getChemicalFormula()
        else:
            if specie.isSetCharge():
                LOGGER.warning("Use of the species charge attribute is "
                               "discouraged, use fbc:charge "
                               "instead: %s", specie)
                met.charge = specie.getCharge()
            else:
                if 'CHARGE' in met.notes:
                    LOGGER.warning("Use of CHARGE in the notes element is "
                                   "discouraged, use fbc:charge "
                                   "instead: %s", specie)
                    try:
                        met.charge = int(met.notes['CHARGE'])
                    except ValueError:
                        # handle nan, na, NA, ...
                        pass

            if 'FORMULA' in met.notes:
                LOGGER.warning("Use of FORMULA in the notes element is "
                               "discouraged, use fbc:chemicalFormula "
                               "instead: %s", specie)
                met.formula = met.notes['FORMULA']

        # Detect boundary metabolites
        if specie.getBoundaryCondition() is True:
            boundary_metabolites.append(met)

        metabolites.append(met)

    cobra_model.add_metabolites(metabolites)

    # Add exchange reactions for boundary metabolites
    ex_reactions = []
    for met in boundary_metabolites:
        ex_rid = "EX_{}".format(met.id)
        ex_reaction = Reaction(ex_rid)
        ex_reaction.name = ex_rid
        ex_reaction.annotation = {
            'sbo': SBO_EXCHANGE_REACTION
        }
        ex_reaction.lower_bound = -float("Inf")
        ex_reaction.upper_bound = float("Inf")
        LOGGER.warning("Adding exchange reaction %s for boundary metabolite: "
                       "%s" % (ex_reaction.id, met.id))
        # species is reactant
        ex_reaction.add_metabolites({met: -1})
        ex_reactions.append(ex_reaction)
    cobra_model.add_reactions(ex_reactions)

    # Genes
    if model_fbc:
        for gp in model_fbc.getListOfGeneProducts():  # noqa: E501 type: libsbml.GeneProduct
            gid = gp.getId()
            if f_replace and F_GENE in f_replace:
                gid = f_replace[F_GENE](gid)
            cobra_gene = Gene(gid)
            cobra_gene.name = gp.getName()
            if cobra_gene.name is None:
                cobra_gene.name = gid
            cobra_gene.annotation = _parse_annotations(gp)
            cobra_gene.notes = _parse_notes_dict(gp)

            cobra_model.genes.append(cobra_gene)
    else:
        for cobra_reaction in model.getListOfReactions():  # noqa: E501 type: libsbml.Reaction
            # fallback to notes information
            notes = _parse_notes_dict(cobra_reaction)
            if "GENE ASSOCIATION" in notes:
                gpr = notes['GENE ASSOCIATION']
            elif "GENE_ASSOCIATION" in notes:
                gpr = notes['GENE_ASSOCIATION']
            else:
                gpr = ''

            if len(gpr) > 0:
                gpr = gpr.replace("(", ";")
                gpr = gpr.replace(")", ";")
                gpr = gpr.replace("or", ";")
                gpr = gpr.replace("and", ";")
                gids = [t.strip() for t in gpr.split(';')]

                # create missing genes
                for gid in gids:
                    if f_replace and F_GENE in f_replace:
                        gid = f_replace[F_GENE](gid)

                    if gid not in cobra_model.genes:
                        cobra_gene = Gene(gid)
                        cobra_gene.name = gid
                        cobra_model.genes.append(cobra_gene)

    # GPR rules
    def process_association(ass):
        """ Recursively convert gpr association to a gpr string.
        Defined as inline functions to not pass the replacement dict around.
        """
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
    missing_bounds = False
    reactions = []
    if model.getNumReactions() == 0:
        LOGGER.warning("No reactions in model")

    for reaction in model.getListOfReactions():  # type: libsbml.Reaction
        rid = _check_required(reaction, reaction.getId(), "id")
        if f_replace and F_REACTION in f_replace:
            rid = f_replace[F_REACTION](rid)
        cobra_reaction = Reaction(rid)
        cobra_reaction.name = reaction.getName()
        cobra_reaction.annotation = _parse_annotations(reaction)
        cobra_reaction.notes = _parse_notes_dict(reaction)

        # set bounds
        p_ub, p_lb = None, None
        r_fbc = reaction.getPlugin("fbc")  # type: libsbml.FbcReactionPlugin
        if r_fbc:
            # bounds in fbc
            lb_id = r_fbc.getLowerFluxBound()
            if lb_id:
                p_lb = model.getParameter(lb_id)  # type: libsbml.Parameter
                if p_lb and p_lb.getConstant() and \
                        (p_lb.getValue() is not None):
                    cobra_reaction.lower_bound = p_lb.getValue()
                else:
                    raise CobraSBMLError("No constant bound '%s' for "
                                         "reaction: %s" % (p_lb, reaction))

            ub_id = r_fbc.getUpperFluxBound()
            if ub_id:
                p_ub = model.getParameter(ub_id)  # type: libsbml.Parameter
                if p_ub and p_ub.getConstant() and \
                        (p_ub.getValue() is not None):
                    cobra_reaction.upper_bound = p_ub.getValue()
                else:
                    raise CobraSBMLError("No constant bound '%s' for "
                                         "reaction: %s" % (p_ub, reaction))

        elif reaction.isSetKineticLaw():
            # some legacy models encode bounds in kinetic laws
            klaw = reaction.getKineticLaw()  # type: libsbml.KineticLaw
            p_lb = klaw.getParameter("LOWER_BOUND")  # noqa: E501 type: libsbml.LocalParameter
            if p_lb:
                cobra_reaction.lower_bound = p_lb.getValue()
            p_ub = klaw.getParameter("UPPER_BOUND")  # noqa: E501 type: libsbml.LocalParameter
            if p_ub:
                cobra_reaction.upper_bound = p_ub.getValue()

            if p_ub is not None or p_lb is not None:
                LOGGER.warning("Encoding LOWER_BOUND and UPPER_BOUND in "
                               "KineticLaw is discouraged, "
                               "use fbc:fluxBounds instead: %s", reaction)

        if p_lb is None:
            missing_bounds = True
            if set_missing_bounds:
                lower_bound = config.lower_bound
            else:
                lower_bound = -float("Inf")
            cobra_reaction.lower_bound = lower_bound
            LOGGER.warning("Missing lower flux bound set to '%s' for "
                           " reaction: '%s'", lower_bound, reaction)

        if p_ub is None:
            missing_bounds = True
            if set_missing_bounds:
                upper_bound = config.upper_bound
            else:
                upper_bound = float("Inf")
            cobra_reaction.upper_bound = upper_bound
            LOGGER.warning("Missing upper flux bound set to '%s' for "
                           " reaction: '%s'", upper_bound, reaction)

        # add reaction
        reactions.append(cobra_reaction)

        # parse equation
        stoichiometry = defaultdict(lambda: 0)
        for sref in reaction.getListOfReactants():  # noqa: E501 type: libsbml.SpeciesReference
            sid = sref.getSpecies()
            if f_replace and F_SPECIE in f_replace:
                sid = f_replace[F_SPECIE](sid)
            stoichiometry[sid] -= number(
                _check_required(sref, sref.getStoichiometry(),
                                "stoichiometry"))

        for sref in reaction.getListOfProducts():  # noqa: E501 type: libsbml.SpeciesReference
            sid = sref.getSpecies()
            if f_replace and F_SPECIE in f_replace:
                sid = f_replace[F_SPECIE](sid)
            stoichiometry[sid] += number(
                _check_required(sref, sref.getStoichiometry(),
                                "stoichiometry"))

        # convert to metabolite objects
        object_stoichiometry = {}
        for met_id in stoichiometry:
            metabolite = cobra_model.metabolites.get_by_id(met_id)
            object_stoichiometry[metabolite] = stoichiometry[met_id]
        cobra_reaction.add_metabolites(object_stoichiometry)

        # GPR
        if r_fbc:
            gpr = ''
            gpa = r_fbc.getGeneProductAssociation()  # noqa: E501 type: libsbml.GeneProductAssociation
            if gpa is not None:
                association = gpa.getAssociation()  # noqa: E501 type: libsbml.FbcAssociation
                gpr = process_association(association)
        else:
            # fallback to notes information
            notes = cobra_reaction.notes
            if "GENE ASSOCIATION" in notes:
                gpr = notes['GENE ASSOCIATION']
            elif "GENE_ASSOCIATION" in notes:
                gpr = notes['GENE_ASSOCIATION']
            else:
                gpr = ''

            if len(gpr) > 0:
                LOGGER.warning("Use of GENE ASSOCIATION or GENE_ASSOCIATION "
                               "in the notes element is discouraged, use "
                               "fbc:gpr instead: %s", reaction)
                if f_replace and F_GENE in f_replace:
                    gpr = " ".join(
                        f_replace[F_GENE](t) for t in gpr.split(' ')
                    )

        # remove outside parenthesis, if any
        if gpr.startswith("(") and gpr.endswith(")"):
            gpr = gpr[1:-1].strip()

        cobra_reaction.gene_reaction_rule = gpr

    cobra_model.add_reactions(reactions)

    # Objective
    obj_direction = "max"
    coefficients = {}
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

            for flux_obj in obj.getListOfFluxObjectives():  # noqa: E501 type: libsbml.FluxObjective
                rid = flux_obj.getReaction()
                if f_replace and F_REACTION in f_replace:
                    rid = f_replace[F_REACTION](rid)
                try:
                    objective_reaction = cobra_model.reactions.get_by_id(rid)
                except KeyError:
                    raise CobraSBMLError("Objective reaction '%s' "
                                         "not found" % rid)
                try:
                    coefficients[objective_reaction] = number(
                        flux_obj.getCoefficient()
                    )
                except ValueError as e:
                    LOGGER.warning(str(e))
    else:
        # some legacy models encode objective coefficients in kinetic laws
        for reaction in model.getListOfReactions():  # noqa: E501 type: libsbml.Reaction
            if reaction.isSetKineticLaw():
                klaw = reaction.getKineticLaw()  # type: libsbml.KineticLaw
                p_oc = klaw.getParameter(
                    "OBJECTIVE_COEFFICIENT")  # noqa: E501 type: libsbml.LocalParameter
                if p_oc:
                    rid = reaction.getId()
                    if f_replace and F_REACTION in f_replace:
                        rid = f_replace[F_REACTION](rid)
                    try:
                        objective_reaction = cobra_model.reactions.get_by_id(
                            rid)
                    except KeyError:
                        raise CobraSBMLError("Objective reaction '%s' "
                                             "not found", rid)
                    try:
                        coefficients[objective_reaction] = number(
                            p_oc.getValue())
                    except ValueError as e:
                        LOGGER.warning(str(e))

                    LOGGER.warning("Encoding OBJECTIVE_COEFFICIENT in "
                                   "KineticLaw is discouraged, "
                                   "use fbc:fluxObjective "
                                   "instead: %s", cobra_reaction)

    if len(coefficients) == 0:
        LOGGER.error("No objective coefficients in model. Unclear what should "
                     "be optimized")
    set_objective(cobra_model, coefficients)
    cobra_model.solver.objective.direction = obj_direction

    # parse groups
    model_groups = model.getPlugin("groups")  # type: libsbml.GroupsModelPlugin
    groups = []
    if model_groups:
        # calculate hashmaps to lookup objects in O(1)
        sid_map = {}
        metaid_map = {}
        for obj_list in [model.getListOfCompartments(),
                         model.getListOfSpecies(),
                         model.getListOfReactions(),
                         model_groups.getListOfGroups()]:

            for sbase in obj_list:  # type: libsbml.SBase
                if sbase.isSetId():
                    sid_map[sbase.getId()] = sbase
                if sbase.isSetMetaId():
                    metaid_map[sbase.getMetaId()] = sbase

        # create groups
        for group in model_groups.getListOfGroups():  # type: libsbml.Group
            cobra_group = Group(group.getId())
            cobra_group.name = group.getName()
            if group.isSetKind():
                cobra_group.kind = group.getKindAsString()
            cobra_group.annotation = _parse_annotations(group)
            cobra_group.notes = _parse_notes_dict(group)

            cobra_members = []
            for member in group.getListOfMembers():  # type: libsbml.Member
                if member.isSetIdRef():
                    obj = sid_map[member.getIdRef()]
                    # obj = doc.getElementBySId(member.getIdRef())
                elif member.isSetMetaIdRef():
                    obj = metaid_map[member.getMetaIdRef()]
                    # obj = doc.getElementByMetaId(member.getMetaIdRef())

                typecode = obj.getTypeCode()
                obj_id = obj.getId()

                # id replacements
                cobra_member = None
                if typecode == libsbml.SBML_SPECIES:
                    if f_replace and F_SPECIE in f_replace:
                        obj_id = f_replace[F_SPECIE](obj_id)
                    cobra_member = cobra_model.metabolites.get_by_id(obj_id)
                elif typecode == libsbml.SBML_REACTION:
                    if f_replace and F_REACTION in f_replace:
                        obj_id = f_replace[F_REACTION](obj_id)
                    cobra_member = cobra_model.reactions.get_by_id(obj_id)
                elif typecode == libsbml.SBML_FBC_GENEPRODUCT:
                    if f_replace and F_GENE in f_replace:
                        obj_id = f_replace[F_GENE](obj_id)
                    cobra_member = cobra_model.genes.get_by_id(obj_id)
                else:
                    LOGGER.warning("Member %s could not be added to group %s."
                                   "unsupported type code: "
                                   "%s" % (member, group, typecode))

                if cobra_member:
                    cobra_members.append(cobra_member)

            cobra_group.add_members(cobra_members)
            groups.append(cobra_group)
    else:
        # parse deprecated subsystems on reactions
        groups_dict = {}
        for cobra_reaction in cobra_model.reactions:
            if "SUBSYSTEM" in cobra_reaction.notes:
                g_name = cobra_reaction.notes["SUBSYSTEM"]
                if g_name in groups_dict:
                    groups_dict[g_name].append(cobra_reaction)
                else:
                    groups_dict[g_name] = [cobra_reaction]

        for gid, cobra_members in groups_dict.items():
            cobra_group = Group(gid, name=gid, kind="collection")
            cobra_group.add_members(cobra_members)
            groups.append(cobra_group)

    cobra_model.add_groups(groups)

    # general hint for missing flux bounds
    if missing_bounds and not set_missing_bounds:
        LOGGER.warning("Missing flux bounds on reactions. As best practise "
                       "and to avoid confusion flux bounds should be set "
                       "explicitly on all reactions. "
                       "To set the missing flux bounds to default bounds "
                       "specified in cobra.Configuration use the flag "
                       "`read_sbml_model(..., set_missing_bounds=True)`.")

    return cobra_model


# -----------------------------------------------------------------------------
# Write SBML
# -----------------------------------------------------------------------------
def write_sbml_model(cobra_model, filename, f_replace=F_REPLACE, **kwargs):
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
    doc = _model_to_sbml(cobra_model, f_replace=f_replace, **kwargs)

    if isinstance(filename, string_types):
        # write to path
        libsbml.writeSBMLToFile(doc, filename)

    elif hasattr(filename, "write"):
        # write to file handle
        sbml_str = libsbml.writeSBMLToString(doc)
        filename.write(sbml_str)


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

    sbml_ns = libsbml.SBMLNamespaces(3, 1)  # SBML L3V1
    sbml_ns.addPackageNamespace("fbc", 2)  # fbc-v2

    doc = libsbml.SBMLDocument(sbml_ns)  # noqa: E501 type: libsbml.SBMLDocument
    doc.setPackageRequired("fbc", False)
    doc.setSBOTerm(SBO_FBA_FRAMEWORK)

    model = doc.createModel()  # type: libsbml.Model
    model_fbc = model.getPlugin("fbc")  # type: libsbml.FbcModelPlugin
    model_fbc.setStrict(True)

    if cobra_model.id is not None:
        model.setId(cobra_model.id)
        model.setMetaId("meta_" + cobra_model.id)
    else:
        model.setMetaId("meta_model")
    if cobra_model.name is not None:
        model.setName(cobra_model.name)

    _sbase_annotations(model, cobra_model.annotation)

    # Meta information (ModelHistory)
    if hasattr(cobra_model, "_sbml"):
        meta = cobra_model._sbml
        if "annotation" in meta:
            _sbase_annotations(doc, meta["annotation"])
        if "notes" in meta:
            _sbase_notes_dict(doc, meta["notes"])

        history = libsbml.ModelHistory()  # type: libsbml.ModelHistory
        if "created" in meta and meta["created"]:
            history.setCreatedDate(meta["created"])
        else:
            time = datetime.datetime.now()
            timestr = time.strftime('%Y-%m-%dT%H:%M:%S')
            date = libsbml.Date(timestr)
            _check(history.setCreatedDate(date), 'set creation date')
            _check(history.setModifiedDate(date), 'set modified date')

        if "creators" in meta:
            for cobra_creator in meta["creators"]:
                creator = libsbml.ModelCreator()  # noqa: E501 type: libsbml.ModelCreator
                if cobra_creator.get("familyName", None):
                    creator.setFamilyName(cobra_creator["familyName"])
                if cobra_creator.get("givenName", None):
                    creator.setGivenName(cobra_creator["givenName"])
                if cobra_creator.get("organisation", None):
                    creator.setOrganisation(cobra_creator["organisation"])
                if cobra_creator.get("email", None):
                    creator.setEmail(cobra_creator["email"])

                _check(history.addCreator(creator),
                       "adding creator to ModelHistory.")
        _check(model.setModelHistory(history), 'set model history')

    # Units
    if units:
        flux_udef = model.createUnitDefinition()  # noqa: E501 type: libsbml.UnitDefinition
        flux_udef.setId(UNITS_FLUX[0])
        for u in UNITS_FLUX[1]:
            unit = flux_udef.createUnit()  # type: libsbml.Unit
            unit.setKind(u.kind)
            unit.setExponent(u.exponent)
            unit.setScale(u.scale)
            unit.setMultiplier(u.multiplier)

    # minimum and maximum value from model
    if len(cobra_model.reactions) > 0:
        min_value = min(cobra_model.reactions.list_attr("lower_bound"))
        max_value = max(cobra_model.reactions.list_attr("upper_bound"))
    else:
        min_value = config.lower_bound
        max_value = config.upper_bound

    _create_parameter(model, pid=LOWER_BOUND_ID,
                      value=min_value, sbo=SBO_DEFAULT_FLUX_BOUND)
    _create_parameter(model, pid=UPPER_BOUND_ID,
                      value=max_value, sbo=SBO_DEFAULT_FLUX_BOUND)
    _create_parameter(model, pid=ZERO_BOUND_ID,
                      value=0, sbo=SBO_DEFAULT_FLUX_BOUND)
    _create_parameter(model, pid=BOUND_MINUS_INF,
                      value=-float("Inf"), sbo=SBO_FLUX_BOUND)
    _create_parameter(model, pid=BOUND_PLUS_INF,
                      value=float("Inf"), sbo=SBO_FLUX_BOUND)

    # Compartments
    # FIXME: use first class compartment model (and write notes & annotations)
    #     (https://github.com/opencobra/cobrapy/issues/811)
    for cid, name in iteritems(cobra_model.compartments):
        compartment = model.createCompartment()  # type: libsbml.Compartment
        compartment.setId(cid)
        compartment.setName(name)
        compartment.setConstant(True)

        # FIXME: write annotations and notes
        # _sbase_notes(c, com.notes)
        # _sbase_annotations(c, com.annotation)

    # Species
    for metabolite in cobra_model.metabolites:
        specie = model.createSpecies()  # type: libsbml.Species
        mid = metabolite.id
        if f_replace and F_SPECIE_REV in f_replace:
            mid = f_replace[F_SPECIE_REV](mid)
        specie.setId(mid)
        specie.setConstant(False)
        specie.setBoundaryCondition(False)
        specie.setHasOnlySubstanceUnits(False)
        specie.setName(metabolite.name)
        specie.setCompartment(metabolite.compartment)
        s_fbc = specie.getPlugin("fbc")  # type: libsbml.FbcSpeciesPlugin
        if metabolite.charge is not None:
            s_fbc.setCharge(metabolite.charge)
        if metabolite.formula is not None:
            s_fbc.setChemicalFormula(metabolite.formula)

        _sbase_annotations(specie, metabolite.annotation)
        _sbase_notes_dict(specie, metabolite.notes)

    # Genes
    for cobra_gene in cobra_model.genes:
        gp = model_fbc.createGeneProduct()  # type: libsbml.GeneProduct
        gid = cobra_gene.id
        if f_replace and F_GENE_REV in f_replace:
            gid = f_replace[F_GENE_REV](gid)
        gp.setId(gid)
        gname = cobra_gene.name
        if gname is None or len(gname) == 0:
            gname = gid
        gp.setName(gname)
        gp.setLabel(gid)

        _sbase_annotations(gp, cobra_gene.annotation)
        _sbase_notes_dict(gp, cobra_gene.notes)

    # Objective
    objective = model_fbc.createObjective()  # type: libsbml.Objective
    objective.setId("obj")
    objective.setType(SHORT_LONG_DIRECTION[cobra_model.objective.direction])
    model_fbc.setActiveObjectiveId("obj")

    # Reactions
    reaction_coefficients = linear_reaction_coefficients(cobra_model)
    for cobra_reaction in cobra_model.reactions:
        rid = cobra_reaction.id
        if f_replace and F_REACTION_REV in f_replace:
            rid = f_replace[F_REACTION_REV](rid)
        reaction = model.createReaction()  # type: libsbml.Reaction
        reaction.setId(rid)
        reaction.setName(cobra_reaction.name)
        reaction.setFast(False)
        reaction.setReversible((cobra_reaction.lower_bound < 0))
        _sbase_annotations(reaction, cobra_reaction.annotation)
        _sbase_notes_dict(reaction, cobra_reaction.notes)

        # stoichiometry
        for metabolite, stoichiometry in iteritems(cobra_reaction._metabolites):  # noqa: E501
            sid = metabolite.id
            if f_replace and F_SPECIE_REV in f_replace:
                sid = f_replace[F_SPECIE_REV](sid)
            if stoichiometry < 0:
                sref = reaction.createReactant()  # noqa: E501 type: libsbml.SpeciesReference
                sref.setSpecies(sid)
                sref.setStoichiometry(-stoichiometry)
                sref.setConstant(True)
            else:
                sref = reaction.createProduct()  # noqa: E501 type: libsbml.SpeciesReference
                sref.setSpecies(sid)
                sref.setStoichiometry(stoichiometry)
                sref.setConstant(True)

        # bounds
        r_fbc = reaction.getPlugin("fbc")  # type: libsbml.FbcReactionPlugin
        r_fbc.setLowerFluxBound(_create_bound(model, cobra_reaction,
                                              "lower_bound",
                                              f_replace=f_replace, units=units,
                                              flux_udef=flux_udef))
        r_fbc.setUpperFluxBound(_create_bound(model, cobra_reaction,
                                              "upper_bound",
                                              f_replace=f_replace, units=units,
                                              flux_udef=flux_udef))

        # GPR
        gpr = cobra_reaction.gene_reaction_rule
        if gpr is not None and len(gpr) > 0:

            # replace ids in string
            if f_replace and F_GENE_REV in f_replace:
                gpr = gpr.replace('(', '( ')
                gpr = gpr.replace(')', ' )')
                tokens = gpr.split(' ')
                for k in range(len(tokens)):
                    if tokens[k] not in [' ', 'and', 'or', '(', ')']:
                        tokens[k] = f_replace[F_GENE_REV](tokens[k])
                gpr_new = " ".join(tokens)

            gpa = r_fbc.createGeneProductAssociation()  # noqa: E501 type: libsbml.GeneProductAssociation
            gpa.setAssociation(gpr_new)

        # objective coefficients
        if reaction_coefficients.get(cobra_reaction, 0) != 0:
            flux_obj = objective.createFluxObjective()  # noqa: E501 type: libsbml.FluxObjective
            flux_obj.setReaction(rid)
            flux_obj.setCoefficient(cobra_reaction.objective_coefficient)

    # write groups
    if len(cobra_model.groups) > 0:
        doc.enablePackage(
            "http://www.sbml.org/sbml/level3/version1/groups/version1",
            "groups", True)
        doc.setPackageRequired("groups", False)
        model_group = model.getPlugin("groups")  # noqa: E501 type: libsbml.GroupsModelPlugin
        for cobra_group in cobra_model.groups:
            group = model_group.createGroup()  # type: libsbml.Group
            group.setId(cobra_group.id)
            group.setName(cobra_group.name)
            group.setKind(cobra_group.kind)

            _sbase_notes_dict(group, cobra_group.notes)
            _sbase_annotations(group, cobra_group.annotation)

            for cobra_member in cobra_group.members:
                member = group.createMember()  # type: libsbml.Member
                mid = cobra_member.id
                m_type = str(type(cobra_member))

                # id replacements
                if "Reaction" in m_type:
                    if f_replace and F_REACTION_REV in f_replace:
                        mid = f_replace[F_REACTION_REV](mid)
                if "Metabolite" in m_type:
                    if f_replace and F_SPECIE_REV in f_replace:
                        mid = f_replace[F_SPECIE_REV](mid)
                if "Gene" in m_type:
                    if f_replace and F_GENE_REV in f_replace:
                        mid = f_replace[F_GENE_REV](mid)

                member.setIdRef(mid)
                if cobra_member.name and len(cobra_member.name) > 0:
                    member.setName(cobra_member.name)

    return doc


def _create_bound(model, reaction, bound_type, f_replace, units=None,
                  flux_udef=None):
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
    f_replace : dict of id replacement functions
    units : flux units

    Returns
    -------
    Id of bound parameter.
    """
    value = getattr(reaction, bound_type)
    if value == config.lower_bound:
        return LOWER_BOUND_ID
    elif value == 0:
        return ZERO_BOUND_ID
    elif value == config.upper_bound:
        return UPPER_BOUND_ID
    elif value == -float("Inf"):
        return BOUND_MINUS_INF
    elif value == float("Inf"):
        return BOUND_PLUS_INF
    else:
        # new parameter
        rid = reaction.id
        if f_replace and F_REACTION_REV in f_replace:
            rid = f_replace[F_REACTION_REV](rid)
        pid = rid + "_" + bound_type
        _create_parameter(model, pid=pid, value=value, sbo=SBO_FLUX_BOUND,
                          units=units, flux_udef=flux_udef)
        return pid


def _create_parameter(model, pid, value, sbo=None, constant=True, units=None,
                      flux_udef=None):
    """Create parameter in SBML model."""
    parameter = model.createParameter()  # type: libsbml.Parameter
    parameter.setId(pid)
    parameter.setValue(value)
    parameter.setConstant(constant)
    if sbo:
        parameter.setSBOTerm(sbo)
    if units:
        parameter.setUnits(flux_udef.getId())


def _check_required(sbase, value, attribute):
    """Get required attribute from SBase.

    Parameters
    ----------
    sbase : libsbml.SBase
    value : existing value
    attribute: name of attribute

    Returns
    -------
    attribute value (or value if already set)
    """

    if (value is None) or (value == ""):
        msg = "Required attribute '%s' cannot be found or parsed in '%s'" % \
              (attribute, sbase)
        if hasattr(sbase, "getId") and sbase.getId():
            msg += " with id '%s'" % sbase.getId()
        elif hasattr(sbase, "getName") and sbase.getName():
            msg += " with name '%s'" % sbase.getName()
        elif hasattr(sbase, "getMetaId") and sbase.getMetaId():
            msg += " with metaId '%s'" % sbase.getName()
        raise CobraSBMLError(msg)
    return value


def _check(value, message):
    """
    Checks the libsbml return value and logs error messages.

    If 'value' is None, logs an error message constructed using
      'message' and then exits with status code 1. If 'value' is an integer,
      it assumes it is a libSBML return status code. If the code value is
      LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
      prints an error message constructed using 'message' along with text from
      libSBML explaining the meaning of the code, and exits with status code 1.

    """
    if value is None:
        LOGGER.error('Error: LibSBML returned a null value trying '
                     'to <' + message + '>.')
    elif type(value) is int:
        if value == libsbml.LIBSBML_OPERATION_SUCCESS:
            return
        else:
            LOGGER.error('Error encountered trying to <' + message + '>.')
            LOGGER.error('LibSBML error code {}: {}'.format(str(value),
                         libsbml.OperationReturnValue_toString(value).strip()))
    else:
        return


# -----------------------------------------------------------------------------
# Notes
# -----------------------------------------------------------------------------
def _parse_notes_dict(sbase):
    """ Creates dictionary of COBRA notes.

    Parameters
    ----------
    sbase : libsbml.SBase

    Returns
    -------
    dict of notes
    """
    notes = sbase.getNotesString()
    if notes and len(notes) > 0:
        pattern = r"<p>\s*(\w+\s*\w*)\s*:\s*([\w|\s]+)<"
        matches = re.findall(pattern, notes)
        d = {k.strip(): v.strip() for (k, v) in matches}
        return {k: v for k, v in d.items() if len(v) > 0}
    else:
        return {}


def _sbase_notes_dict(sbase, notes):
    """Set SBase notes based on dictionary.

    Parameters
    ----------
    sbase : libsbml.SBase
        SBML object to set notes on
    notes : notes object
        notes information from cobra object
    """
    if notes and len(notes) > 0:
        tokens = ['<html xmlns = "http://www.w3.org/1999/xhtml" >'] + \
            ["<p>{}: {}</p>".format(k, v) for (k, v) in notes.items()] + \
            ["</html>"]
        _check(
            sbase.setNotes("\n".join(tokens)),
            "Setting notes on sbase: {}".format(sbase)
        )


# -----------------------------------------------------------------------------
# Annotations
# -----------------------------------------------------------------------------
"""
cobra annotations will be dictionaries of the form:
    object.annotation = {
        'provider' : [(qualifier, entity), ...]
    }
A concrete example for a metabolite would look like the following
    metabolite.annotation = {
        'chebi': [(isVersionOf, "CHEBI:17234), (is, "CHEBI:4167),],
        'kegg.compound': [(is, "C00031")]
    }
The providers are hereby MIRIAM registry keys for collections
https://www.ebi.ac.uk/miriam/main/collections
The qualifiers are biomodel qualifiers
https://co.mbine.org/standards/qualifiers

In the current stage the new annotation format is not completely supported yet.
"""
URL_IDENTIFIERS_PATTERN = re.compile(r"^https?://identifiers.org/(.+?)/(.+)")  # noqa: E501
URL_IDENTIFIERS_PREFIX = "https://identifiers.org"
QUALIFIER_TYPES = {
     "is": libsbml.BQB_IS,
     "hasPart": libsbml.BQB_HAS_PART,
     "isPartOf": libsbml.BQB_IS_PART_OF,
     "isVersionOf": libsbml.BQB_IS_VERSION_OF,
     "hasVersion": libsbml.BQB_HAS_VERSION,
     "isHomologTo": libsbml.BQB_IS_HOMOLOG_TO,
     "isDescribedBy": libsbml.BQB_IS_DESCRIBED_BY,
     "isEncodedBy": libsbml.BQB_IS_ENCODED_BY,
     "encodes": libsbml.BQB_ENCODES,
     "occursIn": libsbml.BQB_OCCURS_IN,
     "hasProperty": libsbml.BQB_HAS_PROPERTY,
     "isPropertyOf": libsbml.BQB_IS_PROPERTY_OF,
     "hasTaxon": libsbml.BQB_HAS_TAXON,
     "unknown": libsbml.BQB_UNKNOWN,
     "bqm_is": libsbml.BQM_IS,
     "bqm_isDescribedBy": libsbml.BQM_IS_DESCRIBED_BY,
     "bqm_isDerivedFrom": libsbml.BQM_IS_DERIVED_FROM,
     "bqm_isInstanceOf": libsbml.BQM_IS_INSTANCE_OF,
     "bqm_hasInstance": libsbml.BQM_HAS_INSTANCE,
     "bqm_unknown": libsbml.BQM_UNKNOWN,
}


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

    FIXME: annotation format must be updated (this is a big collection of
          fixes) - see: https://github.com/opencobra/cobrapy/issues/684)
    """
    annotation = {}

    # SBO term
    if sbase.isSetSBOTerm():
        # FIXME: correct handling of annotations
        annotation["sbo"] = sbase.getSBOTermID()

    # RDF annotation
    cvterms = sbase.getCVTerms()
    if cvterms is None:
        return annotation

    for cvterm in cvterms:  # type: libsbml.CVTerm
        for k in range(cvterm.getNumResources()):
            # FIXME: read and store the qualifier

            uri = cvterm.getResourceURI(k)
            match = URL_IDENTIFIERS_PATTERN.match(uri)
            if not match:
                LOGGER.warning("%s does not conform to "
                               "http(s)://identifiers.org/collection/id", uri)
                continue

            provider, identifier = match.group(1), match.group(2)
            if provider in annotation:
                if isinstance(annotation[provider], string_types):
                    annotation[provider] = [annotation[provider]]
                # FIXME: use a list
                if identifier not in annotation[provider]:
                    annotation[provider].append(identifier)
            else:
                # FIXME: always in list
                annotation[provider] = identifier

    return annotation


def _sbase_annotations(sbase, annotation):
    """Set SBase annotations based on cobra annotations.

    Parameters
    ----------
    sbase : libsbml.SBase
        SBML object to annotate
    annotation : cobra annotation structure
        cobra object with annotation information

    FIXME: annotation format must be updated
        (https://github.com/opencobra/cobrapy/issues/684)
    """

    if not annotation or len(annotation) == 0:
        return

    # standardize annotations
    annotation_data = deepcopy(annotation)

    for key, value in annotation_data.items():
        # handling of non-string annotations (e.g. integers)
        if isinstance(value, (float, int)):
            value = str(value)
        if isinstance(value, string_types):
            annotation_data[key] = [("is", value)]

    for key, value in annotation_data.items():
        for idx, item in enumerate(value):
            if isinstance(item, string_types):
                value[idx] = ("is", item)

    # set metaId
    meta_id = "meta_{}".format(sbase.getId())
    sbase.setMetaId(meta_id)

    # rdf_items = []
    for provider, data in iteritems(annotation_data):

        # set SBOTerm
        if provider in ["SBO", "sbo"]:
            if provider == "SBO":
                LOGGER.warning("'SBO' provider is deprecated, "
                               "use 'sbo' provider instead")
            sbo_term = data[0][1]
            _check(sbase.setSBOTerm(sbo_term),
                   "Setting SBOTerm: {}".format(sbo_term))

            # FIXME: sbo should also be written as CVTerm
            continue

        for item in data:
            qualifier_str, entity = item[0], item[1]
            qualifier = QUALIFIER_TYPES.get(qualifier_str, None)
            if qualifier is None:
                qualifier = libsbml.BQB_IS
                LOGGER.error("Qualifier type is not supported on "
                             "annotation: '{}'".format(qualifier_str))

            qualifier_type = libsbml.BIOLOGICAL_QUALIFIER
            if qualifier_str.startswith("bqm_"):
                qualifier_type = libsbml.MODEL_QUALIFIER

            cv = libsbml.CVTerm()  # type: libsbml.CVTerm
            cv.setQualifierType(qualifier_type)
            if qualifier_type == libsbml.BIOLOGICAL_QUALIFIER:
                cv.setBiologicalQualifierType(qualifier)
            elif qualifier_type == libsbml.MODEL_QUALIFIER:
                cv.setModelQualifierType(qualifier)
            else:
                raise CobraSBMLError('Unsupported qualifier: '
                                     '%s' % qualifier)
            resource = "%s/%s/%s" % (URL_IDENTIFIERS_PREFIX, provider, entity)
            cv.addResource(resource)
            _check(sbase.addCVTerm(cv),
                   "Setting cvterm: {}, resource: {}".format(cv, resource))


# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
def validate_sbml_model(filename,
                        check_model=True,
                        internal_consistency=True,
                        check_units_consistency=False,
                        check_modeling_practice=False, **kwargs):
    """Validate SBML model and returns the model along with a list of errors.

    Parameters
    ----------
    filename : str
        The filename (or SBML string) of the SBML model to be validated.
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
    (model, errors)
    model : :class:`~cobra.core.Model.Model` object
        The cobra model if the file could be read successfully or None
        otherwise.
    errors : dict
        Warnings and errors grouped by their respective types.

    Raises
    ------
    CobraSBMLError
    """
    # Errors and warnings are grouped based on their type. SBML_* types are
    # from the libsbml validator. COBRA_* types are from the cobrapy SBML
    # parser.
    keys = (
        "SBML_FATAL",
        "SBML_ERROR",
        "SBML_SCHEMA_ERROR",
        "SBML_WARNING",

        "COBRA_FATAL",
        "COBRA_ERROR",
        "COBRA_WARNING",
        "COBRA_CHECK",
    )
    errors = {key: [] for key in keys}

    # [1] libsbml validation
    doc = _get_doc_from_filename(filename)  # type: libsbml.SBMLDocument

    # set checking of units & modeling practise
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY,
                             check_units_consistency)
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_MODELING_PRACTICE,
                             check_modeling_practice)

    # check internal consistency
    if internal_consistency:
        doc.checkInternalConsistency()
    doc.checkConsistency()

    for k in range(doc.getNumErrors()):
        e = doc.getError(k)  # type: libsbml.SBMLError
        msg = _error_string(e, k=k)
        sev = e.getSeverity()
        if sev == libsbml.LIBSBML_SEV_FATAL:
            errors["SBML_FATAL"].append(msg)
        elif sev == libsbml.LIBSBML_SEV_ERROR:
            errors["SBML_ERROR"].append(msg)
        elif sev == libsbml.LIBSBML_SEV_SCHEMA_ERROR:
            errors["SBML_SCHEMA_ERROR"].append(msg)
        elif sev == libsbml.LIBSBML_SEV_WARNING:
            errors["SBML_WARNING"].append(msg)

    # [2] cobrapy validation (check that SBML can be read into model)
    # all warnings generated while loading will be logged as errors
    log_stream = StringIO()
    stream_handler = logging.StreamHandler(log_stream)
    formatter = logging.Formatter('%(levelname)s:%(message)s')
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)
    LOGGER.addHandler(stream_handler)
    LOGGER.propagate = False

    try:
        # read model and allow additional parser arguments
        model = _sbml_to_model(doc, **kwargs)
    except CobraSBMLError as e:
        errors["COBRA_ERROR"].append(str(e))
        return None, errors
    except Exception as e:
        errors["COBRA_FATAL"].append(str(e))
        return None, errors

    cobra_errors = log_stream.getvalue().split("\n")
    for cobra_error in cobra_errors:
        tokens = cobra_error.split(":")
        error_type = tokens[0]
        error_msg = ":".join(tokens[1:])

        if error_type == "WARNING":
            errors["COBRA_WARNING"].append(error_msg)
        elif error_type == "ERROR":
            errors["COBRA_ERROR"].append(error_msg)

    # remove stream handler
    LOGGER.removeHandler(stream_handler)
    LOGGER.propagate = True

    # [3] additional model tests
    if check_model:
        errors["COBRA_CHECK"].extend(
            check_metabolite_compartment_formula(model)
        )

    for key in ["SBML_FATAL", "SBML_ERROR", "SBML_SCHEMA_ERROR"]:
        if len(errors[key]) > 0:
            LOGGER.error("SBML errors in validation, check error log "
                         "for details.")
            break
    for key in ["SBML_WARNING"]:
        if len(errors[key]) > 0:
            LOGGER.error("SBML warnings in validation, check error log "
                         "for details.")
            break
    for key in ["COBRA_FATAL", "COBRA_ERROR"]:
        if len(errors[key]) > 0:
            LOGGER.error("COBRA errors in validation, check error log "
                         "for details.")
            break
    for key in ["COBRA_WARNING", "COBRA_CHECK"]:
        if len(errors[key]) > 0:
            LOGGER.error("COBRA warnings in validation, check error log "
                         "for details.")
            break

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

    template = 'E{} ({}): {} ({}, L{}); {}; {}'
    error_str = template.format(k, error.getSeverityAsString(),
                                error.getCategoryAsString(), package,
                                error.getLine(), error.getShortMessage(),
                                error.getMessage())
    return error_str
