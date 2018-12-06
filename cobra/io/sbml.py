# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
from math import isinf, isnan
from os.path import isfile
from warnings import warn

from six import iteritems

from cobra.core import Configuration, Metabolite, Model, Reaction
from cobra.util.solver import set_objective


try:
    import libsbml
except ImportError:
    libsbml = None


CONFIGURATION = Configuration()


def parse_legacy_id(the_id, the_compartment=None, the_type='metabolite',
                    use_hyphens=False):
    """Deals with a bunch of problems due to bigg.ucsd.edu not following SBML
    standards

    Parameters
    ----------
    the_id: String.
    the_compartment: String
    the_type: String
        Currently only 'metabolite' is supported
    use_hyphens: Boolean
        If True, double underscores (__) in an SBML ID will be converted to
        hyphens

    Returns
    -------
    string: the identifier
    """
    if use_hyphens:
        the_id = the_id.replace('__', '-')
    if the_type == 'metabolite':
        if the_id.split('_')[-1] == the_compartment:
            # Reformat Ids to match convention in Palsson Lab.
            the_id = the_id[:-len(the_compartment) - 1]
        the_id += '[%s]' % the_compartment
    return the_id


def create_cobra_model_from_sbml_file(sbml_filename, old_sbml=False,
                                      legacy_metabolite=False,
                                      print_time=False, use_hyphens=False):
    """convert an SBML XML file into a cobra.Model object.

    Supports SBML Level 2 Versions 1 and 4.  The function will detect if the
    SBML fbc package is used in the file and run the converter if the fbc
    package is used.

    Parameters
    ----------
    sbml_filename: string
    old_sbml: bool
        Set to True if the XML file has metabolite formula appended to
        metabolite names. This was a poorly designed artifact that persists in
        some models.
    legacy_metabolite: bool
        If True then assume that the metabolite id has the compartment id
         appended after an underscore (e.g. _c for cytosol). This has not been
         implemented but will be soon.
    print_time: bool
         deprecated
    use_hyphens: bool
        If True, double underscores (__) in an SBML ID will be converted to
        hyphens

    Returns
    -------
    Model : The parsed cobra model
    """
    if not libsbml:
        raise ImportError('create_cobra_model_from_sbml_file '
                          'requires python-libsbml')

    __default_objective_coefficient = 0
    # Ensure that the file exists
    if not isfile(sbml_filename):
        raise IOError('Your SBML file is not found: %s' % sbml_filename)
    # Expressions to change SBML Ids to Palsson Lab Ids
    metabolite_re = re.compile('^M_')
    reaction_re = re.compile('^R_')
    compartment_re = re.compile('^C_')
    if print_time:
        warn("print_time is deprecated", DeprecationWarning)
    model_doc = libsbml.readSBML(sbml_filename)
    if model_doc.getPlugin("fbc") is not None:
        from libsbml import ConversionProperties, LIBSBML_OPERATION_SUCCESS
        conversion_properties = ConversionProperties()
        conversion_properties.addOption(
            "convert fbc to cobra", True, "Convert FBC model to Cobra model")
        result = model_doc.convert(conversion_properties)
        if result != LIBSBML_OPERATION_SUCCESS:
            raise Exception("Conversion of SBML+fbc to COBRA failed")
    sbml_model = model_doc.getModel()
    sbml_model_id = sbml_model.getId()
    sbml_species = sbml_model.getListOfSpecies()
    sbml_reactions = sbml_model.getListOfReactions()
    sbml_compartments = sbml_model.getListOfCompartments()
    compartment_dict = dict([(compartment_re.split(x.getId())[-1], x.getName())
                             for x in sbml_compartments])
    if legacy_metabolite:
        # Deal with the palsson lab appending the compartment id to the
        # metabolite id
        new_dict = {}
        for the_id, the_name in compartment_dict.items():
            if the_name == '':
                new_dict[the_id[0].lower()] = the_id
            else:
                new_dict[the_id] = the_name
        compartment_dict = new_dict
        legacy_compartment_converter = dict(
            [(v, k) for k, v in iteritems(compartment_dict)])

    cobra_model = Model(sbml_model_id)
    metabolites = []
    metabolite_dict = {}
    # Convert sbml_metabolites to cobra.Metabolites
    for sbml_metabolite in sbml_species:
        # Skip sbml boundary species
        if sbml_metabolite.getBoundaryCondition():
            continue

        if (old_sbml or legacy_metabolite) and \
                sbml_metabolite.getId().endswith('_b'):
            # Deal with incorrect sbml from bigg.ucsd.edu
            continue
        tmp_metabolite = Metabolite()
        metabolite_id = tmp_metabolite.id = sbml_metabolite.getId()
        tmp_metabolite.compartment = compartment_re.split(
            sbml_metabolite.getCompartment())[-1]
        if legacy_metabolite:
            if tmp_metabolite.compartment not in compartment_dict:
                tmp_metabolite.compartment = legacy_compartment_converter[
                    tmp_metabolite.compartment]
            tmp_metabolite.id = parse_legacy_id(
                tmp_metabolite.id, tmp_metabolite.compartment,
                use_hyphens=use_hyphens)
        if use_hyphens:
            tmp_metabolite.id = metabolite_re.split(
                tmp_metabolite.id)[-1].replace('__', '-')
        else:
            # Just in case the SBML ids are ill-formed and use -
            tmp_metabolite.id = metabolite_re.split(
                tmp_metabolite.id)[-1].replace('-', '__')
        tmp_metabolite.name = sbml_metabolite.getName()
        tmp_formula = ''
        tmp_metabolite.notes = parse_legacy_sbml_notes(
            sbml_metabolite.getNotesString())
        if sbml_metabolite.isSetCharge():
            tmp_metabolite.charge = sbml_metabolite.getCharge()
        if "CHARGE" in tmp_metabolite.notes:
            note_charge = tmp_metabolite.notes["CHARGE"][0]
            try:
                note_charge = float(note_charge)
                if note_charge == int(note_charge):
                    note_charge = int(note_charge)
            except:
                warn("charge of %s is not a number (%s)" %
                     (tmp_metabolite.id, str(note_charge)))
            else:
                if ((tmp_metabolite.charge is None) or
                        (tmp_metabolite.charge == note_charge)):
                    tmp_metabolite.notes.pop("CHARGE")
                    # set charge to the one from notes if not assigend before
                    # the same
                    tmp_metabolite.charge = note_charge
                else:  # tmp_metabolite.charge != note_charge
                    msg = "different charges specified for %s (%d and %d)"
                    msg = msg % (tmp_metabolite.id,
                                 tmp_metabolite.charge, note_charge)
                    warn(msg)
                    # Chances are a 0 note charge was written by mistake. We
                    # will default to the note_charge in this case.
                    if tmp_metabolite.charge == 0:
                        tmp_metabolite.charge = note_charge

        for the_key in tmp_metabolite.notes.keys():
            if the_key.lower() == 'formula':
                tmp_formula = tmp_metabolite.notes.pop(the_key)[0]
                break
        if tmp_formula == '' and old_sbml:
            tmp_formula = tmp_metabolite.name.split('_')[-1]
            tmp_metabolite.name = tmp_metabolite.name[:-len(tmp_formula) - 1]
        tmp_metabolite.formula = tmp_formula
        metabolite_dict.update({metabolite_id: tmp_metabolite})
        metabolites.append(tmp_metabolite)
    cobra_model.add_metabolites(metabolites)

    # Construct the vectors and matrices for holding connectivity and numerical
    # info to feed to the cobra toolbox.
    # Always assume steady state simulations so b is set to 0
    cobra_reaction_list = []
    coefficients = {}
    for sbml_reaction in sbml_reactions:
        if use_hyphens:
            # Change the ids to match conventions used by the Palsson lab.
            reaction = Reaction(reaction_re.split(
                sbml_reaction.getId())[-1].replace('__', '-'))
        else:
            # Just in case the SBML ids are ill-formed and use -
            reaction = Reaction(reaction_re.split(
                sbml_reaction.getId())[-1].replace('-', '__'))
        cobra_reaction_list.append(reaction)
        # reaction.exchange_reaction = 0
        reaction.name = sbml_reaction.getName()
        cobra_metabolites = {}
        # Use the cobra.Metabolite class here
        for sbml_metabolite in sbml_reaction.getListOfReactants():
            tmp_metabolite_id = sbml_metabolite.getSpecies()
            # This deals with boundary metabolites
            if tmp_metabolite_id in metabolite_dict:
                tmp_metabolite = metabolite_dict[tmp_metabolite_id]
                cobra_metabolites[tmp_metabolite] = - \
                    sbml_metabolite.getStoichiometry()
        for sbml_metabolite in sbml_reaction.getListOfProducts():
            tmp_metabolite_id = sbml_metabolite.getSpecies()
            # This deals with boundary metabolites
            if tmp_metabolite_id in metabolite_dict:
                tmp_metabolite = metabolite_dict[tmp_metabolite_id]
                # Handle the case where the metabolite was specified both
                # as a reactant and as a product.
                if tmp_metabolite in cobra_metabolites:
                    warn("%s appears as a reactant and product %s" %
                         (tmp_metabolite_id, reaction.id))
                    cobra_metabolites[
                        tmp_metabolite] += sbml_metabolite.getStoichiometry()
                    # if the combined stoichiometry is 0, remove the metabolite
                    if cobra_metabolites[tmp_metabolite] == 0:
                        cobra_metabolites.pop(tmp_metabolite)
                else:
                    cobra_metabolites[
                        tmp_metabolite] = sbml_metabolite.getStoichiometry()
        # check for nan
        for met, v in iteritems(cobra_metabolites):
            if isnan(v) or isinf(v):
                warn("invalid value %s for metabolite '%s' in reaction '%s'" %
                     (str(v), met.id, reaction.id))
        reaction.add_metabolites(cobra_metabolites)
        # Parse the kinetic law info here.
        parameter_dict = {}
        # If lower and upper bounds are specified in the Kinetic Law then
        # they override the sbml reversible attribute.  If they are not
        # specified then the bounds are determined by getReversible.
        if not sbml_reaction.getKineticLaw():

            if sbml_reaction.getReversible():
                parameter_dict['lower_bound'] = CONFIGURATION.lower_bound
                parameter_dict['upper_bound'] = CONFIGURATION.upper_bound
            else:
                # Assume that irreversible reactions only proceed from left to
                # right.
                parameter_dict['lower_bound'] = 0
                parameter_dict['upper_bound'] = CONFIGURATION.upper_bound

            parameter_dict[
                'objective_coefficient'] = __default_objective_coefficient
        else:
            for sbml_parameter in \
                    sbml_reaction.getKineticLaw().getListOfParameters():
                parameter_dict[
                    sbml_parameter.getId().lower()] = sbml_parameter.getValue()

        if 'lower_bound' in parameter_dict:
            reaction.lower_bound = parameter_dict['lower_bound']
        elif 'lower bound' in parameter_dict:
            reaction.lower_bound = parameter_dict['lower bound']
        elif sbml_reaction.getReversible():
            reaction.lower_bound = CONFIGURATION.lower_bound
        else:
            reaction.lower_bound = 0

        if 'upper_bound' in parameter_dict:
            reaction.upper_bound = parameter_dict['upper_bound']
        elif 'upper bound' in parameter_dict:
            reaction.upper_bound = parameter_dict['upper bound']
        else:
            reaction.upper_bound = CONFIGURATION.upper_bound

        objective_coefficient = parameter_dict.get(
            'objective_coefficient', parameter_dict.get(
                'objective_coefficient', __default_objective_coefficient))
        if objective_coefficient != 0:
            coefficients[reaction] = objective_coefficient

        # ensure values are not set to nan or inf
        if isnan(reaction.lower_bound) or isinf(reaction.lower_bound):
            reaction.lower_bound = CONFIGURATION.lower_bound
        if isnan(reaction.upper_bound) or isinf(reaction.upper_bound):
            reaction.upper_bound = CONFIGURATION.upper_bound

        reaction_note_dict = parse_legacy_sbml_notes(
            sbml_reaction.getNotesString())
        # Parse the reaction notes.
        # POTENTIAL BUG: DEALING WITH LEGACY 'SBML' THAT IS NOT IN A
        # STANDARD FORMAT
        # TODO: READ IN OTHER NOTES AND GIVE THEM A reaction_ prefix.
        # TODO: Make sure genes get added as objects
        if 'GENE ASSOCIATION' in reaction_note_dict:
            rule = reaction_note_dict['GENE ASSOCIATION'][0]
            try:
                rule.encode('ascii')
            except (UnicodeEncodeError, UnicodeDecodeError):
                warn("gene_reaction_rule '%s' is not ascii compliant" % rule)
            if rule.startswith("&quot;") and rule.endswith("&quot;"):
                rule = rule[6:-6]
            reaction.gene_reaction_rule = rule
            if 'GENE LIST' in reaction_note_dict:
                reaction.systematic_names = reaction_note_dict['GENE LIST'][0]
            elif ('GENES' in reaction_note_dict and
                  reaction_note_dict['GENES'] != ['']):
                reaction.systematic_names = reaction_note_dict['GENES'][0]
            elif 'LOCUS' in reaction_note_dict:
                gene_id_to_object = dict([(x.id, x) for x in reaction._genes])
                for the_row in reaction_note_dict['LOCUS']:
                    tmp_row_dict = {}
                    the_row = 'LOCUS:' + the_row.lstrip('_').rstrip('#')
                    for the_item in the_row.split('#'):
                        k, v = the_item.split(':')
                        tmp_row_dict[k] = v
                    tmp_locus_id = tmp_row_dict['LOCUS']
                    if 'TRANSCRIPT' in tmp_row_dict:
                        tmp_locus_id = tmp_locus_id + \
                                       '.' + tmp_row_dict['TRANSCRIPT']

                    if 'ABBREVIATION' in tmp_row_dict:
                        gene_id_to_object[tmp_locus_id].name = tmp_row_dict[
                            'ABBREVIATION']

        if 'SUBSYSTEM' in reaction_note_dict:
            reaction.subsystem = reaction_note_dict.pop('SUBSYSTEM')[0]

        reaction.notes = reaction_note_dict

    # Now, add all of the reactions to the model.
    cobra_model.id = sbml_model.getId()
    # Populate the compartment list - This will be done based on
    # cobra.Metabolites in cobra.Reactions in the future.
    cobra_model.compartments = compartment_dict

    cobra_model.add_reactions(cobra_reaction_list)
    set_objective(cobra_model, coefficients)
    return cobra_model


def parse_legacy_sbml_notes(note_string, note_delimiter=':'):
    """Deal with various legacy SBML format issues.
    """
    note_dict = {}
    start_tag = '<p>'
    end_tag = '</p>'
    if '<html:p>' in note_string:
        start_tag = '<html:p>'
        end_tag = '</html:p>'
    while start_tag in note_string and end_tag in note_string:
        note_start = note_string.index(start_tag)
        note_end = note_string.index(end_tag)
        the_note = note_string[
                   (note_start + len(start_tag)):note_end].lstrip(' ').rstrip(
            ' ')
        if note_delimiter in the_note:
            note_delimiter_index = the_note.index(note_delimiter)
            note_field = the_note[:note_delimiter_index].lstrip(
                ' ').rstrip(' ').replace('_', ' ').upper()
            note_value = the_note[
                         (note_delimiter_index + 1):].lstrip(' ').rstrip(' ')
            if note_field in note_dict:
                note_dict[note_field].append(note_value)
            else:
                note_dict[note_field] = [note_value]
        note_string = note_string[(note_end + len(end_tag)):]

    if ('CHARGE' in note_dict and
            note_dict['CHARGE'][0].lower() in ['none', 'na', 'nan']):
        note_dict.pop('CHARGE')  # Remove non-numeric charges

    if 'CHARGE' in note_dict and note_dict['CHARGE'][0].lower() in ['none',
                                                                    'na',
                                                                    'nan']:
        note_dict.pop('CHARGE')  # Remove non-numeric charges

    return note_dict


def write_cobra_model_to_sbml_file(cobra_model, sbml_filename,
                                   sbml_level=2, sbml_version=1,
                                   print_time=False,
                                   use_fbc_package=True):
    """Write a cobra.Model object to an SBML XML file.

    Parameters
    ----------
    cobra_model : cobra.core.Model.Model
        The model object to write
    sbml_filename : string
        The file to write the SBML XML to.
    sbml_level : int
        2 is the only supported level.
    sbml_version : int
        1 is the only supported version.
    print_time : bool
        deprecated
    use_fbc_package : bool
        Convert the model to the FBC package format to improve portability.
        http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/Flux_Balance_Constraints_(flux)

    Notes
    -----
    TODO: Update the NOTES to match the SBML standard and provide support for
    Level 2 Version 4
    """
    if not libsbml:
        raise ImportError('write_cobra_model_to_sbml_file '
                          'requires python-libsbml')
    sbml_doc = get_libsbml_document(cobra_model,
                                    sbml_level=sbml_level,
                                    sbml_version=sbml_version,
                                    print_time=print_time,
                                    use_fbc_package=use_fbc_package)

    libsbml.writeSBML(sbml_doc, sbml_filename)


def get_libsbml_document(cobra_model,
                         sbml_level=2, sbml_version=1,
                         print_time=False,
                         use_fbc_package=True):
    """ Return a libsbml document object for writing to a file. This function
    is used by write_cobra_model_to_sbml_file(). """
    note_start_tag, note_end_tag = '<p>', '</p>'
    if sbml_level > 2 or (sbml_level == 2 and sbml_version == 4):
        note_start_tag, note_end_tag = '<html:p>', '</html:p>'

    sbml_doc = libsbml.SBMLDocument(sbml_level, sbml_version)
    sbml_model = sbml_doc.createModel(cobra_model.id.split('.')[0])
    # Note need to set units
    reaction_units = 'mmol_per_gDW_per_hr'
    model_units = sbml_model.createUnitDefinition()
    model_units.setId(reaction_units)
    sbml_unit = model_units.createUnit()
    sbml_unit.setKind(libsbml.UNIT_KIND_MOLE)
    sbml_unit.setScale(-3)
    sbml_unit = model_units.createUnit()
    sbml_unit.setKind(libsbml.UNIT_KIND_GRAM)
    sbml_unit.setExponent(-1)
    sbml_unit = model_units.createUnit()
    sbml_unit.setKind(libsbml.UNIT_KIND_SECOND)
    sbml_unit.setMultiplier(1.0 / 60 / 60)
    sbml_unit.setExponent(-1)

    # Add in the common compartment abbreviations.  If there are additional
    # compartments they also need to be added.
    if not cobra_model.compartments:
        cobra_model.compartments = {'c': 'cytosol',
                                    'p': 'periplasm',
                                    'e': 'extracellular'}
    for the_key in cobra_model.compartments.keys():
        sbml_comp = sbml_model.createCompartment()
        sbml_comp.setId(the_key)
        sbml_comp.setName(cobra_model.compartments[the_key])
        sbml_comp.setSize(1)  # Just to get rid of warnings

    if print_time:
        warn("print_time is deprecated", DeprecationWarning)
    # Use this dict to allow for fast look up of species id
    # for references created in the reaction section.
    metabolite_dict = {}

    for cobra_metabolite in cobra_model.metabolites:
        metabolite_dict[cobra_metabolite.id] = add_sbml_species(
            sbml_model, cobra_metabolite, note_start_tag=note_start_tag,
            note_end_tag=note_end_tag)

    for the_reaction in cobra_model.reactions:
        # This is probably the culprit.  Including cobra.Reaction
        # objects explicitly in cobra.Model will speed this up.
        sbml_reaction = sbml_model.createReaction()
        # Need to remove - for proper SBML.  Replace with __
        the_reaction_id = 'R_' + the_reaction.id.replace('-', '__')
        sbml_reaction.setId(the_reaction_id)
        # The reason we are not using the Reaction.reversibility property
        # is because the SBML definition of reversibility does not quite
        # match with the cobra definition. In cobra, reversibility implies
        # that both positive and negative flux values are feasible. However,
        # SBML requires negative-flux-only reactions to still be classified
        # as reversible. To quote from the SBML Level 3 Version 1 Spec:
        # > However, labeling a reaction as irreversible is interpreted as
        # > an assertion that the rate expression will not have negative
        # > values during a simulation.
        # (Page 60 lines 44-45)
        sbml_reaction.setReversible(the_reaction.lower_bound < 0)
        if the_reaction.name:
            sbml_reaction.setName(the_reaction.name)
        else:
            sbml_reaction.setName(the_reaction.id)
        # Add in the reactant/product references
        for the_metabolite, the_coefficient in \
                iteritems(the_reaction._metabolites):
            sbml_stoichiometry = the_coefficient
            metabolite_id = str(metabolite_dict[the_metabolite.id])
            # Each SpeciesReference must have a unique id
            if sbml_stoichiometry < 0:
                species_reference = sbml_reaction.createReactant()
            else:
                species_reference = sbml_reaction.createProduct()
            species_reference.setId(metabolite_id + '_' + the_reaction_id)
            species_reference.setSpecies(metabolite_id)
            species_reference.setStoichiometry(abs(sbml_stoichiometry))
        # Deal with the case where the reaction is a boundary reaction
        if len(the_reaction._metabolites) == 1:
            the_metabolite, the_coefficient = list(
                the_reaction._metabolites.items())[0]
            metabolite_id = add_sbml_species(sbml_model, the_metabolite,
                                             note_start_tag=note_start_tag,
                                             note_end_tag=note_end_tag,
                                             boundary_metabolite=True)
            sbml_stoichiometry = -the_coefficient
            # Each SpeciesReference must have a unique id
            if sbml_stoichiometry < 0:
                species_reference = sbml_reaction.createReactant()
            else:
                species_reference = sbml_reaction.createProduct()
            species_reference.setId(metabolite_id + '_' + the_reaction_id)
            species_reference.setSpecies(metabolite_id)
            species_reference.setStoichiometry(abs(sbml_stoichiometry))

        # Add in the kineticLaw
        sbml_law = libsbml.KineticLaw(sbml_level, sbml_version)
        if hasattr(sbml_law, 'setId'):
            sbml_law.setId('FLUX_VALUE')
        sbml_law.setFormula('FLUX_VALUE')
        reaction_parameter_dict = {
            'LOWER_BOUND': [the_reaction.lower_bound, reaction_units],
            'UPPER_BOUND': [the_reaction.upper_bound, reaction_units],
            'FLUX_VALUE': [0, reaction_units],
            'OBJECTIVE_COEFFICIENT': [the_reaction.objective_coefficient,
                                      'dimensionless']}
        for k, v in reaction_parameter_dict.items():
            sbml_parameter = libsbml.Parameter(sbml_level, sbml_version)
            sbml_parameter.setId(k)
            if hasattr(v, '__iter__'):
                sbml_parameter.setValue(v[0])
                sbml_parameter.setUnits(v[1])
            else:
                sbml_parameter.setValue(v)
            sbml_law.addParameter(sbml_parameter)
        sbml_reaction.setKineticLaw(sbml_law)

        # Checks if GPR and Subsystem annotations are present in the notes
        # section and if they are the same as those in the reaction's
        # gene_reaction_rule/ subsystem attribute. If they are not identical,
        # they are set to be identical
        note_dict = the_reaction.notes.copy()
        if the_reaction.gene_reaction_rule:
            note_dict['GENE ASSOCIATION'] = [
                str(the_reaction.gene_reaction_rule)]
        if the_reaction.subsystem:
            note_dict['SUBSYSTEM'] = [str(the_reaction.subsystem)]

        # In a cobrapy model the notes section is stored as a dictionary. The
        # following section turns the key-value-pairs of the dictionary into a
        # string and replaces recurring symbols so that the string has the i
        # required syntax for an SBML doc.
        note_str = str(list(iteritems(note_dict)))
        note_start_tag, note_end_tag, note_delimiter = '<p>', '</p>', ':'
        note_str = note_str.replace('(\'', note_start_tag)
        note_str = note_str.replace('\']),', note_end_tag)
        note_str = note_str.replace('\',', note_delimiter)
        note_str = note_str.replace('\']', '')
        note_str = note_str.replace('[\'', '')
        note_str = note_str.replace(
            '[', '<html xmlns="http://www.w3.org/1999/xhtml">')
        note_str = note_str.replace(')]', note_end_tag + '</html>')
        sbml_reaction.setNotes(note_str)

    if use_fbc_package:
        try:
            from libsbml import ConversionProperties, LIBSBML_OPERATION_SUCCESS
            conversion_properties = ConversionProperties()
            conversion_properties.addOption(
                "convert cobra", True, "Convert Cobra model")
            result = sbml_doc.convert(conversion_properties)
            if result != LIBSBML_OPERATION_SUCCESS:
                raise Exception("Conversion of COBRA to SBML+fbc failed")
        except Exception as e:
            error_string = 'Error saving as SBML+fbc. %s'
            try:
                # Check whether the FbcExtension is there
                from libsbml import FbcExtension
                error_string = error_string % e
            except ImportError:
                error_string = error_string % "FbcExtension not available in "
                "libsbml. If use_fbc_package == True then libsbml must be "
                "compiled with the fbc extension."
                from libsbml import getLibSBMLDottedVersion
                _sbml_version = getLibSBMLDottedVersion()
                _major, _minor, _patch = map(int, _sbml_version.split('.'))
                if _major < 5 or (_major == 5 and _minor < 8):
                    error_string += "You've got libsbml %s installed. "
                    "You need 5.8.0 or later with the fbc package"

            raise Exception(error_string)
    return sbml_doc


def add_sbml_species(sbml_model, cobra_metabolite, note_start_tag,
                     note_end_tag, boundary_metabolite=False):
    """A helper function for adding cobra metabolites to an sbml model.

    Parameters
    ----------
    sbml_model: sbml_model object

    cobra_metabolite: a cobra.Metabolite object

    note_start_tag: string
       the start tag for parsing cobra notes. this will eventually
       be supplanted when COBRA is worked into sbml.

    note_end_tag: string
       the end tag for parsing cobra notes. this will eventually
       be supplanted when COBRA is worked into sbml.
    boundary_metabolite: bool
       if metabolite boundary condition should be set or not

    Returns
    -------
    string: the created metabolite identifier
    """
    sbml_species = sbml_model.createSpecies()
    the_id = 'M_' + cobra_metabolite.id.replace('-', '__')
    # Deal with legacy naming issues
    the_compartment = cobra_metabolite.compartment
    if the_id.endswith('[%s]' % the_compartment):
        the_id = the_id[:-len('[%s]' % the_compartment)]
    elif not the_id.endswith('_%s' % the_compartment):
        the_id += '_%s' % the_compartment
    if boundary_metabolite:
        the_id += '_boundary'
    sbml_species.setId(the_id)
    metabolite_id = the_id
    if boundary_metabolite:
        sbml_species.setBoundaryCondition(True)
    if cobra_metabolite.name:
        sbml_species.setName(cobra_metabolite.name)
    else:
        sbml_species.setName(cobra_metabolite.id)
    if the_compartment is not None:
        try:
            sbml_species.setCompartment(the_compartment)
        except:
            warn('metabolite failed: ' + the_id)
            return cobra_metabolite
    if cobra_metabolite.charge is not None:
        sbml_species.setCharge(cobra_metabolite.charge)
    # Deal with cases where the formula in the model is not an object but as a
    # string
    if cobra_metabolite.formula or cobra_metabolite.notes:
        tmp_note = '<html xmlns="http://www.w3.org/1999/xhtml">'
        if hasattr(cobra_metabolite.formula, 'id'):
            tmp_note += '%sFORMULA: %s%s' % (note_start_tag,
                                             cobra_metabolite.formula.id,
                                             note_end_tag)
        else:
            tmp_note += '%sFORMULA: %s%s' % (note_start_tag,
                                             cobra_metabolite.formula,
                                             note_end_tag)
        if hasattr(cobra_metabolite.notes, 'items'):
            for the_id_type, the_id in cobra_metabolite.notes.items():
                if the_id_type.lower() == 'charge':
                    # Use of notes['CHARGE'] has been deprecated in favor of
                    # metabolite.charge
                    continue
                if not isinstance(the_id_type, str):
                    the_id_type = str(the_id_type)
                if hasattr(the_id, '__iter__') and len(the_id) == 1:
                    the_id = the_id[0]
                if not isinstance(the_id, str):
                    the_id = str(the_id)
                tmp_note += '%s%s: %s%s' % (note_start_tag,
                                            the_id_type,
                                            the_id, note_end_tag)
        sbml_species.setNotes(tmp_note + '</html>')
    return metabolite_id


def fix_legacy_id(id, use_hyphens=False, fix_compartments=False):
    id = id.replace('_DASH_', '__')
    id = id.replace('_FSLASH_', '/')
    id = id.replace('_BSLASH_', "\\")
    id = id.replace('_LPAREN_', '(')
    id = id.replace('_LSQBKT_', '[')
    id = id.replace('_RSQBKT_', ']')
    id = id.replace('_RPAREN_', ')')
    id = id.replace('_COMMA_', ',')
    id = id.replace('_PERIOD_', '.')
    id = id.replace('_APOS_', "'")
    id = id.replace('&amp;', '&')
    id = id.replace('&lt;', '<')
    id = id.replace('&gt;', '>')
    id = id.replace('&quot;', '"')
    if use_hyphens:
        id = id.replace('__', '-')
    else:
        id = id.replace("-", "__")
    if fix_compartments:
        if len(id) > 2:
            if (id[-3] == "(" and id[-1] == ")") or \
                    (id[-3] == "[" and id[-1] == "]"):
                id = id[:-3] + "_" + id[-2]
    return id


def read_legacy_sbml(filename, use_hyphens=False):
    """read in an sbml file and fix the sbml id's"""
    model = create_cobra_model_from_sbml_file(filename, old_sbml=True)
    for metabolite in model.metabolites:
        metabolite.id = fix_legacy_id(metabolite.id)
    model.metabolites._generate_index()
    for reaction in model.reactions:
        reaction.id = fix_legacy_id(reaction.id)
        if reaction.id.startswith("EX_") and reaction.id.endswith("(e)"):
            reaction.id = reaction.id[:-3] + "_e"
    model.reactions._generate_index()
    # remove boundary metabolites (end in _b and only present in exchanges)
    for metabolite in list(model.metabolites):
        if not metabolite.id.endswith("_b"):
            continue
        if len(metabolite._reaction) == 1:
            if list(metabolite._reaction)[0].id.startswith("EX_"):
                metabolite.remove_from_model()
    model.metabolites._generate_index()
    return model
