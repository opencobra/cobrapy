#cobra/sbml.py: Tools for reading / writing SBML now contained in
#this module
#System modules
from .. import Model, Reaction, Metabolite, Formula
from os.path import isfile
from copy import deepcopy
from numpy import zeros
from scipy.sparse import lil_matrix
from time import time
import re
#Add in the switch for importing the java sbml if this is run in jython
from libsbml import SBMLDocument, SpeciesReference, KineticLaw, Parameter
from libsbml import readSBML, writeSBML
from libsbml import UNIT_KIND_MOLE, UNIT_KIND_GRAM, UNIT_KIND_SECOND, UNIT_KIND_DIMENSIONLESS
def parse_legacy_id(the_id, the_compartment=None, the_type='metabolite',
                    use_hyphens=False):
    """Deals with a bunch of problems due to bigg.ucsd.edu not following SBML standards

    the_id: String.

    the_compartment: String.

    the_type: String.  Currently only 'metabolite' is supported

    use_hyphens:   Boolean.  If True, double underscores (__) in an SBML ID will be converted to hyphens
    
    """
    if use_hyphens:
        the_id = the_id.replace('__','-')
    if the_type == 'metabolite':
        if the_id.split('_')[-1] == the_compartment:
            #Reformat Ids to match convention in Palsson Lab.
            the_id = the_id[:-len(the_compartment)-1]
        the_id += '[%s]'%the_compartment
    return the_id
def create_cobra_model_from_sbml_file(sbml_filename, old_sbml=False, legacy_metabolite=False,
                                      print_time=False, use_hyphens=False):
    """convert an SBML XML file into a cobra.Model object.  Supports
    SBML Level 2 Versions 1 and 4

    sbml_filename: String.
    
    old_sbml:  Boolean. Set to True if the XML file has metabolite
    formula appended to metabolite names.  This was a poorly designed
    artifact that persists in some models.

    legacy_metabolite: Boolean.  If True then assume that the metabolite id has
    the compartment id appended after an underscore (e.g. _c for cytosol).  This
    has not been implemented but will be soon.

    print_time:  Boolean.  Print the time requirements for different sections

    use_hyphens:   Boolean.  If True, double underscores (__) in an SBML ID will be converted to hyphens

	
    """
     # Ensure that the file exists
    if not isfile(sbml_filename):
        raise IOError('Your SBML file is not found: %s'%sbml_filename)
    #Expressions to change SBML Ids to Palsson Lab Ids
    metabolite_re = re.compile('^M_')
    reaction_re = re.compile('^R_')
    compartment_re = re.compile('^C_')
    if print_time:
        start_time = time()
    model_doc = readSBML(sbml_filename)
    if print_time:
       print 'Loading %s took %1.2f seconds'%(sbml_filename,
                                              time()-start_time)
                                             
     
    sbml_model = model_doc.getModel()
    sbml_model_id = sbml_model.getId()
    sbml_species = sbml_model.getListOfSpecies()
    sbml_reactions = sbml_model.getListOfReactions()
    sbml_compartments = sbml_model.getListOfCompartments()
    compartment_dict = dict([(compartment_re.split(x.getId())[-1], x.getName())
                             for x in sbml_compartments])
    if legacy_metabolite:
        #Deal with the palsson lab appending the compartment id to the metabolite id
        new_dict = {}
        for the_id, the_name in compartment_dict.items():
            if the_name == '':
                new_dict[the_id[0].lower()] = the_id
            else:
                new_dict[the_id] = the_name
        compartment_dict = new_dict
        legacy_compartment_converter = dict([(v,k)
                                             for k, v in compartment_dict.items()])
    if print_time:
        start_time = time()
    metabolite_dict = {}
    #Convert sbml_metabolites to cobra.Metabolites
    for sbml_metabolite in sbml_species:
        #Skip sbml boundary species
        if sbml_metabolite.getBoundaryCondition():
            continue

        if (old_sbml or legacy_metabolite) and \
               sbml_metabolite.getId().endswith('_b'):
            #Deal with incorrect sbml from bigg.ucsd.edu
            continue
        tmp_metabolite = Metabolite()
        metabolite_id = tmp_metabolite.id = sbml_metabolite.getId()
        tmp_metabolite.compartment = compartment_re.split(sbml_metabolite.getCompartment())[-1]
        if legacy_metabolite:
            if tmp_metabolite.compartment not in compartment_dict:
                tmp_metabolite.compartment = legacy_compartment_converter[tmp_metabolite.compartment]
            tmp_metabolite.id = parse_legacy_id(tmp_metabolite.id, tmp_metabolite.compartment,
                                                use_hyphens=use_hyphens)
        if use_hyphens:
            tmp_metabolite.id = metabolite_re.split(tmp_metabolite.id)[-1].replace('__','-')
        else:
            #Just in case the SBML ids are ill-formed and use -
            tmp_metabolite.id = metabolite_re.split(tmp_metabolite.id)[-1].replace('-','__')
        tmp_metabolite.name = sbml_metabolite.getName()
        tmp_formula = ''
        tmp_metabolite.charge = sbml_metabolite.getCharge()
        tmp_metabolite.notes = parse_legacy_sbml_notes(sbml_metabolite.getNotesString())
        for the_key in tmp_metabolite.notes.keys():
            if the_key.lower() == 'formula':
                tmp_formula = tmp_metabolite.notes.pop(the_key)[0]
                break
        if tmp_formula == '' and old_sbml:
            tmp_formula = tmp_metabolite.name.split('_')[-1]
            tmp_metabolite.name = tmp_metabolite.name[:-len(tmp_formula)-1]
        tmp_metabolite.formula = Formula(tmp_formula)
        metabolite_dict.update({metabolite_id: tmp_metabolite})
    if print_time:
       print 'Parsing %s took %1.2f seconds'%('metabolites',
                                              time()-start_time)


    if print_time:
        start_time = time()
    #Construct the vectors and matrices for holding connectivity and numerical info
    #to feed to the cobra toolbox.
    #Always assume steady state simulations so b is set to 0
    cobra_reaction_list = []
    for sbml_reaction in sbml_reactions:
        if use_hyphens:
            #Change the ids to match conventions used by the Palsson lab.
            reaction = Reaction(reaction_re.split(sbml_reaction.getId())[-1].replace('__','-'))
        else:
            #Just in case the SBML ids are ill-formed and use -
            reaction = Reaction(reaction_re.split(sbml_reaction.getId())[-1].replace('-','__'))
        cobra_reaction_list.append(reaction)
        reaction.exchange_reaction = 0
        reaction.name = sbml_reaction.getName()
        cobra_metabolites = {}
        #Use the cobra.Metabolite class here
        for sbml_metabolite in sbml_reaction.getListOfReactants():
            tmp_metabolite_id = sbml_metabolite.getSpecies()
            #This deals with boundary metabolites
            if tmp_metabolite_id in metabolite_dict:
                tmp_metabolite = deepcopy(metabolite_dict[tmp_metabolite_id])
                cobra_metabolites[tmp_metabolite] = -sbml_metabolite.getStoichiometry()
            else:
                reaction.boundary = 'system_boundary'
        for sbml_metabolite in sbml_reaction.getListOfProducts():
            tmp_metabolite_id = sbml_metabolite.getSpecies()
            #This deals with boundary metabolites
            if tmp_metabolite_id in metabolite_dict:
                tmp_metabolite = deepcopy(metabolite_dict[tmp_metabolite_id])
                cobra_metabolites[tmp_metabolite] = sbml_metabolite.getStoichiometry()
            else:
                reaction.boundary = 'system_boundary'

        #Parse the kinetic law info here.
        parameter_dict = {}
        #            if isinstance(the_reaction.getKineticLaw(), NoneType):
        if not sbml_reaction.getKineticLaw():
            parameter_dict['lower_bound'] = -1000 
            parameter_dict['upper_bound'] = 1000 
            parameter_dict['objective_coefficient'] = 0 
        else:
            for sbml_parameter in sbml_reaction.getKineticLaw().getListOfParameters():
                parameter_dict[sbml_parameter.getId().lower()] = sbml_parameter.getValue()

        if 'lower_bound' in parameter_dict:
            the_key = 'lower_bound'
        elif 'lower bound' in parameter_dict:
            the_key = 'lower bound'
        reaction.lower_bound = parameter_dict[the_key]
        if 'upper_bound' in parameter_dict:
            the_key = 'upper_bound'
        elif 'upper bound' in parameter_dict:
            the_key = 'upper bound'
        reaction.upper_bound = parameter_dict[the_key]
        if 'objective_coefficient' in parameter_dict:
            the_key = 'objective_coefficient'
        elif 'objective coefficient' in parameter_dict:
            the_key = 'objective coefficient'
        reaction.objective_coefficient = parameter_dict[the_key]
        reaction_note_dict = parse_legacy_sbml_notes(sbml_reaction.getNotesString())
        #Parse the reaction notes.
        #POTENTIAL BUG: DEALING WITH LEGACY 'SBML' THAT IS NOT IN A
        #STANDARD FORMAT
        #TODO: READ IN OTHER NOTES AND GIVE THEM A reaction_ prefix.
        #TODO: Make sure genes get added as objects
        if reaction_note_dict.has_key('GENE ASSOCIATION'):
            reaction.gene_reaction_rule = reaction_note_dict['GENE ASSOCIATION'][0]
            reaction.parse_gene_association() 
            if reaction_note_dict.has_key('GENE LIST'):
                reaction.systematic_names = reaction_note_dict['GENE LIST'][0]
            elif reaction_note_dict.has_key('GENES') and \
                     reaction_note_dict['GENES'] != ['']:
                reaction.systematic_names = reaction_note_dict['GENES'][0]
            elif reaction_note_dict.has_key('LOCUS'):
                gene_id_to_object = dict([(x.id, x) for x in reaction._genes])
                for the_row in reaction_note_dict['LOCUS']:
                    tmp_row_dict = {}
                    the_row = 'LOCUS:' + the_row.lstrip('_').rstrip('#')
                    for the_item in the_row.split('#'):
                        try:
                            k, v = the_item.split(':')
                            tmp_row_dict[k] = v
                        except ValueError, e:
                            print the_item
                            raise e
                    tmp_locus_id = tmp_row_dict['LOCUS']
                    if 'TRANSCRIPT' in tmp_row_dict:
                            tmp_locus_id = tmp_locus_id + '.' + tmp_row_dict['TRANSCRIPT']
                    gene_id_to_object[tmp_locus_id].name = tmp_row_dict['ABBREVIATION']

        if reaction_note_dict.has_key('SUBSYSTEM'):
            reaction.subsystem = reaction_note_dict['SUBSYSTEM'][0]   

        reaction.reversibility = int(sbml_reaction.getReversible())
        #TODO: Use the cobra.metabolite objects here.
        reaction.add_metabolites(cobra_metabolites)

    if print_time:
       print 'Parsing %s took %1.2f seconds'%('reactions',
                                              time()-start_time)

    if print_time:
        start_time = time()
    #Now, add all of the reactions to the model.
    cobra_model = Model(sbml_model_id)
    cobra_model.description = sbml_model.getId()
    #Populate the compartment list - This will be done based on cobra.Metabolites
    #in cobra.Reactions in the future.
    cobra_model.compartments = compartment_dict

    cobra_model.add_reactions(cobra_reaction_list)
    if print_time:
        print '%s took %1.2f seconds'%('Adding reactions',
                                       time()-start_time)
        #cobra_model.update_rules()
    return cobra_model


def parse_legacy_sbml_notes(note_string, note_delimiter = ':'):
    """Deal with legacy SBML format issues arising from the
    COBRA Toolbox for MATLAB and BiGG.ucsd.edu developers.

	
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
        the_note = note_string[(note_start + len(start_tag)):note_end].lstrip(' ').rstrip(' ')
        if note_delimiter in the_note:
            note_delimiter_index = the_note.index(note_delimiter)
            note_field = the_note[:note_delimiter_index].lstrip(' ').rstrip(' ').replace('_',' ').upper()
            note_value = the_note[(note_delimiter_index+1):].lstrip(' ').rstrip(' ')
            if note_dict.has_key(note_field ):
                note_dict[note_field ].append(note_value)
            else:
                note_dict[note_field] = [note_value]
        note_string = note_string[(note_end+len(end_tag)): ]

    
    return note_dict


def write_cobra_model_to_sbml_file(cobra_model, sbml_filename,
                                   sbml_level=2, sbml_version=1,
                                   print_time=False):
    """Write a cobra.Model object to an SBML XML file.

    cobra_model:  A cobra.Model object

    sbml_filename:  The file to write the SBML XML to.

    sbml_level:  2 is the only level supported at the moment.

    sbml_version: 1 is the only version supported at the moment.

    print_time:  Boolean.  Print the time requirements for different sections

    TODO: Update the NOTES to match the SBML standard and provide support for
    Level 2 Version 4
    
    """
    #Add in the common compartment abbreviations.  If there are additional compartments
    #they also need to be added.
    note_start_tag, note_end_tag = '<p>', '</p>'
    if sbml_level > 2 or (sbml_level == 2 and sbml_version == 4):
        note_start_tag, note_end_tag = '<html:p>', '</html:p>'
        
    if not hasattr(cobra_model, 'compartments'):
        cobra_model.compartments = {'c': 'cytosol',
                                    'p': 'periplasm',
                                    'e': 'extracellular'}
    
    sbml_doc = SBMLDocument(sbml_level, sbml_version)
    sbml_model = sbml_doc.createModel(cobra_model.description.split('.')[0])
    #Note need to set units
    reaction_units = 'mmol_per_gDW_per_hr'
    model_units = sbml_model.createUnitDefinition()
    model_units.setId(reaction_units)
    sbml_unit = model_units.createUnit()
    sbml_unit.setKind(UNIT_KIND_MOLE)
    sbml_unit.setScale(-3)
    sbml_unit = model_units.createUnit()
    sbml_unit.setKind(UNIT_KIND_GRAM)
    sbml_unit.setExponent(-1)
    sbml_unit = model_units.createUnit()
    sbml_unit.setKind(UNIT_KIND_SECOND)
    sbml_unit.setMultiplier(1.0/60/60)
    sbml_unit.setExponent(-1)

    
    for the_key in cobra_model.compartments.keys():
        sbml_comp = sbml_model.createCompartment()
        sbml_comp.setId(the_key)
        sbml_comp.setName(cobra_model.compartments[the_key])
        sbml_comp.setSize(1) #Just to get rid of warnings

    if print_time:
        start_time = time()
    #Use this dict to allow for fast look up of species id
    #for references created in the reaction section.
    metabolite_dict = {}

    for cobra_metabolite in cobra_model.metabolites:
        metabolite_dict[cobra_metabolite.id] =  add_sbml_species(sbml_model,
                                                                 cobra_metabolite,
                                                                 note_start_tag=note_start_tag,
                                                                 note_end_tag=note_end_tag)

    if print_time:
        print 'Adding %s took %1.2f seconds'%('metabolites',
                                              time()-start_time)
    if print_time:
        start_time = time()
    for the_reaction in cobra_model.reactions:
        #This is probably the culprit.  Including cobra.Reaction
        #objects explicitly in cobra.Model will speed this up.
        sbml_reaction = sbml_model.createReaction()
        #Need to remove - for proper SBML.  Replace with __
        the_reaction_id = 'R_' + the_reaction.id.replace('-','__' )
        sbml_reaction.setId(the_reaction_id)
        if the_reaction.reversibility == 1:
            sbml_reaction.setReversible(True)
        else:
            sbml_reaction.setReversible(False)
        if the_reaction.name:
            sbml_reaction.setName(the_reaction.name)
        else:
            sbml_reaction.setName(the_reaction.id)
        #Add in the reactant/product references
        for the_metabolite, the_coefficient in the_reaction._metabolites.items():
            sbml_stoichiometry = the_coefficient
            metabolite_id = str(metabolite_dict[the_metabolite.id])
            #Each SpeciesReference must have a unique id
            if sbml_stoichiometry < 0:
                species_reference = sbml_reaction.createReactant()
            else:
                species_reference = sbml_reaction.createProduct()
            species_reference.setId(metabolite_id + '_' + the_reaction_id)
            species_reference.setSpecies(metabolite_id)
            species_reference.setStoichiometry(abs(sbml_stoichiometry))
        #Deal with the case where the reaction is a boundary reaction
        if len(the_reaction._metabolites) == 1:
            the_metabolite, the_coefficient = the_reaction._metabolites.items()[0]
            the_metabolite = the_metabolite.copy()
            metabolite_id = add_sbml_species(sbml_model, the_metabolite,
                                             note_start_tag=note_start_tag,
                                             note_end_tag=note_end_tag,
                                             boundary_metabolite=True)
            the_coefficient *= -1
            #Each SpeciesReference must have a unique id
            if sbml_stoichiometry < 0:
                species_reference = sbml_reaction.createReactant()
            else:
                species_reference = sbml_reaction.createProduct()
            species_reference.setId(metabolite_id + '_' + the_reaction_id)
            species_reference.setSpecies(metabolite_id)
            species_reference.setStoichiometry(abs(sbml_stoichiometry))
            
        #Add in the kineticLaw
        sbml_law = KineticLaw(sbml_level, sbml_version)
        sbml_law.setId('FLUX_VALUE')
        sbml_law.setFormula('FLUX_VALUE')
        reaction_parameter_dict = {'LOWER_BOUND': [the_reaction.lower_bound, reaction_units],
                                   'UPPER_BOUND': [the_reaction.upper_bound, reaction_units],
                                   'FLUX_VALUE': [0, reaction_units],
                                   'OBJECTIVE_COEFFICIENT': [the_reaction.objective_coefficient,
                                                             'dimensionless']}
        for k, v in reaction_parameter_dict.items():
            sbml_parameter = Parameter(sbml_level, sbml_version)
            sbml_parameter.setId(k)
            if hasattr(v, '__iter__'):
                sbml_parameter.setValue(v[0])
                sbml_parameter.setUnits(v[1])
            else:
                sbml_parameter.setValue(v)
            sbml_law.addParameter(sbml_parameter)
        sbml_reaction.setKineticLaw(sbml_law)
        sbml_reaction.setNotes('<html xmlns="http://www.w3.org/1999/xhtml">%sGENE_ASSOCIATION: %s%s%sSUBSYSTEM: %s%s</html>'%(note_start_tag,
                                                                 the_reaction.gene_reaction_rule,
                                                                 note_end_tag,
                                                                 note_start_tag,
                                                                 the_reaction.subsystem,
                                                                 note_end_tag))

    if print_time:
       print 'Adding %s took %1.2f seconds'%('reactions',
                                             time()-start_time)

    writeSBML(sbml_doc, sbml_filename)

def add_sbml_species(sbml_model, cobra_metabolite,                                                                  note_start_tag, note_end_tag, boundary_metabolite=False):
    """A helper function for adding cobra metabolites to an sbml model.

    sbml_model: sbml_model object

    cobra_metabolite: a cobra.Metabolite object

    note_start_tag: the start tag for parsing cobra notes. this will eventually
    be supplanted when COBRA is worked into sbml.

    note_end_tag: the end tag for parsing cobra notes. this will eventually
    be supplanted when COBRA is worked into sbml.

    """
    sbml_species = sbml_model.createSpecies()
    the_id = 'M_' + cobra_metabolite.id.replace('-', '__')
    #Deal with legacy naming issues
    the_compartment = cobra_metabolite.compartment
    if the_id.endswith('[%s]'%the_compartment):
        the_id = the_id[:-len('[%s]'%the_compartment)]
    elif not the_id.endswith('_%s'%the_compartment):
        the_id += '_%s'%the_compartment
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
    try:
        sbml_species.setCompartment(the_compartment)
    except:
        print 'metabolite failed: %s'%the_id
        return cobra_metabolite
    if cobra_metabolite.charge is not None:
        sbml_species.setCharge(cobra_metabolite.charge)
    if hasattr(cobra_metabolite.formula, 'id') or \
       hasattr(cobra_metabolite.notes, 'items'):
        tmp_note =  '<html xmlns="http://www.w3.org/1999/xhtml">'
        if hasattr(cobra_metabolite.formula, 'id'):
            tmp_note += '%sFORMULA: %s%s'%(note_start_tag,
                                              cobra_metabolite.formula.id,
                                              note_end_tag)
        if hasattr(cobra_metabolite.notes, 'items'):
            for the_id_type, the_id in cobra_metabolite.notes.items():
                if not isinstance(the_id_type, str):
                    the_id_type = repr(the_id_type)
                if hasattr(the_id, '__iter__') and len(the_id) == 1:
                    the_id = the_id[0]
                if not isinstance(the_id, str):
                    the_id = repr(the_id)
                tmp_note += '%s%s: %s%s'%(note_start_tag,
                                             the_id_type,
                                             the_id, note_end_tag)
        sbml_species.setNotes(tmp_note + '</html>')
    return metabolite_id


def fix_legacy_id(id, use_hyphens=False):
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
    return id


def read_legacy_sbml(filename, use_hyphens=False):
    """read in an sbml file and fix the sbml id's"""
    model = create_cobra_model_from_sbml_file(filename)
    for metabolite in model.metabolites:
        metabolite.id = fix_legacy_id(metabolite.id)
    model.metabolites._generate_index()
    for reaction in model.reactions:
        reaction.id = fix_legacy_id(reaction.id)
        if reaction.id.startswith("EX_") and reaction.id.endswith("(e)"):
            reaction.id = reaction.id[:-3] + "_e"
        reaction.reconstruct_reaction()
    model.reactions._generate_index()
    return model
