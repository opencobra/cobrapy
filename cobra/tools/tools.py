#cobra.tools.py
#System modules
import numpy, os, re, scipy, sys
from copy import deepcopy
from cPickle import dump
from numpy import apply_along_axis,array, corrcoef, diag, hstack, isnan, isinf
from numpy import  log, mean, std, tril, vstack, where, unique, zeros
from scipy.stats import norm, randint
from collections import defaultdict
#TODO: Move application specific funtions to their own
#files and split other functions out as well
print 'WARNING: cobra.tools contains miscellaneous functions that are under development' +\
      'use of this module IS NOT RECOMMENDED'
def translate_entrez_via_homologene(gene_list,source_taxa='Homo sapiens',
                                    target_taxa='Mus musculus',
                                    homologene_directory='/Users/danie/sysbep/Data/homologene/64/',
                                     default_value='S0001'):
    """Uses the files in homologene_directory to map the entities in gene_list
    from source_taxa to target_taxa.  If there isn't a homolog then the
    default_value is used for the mapping.

    Returns a list with the corresponding target_taxa genes.
    
    """
    taxa_map = {}
    with open(homologene_directory + 'taxid_taxname.txt') as in_file:
        the_data = [x.rstrip('\r\n').split('\t') for x in in_file.readlines()]
        [taxa_map.update({x[1]: x[0]}) for x in the_data]

    homologene_to_target_gi = {}
    source_gi_to_homologene = {}
    source_id = taxa_map[source_taxa]
    target_id = taxa_map[target_taxa]
    with open(homologene_directory + 'homologene.data') as in_file:
        for the_line in in_file.readlines():
            the_line = the_line.rstrip('\n').split('\t')
            taxa_id = the_line[1]
            if taxa_id == source_id:
                source_gi_to_homologene[the_line[2]] = the_line[0]
            elif taxa_id == target_id:
                homologene_to_target_gi[the_line[0]] = the_line[2]

    target_gene_ids = []
    for the_gene in gene_list:
        gene_id = default_value
        if the_gene in source_gi_to_homologene:
            homologene = source_gi_to_homologene[the_gene]
            if homologene in homologene_to_target_gi:
                gene_id = homologene_to_target_gi[homologene]
        target_gene_ids.append(gene_id)
    return(target_gene_ids)

#Get enzyme inhibition data from BRENDA
def get_druggable_targets_from_brenda(organism_genus_species):
    """get_druggable_targets_from_brenda serves as a wrapper to BRENDAs SOAP
    getInhibitors.  That organizes the output and allows one to use paralell
    python. BRENDA is the enzyme database

    organism_genus_species: The organism_genus_species string from BRENDA for
    the organism of interest

    NOTE: This will sometimes fail due to the presence of non-unicode strings
    in BRENDA even though their web-service defines the encoding as unicode.
    
    """
    return(get_inhibitors_from_brenda({'organism': organism_genus_species}))

def get_inhibitors_from_brenda(brenda_dict):
    """get_inhibitors_from_brenda serves as a wrapper to BRENDAs SOAP
    getInhibitors.  That organizes the output and allows one to use
    paralell python.

    brenda_dict: The query parameters to feed to BRENDA's function
    
    NOTE: This will sometimes fail due to the presence of non-unicode strings
    in BRENDA even though their web-service defines the encoding as unicode.

    """
    from SOAPpy import WSDL, faultType
    from xml.sax import SAXParseException as SAXParseException
    wsdl = 'http://www.brenda-enzymes.org/soap/brenda.wsdl'
    serv = WSDL.Proxy(wsdl)
    try:
        tmp_inhibitors = serv.getInhibitors(brenda_dict)
    except faultType or SAXParseException:
        tmp_inhibitor_list = 'SOAP error: %s'%sys.exc_info()[1]
        tmp_druggable_ec = brenda_dict
    else:
        tmp_inhibitor_list = []
        tmp_druggable_ec = []
        for tmp_inhibitor in tmp_inhibitors:
            tmp_inhibitor = tmp_inhibitor['item']
            tmp_dict = {}
            for tmp_field in tmp_inhibitor:
                tmp_dict[tmp_field['key']] = tmp_field['value']
            tmp_inhibitor_list.append(tmp_dict)
            tmp_druggable_ec.append(tmp_dict['ecNumber'])
        tmp_druggable_ec = list(set(tmp_druggable_ec))
        tmp_druggable_ec.sort()
        if len(tmp_druggable_ec) == 1:
            tmp_druggable_ec = tmp_druggable_ec[0]
    return {'list': tmp_inhibitor_list, 'ec':tmp_druggable_ec}

def select_rows_and_columns(the_symmetric_array, the_indices):
    """select_rows_and_columns will select a set of rows and columns from a numpy
    array, and returns them as an array.

    the_symmetric_array: A symmetric numpy array

    the_indices: A list of indices.

    """
    tmp_array = the_symmetric_array[the_indices[0], :]
    for tmp_idx in the_indices[1:]:
        tmp_array = vstack((tmp_array, the_symmetric_array[tmp_idx, :]))
    n_elements = len(the_indices)
    reduced_array = tmp_array[:, the_indices[0]].reshape(n_elements, 1)
    for tmp_idx in the_indices[1:]:
        reduced_array = hstack((reduced_array,
                                tmp_array[:, tmp_idx].reshape(n_elements, 1)))
    return deepcopy(reduced_array)

def remove_essential_elements(the_symmetric_array, the_element_names=None):
    """Used to remove rows and columns from an array where the
    corresponding diagonal is not 0.  Typically used with analyzing
    output of cobra.flux_analysis.double_deletion

    the_symmetric_array: A symmetric numpy array

    the_element_names: None or a list of strings corresponding to the
    rows of the_symmetric_array.
    
    """
    nonessential_elements = list(where(diag(the_symmetric_array) != 0)[0])
    nonessential_array = select_rows_and_columns(the_symmetric_array,
                                                  nonessential_elements)
    elements_in_array = []
    if not the_element_names:
        elements_in_array = nonessential_elements
    else:
        map(lambda x: elements_in_array.append(the_element_names[x]),
            nonessential_elements)
    return {'array': nonessential_array,
            'elements': elements_in_array}

def remove_nonsynthetic_lethal(the_symmetric_array, the_element_names=None):
    """Used to remove rows and columns from an array where the corresponding
    rows and columns only contain values greater than 0.  Typically used
    with analyzing output of cobra.flux_analysis.double_deletion

    the_symmetric_array: A symmetric numpy array

    the_element_names: None or a list of strings corresponding to the
    rows of the_symmetric_array.
    """
    
    n_elements = the_symmetric_array.shape[0]
    synthetic_lethal_indices = list(where(apply_along_axis(sum, 1,
                                                           the_symmetric_array)\
                                          < n_elements)[0])
    synthetic_lethal_array = select_rows_and_columns(the_symmetric_array,
                                                      synthetic_lethal_indices)
    elements_in_array = []
    if the_element_names == []:
        elements_in_array = synthetic_lethal_indices 
    else:
        map(lambda x: elements_in_array.append(the_element_names[x]),
             synthetic_lethal_indices )
    return {'array': synthetic_lethal_array, 'elements':elements_in_array}



def create_sif(the_symmetric_array, the_elements, the_filename,
                the_interaction=1, interaction_type='syn_lethal'):
    """Convert a symmetric array into a SIF-formatted file that can
    be opened with Cytoscape.

    the_symmetric_array:  A symmetric numpy array.

    the_elements: Human understandable names for each row in the array

    the_filename: The file to save the SIF to.

    the_interaction: The value that must exist at a point in the array
    to indicate that row i interacts with column j.

    interaction_type: The name to call the interaction in the SIF file.
    
    """
    out_file = open(the_filename, 'w')
    for the_row_idx in range(the_symmetric_array.shape[0]):
        tmp_string = the_elements[the_row_idx] + '\t' + interaction_type
        #use the_row_idx as an offset to deal with symmetry
        for the_interactor in list(where(the_symmetric_array[the_row_idx,
                                                             the_row_idx:] \
                                         == the_interaction)[0]):
            tmp_string += '\t%s'%the_elements[the_row_idx + the_interactor]

        if not tmp_string.endswith(interaction_type):
            out_file.write(tmp_string + '\n')
    out_file.close()

def gene_to_reaction_value(cobra_model, gene_value_dict=None,
                           value_type='gene_p', gene_prefix='(STM|PSLT){1}',
                           default_value=1, complex_aggregate=min,
                           return_id=None, spontaneous_value=0):
    """Converts gene associated values to reaction associated values
    based on gene_reaction_rules in cobra_model.reactions[:].gene_reaction_rule

    gene_value_dict: A dictionary where the keys are from cobra_model.genes
    and the values are some number such as the p value.

    value_type: 'gene_p', 'gene_expression' indicate that these are
    p-values or expression intensities for genes.  'protein_p' and
    'protein_expression' are not yet implemented.

    gene_prefix: Text.  A regular expression that distinguishes gene ids from
    junk in cobra_model.gene_reaction_ruiles

    default_value: Float.  The value to assign to genes in a gene_reaction_rule
    that are not in gene_value_dict.

    complex_aggregate: Function that aggregates the values for an enzyme complex

    return_id: None or Identifier for the run to help maintain identity when
    running a large number of these functions in parallel.

    spontaneous_value: The value to assign to spontaneous reactions

    returns a dict {cobra_model.reaction[i]: score[i]}
    
    """
    gene_re = re.compile(gene_prefix + '(\w|\.)+')
    #2010-06-15: Cannot remember by the or was atomized, so deleted it
    #because of problems.  Potential BUG.
    #spontaneous_re = re.compile('[sS]0001( +or)?')
    spontaneous_re = re.compile('[sS]0001')
    spontaneous_ids = ['s0001', 'S0001']
    #If a gene value dict was not fed as a parameter
    #extract from the cobra_model
    if not gene_value_dict:
        if value_type.endswith('_p'):
            gene_value_dict = dict(zip(cobra_model.cobra_model.genes,
                                       cobra_model.p_values[value_type.replace( '_p ', '')]))
        elif value_type.endswith('_expression'):
            gene_value_dict = dict(zip(cobra_model.cobra_model.genes,
                                       cobra_model.expression_values[\
                                           value_type.replace('_expression',
                                                              '')]))
    if hasattr(cobra_model, 'cobra_model'):
         cobra_model = cobra_model.cobra_model
    #First collect all the reactions that are in the gene list
    #Remove the spontaneous elements from the gene_value_dict because they
    #shouldn't have any values
    [gene_value_dict.pop(x)
     for x in spontaneous_ids if x in gene_value_dict]
    #remove the rows that have an NAN or inf value
    [gene_value_dict.pop(k)
      for k, v in gene_value_dict.items()
     if isnan(v).any() or isinf(v).any()]

    gene_id_to_value = dict([(k.id, v) for k, v in gene_value_dict.iteritems()])
    #Since reactions and genes point to each other we
    #can aggregate based on these
    the_reactions = set()
    [the_reactions.update(x._reaction) for x in gene_value_dict]

    #This deals with the fact that some genes have the same base but an additional suffix
    #HERE 2010-04-23: Dealing with ids being substrings of other ids.
    new_rules = {}
    #Make copies of the gene_reaction_rules for the reactions and
    #substitute in the spontaneous value

    for the_reaction in the_reactions:
        the_rule =  deepcopy(the_reaction.gene_reaction_rule)
        if the_rule == '':
            new_rules[the_reaction] = spontaneous_value
            continue
        else:
            the_rule = spontaneous_re.sub(repr(spontaneous_value),
                                          the_rule)

        the_rule = the_rule.replace('+', ' and ').replace(',', ' or ')
        the_gene_list = the_rule.replace('(', '').replace(')', '').replace(' and ', '\t' ).replace(' or ', '\t').replace(' ', '').split('\t')

        the_replacements = []
        for the_gene in the_gene_list:
            if the_gene in gene_id_to_value:
                the_replacements.append(repr(gene_id_to_value[the_gene]))
            elif the_gene in spontaneous_ids:
                the_replacements.append(repr(spontaneous_value))
            else:
                the_replacements.append(repr(default_value))
        #ordered_replace is used to deal with some genes being substrings of other genes
        new_rules[the_reaction] = ordered_replace(the_rule, the_gene_list,
                                                  the_replacements)
    #Now calculate the p value for the reaction based on the complex_aggregate
    #flag.
    #TODO:  This needs to be retooled to deal with enzyme complexes
    for the_reaction, the_rule in new_rules.iteritems():
        if the_rule.startswith('(') and (the_rule.find('or') != -1
                                           or the_rule.find('and') != -1):
            new_rules[the_reaction] = complex_aggregate(eval('(%s)'%the_rule.replace('and', ',').replace('or', ',').replace('(', '').replace(')', '')))
        else:
            new_rules[the_reaction] = eval(the_rule.lstrip('(').rstrip(')'))
    if return_id:
        return([return_id, new_rules])
    else:
        return(new_rules)



def ordered_replace(the_string, elements_to_replace, the_replacements):
    """Some times there are strings that need elements to be replaced,
    however, the elements may be substrings of other elements or even the
    replacement elements.  This function will replace each element in order.

    the_string: The string that needs to be modified.

    elements_to_replace:  A list of elements to replace in a string
    in order.  Only the first instance of an element will be replaced

    the_replacements: The corresponding list of replacements.
    
    """
    new_string = ''
    for s, v in zip(elements_to_replace, the_replacements):
        the_prefix, the_string = the_string.split(s, 1)
        new_string += the_prefix + v
    the_string = new_string + the_string
    return the_string



def log_function(x, base=2.):
    return log(x) / log(base)
def exp_function(x, base = 2.):
    return base**x


