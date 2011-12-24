#cobra.io.feature_extraction.py
#Tools for parsing agilent feature extraction files.
import pdb
from time import time
from rpy2 import robjects
r = robjects.r
import rpy2.robjects.numpy2ri
robjects.numpy2ri.activate()

from copy import deepcopy
import os, sys
from collections import defaultdict
from numpy import mean, array,std, log10
from cobra.tools import log_function
from cobra.stats.stats import combine_p_values, error_weighted
from warnings import warn
warn("WARNING: cobra.io.feature_extraction is not ready for general use ")

def parse_file(in_file, polarity=1, quality_control=True, return_id='accession',
               normalization=None, lowess_parameter=0.33, log_base=10, single_channel=False,
               single_field='gProcessedSignal'):
    """Extract data from feature extraction >=9.5 files.  Returns
    the average log ratios for each return_id.

    in_file: String.  Name of input file.

    polarity: 1 or -1.  Indicates whether to do red over green (1) or
    green over red.  If normalization isn't performed then the polarity
    is multiplied by the log ratio.

    return_id: 'accession' or 'probe'

    normalization: 'lowess' or None

    lowess_parameter: Float.  Smoothing parameter for lowess normalization

    """
    ## in_file = '/Users/danie/tg/Work/Li_et_al/Data/Feature_Extraction/for GEO/Oxam+.txt'
    ## polarity = 1
    ## quality_control = True
    ## return_id = 'accession'
    ## normalization = 'lowess'
    ## lowess_parameter = 0.33
    ## log_base = 10
    channel_prefixes = ['r','g']
    if single_channel:
        channel_prefixes = [x for x in channel_prefixes
                            if single_field.startswith(x)]
    with open(in_file) as in_file_handle:
        the_header = in_file_handle.readline().rstrip('\r\n').split('\t')
        while not the_header[0] == 'FEATURES':
            the_header = in_file_handle.readline().rstrip('\r\n').split('\t')
        #parse the rows.  skip the non-data spots
        the_data = [x.rstrip('\r\n').split('\t')
                     for x in in_file_handle.readlines()
                    if 'GE_BrightCorner' not in x and 'DarkCorner' not in x]
 
    if return_id.lower() == 'accession':
        gene_index = the_header.index('SystematicName')
    elif return_id.lower() == 'probe':
        gene_index = the_header.index('ProbeName')
    else:
        raise Exception("return_id must be 'accession' or 'probe' not '%s'"%return_id)

    
    if quality_control:
        #Get the column indices for the quality control statistics and
        #assign the value that must be met to pass the test.
        quality_control_indices = {}
        [quality_control_indices.update({the_header.index(k): '1'})
         for k in map(lambda x: x + 'IsFound', channel_prefixes)
         if k in the_header] #1 means the spot is found
        [[quality_control_indices.update({the_header.index(k): '0'})
         for k in map(lambda x: x + y, channel_prefixes)
          if k in the_header]
         for y in ['IsSaturated', 'IsFeatNonUnifOL']] #0 means the spot is not saturatd and not a nonuniform outlier
        
        quality_control_indices[the_header.index('IsManualFlag')] = '0' #0 means it wasn't flagged by the user.
        filtered_data = []
        flagged_rows = []
        for the_row in the_data:
            the_row_is_good = True
            for the_index, passing_flag in quality_control_indices.items():
                if the_row[the_index] != passing_flag:
                    the_row_is_good = False
                    break
            if the_row_is_good:
                filtered_data.append(the_row)
            else:
                flagged_rows.append(the_row)
        the_data = filtered_data #Use the filtered data from here on out


    if single_channel: #Single channel experiments don't have ratios or polarities
        column_to_index = {'intensity_1': the_header.index(single_field)}
    else:

        column_to_index = {'log_ratio': the_header.index('LogRatio'),
                           'log_error': the_header.index('LogRatioError'),
                           'p_value': the_header.index('PValueLogRatio'),
                           'intensity_1': the_header.index('gProcessedSignal'),
                           'intensity_2': the_header.index('rProcessedSignal'),
                           #The last two are in case lowess normalization needs to be
                           #performed.
                           'background_subtracted_1': the_header.index('gBGSubSignal'),
                           'background_subtracted_2': the_header.index('rBGSubSignal')}
        for the_row in the_data: #Speed this up
            the_row[column_to_index['log_ratio']] = polarity*float(the_row[column_to_index['log_ratio']])
        if polarity == -1: #Change the polarity if requested.  Doing so here makes it easier
            #to run the calculations downstream.
            column_to_index.update({'intensity_1': the_header.index('rProcessedSignal'),
                                    'intensity_2': the_header.index('gProcessedSignal'),
                                    'background_subtracted_1': the_header.index('rBGSubSignal'),
                                    'background_subtracted_2': the_header.index('gBGSubSignal')})

        #These need to be popped because they're only used during lowess normalization
        channel_2_index = column_to_index.pop('background_subtracted_2')
        channel_1_index = column_to_index.pop('background_subtracted_1')
        
        #Apply lowess normalization to double channel data if requested.
        if normalization is not None and normalization.lower() == 'lowess':
            if single_channel:
                raise Exception('Lowess normalization does not work with single channel arrays')
            warn('Lowess Normalization looks off, please correct')
            #Polarity is already adjusted above
            log_ratio_index = column_to_index['log_ratio']
            #Now that we're not using the processed values we need to do change intensity
            #indices to use the unprocessed intensities.
            column_to_index['intensity_1'] = channel_1_index
            column_to_index['intensity_2'] = channel_2_index
            channel_2_signal = []
            channel_1_signal = []
            data_to_normalize = []
            #We can't take log ratios of 0 or divide by 0
            for the_row in the_data:
                channel_1_value = float(the_row[channel_1_index])
                channel_2_value = float(the_row[channel_2_index])
                if channel_1_value > 0. and channel_2_value > 0.:
                    data_to_normalize.append(the_row)
                    channel_1_signal.append(channel_1_value)
                    channel_2_signal.append(channel_2_value)
            the_data = data_to_normalize
            #Now perform the lowess normalization
            channel_2_signal = array(channel_2_signal)
            channel_1_signal = array(channel_1_signal)
            a = log10(channel_2_signal*channel_1_signal)/2. 
            m = log10(channel_2_signal/channel_1_signal)
            lowess_fit = array(r.lowess(a, m, lowess_parameter)).T
            m = m - lowess_fit[:,1]
            #Now update the_data list to use the normalized values
            for the_row, the_log_ratio in zip(the_data, list(m)):
                the_row[log_ratio_index] = float(the_log_ratio)


    #Make a dictionary to return the requested values
    gene_dict = dict([(x[gene_index], defaultdict(list))
                      for x in the_data])
    for the_row in the_data:
        the_dict = gene_dict[the_row[gene_index]]
        [the_dict[k].append(float(the_row[v]))
         for k, v in column_to_index.items()]
    #Now combined the items       SPEED THIS UP 
    [collapse_fields(v) for v in gene_dict.values()]
    return gene_dict
    



def parse_annotation(annotation_file, key_type='PrimaryAccession',
                     value_type='EntrezGeneID'):
    """Parse an agilent array annotation file into a dict.

    annotation_file: String.  Name of the annotation file.

    key_type: The column from the annotation_file to use as the key
    for the dict.

    value_type: The column from the annotation_file to use as the value
    for the dict.
    
    """
    annotation_dict = {}

    with open(annotation_file) as in_file:
        the_header = in_file.readline().rstrip('\r\n').split('\t')
        key_index = the_header.index(key_type)
        value_index = the_header.index(value_type)
        the_data = [x.rstrip('\r\n').split('\t') for x in in_file.readlines()]
        [annotation_dict.update({x[key_index]: int(x[value_index])})
         for x in the_data if len(x) >= value_index and x[value_index] != '']
    return annotation_dict


def collapse_fields(data_dict, quantitative_fields=['intensity_1',
                                                    'intensity_2'],
                    log_fields=['log_ratio', 'log_error'],
                    p_fields=['p_value'], log_base=10, error_weighting=True):
    """Collapses the rows for each field from a feature extraction element
    based on the data type.  Quantative fields are averaged, log_fields are
    converted to linear scale, averaged, and then log10 is taken. p-value
    fields are combined using the default method from cobra.stats.stats.
    combine_p_values.

    data_dict: A dictionary of numerical data values.

    quantitative_fields:  A list of the fields in data_dict that can be
    directly averaged.

    log_fields: A list of the fields in data_dict that are in log form
    and must be transformed before averaging.

    p_fields:  A list of the fields in data_dict that are p-values.

    log_base:  The base for the log transform.

    error_weighting:  Boolean.  If True return the error weighted-mean and error.
    Assumes the mean and std devs are the 1st two fields of log_fields.

    NOTE: This happens in place so the data_dict will be modified.

    """
    [data_dict.update({k: mean(data_dict[k])})
     for k in quantitative_fields
     if k in data_dict]
    if log_fields[0] in data_dict:
        if error_weighting:
            log_fields = deepcopy(log_fields)
            mean_field = log_fields.pop(0)
            std_field = log_fields.pop(0)
            the_means = map(lambda x: log_base**x, data_dict[mean_field])
            the_stds = map(lambda x: log_base**x, data_dict[std_field])
            weighted_mean, weighted_std = error_weighted(the_means, the_stds)
            data_dict[mean_field] = log_function(weighted_mean, log_base)
            data_dict[std_field] = log_function(weighted_std, log_base)
        #TODO:  This needs to be updated for log_error
        else:
            [data_dict.update({k: log_function(mean(map(lambda x: log_base**x,
                                                        data_dict[k])),
                                               log_base)})
             for k in log_fields]
    if p_fields[0] in data_dict:
        [data_dict.update({k: combine_p_values(data_dict[k])})
         for k in p_fields]

def combine_files(file_list, annotation_file=None, polarity_list=None,
                  print_time=False, quality_control=False, return_id='entrez',
                  normalization=None):
    """Parse feature extraction files.  This function
    combines multiple technical replicates at the RNA level into a single
    experiment.  Typically multiple replicates will be used when dye-swapping
    is employed.  The combined values are collapsed to entrez gene ids

    file_list: A list of feature extraction file names.  Or a single file.
    A list is only provided if the values are to be combined across all files.

    annotation_file: The agilent annotation file corresponding to the feature
    extraction files in file_list.

    polarity_list:  A list of integers (1 or -1) indicating the polarity of
    the experiment.  1 indicates that the red channel is divided by the green
    channel and -1 indicates to divide the green channel by the red channel

    return_id: String.  Either 'entrez' or 'accession'.  If 'entrez' then
    an annotation file must be supplied to convert from the accession
    ids on the chip to the entrez ids.

    normalization: None or 'lowess'
    
    """
    if print_time:
        start_time = time()
    if return_id.lower() == 'entrez':
        if annotation_file is None:
            raise Exception("An annotation file must be supplied if you want " +\
                            "entrez ids")
        accession_to_entrez = parse_annotation(annotation_file)
    if print_time:
        print '%s %f'%('annotation time',
                       time() - start_time)
        start_time = time()

    parsed_files = []
    if not hasattr(file_list, '__iter__'):
        file_list = [file_list]
    if not hasattr(polarity_list, '__iter__'):
        if polarity_list is None:
            polarity_list = [1]*len(file_list)
        else:
            polarity_list = list(polarity_list)
    #Extract the values to load into the database
    [parsed_files.append(parse_file(the_file,
                                     the_polarity,
                                      quality_control=quality_control,
                                    normalization=normalization))
     for the_file, the_polarity in zip(file_list,
                                       polarity_list)]
    if print_time:
        print '%s %f'%('parse time',
                       time() - start_time)
        start_time = time()


    #If quality filtering is on then some genes may not exist in
    #both files.
    the_keys = set(parsed_files[0]).union(parsed_files[1])
    if return_id.lower() == 'entrez':
        #Only look at items that are in the annotation file
        the_keys = the_keys.intersection(accession_to_entrez)
    combined_data = {}
    the_data_fields = parsed_files[0].values()[0].keys()
    #Merge the results for each field across files into a list
    for the_key in the_keys:
        tmp_dict = combined_data[the_key] = defaultdict(list)
        for data_dict in parsed_files:
            #If quality filtering is on then some genes may not exist in
            #both files.
            try:
                data_dict = data_dict[the_key] 
            except:
                continue
            [tmp_dict[the_field].append(data_dict[the_field])
             for the_field in the_data_fields]
    if print_time:
        print '%s %f'%('combine time',
                       time() - start_time)
        start_time = time()

    #Collapse the lists for each field for each probe.
    [collapse_fields(v) for v in combined_data.values()]
    if print_time:
        print '%s %f'%('collapse time',
                       time() - start_time)
        start_time = time()
    if return_id.lower() == 'entrez':
        #Now get the entrez ids for inserting into the database.
        the_entrez_ids = set([accession_to_entrez[the_key] for the_key in the_keys])
        #Now need to collapse to entrez gene ids
        entrez_data = dict([(k, defaultdict(list))
                            for k in the_entrez_ids])
        for the_id, the_dict in combined_data.items():
            tmp_dict = entrez_data[accession_to_entrez[the_id]]
            [tmp_dict[k].append(v)
             for k, v in the_dict.items()]
        #Now collapse the instances where the ids are have repeats.
        [collapse_fields(v) for v in entrez_data.values()]
        if print_time:
            print '%s %f'%('collapse to entrez and combine time',
                           time() - start_time)
        combined_data = entrez_data
    #
        
    return combined_data
