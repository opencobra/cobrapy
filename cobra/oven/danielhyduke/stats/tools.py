from copy import deepcopy
from numpy import mean
from cobra.oven.danielhyduke.tools import log_function
from cobra.stats.stats import combine_p_values, error_weighted
def collapse_fields(data_dict, quantitative_fields=['intensity_1',
                                                    'intensity_2',
                                                    'intensity_1_error',
                                                    'intensity_2_error'],
                    log_fields=['log_ratio', 'log_error'],
                    p_fields=['p_value'], log_base=10, error_weighting=False):
    """Collapses the rows for each field from a feature extraction element
    based on the data type.  Quantative fields are averaged, log_fields are
    converted to linear scale, averaged, and then log10 is taken. p-value
    fields are combined using the default method from cobra.stats.stats.
    combine_p_values.

    data_dict: A dictionary of  numerical data values associated with an element
    id.  string: dict where string is the id and dict has the *fields as keys
    and a list or vector as the values.

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
