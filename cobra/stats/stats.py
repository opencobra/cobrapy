#cobra.stats.py: A place to warehouse generic statistical functions and
#to abstract away more of the interface to R and other statistical tools
from numpy import mean, array, log10, isinf, sign, sum
from scipy import log
from scipy.stats import t, chi2, norm

def p_adjust(the_p_values, correction_method='bh'):
    """Adjusts the p-values in a list for multiple hypothesis testing.

    the_p_values: a list of p-values

    correction_method: String.  'bh'|'by'|'bonferroni'.  'bh' for
    Benjamini & Hochberg, 'by' for Benjamini & Yekuieli, and
    'bonferroni' for Bonferroni.

    NOTE: 'bh' and 'by' require R and rpy2 installed
    
    """
    correction_method = correction_method.lower()
    if correction_method == 'bh' or correction_method == 'by':
        correction_method = correction_method.upper()
    if correction_method == 'bonferroni':
        the_p_values = [min(1.0, x * float(len(the_p_values))) for x in the_p_values]
    else:
        #BUG: This may need to be reworked for the jython set up.
        try:
            from cobra.stats import r
        except Exception as e:
            print "Couldn't load cobra.r perhaps no rpy2 modules:%s'"%e
        the_p_values = [float(x)
                   for x in list(array(r.adjust_p(array(the_p_values),
                                                  correction_method)))]
    return the_p_values

def combine_p_values(the_p_values, method='z'):
    """Combines p-values from repeat measurements into a single
    p-value.

    the_p_values: a list of p-values.

    method: String. 'z'|'fisher'.  'z' for using the weighted z-score.
    'fisher' for using fisher's combined probability test.
    """
    if len(the_p_values) == 1 or sum(the_p_values) == 0:
        combined_p_value = the_p_values[0]
    elif method.lower() == 'z':
        #combine p-values using weighted z-score.  To not deal with inifinite
        #values replace 
        new_p = []
        for the_p in the_p_values:
            the_p = norm.ppf(1.-the_p)
            if isinf(the_p):
                if sign(the_p) < 0:
                    the_p = -7.
                else:
                    the_p = 7.
        new_p.append(the_p)
        combined_p_value = norm.sf(sum(new_p) / len(new_p)**0.5)
    elif method.lower() == 'fisher':
        combined_p_value = 1-chi2.cdf(-2*sum(map(log,
                                                    the_p_values)),
                                         2*len(the_p_values))


    return combined_p_value

def error_weighted(the_means, the_stds):
    """Calculate the error-weighted mean and standard deviation.

    the_means: A list or numpy array of floats.

    the_stds: A list or numpy array of floats.
    
    """
    #Allow the input to be a list or an array
    mean_list = array(the_means).tolist()
    std_list = array(the_stds).tolist()
    if len(mean_list) == 1:
        return (mean_list[0], std_list[0])
    the_variances = array(map(float, the_stds))**2
    weighted_std = (1/sum(1/the_variances))**0.5
    weighted_mean = sum(array(the_means)/the_variances) / sum(1/the_variances)
    return (weighted_mean, weighted_std)
