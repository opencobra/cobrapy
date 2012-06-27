#tm4.py: A set of modules for parsing data into / out of tm4 (tm4.org) file formats
#(such as mev) provided by JCVI.
#
#We assume that for the signal the statistics are calculated from the number of
#pixel (SA)s and for the background the signal is calculated by the pixels surrounding
#the spot drawn from an annulus with area =  pi*radius**2 - SA.  Assume spot
#radius is 10.
#
#BUG: FUNCTIONS THAT RUN THROUGH R ARE NOT COMPATIBLE WITH JYTHON
#YET.  LOOK FOR A JAVA/R CONNECTOR TO ALLEVIATE THIS.
#
#TODO: CREATE OBJECTS FOR THE DATA BECAUSE THE DICTS WILL BREAK SOON.
#
#Spotfinder Flags (This should be the default for the arrays)
#  C: Spot Area > 50 pixels
#  B: 30 pixels =< Spot Area =< 50 pixels
#  A: Spot Area < 30 pixels 
#BAD FLAGS

#  S: Saturated spot
#  U: User flagged
#  X: QC says bad
#  Y: Background > Signal
#  Z: Spot not detected
######
#mev column headers
#  UID	Unique identifier for this spot
#  IA	Intensity value in channel A (Background Subtracted: Integrated Intensity A - BkgA)
#  IB	Intensity value in channel B (Background Subtracted: Integrated Intensity B - BkgB)
#  R	Row (slide row)
#  C	Column (slide column)
#  MR	Meta-row (block row)
#  MC	Meta-column (block column)
#  SR	Sub-row
#  SC	Sub-column
#  FlagA	TIGR Spotfinder flag value in channel A
#  FlagB	TIGR Spotfinder flag value in channel B
#  SA	Actual spot area (in pixels)
#  SF	Saturation factor
#  QC	Cumulative quality control score
#  QCA	Quality control score in channel A
#  QCB	Quality control score in channel B
#  BkgA	Background value in channel A (Total intensity) = MedBkgA*SA
#  BkgB	Background value in channel B (Total intensity) = MedBkgB*SA
#  SDA	Standard deviation for spot pixels in channel A
#  SDB	Standard deviation for spot pixels in channel B
#  SDBkgA Standard deviation of the background value in channel A
#  SDBkgB Standard deviation of the background value in channel B
#  MedA	Median intensity value in channel A
#  MedB	Median intensity value in channel B
#  MNA	Mean intensity value in channel A ( IA / ( SA * SF ) )
#  MNB	Mean intensity value in channel B ( IB / ( SA * SF ) )
#  MedBkgA Median background in channel A
#  MedBkgB Median background in channel B.
#  X	X coordinate of the spot cell rectangle
#  Y	Y coordinate of the spot cell rectangle
######NOTE THAT THE P-values are not accurate.  They are calculated by
####### gsl_ran_tdist_pdf (t, df) where t-statistics is
####### t = (spotmean-bkgmean)/sqrt(spotmean*bkgmean)
####### df = spotmean+bkgmean-2.
#  PValueA P-value in channel A
#  PValueB P-value in channel B
#  OligoID Probe ID
#end mev column headers
######
#Use this threshold and the spot finder flags for automatically marking a
#spot as not detected.
#
from warnings import warn
warn('DO NOT USE cobra.io.tm4, it is not ready for general use')
#raise Exception('cobra.io.tm4 is not ready for general use')
bad_flags = {'S': 'Saturated', 'U': 'User Flagged', 'X': 'Failed QC',
             'Y': 'Background > Signal' ,'Z': 'Spot Not Detected'}
import re
from scipy.stats import t
from copy import deepcopy
from numpy import array
from math import pi
from cobra.tools import log_function, exp_function
from time import time
from warnings import warn
try:
    from cobra.db_tools.tm4 import mev_pid_to_locus_id, mev_uid_to_locus_id, \
         mev_uid_to_dbid
         
    #Note the dictionaries for this conversion should be supplied to the function
    #instead of calling the database here.
except:
    warn("cobra.io.tm4 is not yet ready for usage.")
def parse_mev_file(file_name, file_data, file_name_delimiter='_', file_suffix='mev',
                   calculate_p_values=True, qc_threshold=0.0, maximum_radius=10.0, db_cursor=None):
    """Returns a dict of a parsed mev file.  file_name is the name of the file, excluding the directory;
    file_data is a list of the data in the file obtained by open(file_name).readlines()
    POTENTIAL BUG: assumes that file_name follows a specific format
    BUG: Only set to parse LPM medium
    """
    vs_test = file_name_delimiter + 'vs' + file_name_delimiter
    if vs_test not in file_name:
        raise Exception('%s does not conform to naming standards'%file_name)
    data_dict = {}
    the_delimiter = '\t'
    the_newline = '\r\n'
    #In the title, the sample that comes first in the title is labeled with cy5.  And cy5 is usually channel B
    the_channel_ids = ['A', 'B']
    the_media = ['MgM', 'LB', 'LPM']
    #Parse the relevant experimental configuration data from the file name
    #TODO: Describe the filename structure.
    channel2_strain, channel2_condition, tmp_crap, channel1_strain, channel1_condition, experiment_id = file_name.rstrip('.'+file_suffix).split(file_name_delimiter)[:6]
    channel1_medium = channel2_medium = 'LPM'
    for the_medium in the_media:
        if the_medium in channel1_condition:
            channel1_medium = the_medium
            break
    for the_medium in the_media:
        if the_medium in channel2_condition:
            channel2_medium = the_medium
            break
    channel1_time = channel1_condition.lstrip(channel1_medium).rstrip('hr')
    channel2_time = channel2_condition.lstrip(channel2_medium).rstrip('hr')
    file_dict = {}
    file_dict['A'] = {'Strain': channel1_strain, 'Medium': channel1_medium, 'Time': channel1_time}
    file_dict['B'] = {'Strain': channel2_strain, 'Medium': channel2_medium, 'Time': channel2_time}
    #Strip out the comment lines and strip newlines and split
    file_data = [x.rstrip(the_newline).split(the_delimiter)
                 for x in file_data if not x.startswith('#') and not x.startswith('"')]
    #TODO:  The remainder can be moved to a separate function that can be used with any
    #mev file
    file_dict['Data'] = parse_mev_data(file_data, the_channel_ids,
                                       calculate_p_values=calculate_p_values,
                                       qc_threshold=qc_threshold, maximum_radius=maximum_radius,
                                       db_cursor=db_cursor)
    return(file_dict)

def parse_mev_data(file_data, the_channel_ids,  calculate_p_values=True,
                   qc_threshold=0.0, maximum_radius=10.0, db_cursor=None,
                   print_time=False):
    """Returns a dict of the parsed data rows from an mev file.
    qc_threshold is useful when wanting to select a subset of measurements passing the threshold, however,
    it should be 0 when loading into data warehousing software.
    
    TODO: Consider accessing a function that translates pid into locus id and appends that into the data dict for
    the spot.

    """
    #Some of the 'raw' stats in SpotFinder are actually calculated and this was not stated in the manual.
    #See correspondence with SpotFinder Engineer (Sharov, Vasily A. <VSharov@jcvi.org>):
    the_header = file_data.pop(0)
    #Some of the older mev files used OligoID instead of locus.
    if 'OligoID' in the_header:
        pid_type = 'oligo'
        if 'Locus' in the_header:
            the_header[the_header.index('Locus')] = 'Gene Name'
        probe_index = the_header.index('OligoID')
    elif 'Locus' in the_header:
        pid_type = 'locus'
        probe_index = the_header.index('Locus')
    else:
        pid_type = 'uid'
        probe_index = the_header.index('UID')

    cell_area = pi * (maximum_radius)**2
    the_data = {}
    the_pids = []
    if print_time:
        start_time = time()
    for the_line in file_data:
        the_spot = {}
        the_spot['pid'] = the_line[probe_index]
        spot_pixels = the_spot['pixels'] = float(the_line[the_header.index('SA')])
        saturation_factor = the_spot['saturation factor'] = float(the_line[the_header.index('SF')])
        the_spot['QC'] = float(the_line[ the_header.index('QC')])
        #The way that the stats were calculated involved setting the background area to the
        #spot area.
        #basically, there's not enough information to make an informed call on statistics
        #calculations for a single channel as SpotFinder doesn't provide enough
        #raw data.
        background_pixels = spot_pixels = spot_pixels*the_spot['saturation factor']


        for the_channel in the_channel_ids:
            #THE P-values should be recalculated here because the ones
            #presented in spotfinder are not accurate.
            #Recalculating the p-values here will reduce overhead in downstream
            #calculations
            the_measurement = {}
            the_measurement['Total Intensity'] = float(the_line[the_header.index('I' + the_channel)])
            the_mean = the_measurement['Mean'] = float(the_line[the_header.index('MN' + the_channel)])
            the_measurement['Median'] = float(the_line[the_header.index('Med'+the_channel)])
            the_std_dev = the_measurement['Std Dev'] = float(the_line[the_header.index('SD' + the_channel)])
            the_measurement['Background Total'] = float(the_line[the_header.index('Bkg' + the_channel)])
            background_mean = the_measurement['Background Mean'] = float(the_line[the_header.index('MedBkg' + the_channel)])
            background_std_dev = the_measurement['Background Std Dev'] = float(the_line[the_header.index('SDBkg' + the_channel)])            
            the_qc = the_measurement['QC'] = float(the_line[the_header.index('QC' + the_channel)])
            the_flag = the_measurement['Flag'] = the_line[the_header.index('Flag' + the_channel)]
           
            if the_flag in bad_flags.keys() or the_qc < qc_threshold:
                the_p = 1 #For spots that have a bad flag set p = 1
            else:
                #Since the background area is adjusted to match the signal area
                #we can change the t-statistic to the one for equal variance
                #pooled_std_dev = ((the_std_dev**2 + background_std_dev**2)/2.)**0.5
                #the_t = (the_mean - background_mean)/(pooled_std_dev*(2/spot_pixels) **0.5)
                #degrees_of_freedom = background_pixels + spot_pixels - 2
                #warn("Doing unequal variance test")
                #However, I don't really know if we can assume that the signal and
                #background should have equal variances as background should be dominated
                #by nonbiological.
                deviation_estimate = (the_std_dev**2/spot_pixels + background_std_dev**2/background_pixels)**0.5
                the_t = (the_mean - background_mean) / deviation_estimate
                degrees_of_freedom = (the_std_dev**2/spot_pixels + background_std_dev**2/background_pixels)**2/\
                                     ((the_std_dev**2/spot_pixels)**2/(spot_pixels-1) +\
                                      (background_std_dev**2/background_pixels)**2/background_pixels-1)
                the_p = 1 - t.cdf(the_t, degrees_of_freedom)
            the_measurement['p'] = the_p
            the_spot[the_channel] = the_measurement
        the_uid = int(the_line[the_header.index('UID')])
        the_data[the_uid] = the_spot
        the_pids.append(the_spot['pid'])

    if print_time:
        print 'Parsing data took %1.2f minutes'%((time()-start_time)/60)
        start_time = time()
    if db_cursor:
        #Only convert the pid if it is an oligo id.
         if pid_type == 'oligo': 
             pid_to_locus = mev_pid_to_locus_id(db_cursor, set(the_pids))
             [v.update({'locus': pid_to_locus[v['pid']]})
              for v in the_data.itervalues()]
         elif pid_type == 'locus':
             [v.update({'locus': v['pid']}) for v in the_data.itervalues()]
         elif pid_type == 'uid':
             uid_to_locus = mev_uid_to_locus_id(db_cursor, set(map(int, the_pids)))
             [v.update({'locus': uid_to_locus[int(v['pid'])]})
              for v in the_data.itervalues()]
         else:
             raise Exception("Unkown pid_type %s"%repr(pid_type))
         #Get the database unique id for the spot to aid in loading.
         #TODO: The repeated calls makes this slow.  Change this to a single
         #call and allow mev_uid_pid_to_dbid to use prepared statements
         uid_to_dbid = mev_uid_to_dbid(db_cursor, the_data.keys())
         [v.update({'id': uid_to_dbid[k]})
          for k,v in the_data.iteritems()]
         #[v.update({'id': mev_uid_pid_to_dbid(k, v['pid'], db_cursor)})
         #for k, v in the_data.items()]
    if print_time:
        print 'Getting ids took %1.2f minutes'%((time()-start_time)/60)

    return(the_data)
#

