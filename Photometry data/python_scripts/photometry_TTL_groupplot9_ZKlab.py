# -*- coding: utf-8 -*-
"""
Created on Tue May 26 23:03:32 2020


Processing order for TTL lick fiber photometry data
1. photometry_TTL_raw_PSTH.py (or photometry_TTL_raw_PSTH_ds.py for downsampling dFF)
    1-1. photometry_TTL_groupplot_PSTH.py 
2. photometry_TTL_findbdecay.py  
3. photometry_TTL_groupplot.py(or photometry_TTL_groupplot_ds.py for downsampling dFF)
      ********** this script  


Input:(dFFfall added) regular .pkl file saved from photometry_TTL_findbdecay.py 
        
Process
    1)   Find time spent to return to baseline and compare between groups: 
            a) duration 
            b) dFF change
            c) number of licks
            d) dFF change per lick
    
        i) Total, ii) first bout, iii) excluding first bout

    2)   For first 3 mins from first bout, compare between groups: 
            a) dFF change
            b) number of licks
            c) dFF change per lick

        i) Total, ii) first bout, iii) excluding first bout

    3)   mean lick rate + dFF, compare between groups
    
    4)   mean cumulative licks + dFF, compare between groups
    
    5)   dFF change over time 
    
    6)   dFF changes over bout

    
Output: Plots, statistical values 


@author: heeun
"""

import os
import sys
import glob
import pickle as pk
import numpy as np
import pandas as pd
import math
#from datetime import date
import scipy.stats as ss 
from scipy.integrate import simps
from scipy import signal
from sklearn.neighbors import KDTree
import seaborn as sns
# import scikit_posthocs as sp
import matplotlib.pyplot as plt  # standard Python plotting library
import matplotlib
matplotlib.rcParams['font.size'] = 18 

print("Enter the full path for directory with .pkl files: ")
path      = input()

   
analysis_path = path + "/analysis" 
if not os.path.exists(analysis_path):
    os.makedirs(analysis_path)   

###  Smoothing 
def smooth (x, window_len, window='hanning'):
        s = np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
            w = np.ones(window_len,'d')
        else:  
            w = eval('np.'+window+'(window_len)')
        y = np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]

def lick_rate_by_sec (lick_on, time_sec):

    lick_on_ceil = np.ceil(lick_on)
    lick_on_ceil = lick_on_ceil.astype(int)
        
    lick_rate_sec_ = np.histogram(lick_on_ceil, bins = time_sec)   
    lick_rate_sec  = np.concatenate([[0], lick_rate_sec_[0]])
    lick_rate_cum  = np.cumsum(lick_rate_sec)

    return (lick_rate_sec, lick_rate_cum)

palette1       = ["black", "red"]
palette2       = ["black", "dimgray", "red", "maroon"]

# palette_mean      = ["dimgray", "maroon"]
palette_mean_GCaMP      = ["darkgreen", "maroon"]
palette_ste_GCaMP       = ["forestgreen", "pink"]

palette_mean_lick       = ["dimgrey", "dimgrey"]
palette_ste_lick        = ["darkgrey", "darkgrey"]


frame_rate      = 1017.2527
pre_length      = 10     ## in minutes
length_cutoff   = 90     ## in minutes 

lickvol         = 1.205697403329695/1000  # uL/lick value converted to mL/lick


dsfactor            = 50 
new_length_cutoff   = 90 ## in minutes


all_files                = glob.glob(path + "/*.pkl")
file_list                = []
minlendFF_firstlick_list = []
minlendFF_firstbout_list = []

data = pd.DataFrame(columns = ["expID", "mouseID", "expdate", 
                               "dFF_ds_sm", 
                               "dFF_ds_sm_firstlick", 
                               "dFF_ds_sm_firstbout", 
                               
                               "lickrate_sec", "lickrate_min", "lickrate_sec_sm60",
                               "lickon_firstlick", "lickrate_sec_firstlick", 
                               "lickrate_min_firstlick", "lickrate_sec_sm60_firstlick", 
                               
                               "lickon_firstbout", "lickrate_sec_firstbout", 
                               "lickrate_min_firstbout", "lickrate_sec_sm60_firstbout",

                               "time_sec_ds", "time_sec_bin", "time_min_bin",
                               "time_sec_ds_firstlick", "time_sec_ds_firstbout",
                                          
                               "lick_on", "lick_off", "lick_x", "lick_y",
                               "bout_on", "bout_off", "bout_x", "bout_y",
                               "bout_on_ind", "bout_off_ind", "bout_on_dFF", "bout_off_dFF", 
                               "bout_size", "bout_length", "bout_number",
                               "bout_AUC", "bout_AUCperlick", "bout_meandFF", "all_lick_meandFF",
                               "all_bout_meandFF", "all_interbout_meandFF",
                               "interbout_size", "interbout_AUC", "interbout_AUCperlick", "interbout_meandFF",
                                
                               "cumlick_sec", "cumbout_sec", 
                               "cumlick_sec_firstlick", "cumbout_sec_firstlick", 
                               "cumlick_sec_firstbout", "cumbout_sec_firstbout", 
                               "bout_on_PSTH"])
    
data["dFF_ds_sm"]                   = data["dFF_ds_sm"].astype(object) 
data["dFF_ds_sm_firstlick"]         = data["dFF_ds_sm_firstlick"].astype(object)

data["dFF_ds_sm_firstbout"]         = data["dFF_ds_sm_firstbout"].astype(object)

data["lickrate_sec"]                = data["lickrate_sec"].astype(object)
data["lickrate_min"]                = data["lickrate_min"].astype(object)
data["lickrate_sec_sm60"]           = data["lickrate_sec_sm60"].astype(object)

data["lickon_firstlick"]            = data["lickon_firstlick"].astype(object)
data["lickrate_min_firstlick"]      = data["lickrate_min_firstlick"].astype(object)
data["lickrate_sec_sm60_firstlick"] = data["lickrate_sec_sm60_firstlick"].astype(object)
   
data["lickon_firstbout"]            = data["lickon_firstbout"].astype(object)
data["lickrate_min_firstbout"]      = data["lickrate_min_firstbout"].astype(object)
data["lickrate_sec_sm60_firstbout"] = data["lickrate_sec_sm60_firstbout"].astype(object)

data["time_sec_ds"]           = data["time_sec_ds"].astype(object) 
data["time_sec_ds_firstlick"] = data["time_sec_ds_firstlick"].astype(object)
data["time_sec_ds_firstbout"] = data["time_sec_ds_firstbout"].astype(object)

data["lick_on"]           = data["lick_on"].astype(object) 
data["lick_off"]          = data["lick_off"].astype(object) 
data["lick_x"]            = data["lick_x"].astype(object) 
data["lick_y"]            = data["lick_y"].astype(object) 
data["bout_on"]           = data["bout_on"].astype(object) 
data["bout_off"]          = data["bout_off"].astype(object) 
data["bout_x"]            = data["bout_x"].astype(object) 
data["bout_y"]            = data["bout_y"].astype(object) 
data["bout_size"]         = data["bout_size"].astype(object)
data["bout_length"]       = data["bout_length"].astype(object)

data["bout_number"]       = data["bout_number"].astype(object)
data["bout_AUC"]          = data["bout_AUC"].astype(object)
data["bout_AUCperlick"]   = data["bout_AUCperlick"].astype(object)
data["bout_meandFF"]      = data["bout_meandFF"].astype(object)

data["interbout_AUC"]          = data["interbout_AUC"].astype(object)
data["interbout_AUCperlick"]   = data["interbout_AUCperlick"].astype(object)
data["interbout_meandFF"]      = data["interbout_meandFF"].astype(object)

data["bout_on_ind"]         = data["bout_on_ind"].astype(object)  
data["bout_off_ind"]        = data["bout_off_ind"].astype(object)       
data["bout_on_dFF"]         = data["bout_on_dFF"].astype(object)     
data["bout_off_dFF"]        = data["bout_off_dFF"].astype(object)   

data["cumlick_sec"]       = data["cumlick_sec"].astype(object)
data["cumbout_sec"]       = data["cumbout_sec"].astype(object)

data["cumlick_sec_firstlick"]       = data["cumlick_sec_firstlick"].astype(object)
data["cumbout_sec_firstlick"]       = data["cumbout_sec_firstlick"].astype(object)

data["cumlick_sec_firstbout"]       = data["cumlick_sec_firstbout"].astype(object)
data["cumbout_sec_firstbout"]       = data["cumbout_sec_firstbout"].astype(object)

data["bout_on_PSTH"]        = data["bout_on_PSTH"].astype(object)

time_min_bin              = np.arange(-pre_length*60, length_cutoff*60+1, 60)
time_sec_bin              = np.arange(-pre_length*60, length_cutoff*60+1, 1)

#%% Go through each file in the path to save the data into a single dataframe
tree           = {}
tree_firstlick = {}
tree_firstbout = {}
expID_list     = []

for i, file_path in enumerate (all_files):
    expID            = file_path[(len(path)+1):-4] ## previous called as "file_name_cut
    file_list.append(expID)  
    inputfile        = open(file_path, 'rb')
    tempdata         = pk.load(inputfile)
    inputfile.close()
    
    if i == 0:
        stimulus  = tempdata["stimulus"]
        fluid     = tempdata["LICK_fluid"]
    data.at[i, "expID"]         = expID
    expID_list.append(expID)
    
    data.at[i, "mouseID"]       = tempdata["mouseID"]
    data.at[i, "expdate"]       = tempdata["date"]
          
    data.at[i, "lick_on"]       = tempdata["LICK_on"]
    data.at[i, "lick_off"]      = tempdata["LICK_off"]
    data.at[i, "lick_x"]        = tempdata["LICK_x"]
    data.at[i, "lick_y"]        = tempdata["LICK_y"]    
    data.at[i, "bout_on"]       = tempdata["BOUT_on"]
    data.at[i, "bout_off"]      = tempdata["BOUT_off"]
    data.at[i, "bout_x"]        = tempdata["BOUT_x"]
    data.at[i, "bout_y"]        = tempdata["BOUT_y"] 
    
    
    temp_bout_size_nan          = [np.nan]*(50 - len(tempdata["BOUT_size"]))
    data.at[i, "bout_size"]     = np.array(list(tempdata["BOUT_size"]) + temp_bout_size_nan)
    
    temp_bout_length            = list(np.array(data.at[i, "bout_off"])  - np.array(data.at[i, "bout_on"]))
    temp_bout_length_nan        = [np.nan]*(50 - len(tempdata["BOUT_size"]))
    data.at[i, "bout_length"]   = np.array(temp_bout_length + temp_bout_length_nan)

    data.at[i, "bout_AUC"]              = tempdata["df_bout_dFF"]["bout_AUC"]
    data.at[i, "bout_size"]             = tempdata["df_bout_dFF"]["bout_size"]
    data.at[i, "bout_AUCperlick"]       = tempdata["df_bout_dFF"]["bout_AUCperlick"]
    data.at[i, "bout_meandFF"]          = tempdata["df_bout_dFF"]["bout_meandFF"]
    
    data.at[i, "interbout_AUC"]         = tempdata["df_bout_dFF"]["interbout_AUC"]
    data.at[i, "interbout_size"]        = tempdata["df_bout_dFF"]["interbout_size"]
    data.at[i, "interbout_AUCperlick"]  = tempdata["df_bout_dFF"]["interbout_AUCperlick"]
    data.at[i, "interbout_meandFF"]     = tempdata["df_bout_dFF"]["interbout_meandFF"]

    lickrate_sec, lickrate_cum          = lick_rate_by_sec(tempdata["LICK_on"], time_sec_bin)      
    lickrate_min, dump                  = lick_rate_by_sec(tempdata["LICK_on"], time_min_bin) 
    
    data.at[i, "lickrate_sec"]          = lickrate_sec
    data.at[i, "cumlick_sec"]           = lickrate_cum+1
    dump, data.at[i, "cumbout_sec"]     = lick_rate_by_sec(tempdata["BOUT_on"], time_sec_bin)
         
    data.at[i, "lickrate_min"]          = lickrate_min
    data.at[i, "lickrate_sec_sm60"]     = smooth(lickrate_sec, window_len = 60)
    
    ## downsample all time series by a factor of 50 (1000Hz >> 20Hz. 1ms >> 50ms exposure time )
    data.at[i, "dFF_ds_sm"]        = tempdata["dFF_ds_sm"]
    data.at[i, "time_sec_ds"]      = tempdata["time_sec_ds"]
    
    time_sec_dims                  = np.expand_dims(data.at[i, "time_sec_ds"], axis = 1)
    tree[i]                        = KDTree(time_sec_dims, leaf_size = 2)  
 

    
    ############################################################################
    ### lick rate: time_min. aligned to first lick
    ############################################################################ 
    firstlick                           = data.at[i, "lick_on"][0]
    data.at[i, "lickon_firstlick"]      = data.at[i, "lick_on"] - firstlick
    data.at[i, "time_sec_ds_firstlick"] = data.at[i, "time_sec_ds"] - firstlick

         
    firstlick                                      = firstlick.reshape(1, -1)
    firstlick_ind                                  = tree[i].query(firstlick)[1][0][0]
    
    timezero_ind                                   = tree[i].query([[0]])[1][0][0]
    diff_ind                                       = firstlick_ind - timezero_ind 
    data.at[i, "dFF_ds_sm_firstlick"]              = data.at[i, "dFF_ds_sm"][diff_ind:] 
    data.at[i, "time_sec_ds_firstlick"]            = data.at[i, "time_sec_ds_firstlick"][diff_ind:]

    minlendFF_temp   = len(data.at[i, "time_sec_ds_firstlick"])
    minlendFF_firstlick_list.append(minlendFF_temp)
    
    lickrate_sec_firstlick, lickrate_cum_firstlick = lick_rate_by_sec(data.at[i, "lickon_firstlick"], time_sec_bin)
    data.at[i, "lickrate_sec_firstlick"]           = lickrate_sec_firstlick
    data.at[i, "cumlick_sec_firstlick"]            = lickrate_cum_firstlick + 1
    dump, cumbout_sec_firstlick                    = lick_rate_by_sec(data.at[i, "bout_on"] - firstlick, time_sec_bin)
    data.at[i, "cumbout_sec_firstlick"]            = cumbout_sec_firstlick + 1

    lickrate_min_firstlick, dump     = lick_rate_by_sec(data.at[i, "lickon_firstlick"], time_min_bin)    
 
    data.at[i, "lickrate_min_firstlick"] = lickrate_min_firstlick
 
    ############################################################################     
    ### lick rate: time_sec (moving average of 60 sec). aligned to first lick
    ############################################################################       
    lickrate_sec_firstlick, dump                     = lick_rate_by_sec(data.at[i, "lickon_firstlick"], time_sec_bin)      
    data.at[i, "lickrate_sec_sm60_firstlick"]        = smooth(lickrate_sec_firstlick, window_len= 60)
    lickrate_sec_firstlick_nan                       = [np.nan if x== 0 else x for x in data.at[i, "lickrate_sec_sm60_firstlick"]]      


    ############################################################################
    ### lick rate: time_min. aligned to first bout
    ############################################################################   
    firstbout                           = data.at[i, "bout_on"][0]
    data.at[i, "lickon_firstbout"]      = data.at[i, "lick_on"] - firstbout
    data.at[i, "time_sec_ds_firstbout"] = data.at[i, "time_sec_ds"] - firstbout
    
    firstbout                           = firstbout.reshape(1, -1)
    firstbout_ind                       = tree[i].query(firstbout)[1][0][0]
    
    timezero_ind                        = tree[i].query([[0]])[1][0][0]
    diff_ind                            = firstbout_ind - timezero_ind 
    data.at[i, "dFF_ds_sm_firstbout"]   = data.at[i, "dFF_ds_sm"][diff_ind:] 
    data.at[i, "time_sec_ds_firstbout"] = data.at[i, "time_sec_ds_firstbout"][diff_ind:]

    minlendFF_temp   = len(data.at[i, "time_sec_ds_firstbout"])
    minlendFF_firstbout_list.append(minlendFF_temp)
    
    lickrate_sec_firstbout, lickrate_cum_firstbout = lick_rate_by_sec(data.at[i, "lickon_firstbout"], time_sec_bin)
    data.at[i, "lickrate_sec_firstbout"]           = lickrate_sec_firstbout
    data.at[i, "cumlick_sec_firstbout"]            = lickrate_cum_firstbout + 1
    dump, cumbout_sec_firstbout                    = lick_rate_by_sec(data.at[i, "bout_on"] - firstbout, time_sec_bin)
    data.at[i, "cumbout_sec_firstbout"]            = cumbout_sec_firstbout + 1
    data.at[i, "lickrate_min_firstbout"], dump     = lick_rate_by_sec(data.at[i, "lickon_firstbout"], time_min_bin)
     
    lickrate_min_firstbout, dump     = lick_rate_by_sec(data.at[i, "lickon_firstbout"], time_min_bin)      
    lickrate_min_firstbout_nan       = [np.nan if x== 0 else x for x in lickrate_min_firstbout]
    data.at[i, "lickrate_min_firstbout"] = lickrate_min_firstbout

    ############################################################################     
    ### lick rate: time_sec (moving average of 60 sec). aligned to first lick
    ############################################################################  
    lickrate_sec_firstbout, dump                     = lick_rate_by_sec(data.at[i, "lickon_firstbout"], time_sec_bin)      
    data.at[i, "lickrate_sec_sm60_firstbout"]        = smooth(lickrate_sec_firstbout, window_len= 60)
    lickrate_sec_firstbout_nan                       = [np.nan if x== 0 else x for x in data.at[i, "lickrate_sec_sm60_firstbout"]]      

    ##############################################################
    ### mean dFF for individual lick events 
    ##############################################################
    df_lick_dFF     = pd.DataFrame(columns = ["lick_meandFF"])
    for k, lick_on in enumerate (data.at[i, "lick_on"]):
        lick_on                                = lick_on.reshape(1, -1)       				
        ind_lick_on                            = tree[i].query(lick_on)[1][0][0]
        df_lick_dFF.at[k, "lick_meandFF"]      = data.at[i, "dFF_ds_sm"][ind_lick_on]
    data.at[i, "all_lick_meandFF"]             = df_lick_dFF["lick_meandFF"].mean()
    
    #%%#############################################################
    ### PSTH for first lick of individual bouts 
    ##############################################################                   
for i, file_path in enumerate (all_files):
    psth_width_sec  = 120                                         
    psth_width      = math.floor(frame_rate*psth_width_sec/dsfactor)
    psth_num_points = 2*psth_width+1
    psth_plot_times = np.linspace(-psth_width_sec, psth_width_sec, psth_num_points, endpoint = True, dtype = None)
            
    psth            = np.empty([len(data.at[i,"bout_size"]), psth_num_points])
    psth[:]         = np.nan
    psth_y_min_list = []
    psth_y_max_list = []
    
    for k, bout_on in enumerate (data.at[i, "bout_on"]):
        bout_on     = bout_on.reshape(1, -1)
        ind         = tree[i].query(bout_on)[1][0][0] 

        ## Process without error if PSTH post-width exceeds the recording, fill with nan.
        if ind + psth_width + 1 >= len(data.at[i, "dFF_ds_sm"]):  
            last_post_psth                          = len(data.at[i, "dFF_ds_sm"]) - ind
            psth[k, : psth_width + last_post_psth]  = data.at[i, "dFF_ds_sm"][ind - psth_width : ] - np.nanmedian(data.at[i, "dFF_ds_sm"][ind - psth_width : ind]) 
        else:    
            psth[k, :]  = data.at[i, "dFF_ds_sm"][ind - psth_width : ind + psth_width +1] - np.nanmedian(data.at[i, "dFF_ds_sm"][ind - psth_width : ind])

        psth_y_min_list.append(np.nanmin(psth[k, :]))
        psth_y_max_list.append(np.nanmax(psth[k, :]))
    data.at[i, "bout_on_PSTH"] = pd.DataFrame(psth)
    
#%%#######################################################################
#   if trace length is off after firstbout adjustment, trim the long one
##########################################################################
minlendFF_firstlick  = np.nanmin(np.array(minlendFF_firstlick_list))
minlendFF_firstbout  = np.nanmin(np.array(minlendFF_firstbout_list))

for i in np.arange(len(data)):
    data.at[i, "time_sec_ds_firstlick"]       = data.at[i, "time_sec_ds_firstlick"][:minlendFF_firstlick]
    data.at[i, "dFF_ds_sm_firstlick"]         = data.at[i, "dFF_ds_sm_firstlick"][:minlendFF_firstlick]
    
    data.at[i, "time_sec_ds_firstbout"]       = data.at[i, "time_sec_ds_firstbout"][:minlendFF_firstbout]
    data.at[i, "dFF_ds_sm_firstbout"]         = data.at[i, "dFF_ds_sm_firstbout"][:minlendFF_firstbout]

time_sec              = data.at[i, "time_sec_ds"]  
time_sec_firstlick    = data.at[i, "time_sec_ds_firstlick"][:minlendFF_firstlick] 
time_sec_firstlick    = time_sec_firstlick[time_sec_firstlick/60 <= length_cutoff]
time_sec_firstbout    = data.at[i, "time_sec_ds_firstbout"][:minlendFF_firstbout]  
print ("files analyzed are: %s" %file_list)       
    

#%%###################################################################
#  get the smallest number of bouts  
######################################################################
bout_number_list = []
for i in range(len(data)):
    bout_number_list.append(len(data.at[i, "bout_on"]))
         
min_num_bouts    = np.nanmin(bout_number_list)
bout_number      = np.arange(1, min_num_bouts+1, 1)
interbout_number = bout_number[:-1]

for i in np.arange(len(data)):
    data.at[i, "bout_AUC"]              = np.array(data.at[i, "bout_AUC"][:min_num_bouts])
    data.at[i, "bout_size"]             = np.array(data.at[i, "bout_size"][:min_num_bouts])
    data.at[i, "bout_AUCperlick"]       = np.array(data.at[i, "bout_AUCperlick"][:min_num_bouts])
    data.at[i, "all_bout_meandFF"]      = data.at[i, "bout_meandFF"].mean()
    data.at[i, "bout_meandFF"]          = np.array(data.at[i, "bout_meandFF"][:min_num_bouts])
    
    data.at[i, "interbout_AUC"]         = np.array(data.at[i, "interbout_AUC"][:min_num_bouts-1])
    data.at[i, "interbout_size"]        = np.array(data.at[i, "interbout_size"][:min_num_bouts-1])
    data.at[i, "interbout_AUCperlick"]  = np.array(data.at[i, "interbout_AUCperlick"][:min_num_bouts-1])
    data.at[i, "all_interbout_meandFF"] = data.at[i, "interbout_meandFF"].mean()
    data.at[i, "interbout_meandFF"]     = np.array(data.at[i, "interbout_meandFF"][:min_num_bouts-1])
    
#%%#############################################################################################
#   PROCESS DATA: "Groupby" by "mouseID", save the mean dFF into the new dataframe
################################################################################################
meandata   = pd.DataFrame({"dFF_ds_sm_firstlick": data.groupby("mouseID")["dFF_ds_sm_firstlick"].apply(np.mean), 
                           "dFF_ds_sm_firstbout": data.groupby("mouseID")["dFF_ds_sm_firstbout"].apply(np.mean), 

                           "lickrate_min_firstlick": data.groupby("mouseID")["lickrate_min_firstlick"].apply(np.mean), 
                           "lickrate_sec_sm60_firstlick": data.groupby("mouseID")["lickrate_sec_sm60_firstlick"].apply(np.mean), 
                           "cumlick_sec_firstlick": data.groupby("mouseID")["cumlick_sec_firstlick"].apply(np.mean), 
                           "cumbout_sec_firstlick": data.groupby("mouseID")["cumbout_sec_firstlick"].apply(np.mean), 

                           "lickrate_min_firstbout": data.groupby("mouseID")["lickrate_min_firstbout"].apply(np.mean), 
                           "lickrate_sec_sm60_firstbout": data.groupby("mouseID")["lickrate_sec_sm60_firstbout"].apply(np.mean), 
                           "cumlick_sec_firstbout": data.groupby("mouseID")["cumlick_sec_firstbout"].apply(np.mean), 
                           "cumbout_sec_firstbout": data.groupby("mouseID")["cumbout_sec_firstbout"].apply(np.mean), 

                            "bout_AUC": data.groupby("mouseID")["bout_AUC"].apply(np.mean),
                            "bout_size": data.groupby("mouseID")["bout_size"].apply(np.mean),
                            "bout_AUCperlick": data.groupby("mouseID")["bout_AUCperlick"].apply(np.mean),
                            "bout_meandFF": data.groupby("mouseID")["bout_meandFF"].apply(np.mean),

                            "interbout_AUC": data.groupby("mouseID")["interbout_AUC"].apply(np.mean),
                            "interbout_size": data.groupby("mouseID")["interbout_size"].apply(np.mean),
                            "interbout_AUCperlick": data.groupby("mouseID")["interbout_AUCperlick"].apply(np.mean),
                            "interbout_meandFF": data.groupby("mouseID")["interbout_meandFF"].apply(np.mean),
                           
                            "all_lick_meandFF": data.groupby("mouseID")["all_lick_meandFF"].apply(np.mean),
                            "all_bout_meandFF": data.groupby("mouseID")["all_bout_meandFF"].apply(np.mean),
                            "all_interbout_meandFF": data.groupby("mouseID")["all_interbout_meandFF"].apply(np.mean),
                           
                           }).reset_index()   

time_sec_firstlick     = data.at[0, "time_sec_ds_firstlick"]
time_sec_firstbout     = data.at[0, "time_sec_ds_firstbout"]


#%%##################################################################################
#   Figure 1. lick rate: averaged dFF trace with lick rate
#####################################################################################

fig, ax     = plt.subplots(nrows = 2, ncols = 1, figsize = (8*2, 12))

###############################################################
#   subplot 1-1: mean dFF trace +  lick rate (min)
###############################################################

ax1             = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)      
dFF_values      = meandata["dFF_ds_sm_firstlick"].values
size            = len(meandata)
title           = "lick rate (n=%d)" %(size)        

if size > 1:
    lns1 = ax1.plot(time_sec_firstlick/60, dFF_values.mean(), linewidth = 2, color = palette_mean_GCaMP[0], label = "GCaMP") 
    yerr = dFF_values.std()/math.sqrt(size-1)
    ax1.fill_between(time_sec_firstlick/60, dFF_values.mean() - yerr, dFF_values.mean() + yerr, color = palette_ste_GCaMP[0], alpha = 0.5) # 
else:
    lns1 = ax1.plot(time_sec_firstlick/60, dFF_values.mean(), linewidth = 2, color = palette_mean_GCaMP[0], label = "GCaMP")          

ax1.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray") #
ax1.set_ylabel("dFF (%)", fontsize = 20)
ax1.set_xlabel("time (min)")

ax2                  = ax1.twinx()    
lickrate_values      = meandata["lickrate_min_firstlick"].values
    
if size > 1:
    lns2 = ax2.plot(time_min_bin/60, lickrate_values.mean(), linewidth = 2, color = "black", label = "lick rate") #
    yerr = lickrate_values.std()/math.sqrt(size-1)
    ax2.fill_between(time_min_bin/60, lickrate_values.mean() - yerr, lickrate_values.mean() + yerr, color = "lightgrey", alpha = 0.5) #

else:
    lns2 = ax2.plot(time_min_bin/60, lickrate_values.mean(), linewidth = 2, color = "black", label = "lick rate (min)") #           

ax2.set_ylabel("lick rate (/min)", fontsize = 20)
   
lns     = lns1 + lns2
labels  = [l.get_label() for l in lns]        
#        ax1.legend(lns, labels, loc="lower left", mode = "expand", ncol = 2, bbox_to_anchor = (0, 1.02, 1, 0.2), fontsize = 16)
ax1.legend(lns, labels, loc="upper right", fontsize = 16)
ax1.set_title(title, fontsize = 20)

###############################################################
#   subplot 1-2: mean dFF trace +  lick rate (sec_sm60)
###############################################################   
ax3         = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)          

title       = "lick rate (n=%d)" %(size)      

if size > 1:
    lns3 = ax3.plot(time_sec_firstlick/60, dFF_values.mean(),  linewidth = 2, color = palette_mean_GCaMP[0], label = "GCaMP")
    yerr = dFF_values.std()/math.sqrt(size-1)
    ax3.fill_between(time_sec_firstlick/60, dFF_values.mean() - yerr, dFF_values.mean() + yerr, color = palette_ste_GCaMP[0], alpha = 0.5) ,
else:
    lns3 = ax3.plot(time_sec_firstlick/60, dFF_values.mean(), linewidth = 2, color = palette_mean_GCaMP[0], label = "GCaMP")   

ax3.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
ax3.set_ylabel("dFF (%)", fontsize = 20)
ax3.set_xlabel("time (min)")

ax4                  = ax3.twinx()       
lickrate_values      = meandata["lickrate_sec_sm60_firstlick"].values

if size > 1:
    lns4 = ax4.plot(time_sec_bin/60, lickrate_values.mean(), linewidth = 2, color ="black", label = "lick rate (sm60)") #
    yerr = lickrate_values.std()/math.sqrt(size-1)
    ax4.fill_between(time_sec_bin/60, lickrate_values.mean() - yerr, lickrate_values.mean() + yerr, color = "lightgrey", alpha = 0.5) #
else:
    lns4 = ax4.plot(time_sec_bin/60, lickrate_values.mean(), linewidth = 2, color ="black", label = "lick rate (sm60)") #color = "black",       
ax4.set_ylabel("lick rate sm60 (/sec)", fontsize = 20)

lns     = lns3 + lns4
labels  = [l.get_label() for l in lns]  
ax3.legend(lns, labels, loc="upper right", fontsize = 16)
ax3.set_title(title, fontsize = 20)

suptitle = stimulus + "_" + fluid + "_dFF_lickrate"   
plt.suptitle("%s, first lick = 0" %suptitle, fontsize = 18) 
plt.tight_layout(rect=[0, 0.03, 1, 0.95])           
save_file_path = analysis_path + "/" + suptitle + ".png"
plt.savefig(save_file_path)
save_file_path = analysis_path + "/" + suptitle + ".pdf"
plt.savefig(save_file_path)
plt.show()

#%%##################################################################################
#   Figure 2. cum.licks : averaged dFF trace with cumulative licks and cumulative bouts
#####################################################################################

###############################################################
#   subplot 2-1: mean dFF trace +  cumulative licks
###############################################################
fig, ax     = plt.subplots(nrows = 2, ncols = 1, figsize = (8*2, 12))
        
ax1             = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)      
dFF_values      = meandata["dFF_ds_sm_firstlick"].values
size            = len(meandata)
title           = "cumulative licks (n=%d)" %(size)   

if size > 1:
    lns1 = ax1.plot(time_sec_firstlick/60, dFF_values.mean(), linewidth = 2, color = palette_mean_GCaMP[0], label = "GCaMP")  
    yerr = dFF_values.std()/math.sqrt(size-1)
    ax1.fill_between(time_sec_firstlick/60, dFF_values.mean() - yerr, dFF_values.mean() + yerr, color = palette_ste_GCaMP[0], alpha = 0.5) #
else:
    lns1 = ax1.plot(time_sec_firstlick/60, dFF_values.mean(), linewidth = 2, color = palette_mean_GCaMP[0],label = "GCaMP")        

ax1.axvline (x = 0, linewidth = 2, linestyle = "--")#, color = "gray"
ax1.set_ylabel("dFF (%)", fontsize = 20)
ax1.set_xlabel("time (min)")
   
ax2                 = ax1.twinx()       
cumlick_values      = meandata["cumlick_sec_firstlick"].values

if size > 1:
    lns2 = ax2.plot(time_sec_bin/60, cumlick_values.mean(), linewidth = 2, color = "black", label = "no. licks") # 
    yerr = cumlick_values.std()/math.sqrt(size-1)
    ax2.fill_between(time_sec_bin/60, cumlick_values.mean() - yerr, cumlick_values.mean() + yerr, color = "lightgrey", alpha = 0.5) # 
else:
    lns2 = ax2.plot(time_sec_bin/60, cumlick_values.mean(), linewidth = 2, color = "black", label = "no. licks") #             
ax2.set_ylabel("total # licks", fontsize = 20)
 
lns     = lns1 + lns2
labels  = [l.get_label() for l in lns]        
ax1.set_title(title, fontsize = 20)    

###############################################################
#   subplot 2-2: mean dFF trace +  cumulative bouts
###############################################################    

ax3             = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)      

title           = "cumulative licks (n=%d)" %(size)

if size > 1:
    lns1 = ax3.plot(time_sec_firstlick/60, dFF_values.mean(), color = palette_mean_GCaMP[0], linewidth = 2, label = "GCaMP")
    yerr = dFF_values.std()/math.sqrt(size-1)
    ax3.fill_between(time_sec_firstlick/60, dFF_values.mean() - yerr, dFF_values.mean() + yerr, color = palette_ste_GCaMP[0], alpha = 0.5)
else:
    lns1 = ax3.plot(time_sec_firstlick/60, dFF_values.mean(), color = palette_mean_GCaMP[0], linewidth = 2, label = "GCaMP")            

ax3.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
ax3.set_ylabel("dFF (%)", fontsize = 20)
ax3.set_xlabel("time (min)")
     
ax4                 = ax3.twinx()       
cumbout_values      = meandata["cumbout_sec_firstlick"].values

if size > 1:
    lns2 = ax4.plot(time_sec_bin/60, cumbout_values.mean(), color = "black", linewidth = 2, label = "no. bouts")
    yerr = cumbout_values.std()/math.sqrt(size-1)
    ax4.fill_between(time_sec_bin/60, cumbout_values.mean() - yerr, cumbout_values.mean() + yerr, color = "lightgrey", alpha = 0.5)

else:
    lns2 = ax4.plot(time_sec_bin/60, cumbout_values.mean(), color = "black", linewidth = 2, label = "no. bouts")            

ax4.set_ylabel("total # bouts", fontsize = 20)
 
ax4.set_title(title, fontsize = 20) 
lns     = lns1 + lns2
labels  = [l.get_label() for l in lns]        


suptitle = stimulus + "_" + fluid + "_dFF_totlick"     
plt.suptitle("%s, first lick = 0" %suptitle, fontsize = 18)       
plt.tight_layout(rect=[0, 0.03, 1, 0.95])           
save_file_path = analysis_path + "/" + suptitle + ".png"
plt.savefig(save_file_path)
save_file_path = analysis_path + "/" + suptitle + ".pdf"
plt.savefig(save_file_path)
plt.show()
    




#%%##############################################################
#       Processing for mean PSTH plot 
#################################################################

psth_temp            = {}
psth_mean            = {}
psth_y_min_list      = []
psth_y_max_list      = []


#########################################################################
#   Original "data[bout_on_PSTH]" is a dataframe with 
#    x axis (rows): time series for PSTH, y axis (columns): bout number
#   Now convert it to a dictionary "psth_temp" where each key is a bout number (column = time(sec))
#    and a value is a dataframe with x axis: exp ID, y axis: time series for PSTH. (column = expID)
######################################################################### 
for k in range(min_num_bouts):
    psth_temp[k] = pd.DataFrame([np.nan]*len(psth_plot_times), columns = ["dummy"])
    for i in np.arange(len(data)):
        expID    = data["expID"][i]
        try:
            psth_temp[k][expID]  = pd.Series(data.at[i, "bout_on_PSTH"].loc[k, :].values)
        except:                 
            pass
 
#%%###################################################################################
#   Take the mean of technical replicates from "psth_temp" and save it to a new dataframe "psth_mean"
#    with x axis: mouse ID, y axis: time series for PSTH. (column = mouseID)
####################################################################################          

for k in range(min_num_bouts):
    psth_mean[k] = pd.DataFrame([np.nan]*len(psth_plot_times), columns = ["dummy"])
    mouseID_list = []    
    for expID in expID_list:
        mouseID = expID.split("_")[0]
        if mouseID in mouseID_list:
            continue   # skip to the next loop if this mouse has been updated with replicates already
        mouseID_list.append(mouseID) 
        psth_mousemean_temp = pd.DataFrame(psth_temp[k][expID], columns = [expID])   
        for expID2 in expID_list:
            mouseID2 = expID2.split("_")[0]       
            if (mouseID == mouseID2) and (expID != expID2):
                psth_mousemean_temp[expID2] = pd.Series(psth_temp[k][expID2])    
        psth_mean[k][mouseID]      = pd.Series(psth_mousemean_temp.mean(axis=1))


#%%#############################################################################################
#   Remove the "dummy" column from psth_mean and calculate and add columns for mean and stder
###############################################################################################
for k in range(min_num_bouts):
    try:
        psth_mean[k] = psth_mean[k].drop(columns= ["dummy"])
    except:
        pass    
    psth_mean[k]["mean"]           = np.nanmean(psth_mean[k], axis = 1)      
    if (len(psth_mean[k].columns) <= 2) :
        psth_mean[k]["stderr"] = [0]*len(psth_mean[k]["mean"])     
    else:
        temp_std               = np.nanstd(psth_mean[k], axis = 1)
        psth_mean[k]["stderr"] = temp_std/math.sqrt(len(psth_mean[k].columns)-2)              
    psth_y_min_list.append(np.nanmin(psth_mean[k]["mean"] - psth_mean[k]["stderr"]))  
    psth_y_max_list.append(np.nanmax(psth_mean[k]["mean"] + psth_mean[k]["stderr"]))  

psth_y_min        = np.nanmin(psth_y_min_list)
psth_y_max        = np.nanmax(psth_y_max_list)


# #################################################################################################
# #    SAVE DATA: filter selected columns of "data" dataframe into a new dataframe and save into a .csv file
# #################################################################################################
bout_AUC              = pd.DataFrame()
bout_size             = pd.DataFrame()
bout_AUCperlick       = pd.DataFrame()
bout_meandFF          = pd.DataFrame()

interbout_AUC         = pd.DataFrame()
interbout_size        = pd.DataFrame()
interbout_AUCperlick  = pd.DataFrame()
interbout_meandFF     = pd.DataFrame()

all_lick_meandFF      = pd.DataFrame()

for i in np.arange(len(meandata)):
    mouseID                  = meandata["mouseID"][i]
    bout_AUC[mouseID]        = meandata["bout_AUC"][i]
    bout_size[mouseID]       = meandata["bout_size"][i]
    bout_AUCperlick[mouseID] = meandata["bout_AUCperlick"][i]    
    bout_meandFF[mouseID]    = meandata["bout_meandFF"][i]      

    interbout_AUC[mouseID]        = meandata["interbout_AUC"][i]
    interbout_size[mouseID]       = meandata["interbout_size"][i]
    interbout_AUCperlick[mouseID] = meandata["interbout_AUCperlick"][i]    
    interbout_meandFF[mouseID]    = meandata["interbout_meandFF"][i]  

for i in np.arange(len(meandata)):
    all_lick_meandFF["mouseID"]               = meandata["mouseID"] 
    all_lick_meandFF["all_lick_meandFF"]      = meandata["all_lick_meandFF"] 
    all_lick_meandFF["all_bout_meandFF"]      = meandata["all_bout_meandFF"] 
    all_lick_meandFF["all_interbout_meandFF"] = meandata["all_interbout_meandFF"] 

suptitle       = stimulus + "_" + fluid 
save_file_path = analysis_path + "/" + suptitle + "_bout_AUC.csv"
bout_AUC.to_csv(save_file_path)
save_file_path = analysis_path + "/" + suptitle + "_bout_size.csv"
bout_size.to_csv(save_file_path)
save_file_path = analysis_path + "/" + suptitle + "_bout_AUCperlick.csv"
bout_AUCperlick.to_csv(save_file_path)
save_file_path = analysis_path + "/" + suptitle + "_bout_meandFF.csv"
bout_meandFF.to_csv(save_file_path)
save_file_path = analysis_path + "/" + suptitle + "_interbout_AUC.csv"
interbout_AUC.to_csv(save_file_path)
save_file_path = analysis_path + "/" + suptitle + "_interbout_size.csv"
interbout_size.to_csv(save_file_path)
save_file_path = analysis_path + "/" + suptitle + "_interbout_AUCperlick.csv"
interbout_AUCperlick.to_csv(save_file_path)
save_file_path = analysis_path + "/" + suptitle + "_interbout_meandFF.csv"
interbout_meandFF.to_csv(save_file_path)
save_file_path = analysis_path + "/" + suptitle + "_all_lick_meandFF.csv"
all_lick_meandFF.to_csv(save_file_path)


#%%########################################################################################
#   FIGURES 
###########################################################################################

# #%%################################################
# #   ADDITIONAL PLOTS: mean traces
# #################################################

###################################################################################
## Additional Figure 1: bar graph of AUC and AUC per lick for individual bout      
#####################################################################################
fig, axes     = plt.subplots(nrows = 4, ncols = 1, figsize = (2.5*len(bout_number), 4*4)) ## 4*of rows

bar_width     = 0.75
size          = len(data)
labels        = list(bout_number)

values        = list(np.array(meandata["bout_AUC"]).mean())
if size > 1 :
    errorbars     = list(np.array(meandata["bout_AUC"]).std()/math.sqrt(size-1))     
else:
    errorbars = [np.nan]*len(bout_number)
axes[0].bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 2, alpha = 0.5)
axes[0].set_title("dFF AUC")
axes[0].set_xlabel("bout number")
axes[0].set_ylabel("AUC (a.u.)")
 			
values        = list(np.array(meandata["bout_size"]).mean())
if size > 1 :
    errorbars = list(np.array(meandata["bout_size"]).std()/math.sqrt(size-1))     
else:
    errorbars = [np.nan]*len(bout_number)
axes[1].bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)
axes[1].set_title("bout size")
axes[1].set_xlabel("bout number")
axes[1].set_ylabel("# of licks")

values        = list(np.array(meandata["bout_AUCperlick"]).mean())
if size > 1 :
    errorbars = list(np.array(meandata["bout_AUCperlick"]).std()/math.sqrt(size-1))     
else:
    errorbars = [np.nan]*len(bout_number)
axes[2].bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)
axes[2].set_title("dFF AUC per lick")
axes[2].set_xlabel("bout number")
axes[2].set_ylabel("AUC(a.u.)/# licks")

values        = list(meandata["bout_meandFF"].mean())
if size > 1 :
    errorbars = list(np.array(meandata["bout_meandFF"]).std()/math.sqrt(size-1))     
else:	   
    errorbars = [np.nan]*len(bout_number)   
axes[3].bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)
axes[3].set_title("dFF mean")
axes[3].set_xlabel("bout number")
axes[3].set_ylabel("dFF(%)")

suptitle = stimulus + "_" + fluid + "_boutAUC"                   
plt.suptitle("%s" %suptitle, size = 18)
plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
save_file_path = analysis_path + "/" + suptitle + ".png"
plt.savefig(save_file_path)
save_file_path = analysis_path + "/" + suptitle + ".pdf"
plt.savefig(save_file_path)
plt.show()    
 		 
#%%##################################################################
## Additional Figure 2: bar graph of AUC and AUC per lick for individual interbout interval
####################################################################
if interbout_number > 0:
    
    fig, axes     = plt.subplots(nrows = 4, ncols = 1, figsize = (2.5*len(interbout_number), 4*4)) ## 4*of rows
    
    bar_width     = 0.75
    size          = len(data)
    labels        = list(interbout_number)
    
    values        = list(np.array(meandata["interbout_AUC"]).mean())
    if size > 1 :
        errorbars     = list(np.array(meandata["interbout_AUC"]).std()/math.sqrt(size-1))     
    else:
        errorbars = [np.nan]*len(interbout_number)
    axes[0].bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 2, alpha = 0.5)
    axes[0].set_title("dFF AUC")
    axes[0].set_xlabel("interbout number")
    axes[0].set_ylabel("AUC (a.u.)")
     			
    values        = list(np.array(meandata["interbout_size"]).mean())
    if size > 1 :
        errorbars = list(np.array(meandata["interbout_size"]).std()/math.sqrt(size-1))     
    else:
        errorbars = [np.nan]*len(interbout_number)
    axes[1].bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)
    axes[1].set_title("interbout size")
    axes[1].set_xlabel("interbout number")
    axes[1].set_ylabel("# of licks")
    
    values        = list(np.array(meandata["interbout_AUCperlick"]).mean())
    if size > 1 :
        errorbars = list(np.array(meandata["interbout_AUCperlick"]).std()/math.sqrt(size-1))     
    else:
        errorbars = [np.nan]*len(interbout_number)
    axes[2].bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)
    axes[2].set_title("dFF AUC per lick")
    axes[2].set_xlabel("interbout number")
    axes[2].set_ylabel("AUC(a.u.)/# licks")
    
    values        = list(meandata["interbout_meandFF"].mean())
    if size > 1 :
        errorbars = list(np.array(meandata["interbout_meandFF"]).std()/math.sqrt(size-1))     
    else:	   
        errorbars = [np.nan]*len(interbout_number)   
    axes[3].bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)
    axes[3].set_title("dFF mean")
    axes[3].set_xlabel("interbout number")
    axes[3].set_ylabel("dFF(%)")
    
    suptitle = stimulus + "_" + fluid + "_interboutAUC"                   
    plt.suptitle("%s" %suptitle, size = 18)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
    save_file_path = analysis_path + "/" + suptitle + ".png"
    plt.savefig(save_file_path)
    save_file_path = analysis_path + "/" + suptitle + ".pdf"
    plt.savefig(save_file_path)
    plt.show()    

 			  
#%% ###################################################################################
# ## Additional Figure 3: bar graph of (i) AUC dFF during lick vs during not lick    
#                                      (ii) mean dFF during lick vs during not lick   
# ###################################################################################
fig, ax       = plt.subplots(nrows = 1, ncols = 2, figsize = (8, 6)) ## 4*of rows
bar_width     = 0.75

## 1. mean dFF during all licks vs mean dFF during interbout interval 
values        = [meandata["all_lick_meandFF"].mean(), meandata["all_interbout_meandFF"].mean()]
labels        = ['all licks', 'no lick']
errorbars     = [ss.sem(meandata["all_lick_meandFF"]), ss.sem(meandata["all_interbout_meandFF"])]      
ax1           = plt.subplot2grid((1, 2), (0, 0), rowspan = 1, colspan = 1)
ax1.bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)

ax1.set_ylabel("dFF(%)")
ax1.set_title("mean dF/F0(%)")
 			
## 2. mean dFF during all licks that are part of bouts vs mean dFF during interbout interval 
values        = [meandata["all_bout_meandFF"].mean(), meandata["all_interbout_meandFF"].mean()]
labels        = ['all bouts', 'interbout']
errorbars     = [ss.sem(meandata["all_bout_meandFF"]), ss.sem(meandata["all_interbout_meandFF"])]   
ax2           = plt.subplot2grid((1, 2), (0, 1), rowspan = 1, colspan = 1)            
ax2.bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)

ax2.set_ylabel("dFF(%)")
ax2.set_title("mean dF/F0(%)")

suptitle = stimulus + "_" + fluid + "_lick_meandFF"  
plt.suptitle("%s" %suptitle, size = 18)
plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
save_file_path = analysis_path + "/" + suptitle + ".png"
plt.savefig(save_file_path)
save_file_path = analysis_path + "/" + suptitle + ".pdf"
plt.savefig(save_file_path)
plt.show()  


 			
#%%####################################################################################        
##      FIGURE 1. mean PSTH of first lick of all individual bouts in the order 
#####################################################################################      
fig, ax     = plt.subplots(min_num_bouts, 1, figsize = (6, 2*min_num_bouts))       

for k in range(min_num_bouts):  
    ax      = plt.subplot2grid((min_num_bouts, 1), (k, 0), rowspan=1, colspan=1)       
    legend  = "bout # %d" %(k+1)
    if k == 0:
        plt.errorbar(psth_plot_times, psth_mean[k]["mean"], psth_mean[k]["stderr"], linewidth = 2, color = "red", ecolor = "lightcoral")
        plt.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18, color = "red")            
    else:
        plt.errorbar(psth_plot_times, psth_mean[k]["mean"], psth_mean[k]["stderr"], linewidth = 2, color = "black", ecolor = "lightgray")
        plt.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18)           
    plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
    plt.xlim(-psth_width_sec, psth_width_sec)
    plt.ylim(psth_y_min*1.05, psth_y_max*1.05)   
    plt.yticks(size = 18)
    plt.ylabel("dFF (%)", size = 18)
    if k == (min_num_bouts - 1) :
        plt.xticks(size = 18)
        plt.xlabel("Time (sec)", size = 18)
    else: 
        ax.axes.get_xaxis().set_visible(False)               
suptitle = stimulus + "_" + fluid + "_meanPSTH"              
plt.suptitle("%s\nfirst lick in each bout (n= %d)" %(suptitle, len(psth_mean[0].columns)-2), size = 18)
plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
save_file_path = analysis_path + "/" + suptitle + ".png"
plt.savefig(save_file_path)
save_file_path = analysis_path + "/" + suptitle + ".pdf"
plt.savefig(save_file_path)
plt.show()           


#%%##################################################################################        
##      FIGURE 2. mean PSTH of first lick in first bout VS first lick in other bouts
#####################################################################################       
try:
    fig, ax     = plt.subplots(1, 1, figsize = (8, 4))       
             
    psth_mean_list          = [psth_mean[k]["mean"] for k in np.arange(1, len(psth_mean), 1)]
    psth_mean_otherbouts    = np.nanmean(psth_mean_list, axis = 0)
    psth_stderr_otherbouts  = np.nanstd(psth_mean_list, axis = 0)/math.sqrt(len(psth_mean_list)-1)
    
    ax      = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)       
    plt.errorbar(psth_plot_times, psth_mean[0]["mean"], psth_mean[0]["stderr"], linewidth = 3, color = "red", ecolor = "lightcoral", label = "first bout")
    plt.errorbar(psth_plot_times, psth_mean_otherbouts, psth_stderr_otherbouts, linewidth = 3, color = "black", ecolor = "lightgray", label = "all other bouts")
    plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
    plt.xlim(-psth_width_sec, psth_width_sec)
    plt.ylim(psth_y_min*1.05, psth_y_max*1.05)   
    
    plt.xticks(size = 18)
    plt.yticks(size = 18)
    plt.xlabel("Time (sec)", size = 18)
    plt.ylabel("dFF (%)", size = 18)
    plt.legend(loc="upper left", bbox_to_anchor = (1, 1), fontsize = 14)
    
    suptitle = stimulus + "_" + fluid + "_allmeanPSTH"      
    plt.suptitle("%s (n = %d)" %(suptitle, len(psth_mean[0].columns)-2), size = 18)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
    save_file_path = analysis_path + "/" + suptitle + ".png"
    plt.savefig(save_file_path)
    save_file_path = analysis_path + "/" + suptitle + ".pdf"
    plt.savefig(save_file_path)
    plt.show()

except:
    pass

#%%##################################################################################
#   HEATMAP   LICK AND BOUT Aligned to the water access. Mark first lick with a dotted line
#####################################################################################
df_heatmap                = pd.DataFrame(data.at[0, "dFF_ds_sm"], columns = [data.at[0, "mouseID"]])
for i in range (1, len(data)):
    df_heatmap[data.at[i, "expID"]] = pd.Series(data.at[i, "dFF_ds_sm"])
    
ax        = plt.axes()
sns.heatmap(df_heatmap.T, cmap = "jet", yticklabels = True, ax = ax)
# sns.heatmap(df_heatmap.T, cmap = "PiYG", yticklabels = True, vmax = 5, vmin = -20, ax = ax)
x_labels  = np.arange(-10, 50, 5)
x_ticks   = np.linspace(0, len(df_heatmap), len(x_labels))
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels)
ax.set_xlabel ("Time (min)")

for i in range(0, len(data)):
    ax.scatter(x = (pre_length*60 + data.at[i, "lick_on"][0])*frame_rate/dsfactor, y = 0.5 + i, marker = ".", s = 150, color = "black")
    # ax.axvline(x = data.at[i, @@@@] , ymin = 0.5 + i, ymax = 1.5 + i, linewidth = 5, linestyle = "--", color = "black")
ax.axvline(pre_length*60*frame_rate/dsfactor, linewidth = 3, linestyle = "--", color = "dimgrey") 

for i in range(0, len(data)):
    ax.scatter(x = (pre_length*60 + data.at[i, "bout_on"][0])*frame_rate/dsfactor, y = 0.5 + i, marker = "|", s = 500, color = "blue")
    # ax.axvline(x = data.at[i, @@@@] , ymin = 0.5 + i, ymax = 1.5 + i, linewidth = 5, linestyle = "--", color = "black")
ax.axvline(pre_length*60*frame_rate/dsfactor, linewidth = 3, linestyle = "--", color = "dimgrey") 

save_file_path = analysis_path + "/" + "heatmap_firstlickbout.png"
figure = ax.get_figure() 
figure.savefig(save_file_path, dpi=1200)
save_file_path = analysis_path + "/" + "heatmap_firstlickbout.pdf"
figure = ax.get_figure() 
figure.savefig(save_file_path, dpi=1200)
plt.show()


# ##################################################
# #   POP open the analysis folder
# ##################################################
os.startfile(analysis_path)

