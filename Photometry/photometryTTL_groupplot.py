# -*- coding: utf-8 -*-
"""
Created on Tue May 26 23:03:32 2020

Fiber photometry and licking behavior: TTL-aligned group analysis.

Loads processed photometry .pkl files and quantifies neural activity
and licking behavior aligned to stimulus, first lick, and first bout.

Calculates z-score responses, AUC, response kinetics, lick rates,
cumulative licking, and bout/interbout measures. Generates individual
and mouse-averaged plots, PSTHs, and CSV summaries.

Repeated recordings are averaged per mouse for group analysis.

CAUTION: Include only recordings from the same treatment/stimulus
in the input directory.


"""

import os
import glob
import pickle as pk
import numpy as np
import pandas as pd
import math
from datetime import datetime
from scipy.integrate import simps
from scipy.optimize import curve_fit
from sklearn.neighbors import KDTree
# import scikit_posthocs as sp
import matplotlib.pyplot as plt  # standard Python plotting library
import matplotlib

matplotlib.rcParams['font.size'] = 18 

    
def smooth (x, window_len, window='hanning'):
    if x.ndim != 1:
            raise ValueError ("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
            raise ValueError ("Input vector needs to be bigger than window size.")
    if window_len<3:
            return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError ("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
            w=np.ones(window_len,'d')
    else:  
            w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]

def lick_rate_by_sec (lick_on, time_sec):

    lick_on_ceil = np.ceil(lick_on)
    lick_on_ceil = lick_on_ceil.astype(int)
        
    lick_rate_sec_ = np.histogram(lick_on_ceil, bins = time_sec)   
    lick_rate_sec  = np.concatenate([[0], lick_rate_sec_[0]])
    lick_rate_cum  = np.cumsum(lick_rate_sec)

    return (lick_rate_sec, lick_rate_cum)

def two_lines_max(x, a, b, c, d):
    one = a*x + b
    two = c*x + d
    return np.maximum(one, two)

def two_lines_min(x, a, b, c, d):
    one = a*x + b
    two = c*x + d
    return np.minimum(one, two)

def single_linear(x, a, b):
    y = a*x + b
    return (y)


def three_lines_min(x, a, b, c, d, e, f):
    one = a*x + b
    two = c*x + d
    three = e*x + f
    return np.minimum(one, two, three)

print("Enter the full path for directory with .pkl files: ")
path      = input()


# %%###
if 1== 1:

    add_to_csv = path.split("\\")[-1]

    TTLpath = r"D:\Data_Dropbox_Sync\Mouse_metadata\photometrymiceTTL__.csv" # __.csv for entire timecourse (early stim time)
    # TTLpath = r"D:\Data_Dropbox_Sync\Mouse_metadata\photometrymiceTTL.csv" # .csv for water alignment (late stim time)
    TTLdata = pd.read_csv(TTLpath)

    lick_g          = 1.21 
    
    today = datetime.today().strftime('%Y%m%d')
    analysis_path = path + "/groupAnalysis_"  + today + "_3minPSTH_25minlick_sm501_script17"
    if not os.path.exists(analysis_path):
        os.makedirs(analysis_path)          
        
    if "young" in path:
        palette_mean_GCaMP      = ["black"] 
        palette_ste_GCaMP       = ["silver"]
            
    
    else:
        # palette_mean_GCaMP      = ["black"] 
        # palette_ste_GCaMP       = ["silver"]
        palette_mean_GCaMP      = ["maroon"] 
        palette_ste_GCaMP       = ["pink"]
        
        
    palette_mean_lick       = ["dimgrey", "dimgrey"]
    palette_ste_lick        = ["darkgrey", "darkgrey"]
    
    PSTH_fixed_y        = 1    
    fixed_ylim_mean     = 1
    if fixed_ylim_mean  == 1:
        # fixed_ylim_mean_range_zs               = [-2.5, 0.5]  ## for NAc GRAB-DA3m, FD chow water
        # fixed_ylim_mean_range_zs               = [-1, 4]  
        fixed_ylim_mean_range_lickrate         = [-0.2, 10]
        fixed_ylim_mean_range_cumlick_mincum   = [ -20, 1800]
        fixed_ylim_mean_range_cumwaterBSA      = [ -10, 200]
        
        # min_zs                = fixed_ylim_mean_range_zs[0]
        # max_zs                = fixed_ylim_mean_range_zs[1]
        min_lickrate         = fixed_ylim_mean_range_lickrate[0]
        max_lickrate         = fixed_ylim_mean_range_lickrate[1]
        min_cumlick          = fixed_ylim_mean_range_cumlick_mincum[0]
        max_cumlick          = fixed_ylim_mean_range_cumlick_mincum[1]
        min_cumwaterBW       = fixed_ylim_mean_range_cumwaterBSA[0]
        max_cumwaterBW       = fixed_ylim_mean_range_cumwaterBSA[1]   
      
    max_lickrate_      = 14
    
    sampling_rate      = 1017.2527
    pre_min            = 10     ## in minutes  ## 20 for PEG inj TTL water !!!!!!!!
    post_min           = 25    ## in minutes 

    pre_min_heatmap    = 5
    
    if "GRAB-DA" in path:
        pre_min            = 10  
        post_min           = 25

        pre_min_plot       = 10
        post_min_plot      = 25
        
        water_delay_min    = 0
        pre_min_heatmap    = 5
        post_min_heatmap   = 25         
        heatmap_min        = -2
        heatmap_max        = 10         
        min_zs_            = -2
        max_zs_            = 10      
        min_zs             = -2
        max_zs             = 10 
        max_cumlick        = 1000 
        if "ensure" in path:
            max_cumlick    = 1800
        max_cumlever       = 200
        min_cumlever       = -10
       
        if PSTH_fixed_y == 1:
            PSTH_ymin              = -1
            PSTH_ymax              = 10    
    

    elif "AVP"in path:
        water_delay_min    = 15
        post_min_heatmap   = 45
        min_zs_      = -1
        max_zs_      = 5
        min_zs       = -1
        max_zs       = 5
        heatmap_min        = -1
        heatmap_max        = 5 
        max_cumlick        = 1200
        pre_min_           = 10
        post_min_          = 25        

    
    elif "chow" in path:
        water_delay_min    = 30
        post_min_heatmap   = 60         
        heatmap_min        = -2
        heatmap_max        = 2          
        min_zs_            = -2
        max_zs_            = 2      
        min_zs             = -2
        max_zs             = 2 
        max_cumlick        = 1200
        pre_min_           = 10
        post_min_          = 25


        
    elif ("NaCl" in path) or ("PEG" in path):
        water_delay_min    = 15    ## 15 mins after NaCl or PEG injection. 30min after chow. 
        post_min_heatmap   = 40
        heatmap_min        = -1
        heatmap_max        = 4
        min_zs_            = -1
        max_zs_            = 4  
        min_zs             = -1
        max_zs             = 4 
        if "PEG" in path:
            max_cumlick    = 500
        
        if PSTH_fixed_y == 1:
            PSTH_ymin              = -1
            PSTH_ymax              = 5
    
    else:
        water_delay_min    = 0
        post_min_heatmap   = 30      ## for entire timecourse    
        heatmap_min        = -3
        heatmap_max        = 1    
        min_zs_            = -3
        max_zs_            = 1 
        min_zs             = -3
        max_zs             = 1     

        if PSTH_fixed_y == 1:
            PSTH_ymin              = -1
            PSTH_ymax              = 10

    setnumlicks           = 1
    
    ds_factor1            = 5
    ds_factor2            = 5
    ds_factor3            = 5   
    ds_factor             = int(ds_factor1*ds_factor2*ds_factor3)

    if ("PEG" in path) and ("SFO" in path):
        sm_factor = 51
    else:
        # sm_factor = 51
        # sm_factor = 1001
        sm_factor = 501
        
        
    psth_width_sec        = 180  
    psth_sec_for_AUC_same = 1   ## 0 if different
    if psth_sec_for_AUC_same == 0:
        psth_sec_for_AUC     = 30 ## in seconds, edit here (for GRAB-DA3m, use 4 sec)
    else:
        psth_sec_for_AUC    = psth_width_sec
    
     
    lickrate_min_for_AUC  = post_min ## in minutes, edit here 
    zs_end_min_for_AUC    = post_min
    

       
    frame_start_for_zs_AUC = int(pre_min*60*sampling_rate/ds_factor) 
    frame_end_for_zs_AUC   = int((pre_min+zs_end_min_for_AUC)*60*sampling_rate/ds_factor) 

    
    if ("WD" in path) or ("GRAB-DA" in path) or ("300mMNaCl" in path):
        frame_end_for_zs_rise_AUC = int((pre_min + post_min_heatmap)*60*sampling_rate/ds_factor)
    else:
        frame_end_for_zs_rise_AUC = int((pre_min + water_delay_min)*60*sampling_rate/ds_factor)   
    
    
    all_files                = glob.glob(path + "/*.pkl")
    file_list                = []
    minlenzs_raw_list        = []
    minlenzs_firstlick_list  = []
    minlenzs_firstbout_list  = []
    
    
    data = pd.DataFrame(columns = ["expID", "mouseID", "expdate", "BW", "age_mo", "tau_half_firstlick", "tau_half_firstbout", 
                                   "zscore", "zscore_ds", "zscore_firstlick","zscore_firstbout", "zscore_heatmap", "zscore_peak", "zscore_firstbout_sm", "zscore_heatmap_sm",
                                   
                                   "lickrate_sec_raw", "lickrate_min", "lickrate_sec_heatmap",
                                   "lickon_firstlick", "lickrate_sec_firstlick", 
                                   "lickrate_min_firstlick", #"lickrate_sec_sm60_firstlick", 
                                   "lickrate_min_firstlick_BW", "lickrate_sec_firstlick_BW", 
                                     
                                   "lickon_firstbout", "lickrate_sec_firstbout", 
                                   "lickrate_min_firstbout", "lickrate_sec_sm60_firstbout",
                                   "taulick_half_firstlick", "taulick_half_firstbout",
    
                                   "time_sec_ds_zs", "time_sec_bin", "time_min_bin",
                                   "time_sec_firstlick", "time_sec_firstbout",
                                   
                                   "T_halfdrop_sm_firstbout", "T_halfdrop_sm", "T_halfrise_sm", "T_halfdrop_cumlick_firstbout",
                                   "bout_start_zscore", "T_halfdrop_zscoreDiff_firstbout", "T_halfdrop_zscoreDiff_over_licks",
                                   "tau_half_cumlick_firstbout","tau_half_zscore_over_lick",
                                              
                                   "lick_on", "lick_off", "lick_x", "lick_y",
                                   "bout_on", "bout_off", "bout_x", "bout_y",
                                   "bout_on_ind", "bout_off_ind", "bout_on_zs", "bout_off_zs", 
                                   "bout_size", "bout_length", "bout_number",
                                   "bout_AUC", "bout_AUCperlick", "bout_meanzs", "all_lick_meanzs",
                                   "all_bout_meanzs", "all_interbout_meanzs",
                                   "interbout_number", "interbout_size", "interbout_AUC", "interbout_AUCperlick", "interbout_meanzs",
                                    
                                   "cumlick_sec", "cumbout_sec", 
                                   "cumlick_sec_firstlick", "cumbout_sec_firstlick", "cumwater_sec_firstlick_BW",
                                   "cumlick_sec_firstbout", "cumbout_sec_firstbout", "cumwater_sec_firstbout_BW",
                                   "bout_on_PSTH", "bout_off_PSTH", 
                                           
                                   "cumlick_firstlick", "cumlick_firstbout", "tau_cumlick_firstlick", "tau_cumlick_firstbout",
                                   "nontau_cumlick_firstlick", "nontau_cumlick_firstbout",
                                   "cumwaterBW_firstlick", "cumwaterBW_firstbout",
                                   "cumwater_firstlick_BW_AUC", "cumwater_firstbout_BW_AUC",
                                   "cumwaterBW_firstlick_tauAUC", "cumwaterBW_firstbout_tauAUC",
                                    "zscore_firstlick_setnumlickAUC", "zscore_firstbout_setnumlickAUC",
                                   "zscore_firstlick_setnumlickDiff", "zscore_firstbout_setnumlickDiff",
                                                                      
                                   "zscore_rise_AUC", "zscore_rise_meanz",
                                   "zscore_firstlick_AUC", "zscore_firstbout_AUC", "zscore_firstbout_meanz",
                                   "zscore_firstlick_tauAUC", "zscore_firstbout_tauAUC", 
                                   "zscore_firstlick_tauDiff", "zscore_firstbout_tauDiff",
                                   "lickrate_min_firstlick_AUC","lickrate_sec_firstlick_AUC","cumlick_sec_firstlick_AUC",
                                   "lickrate_min_firstlick_tauAUC","lickrate_sec_firstlick_tauAUC","cumlick_sec_firstlick_tauAUC",
                                   "lickrate_min_firstlick_BW_AUC","lickrate_sec_firstlick_BW_AUC",
                                   "lickrate_sec_firstbout_AUC","cumlick_sec_firstbout_AUC", #"cumlick_sec_firstbout_BW_AUC",
                                   "lickrate_sec_firstbout_tauAUC","cumlick_sec_firstbout_tauAUC"
                                   ])



        
    data["zscore"]                   = data["zscore"].astype(object) 
    data["zscore_ds"]                = data["zscore_ds"].astype(object)
    data["zscore_firstlick"]         = data["zscore_firstlick"].astype(object)   
    data["zscore_firstbout"]         = data["zscore_firstbout"].astype(object)
    data["zscore_heatmap"]           = data["zscore_heatmap"].astype(object)
    data["zscore_heatmap_sm"]        = data["zscore_heatmap_sm"].astype(object)
    data["zscore_firstbout_sm"]      = data["zscore_firstbout_sm"].astype(object)
    
    data["lickrate_sec_raw"]            = data["lickrate_sec_raw"].astype(object)
    data["lickrate_min"]                = data["lickrate_min"].astype(object)
    # data["lickrate_sec_sm60"]           = data["lickrate_sec_sm60"].astype(object)
    data["lickrate_sec_heatmap"]        = data["lickrate_sec_heatmap"].astype(object)  
    
    data["lickon_firstlick"]                = data["lickon_firstlick"].astype(object)
    data["lickrate_min_firstlick"]          = data["lickrate_min_firstlick"].astype(object)
    data["lickrate_min_firstlick_BW"]       = data["lickrate_min_firstlick_BW"].astype(object)
    data["lickrate_sec_firstlick"]          = data["lickrate_sec_firstlick"].astype(object)
    data["lickrate_sec_firstlick_BW"]       = data["lickrate_sec_firstlick_BW"].astype(object)
       
    data["lickon_firstbout"]                = data["lickon_firstbout"].astype(object)
    data["lickrate_min_firstbout"]          = data["lickrate_min_firstbout"].astype(object)
    data["lickrate_sec_firstbout"]          = data["lickrate_sec_firstbout"].astype(object)
    
    data["time_sec_ds_zs"]                  = data["time_sec_ds_zs"].astype(object) 
    data["time_sec_firstlick"]              = data["time_sec_firstlick"].astype(object)
    data["time_sec_firstbout"]              = data["time_sec_firstbout"].astype(object)
    
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
    data["bout_meanzs"]       = data["bout_meanzs"].astype(object)

    data["interbout_number"]       = data["interbout_number"].astype(object)    
    data["interbout_AUC"]          = data["interbout_AUC"].astype(object)
    data["interbout_size"]         = data["interbout_size"].astype(object)
    data["interbout_AUCperlick"]   = data["interbout_AUCperlick"].astype(object)
    data["interbout_meanzs"]       = data["interbout_meanzs"].astype(object)
    
    data["bout_on_ind"]         = data["bout_on_ind"].astype(object)  
    data["bout_off_ind"]        = data["bout_off_ind"].astype(object)       
    data["bout_on_zs"]         = data["bout_on_zs"].astype(object)     
    data["bout_off_zs"]        = data["bout_off_zs"].astype(object)   
    
    data["cumlick_sec"]       = data["cumlick_sec"].astype(object)
    data["cumbout_sec"]       = data["cumbout_sec"].astype(object)
    
    data["cumlick_sec_firstlick"]       = data["cumlick_sec_firstlick"].astype(object)
    data["cumwater_sec_firstlick_BW"]   = data["cumwater_sec_firstlick_BW"].astype(object)
    data["cumbout_sec_firstlick"]       = data["cumbout_sec_firstlick"].astype(object)
    
    data["cumlick_sec_firstbout"]       = data["cumlick_sec_firstbout"].astype(object)
    data["cumwater_sec_firstbout_BW"]   = data["cumwater_sec_firstbout_BW"].astype(object)
    data["cumbout_sec_firstbout"]       = data["cumbout_sec_firstbout"].astype(object)
    
    data["bout_on_PSTH"]        = data["bout_on_PSTH"].astype(object)
    data["bout_off_PSTH"]       = data["bout_off_PSTH"].astype(object)


    time_min_bin              = np.arange(-pre_min*60, post_min*60+1, 60)
    time_sec_bin              = np.arange(-pre_min*60, post_min*60+1, 1)
    time_sec_bin_heatmap      = np.arange(-pre_min_heatmap*60, post_min_heatmap*60+1, 1)
    
    if "GRAB" in path:
        time_min_bin              = np.arange(-lickrate_min_for_AUC*60, lickrate_min_for_AUC*60+1, 60)
        time_sec_bin              = np.arange(-lickrate_min_for_AUC*60, lickrate_min_for_AUC*60+1, 1)

        
    
    #% Go through each file in the path to save the data into a single dataframe
    expID_list     = []
    df_lever_zs    = {}
    
    for i, file_path in enumerate (all_files):
        expID            = file_path[(len(path)+1):-4] ## previous called as "file_name_cut
        file_list.append(expID)  
        inputfile        = open(file_path, 'rb')
        tempdata         = pk.load(inputfile)
        inputfile.close()
        
        if i == 0:
            stimulus  = tempdata["stimulus"]
            fluid     = tempdata["LICK_fluid"]
        data.at[i, "expID"]         = str(expID)
        expID_list.append(expID)
        
        data.at[i, "mouseID"]       = tempdata["mouseID"]
        data.at[i, "expdate"]       = tempdata["date"]
        data.at[i, "age_mo"]        = int(tempdata["age_mo"])        
        try:
            data.at[i, "BW"]        = float(tempdata["BW"])
        except:
            for j in range(len(TTLdata)):
                if (str(int(TTLdata.at[j, "ExpDate"])) == str(data.at[i, "expdate"])) & (TTLdata.at[j, "mouseID"] == data.at[i, "mouseID"]):
                    data.at[i, "BW"]          = float(TTLdata.at[j, "BW"])
                    break    
        
        if data.at[i, "BW"] > 0:
            pass
        else:
            print ("%s has no BW data" %expID)
                    
        data.at[i, "lick_on"]       = tempdata["LICK_on"][tempdata["LICK_on"]>water_delay_min*60]
        data.at[i, "lick_off"]      = tempdata["LICK_off"][tempdata["LICK_off"]>water_delay_min*60]   
        data.at[i, "bout_on"]       = tempdata["BOUT_on"][tempdata["BOUT_on"]>water_delay_min*60]
        data.at[i, "bout_off"]      = tempdata["BOUT_off"][tempdata["BOUT_off"]>water_delay_min*60]
    
        data.at[i, "lick_x"]    = np.append(np.append(tempdata["time_sec"][0], np.reshape(np.kron([data.at[i, "lick_on"] , data.at[i, "lick_off"] ],
                                                       np.array([[1], [1]])).T, [1,-1])[0]), tempdata["time_sec"][-1])        
        LICK_size               = len(data.at[i, "lick_on"])
        d                       = np.ones(LICK_size)
        data.at[i, "lick_y"]    = np.append(np.append(0, np.reshape(np.vstack([np.zeros(LICK_size), d, d, np.zeros(LICK_size)]).T, [1, -1])[0]),0)
        
  
        # LICK_EVENT_x              = np.append(tempdata["time_sec"][0], np.append(np.reshape(np.kron([data.at[i, "bout_on"], data.at[i, "bout_off"]], np.array([[1], [1]])).T, [1,-1])[0], tempdata["time_sec"][-1]))          
        # LICK_EVENT_sz             = len(data.at[i, "bout_on"])
        # LICK_EVENT_y              = np.append(np.append(0, np.reshape(np.vstack([np.zeros(LICK_EVENT_sz), bout_data, bout_data, np.zeros(LICK_EVENT_sz)]).T, [1 ,-1])[0]), 0)
        # data.at[i, "bout_x"]      = LICK_EVENT_x
        # data.at[i, "bout_y"]      = LICK_EVENT_y
        data.at[i, "bout_x"]    = tempdata["BOUT_x"]
        data.at[i, "bout_y"]    = tempdata["BOUT_y"]
                    
        
        
        temp_bout_length            = list(np.array(data.at[i, "bout_off"])  - np.array(data.at[i, "bout_on"]))
        temp_bout_length_nan        = [np.nan]*(50 - len(tempdata["BOUT_size"]))
        data.at[i, "bout_length"]   = np.array(temp_bout_length + temp_bout_length_nan)

        data.at[i, "bout_size"]             = tempdata["df_bout_zs"]["bout_size"]

        interbout_lick_size    = np.asarray([])
        interbout_lick_length  = np.asarray([])
        		
        #save interbout lick on times 
        for k in range(len(data.at[i, "bout_on"])-1):
            interbout_lickons = data.at[i, "lick_on"][data.at[i, "lick_on"] > data.at[i, "bout_off"][k]]
            interbout_lickons = interbout_lickons[interbout_lickons < data.at[i, "bout_on"][k+1]]
            try: 
                lick_length       = interbout_lickons[-1] - interbout_lickons[0]
            except:
                lick_length       = 0
            interbout_lick_size   = np.append(interbout_lick_size, len(interbout_lickons))
            interbout_lick_length = np.append(interbout_lick_length, lick_length)
        
    
        lickrate_sec, lickrate_cum          = lick_rate_by_sec(data.at[i, "lick_on"], time_sec_bin)      
        lickrate_min, dump                  = lick_rate_by_sec(data.at[i, "lick_on"], time_min_bin) 
        
        data.at[i, "lickrate_sec_raw"]      = lickrate_sec
        data.at[i, "cumlick_sec"]           = lickrate_cum+1
        dump, data.at[i, "cumbout_sec"]     = lick_rate_by_sec(data.at[i, "bout_on"], time_sec_bin)
             
        
        data.at[i, "lickrate_sec_heatmap"], dump = lick_rate_by_sec(data.at[i, "lick_on"], time_sec_bin_heatmap)
        
        data.at[i, "lickrate_min"]               = lickrate_min
    
        # data.at[i, "zscore"]                = tempdata["zscore"]
        data.at[i, "zscore_ds"]               = tempdata["zscore_ds"]
        data.at[i, "zscore_heatmap"]          = tempdata["zscore_ds"]
        data.at[i, "time_sec_ds_zs"]          = tempdata["time_sec_ds_zs"]
     
        minlenzs_temp                      = len(data.at[i, "time_sec_ds_zs"])
        minlenzs_raw_list.append(minlenzs_temp)
        
        time_sec_dims                      = np.expand_dims(data.at[i, "time_sec_ds_zs"], axis = 1)
        tree                               = KDTree(time_sec_dims, leaf_size = 2)  
      
        ############################################################################
        ### lick rate: time_sec, time_min. aligned to first lick, trim
        ############################################################################ 
        firstlick                              = data.at[i, "lick_on"][0]
        data.at[i, "lickon_firstlick"]         = data.at[i, "lick_on"] - firstlick
        data.at[i, "time_sec_firstlick"]       = data.at[i, "time_sec_ds_zs"] - firstlick
            
        firstlick_reshaped                             = firstlick.reshape(1, -1)
        firstlick_ind                                  = tree.query(firstlick_reshaped)[1][0][0]
        
        timezero_ind                                   = tree.query([[0]])[1][0][0]
        diff_ind                                       = firstlick_ind - timezero_ind 
        data.at[i, "zscore_firstlick"]                 = data.at[i, "zscore_ds"][diff_ind:] 
        data.at[i, "time_sec_firstlick"]               = data.at[i, "time_sec_firstlick"][diff_ind:]
    
        minlenzs_temp                                  = len(data.at[i, "time_sec_firstlick"])
        minlenzs_firstlick_list.append(minlenzs_temp)
        
        lickrate_sec_firstlick, lickrate_cum_firstlick = lick_rate_by_sec(data.at[i, "lickon_firstlick"], time_sec_bin)
        data.at[i, "lickrate_sec_firstlick"]           = lickrate_sec_firstlick
        data.at[i, "cumlick_sec_firstlick"]            = lickrate_cum_firstlick + 1
      
        dump, cumbout_sec_firstlick                    = lick_rate_by_sec(data.at[i, "bout_on"] - firstlick, time_sec_bin)
        data.at[i, "cumbout_sec_firstlick"]            = cumbout_sec_firstlick + 1
    
        lickrate_min_firstlick, dump                   = lick_rate_by_sec(data.at[i, "lickon_firstlick"], time_min_bin)    
        data.at[i, "lickrate_min_firstlick"]           = lickrate_min_firstlick

        ############################################################################
        ### lick rate: time_min. aligned to first bout, trim
        ############################################################################   
        firstbout                           = data.at[i, "bout_on"][0]
                

        data.at[i, "lickon_firstbout"]      = data.at[i, "lick_on"] - firstbout
        data.at[i, "time_sec_firstbout"]    = data.at[i, "time_sec_ds_zs"] - firstbout
        
        firstbout_reshaped                  = firstbout.reshape(1, -1)
        firstbout_ind                       = tree.query(firstbout_reshaped)[1][0][0]
        
        timezero_ind                        = tree.query([[0]])[1][0][0]
        diff_ind                            = firstbout_ind - timezero_ind 
        data.at[i, "zscore_firstbout"]      = data.at[i, "zscore_ds"][diff_ind:] 
        data.at[i, "time_sec_firstbout"]    = data.at[i, "time_sec_firstbout"][diff_ind:]
    
        minlenzs_temp                       = len(data.at[i, "time_sec_firstbout"])
        minlenzs_firstbout_list.append(minlenzs_temp)
        
        lickrate_sec_firstbout, lickrate_cum_firstbout = lick_rate_by_sec(data.at[i, "lickon_firstbout"], time_sec_bin)
        data.at[i, "lickrate_sec_firstbout"]           = lickrate_sec_firstbout
        data.at[i, "cumlick_sec_firstbout"]            = lickrate_cum_firstbout + 1
        
        dump, cumbout_sec_firstbout                    = lick_rate_by_sec(data.at[i, "bout_on"] - firstbout, time_sec_bin)
        data.at[i, "cumbout_sec_firstbout"]            = cumbout_sec_firstbout + 1
         
        lickrate_min_firstbout, dump         = lick_rate_by_sec(data.at[i, "lickon_firstbout"], time_min_bin)      
        lickrate_min_firstbout_nan           = [np.nan if x== 0 else x for x in lickrate_min_firstbout]
        data.at[i, "lickrate_min_firstbout"] = lickrate_min_firstbout
 

        ############################################################################
        ### cumulative water intake (g) normalized by BSA or BS
        ############################################################################   

        if "WD" in stimulus:     # for WD (or prandial thirst) - water or their controls
            data.at[i, "cumwater_sec_firstlick_BW"]     = data.at[i, "cumlick_sec_firstlick"]*lick_g/data.at[i, "BW"]**0.667
            data.at[i, "cumwater_sec_firstbout_BW"]     = data.at[i, "cumlick_sec_firstbout"]*lick_g/data.at[i, "BW"]**0.667
            # print (data.at[i, "cumwater_sec_firstbout_BW"])
            
        else:
            data.at[i, "cumwater_sec_firstlick_BW"]     = data.at[i, "cumlick_sec_firstlick"]*lick_g/data.at[i, "BW"]
            data.at[i, "cumwater_sec_firstbout_BW"]     = data.at[i, "cumlick_sec_firstbout"]*lick_g/data.at[i, "BW"]
 
       
            
        #%%#############################################################

          
        ##############################################################
        ### mean z score for individual lick events 
        ##############################################################
        lick_ind_list  = np.array([])
        df_lick_zs     = pd.DataFrame(columns = ["lick_meanzs"])
        for k, lick_on in enumerate (data.at[i, "lick_on"]):
            lick_on                                = lick_on.reshape(1, -1)       				
            ind_lick_on                            = tree.query(lick_on)[1][0][0]
            df_lick_zs.at[k, "lick_meanzs"]        = data.at[i, "zscore_ds"][ind_lick_on]
            lick_ind_list                          = np.append(lick_ind_list, ind_lick_on)
        data.at[i, "all_lick_meanzs"]              = df_lick_zs["lick_meanzs"].mean()
        
        ##############################################################
        ### bout and interbout analysis 
        ##############################################################       
        # data.at[i, "bout"]
        for k, (ind_on, ind_off) in enumerate (zip(tempdata["df_bout_zs"]["bout_on_ind"], tempdata["df_bout_zs"]["bout_off_ind"])):
            a                                = simps(tempdata["zscore_ds"][int(ind_on):int(ind_off+1)])    
            data.at[i, "bout_AUC"]           = np.append(data.at[i, "bout_AUC"], a)
            data.at[i, "bout_AUC"][k]        = a
            data.at[i, "bout_AUCperlick"]    = np.append(data.at[i, "bout_AUCperlick"], a/tempdata["df_bout_zs"]["bout_size"][k])        				 
            data.at[i, "bout_meanzs"]        = np.append(data.at[i, "bout_meanzs"], np.mean(tempdata["zscore_ds"][int(ind_on):int(ind_off+1)]))
            # data.at[i, "bout_AUCperlick"][k] = data.at[i, "bout_AUC"][k]/data.at[i, "bout_size"][k]        				 
            # data.at[i, "bout_meanzs"][k]     = np.mean(data.at[i, "zscore_ds"][int(ind_on):int(ind_off+1)])
    

        for k in range(len(tempdata["df_bout_zs"])-1):       
            data.at[i, "interbout_number"]      = np.append(data.at[i, "interbout_number"], int(k+1))
            # data.at[i, "interbout_size"]        = np.append(data.at[i, "interbout_size"], interbout_lick_size[k])
            data.at[i, "interbout_number"][k]   = int(k+1)	 
            interbout_on_ind                    = int(tempdata["df_bout_zs"].at[k, "bout_off_ind"])
            interbout_off_ind                   = int(tempdata["df_bout_zs"].at[k+1, "bout_on_ind"])
        
            try:
                b                                   = simps(tempdata["zscore_ds"][int(interbout_on_ind):int(interbout_off_ind)])    
                data.at[i, "interbout_AUC"]         = np.append(data.at[i, "interbout_AUC"], b)
                data.at[i, "interbout_AUCperlick"]  = np.append(data.at[i, "interbout_AUCperlick"], b/interbout_lick_size[k])
            except:
                data.at[i, "interbout_AUC"]         = np.nan
                data.at[i, "interbout_AUCperlick"]  = np.nan
                
            try:
                data.at[i, "interbout_meanzs"]      = np.append(data.at[i, "interbout_meanzs"], np.mean(tempdata["zscore_ds"][int(interbout_on_ind):int(interbout_off_ind+1)]))
            except:
    
                data.at[i, "interbout_meanzs"]      = np.nan
            # data.at[i, "interbout_size"]       = interbout_lick_size
            data.at[i, "all_bout_meanzs"]      = np.nanmean(data.at[i, "bout_meanzs"])
            data.at[i, "all_interbout_meanzs"] = np.nanmean(data.at[i, "interbout_meanzs"])
    

            
        #%##################################################################################
    #   #   z score AUC by pre-defined minutes 
    #   #####################################################################################

        data.at[i, "zscore_firstlick_AUC"]   = simps(data.at[i, "zscore_firstlick"][frame_start_for_zs_AUC:frame_end_for_zs_AUC], dx=1)
        data.at[i, "zscore_firstbout_AUC"]   = simps(data.at[i, "zscore_firstbout"][frame_start_for_zs_AUC:frame_end_for_zs_AUC], dx=1)
        data.at[i, "zscore_firstbout_meanz"] = np.mean(data.at[i, "zscore_firstbout"][frame_start_for_zs_AUC:frame_end_for_zs_AUC])
        
        data.at[i, "zscore_rise_meanz"] = np.mean(data.at[i, "zscore_ds"][frame_start_for_zs_AUC:frame_end_for_zs_rise_AUC])
        data.at[i, "zscore_rise_AUC"]   = simps(data.at[i, "zscore_ds"][frame_start_for_zs_AUC:frame_end_for_zs_rise_AUC], dx=1)
                
        
    #   #%%##################################################################################
    #   #   lickrate and cumulative lick AUC by pre-defined minutes 
    #   #####################################################################################
        
        # data.at[i, "lickrate_min_firstlick_AUC"] = simps(data.at[i,"lickrate_min_firstlick"][(pre_min+1):int(pre_min+1+lickrate_min_for_AUC)])
        data.at[i, "lickrate_sec_firstlick_AUC"] = simps(data.at[i,"lickrate_sec_firstlick"][pre_min*60:int((pre_min+lickrate_min_for_AUC)*60)])
        data.at[i, "cumlick_sec_firstlick_AUC"]  = simps(data.at[i,"cumlick_sec_firstlick"][pre_min*60:int((pre_min+lickrate_min_for_AUC)*60)])
    
        # data.at[i, "lickrate_min_firstlick_BW_AUC"] = simps(data.at[i,"lickrate_min_firstlick_BW"][(pre_min+1):(pre_min+1+lickrate_min_for_AUC)])
        # data.at[i, "lickrate_sec_firstlick_BW_AUC"] = simps(data.at[i,"lickrate_sec_firstlick_BW"][pre_min*60:(pre_min+lickrate_min_for_AUC)*60])
        data.at[i, "cumwater_firstlick_BW_AUC"]  = simps(data.at[i,"cumwater_sec_firstlick_BW"][pre_min*60:int((pre_min+lickrate_min_for_AUC)*60)])
        data.at[i, "cumwater_firstbout_BW_AUC"]  = simps(data.at[i,"cumwater_sec_firstbout_BW"][pre_min*60:int((pre_min+lickrate_min_for_AUC)*60)])
        
        data.at[i, "lickrate_sec_firstbout_AUC"]    = simps(data.at[i,"lickrate_sec_firstbout"][pre_min*60:int((pre_min+lickrate_min_for_AUC)*60)])
        data.at[i, "cumlick_sec_firstbout_AUC"]     = simps(data.at[i,"cumlick_sec_firstbout"][pre_min*60:int((pre_min+lickrate_min_for_AUC)*60)])
    
        data.at[i, "cumlick_firstlick"]     = data.at[i, "cumlick_sec_firstlick"][-1]
        data.at[i, "cumlick_firstbout"]     = data.at[i, "cumlick_sec_firstbout"][-1]

        data.at[i, "cumwaterBW_firstlick"]  =  data.at[i, "cumwater_sec_firstlick_BW"][-1]        
        data.at[i, "cumwaterBW_firstbout"]  =  data.at[i, "cumwater_sec_firstbout_BW"][-1]
        
        
        #%%##################################################################################
    #   #   trim z score trace by pre-min and post-min 
    #   ####################        if "GRAB" in path:
        data.at[i, "zscore_ds"]          = data.at[i, "zscore_ds"][:int((pre_min+post_min)*60*sampling_rate/ds_factor)] 
        data.at[i, "zscore_firstlick"]   = data.at[i, "zscore_firstlick"][:int((pre_min+post_min)*60*sampling_rate/ds_factor)] 
        data.at[i, "zscore_firstbout"]   = data.at[i, "zscore_firstbout"][:int((pre_min+post_min)*60*sampling_rate/ds_factor)] 
        if "GRAB" in path:
            data.at[i, "zscore_heatmap"]     = data.at[i, "zscore_heatmap"][:int((pre_min+post_min)*60*sampling_rate/ds_factor)]          
        else:
            data.at[i, "zscore_heatmap"]     = data.at[i, "zscore_heatmap"][int((10-pre_min_heatmap)*60*sampling_rate/ds_factor):int((10+post_min_heatmap)*60*sampling_rate/ds_factor)]
        data.at[i, "zscore_peak"]            = np.max(data.at[i, "zscore_heatmap"])

        time_sec_raw          = data.at[i, "time_sec_ds_zs"][int((10-pre_min_heatmap)*60*sampling_rate/ds_factor):int((10+post_min_heatmap)*60*sampling_rate/ds_factor)]
        time_sec_firstlick    = data.at[i, "time_sec_firstlick"][:int((pre_min+post_min)*60*sampling_rate/ds_factor)] 
        time_sec_firstbout    = data.at[i, "time_sec_firstbout"][:int((pre_min+post_min)*60*sampling_rate/ds_factor)] 



        #%#######################################################################
        #   curve_fit cumulative lick plot by two lines, get the crossover time  
        ##########################################################################         
        xlick_firstlick = time_sec_bin[600:]
        ylick_firstlick = data.at[i, "cumlick_sec_firstlick"][600:]
            
        pw0 = (.67, 0, 0.1, 150) # a guess for slope, intercept, slope, intercept
        pwlick_firstlick, cov   = curve_fit(two_lines_min, xlick_firstlick, ylick_firstlick, pw0)
        crossover_firstlick     = (pwlick_firstlick[3] - pwlick_firstlick[1]) / (pwlick_firstlick[0] - pwlick_firstlick[2]) 
        data.at[i, "taulick_half_firstlick"] = round(crossover_firstlick/2, 3)

        xlick_firstbout = time_sec_bin[600:]
        # ylick_firstbout = smooth(data.at[i, "cumlick_sec_firstbout"][600:], window_len = 61)
        ylick_firstbout = data.at[i, "cumlick_sec_firstbout"][600:]
            
        pw0 = (5, 0, 0.5, 200) # a guess for slope, intercept, slope, intercept
        pwlick_firstbout, cov  = curve_fit(two_lines_min, xlick_firstbout, ylick_firstbout, pw0)
        crossover_firstbout    = (pwlick_firstbout[3] - pwlick_firstbout[1]) / (pwlick_firstbout[0] - pwlick_firstbout[2]) 
        data.at[i, "taulick_half_firstbout"] = round(crossover_firstbout/2, 3)

        xlick_firstbout = time_sec_bin[600:]
        # ylick_firstbout = smooth(data.at[i, "cumlick_sec_firstbout"][600:], window_len = 61)
        ylick_firstbout = data.at[i, "cumlick_sec_firstbout"][600:]
            
        pw0 = (5, 0, 2.5, 100, 1, 200) # a guess for slope, intercept, slope, intercept
        pwlick_firstbout, cov  = curve_fit(three_lines_min, xlick_firstbout, ylick_firstbout, pw0)
        crossover_firstbout    = (pwlick_firstbout[3] - pwlick_firstbout[1]) / (pwlick_firstbout[0] - pwlick_firstbout[2]) 
        data.at[i, "taulick_half_firstbout"] = round(crossover_firstbout/2, 3)

    
       

         #%%##################################################################################
    #   #   Find time to peak value (rise), half time to baseline (drop)
    #   #####################################################################################    

            
        data.at[i, "zscore_firstbout_sm"] = smooth(data.at[i, "zscore_firstbout"], window_len = sm_factor)
        data.at[i, "bout_start_zscore"]   = data.at[i, "zscore_firstbout_sm"][int(pre_min*60*sampling_rate/ds_factor)]

        if "PEG" in path:
            diff                          = np.absolute(data.at[i, "zscore_firstbout_sm"][int(pre_min*60*sampling_rate/ds_factor):int(pre_min*3*60*sampling_rate/ds_factor)] - data.at[i, "bout_start_zscore"]/3*2)
        # elif "WD" in path:
            # diff                          = np.absolute(data.at[i, "zscore_firstbout_sm"][int(pre_min*60*sampling_rate/ds_factor):int(pre_min*3*60*sampling_rate/ds_factor)] - data.at[i, "bout_start_zscore"]/2)
        else:             
            diff                          = np.absolute(data.at[i, "zscore_firstbout_sm"][int(pre_min*60*sampling_rate/ds_factor):int(pre_min*1.5*60*sampling_rate/ds_factor)] - data.at[i, "bout_start_zscore"]/2)

        ind_T_halfdrop_sm                                   = diff.argmin() + int(pre_min*60*sampling_rate/ds_factor)
        data.at[i, "T_halfdrop_sm_firstbout"]               = data.at[i, "time_sec_firstbout"][ind_T_halfdrop_sm] 
        data.at[i, "T_halfdrop_sm"]                         = data.at[i, "time_sec_firstbout"][ind_T_halfdrop_sm] + firstbout
        data.at[i, "T_halfdrop_cumlick_firstbout"]          = data.at[i, "cumlick_sec_firstbout"][int(round(data.at[i, "T_halfdrop_sm_firstbout"])+pre_min*60)]
        data.at[i, "T_halfdrop_zscoreDiff_firstbout"]       = data.at[i, "bout_start_zscore"] - data.at[i, "zscore_firstbout_sm"][int(round(ind_T_halfdrop_sm))]
        data.at[i, "T_halfdrop_zscoreDiff_over_licks"]      = data.at[i, "T_halfdrop_zscoreDiff_firstbout"]/data.at[i, "T_halfdrop_cumlick_firstbout"] 

        
        zscore_heatmap_sm                 = smooth(data.at[i, "zscore_heatmap"], window_len = sm_factor)      
        peak_zs                           = np.max(zscore_heatmap_sm)
        if "PEG" in path:
            diff                              = np.absolute(zscore_heatmap_sm[int(pre_min_heatmap*1.75*60*sampling_rate/ds_factor):int(pre_min_heatmap*3.25*60*sampling_rate/ds_factor)]- peak_zs/3)
            ind_T_halfrise_sm                 = diff.argmin() +int(pre_min_heatmap*1.75*60*sampling_rate/ds_factor)
        elif (("AVP" in path) and ("chow" in path)):
            diff                              = np.absolute(zscore_heatmap_sm[int(pre_min_heatmap*60*sampling_rate/ds_factor):int(pre_min_heatmap*5*60*sampling_rate/ds_factor)]- peak_zs/2)            
            ind_T_halfrise_sm                 = diff.argmin() +int(pre_min_heatmap*60*sampling_rate/ds_factor)            

        elif ("GRAB" in path):
            diff                              = np.absolute(zscore_heatmap_sm[int(pre_min_heatmap*60*sampling_rate/ds_factor):int(pre_min_heatmap*2*60*sampling_rate/ds_factor)]- peak_zs/2)            
            ind_T_halfrise_sm                 = diff.argmin() +int(pre_min_heatmap*60*sampling_rate/ds_factor)
        else:
            diff                              = np.absolute(zscore_heatmap_sm[int(pre_min_heatmap*60*sampling_rate/ds_factor):int(pre_min_heatmap*3*60*sampling_rate/ds_factor)]- peak_zs/2)            
            ind_T_halfrise_sm                 = diff.argmin() +int(pre_min_heatmap*60*sampling_rate/ds_factor)
        data.at[i, "T_halfrise_sm"]           = time_sec_raw[ind_T_halfrise_sm]
        data.at[i, "zscore_heatmap_sm"]       = zscore_heatmap_sm
      
        if "GRAB-DA" in path:
            continue


        #%%#######################################################################
        #   curve_fit the trace by two lines, get the half time to reach baseline 
        ##########################################################################         
        x_firstlick = data.at[i, "time_sec_firstlick"][int(pre_min*60*sampling_rate/ds_factor):int((pre_min+post_min)*60*sampling_rate/ds_factor)]
        y_firstlick = data.at[i, "zscore_firstlick"][int(pre_min*60*sampling_rate/ds_factor):int((pre_min+post_min)*60*sampling_rate/ds_factor)]


        if "PEG" in path:
            y_firstlick = smooth(y_firstlick, 1001)
            
        pw0 = (.0007, 0.5, -0.007, 1.5) # a guess for slope, intercept, slope, intercept            
        # pw0 = (.02, -2, -0.1, 3) # a guess for slope, intercept, slope, intercept
        pw_firstlick, cov = curve_fit(two_lines_max, x_firstlick, y_firstlick, pw0)
        crossover_firstlick = (pw_firstlick[3] - pw_firstlick[1]) / (pw_firstlick[0] - pw_firstlick[2]) 
        data.at[i, "tau_half_firstlick"] = round(crossover_firstlick/2, 3)

        x_firstbout = data.at[i, "time_sec_firstbout"][int(pre_min*60*sampling_rate/ds_factor):int((pre_min+post_min)*60*sampling_rate/ds_factor)]
        y_firstbout = data.at[i, "zscore_firstbout"][int(pre_min*60*sampling_rate/ds_factor):int((pre_min+post_min)*60*sampling_rate/ds_factor)]

       
        if "PEG" in path:
            y_firstbout = smooth(y_firstbout, 1001)
           
        try:
            pw_firstbout, cov = curve_fit(two_lines_max, x_firstbout, y_firstbout, pw0)
            crossover_firstbout = (pw_firstbout[3] - pw_firstbout[1]) / (pw_firstbout[0] - pw_firstbout[2]) 
        except:
            crossover_firstbout = crossover_firstlick
            
        data.at[i, "tau_half_firstbout"] = round(crossover_firstbout/2, 3)    
        #############################################################################
        ### cumulative # licks during post_min and during tau (crossover) after curve fit
        ############################################################################3

        try:    crossover_ind_firstlick = int(round(crossover_firstlick) + 1 + pre_min*60)
        except: crossover_ind_firstlick = np.nan
            
        try:
            data.at[i, "tau_cumlick_firstlick"]     = data.at[i, "cumlick_sec_firstlick"][crossover_ind_firstlick]
            data.at[i, "nontau_cumlick_firstlick"]  = data.at[i, "cumlick_firstlick"] - data.at[i, "tau_cumlick_firstlick"]
        
        except:
            data.at[i, "tau_cumlick_firstlick"]     = np.nan
            data.at[i, "nontau_cumlick_firstlick"]  = np.nan
                    
            
        crossover_ind_firstbout = int(round(crossover_firstbout) + 1 + pre_min*60)
        try:
            data.at[i, "tau_cumlick_firstbout"]     = data.at[i, "cumlick_sec_firstbout"][crossover_ind_firstbout]
            data.at[i, "nontau_cumlick_firstbout"]  = data.at[i, "cumlick_firstbout"] - data.at[i, "tau_cumlick_firstbout"]
        except:
            data.at[i, "tau_cumlick_firstbout"]     = np.nan
            data.at[i, "nontau_cumlick_firstbout"]  = np.nan

        #%##################################################################################
    #   #   z score AUC by tau (decay time)
    #   #####################################################################################

        try:    data.at[i, "zscore_firstlick_tauAUC"]   = simps(data.at[i, "zscore_firstlick"][frame_start_for_zs_AUC:int(crossover_ind_firstlick*sampling_rate/ds_factor)], dx=1) 
        except: data.at[i, "zscore_firstlick_tauAUC"]   = np.nan

        try:    data.at[i, "zscore_firstbout_tauAUC"]   = simps(data.at[i, "zscore_firstbout"][frame_start_for_zs_AUC:int(crossover_ind_firstbout*sampling_rate/ds_factor)], dx=1)
        except: data.at[i, "zscore_firstbout_tauAUC"]   = np.nan

        try:    data.at[i, "zscore_firstlick_tauDiff"]  = data.at[i, "zscore_firstlick"][int(crossover_ind_firstlick*sampling_rate/ds_factor)]-data.at[i, "zscore_firstlick"][frame_start_for_zs_AUC]
        except: data.at[i, "zscore_firstlick_tauDiff"]  = np.nan

        try:    data.at[i, "zscore_firstbout_tauDiff"]  = data.at[i, "zscore_firstbout"][int(crossover_ind_firstbout*sampling_rate/ds_factor)]-data.at[i, "zscore_firstbout"][frame_start_for_zs_AUC]
        except: data.at[i, "zscore_firstbout_tauDiff"]  = np.nan


        try:
            data.at[i, "tau_half_cumlick_firstbout"]  = data.at[i, "cumlick_sec_firstbout"][int(round(crossover_firstlick/2) + 1 + pre_min*60)]
            data.at[i, "tau_half_zscore_over_lick"]   = data.at[i, "zscore_firstlick_tauDiff"]/2/data.at[i, "tau_half_cumlick_firstbout"]
        except:
            data.at[i, "tau_half_cumlick_firstbout"]  = np.nan
            data.at[i, "tau_half_zscore_over_lick"]   = np.nan
        
    #   #%%##################################################################################
    #   #   lickrate and cumulative lick AUC by tau (decay time)
    #   #####################################################################################
        
        # data.at[i, "lickrate_min_firstlick_AUC"] = simps(data.at[i,"lickrate_min_firstlick"][(pre_min+1):(pre_min+1+lickrate_min_for_AUC)])
        try:
            data.at[i, "lickrate_sec_firstlick_tauAUC"] = simps(data.at[i,"lickrate_sec_firstlick"][pre_min*60:crossover_ind_firstlick])
            data.at[i, "cumlick_sec_firstlick_tauAUC"]  = simps(data.at[i,"cumlick_sec_firstlick"][pre_min*60:crossover_ind_firstlick])
            data.at[i, "cumwaterBW_firstlick_tauAUC"]  = simps(data.at[i, "cumwater_sec_firstlick_BW"][pre_min*60:crossover_ind_firstlick]) 

        except:
            data.at[i, "lickrate_sec_firstlick_tauAUC"] = np.nan 
            data.at[i, "cumlick_sec_firstlick_tauAUC"]  = np.nan
            data.at[i, "cumwaterBW_firstlick_tauAUC"]  = np.nan
    
        try:
            data.at[i, "lickrate_sec_firstbout_tauAUC"] = simps(data.at[i,"lickrate_sec_firstbout"][pre_min*60:crossover_ind_firstbout])
            data.at[i, "cumlick_sec_firstbout_tauAUC"]  = simps(data.at[i,"cumlick_sec_firstbout"][pre_min*60:crossover_ind_firstbout])
            data.at[i, "cumwaterBW_firstbout_tauAUC"]  = simps(data.at[i, "cumwater_sec_firstbout_BW"][pre_min*60:crossover_ind_firstbout])

        except:
            data.at[i, "lickrate_sec_firstbout_tauAUC"] = np.nan 
            data.at[i, "cumlick_sec_firstbout_tauAUC"]  = np.nan
            data.at[i, "cumwaterBW_firstbout_tauAUC"]  = np.nan            


        data.at[i, "zscore_firstbout_sm"] = smooth(data.at[i, "zscore_firstbout"], window_len = 101)    
        #%%##############################################################################
          #       FIGURE 1.  indivudal plot: z score, curve_fit lines, lick rate (Hz)
          ##############################################################################
          # Cosmetics for plotting
        if "WD" in stimulus: 
              y_scale    = -np.min(data.at[i,"zscore_ds"])/3
              y_shift    = np.min(data.at[i,"zscore_ds"])/3*4
        else:
              y_scale    = np.max(data.at[i,"zscore_ds"])/3
              y_shift    = -np.max(data.at[i,"zscore_ds"])  
        
        fig, ax = plt.subplots(2, 1, figsize = (10, 10))
    
    
        ax0   = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)   
        lns1    = ax0.plot(time_sec_firstlick/60, data.at[i,"zscore_firstlick"], linewidth = 1, color='grey', label='GCaMP')
        ax0.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
        ax0.set_ylabel("z score", fontsize = 20)
        ax0.set_xlabel("time (min)")
        ax0.set_ylim (min_zs_, max_zs_)
        
        ax1   = ax0.twinx()    
        lns2    = ax1.plot(x_firstlick/60, two_lines_max(x_firstlick, *pw_firstlick), "-", linewidth = 6, color = "red", label = "fit") 
        ax1.set_ylim (min_zs_, max_zs_)
        ax1.axes.get_yaxis().set_visible(False)
        
        ax2   = ax0.twinx()
        lns3    = ax2.plot(time_sec_bin/60, data.at[i, "lickrate_sec_firstlick"], linewidth = 1, color ="maroon", label = "lick")
        ax2.set_ylabel("lick rate (Hz)", fontsize = 20)
        ax2.set_ylim (min_lickrate, max_lickrate_)
        ax2.set_xlim (-pre_min, post_min)
        # ax2.axes.get_yaxis().set_visible(False)
        
        lns     = lns1 + lns2 + lns3
        labels  = [l.get_label() for l in lns]        
        ax0.legend(lns, labels, loc="upper right", bbox_to_anchor = (1, 1), fontsize = 18)
        ax0.set_title("aligned by time 0 = first lick", fontsize = 20)
    
    
        ax3   = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)   
        lns4    = ax3.plot(time_sec_firstbout/60, data.at[i,"zscore_firstbout"], linewidth = 1, color='grey', label='GCaMP')
        ax3.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
        ax3.set_ylabel("z score", fontsize = 20)
        ax3.set_xlabel("time (min)")
        ax3.set_ylim (min_zs_, max_zs_)
        
        ax4   = ax3.twinx()
        lns5    = ax4.plot(x_firstbout/60, two_lines_max(x_firstbout, *pw_firstbout), "-", linewidth = 6, color = "red", label = "fit") 
        ax4.set_ylim (min_zs_, max_zs_)
        ax4.axes.get_yaxis().set_visible(False)
        
        ax5   = ax3.twinx()
        lns6    = ax5.plot(time_sec_bin/60, data.at[i, "lickrate_sec_firstbout"], linewidth = 1, color ="maroon", label = "lick")
        ax5.set_ylabel("lick rate (Hz)", fontsize = 20)
        ax5.set_ylim (min_lickrate, max_lickrate_)
        ax5.set_xlim (-pre_min, post_min)
          # ax3.axes.get_yaxis().set_visible(False)
        
        lns     = lns4 + lns5 + lns6
        labels  = [l.get_label() for l in lns]        
        ax3.legend(lns, labels, loc="upper right", bbox_to_anchor = (1, 1), fontsize = 18)
        ax3.set_title("aligned by time 0 = first bout", fontsize = 20)        
        
        suptitle = data.at[i,"mouseID"] + "_" + data.at[i,"expdate"] + "_" + stimulus + "_" + fluid + "_" + str(int(data.at[i,"age_mo"])) + "mo" + "_zs_fit"        
        plt.suptitle("%s" %suptitle, fontsize = 18) 
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])           
        save_file_path = analysis_path + "/" + suptitle + ".png"
        plt.savefig(save_file_path)
          # save_file_path = analysis_path + "/" + suptitle + ".pdf"
          # plt.savefig(save_file_path)
          # plt.show()
        plt.close()
    
        #%##############################################################################
          #       FIGURE 2.  indivudal plot: z score, cumulative licks and lick curve fit 
          ##############################################################################
          # Cosmetics for plotting
        if "WD" in stimulus: 
              y_scale    = -np.min(data.at[i,"zscore_ds"])/3
              y_shift    = np.min(data.at[i,"zscore_ds"])/3*4
        else:
              y_scale    = np.max(data.at[i,"zscore_ds"])/3
              y_shift    = -np.max(data.at[i,"zscore_ds"])  
        
        fig, ax = plt.subplots(2, 1, figsize = (10, 10))
    
        ax0   = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)   
        lns1  = ax0.plot(time_sec_firstlick/60, data.at[i,"zscore_firstlick"], linewidth = 2, color='grey', label='GCaMP')
        ax0.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
        ax0.set_ylabel("z score", fontsize = 20)
        ax0.set_xlabel("time (min)")
        ax0.set_ylim (min_zs, max_zs)
    
        
        ax6   = ax0.twinx()    
        lns8  = ax6.plot(x_firstlick/60, two_lines_max(x_firstlick, *pw_firstlick), "-", linewidth = 6, color = "red", label = "fit") 
        ax6.set_ylim (min_zs, max_zs)
        ax6.axes.get_yaxis().set_visible(False)
                
        ax2   = ax0.twinx()
        lns3  = ax2.plot(time_sec_bin/60, data.at[i, "cumlick_sec_firstlick"], linewidth = 2, color ="maroon", label = "lick")
        ax2.set_ylabel("# licks", fontsize = 20)
        ax2.set_ylim (min_cumlick, max_cumlick)
        ax2.set_xlim (-pre_min, post_min)
          # ax2.axes.get_yaxis().set_visible(False)
    
        ax1   = ax0.twinx()    
        ax1.set_ylim (min_cumlick, max_cumlick)
        ax1.axes.get_yaxis().set_visible(False)
        
        lns     = lns1 + lns3 + lns8
        labels  = [l.get_label() for l in lns]     
        ax0.set_title("aligned by time 0 = first lick", fontsize = 20)
    
    
        ax3   = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)   
        lns4  = ax3.plot(time_sec_firstbout/60, data.at[i,"zscore_firstbout"], linewidth = 2, color='grey', label='GCaMP')
        ax3.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
        ax3.set_ylabel("z score", fontsize = 20)
        ax3.set_xlabel("time (min)")
        ax3.set_ylim (min_zs, max_zs)
        
                
        ax7   = ax3.twinx()    
        lns7  = ax3.plot(x_firstbout/60, two_lines_max(x_firstbout, *pw_firstbout), "-", linewidth = 6, color = "red", label = "fit") 
        ax7.set_ylim (min_zs_, max_zs_)
        ax7.axes.get_yaxis().set_visible(False)
        
        ax5   = ax3.twinx()
        lns6  = ax5.plot(time_sec_bin/60, data.at[i, "cumlick_sec_firstbout"], linewidth = 2, color ="maroon", label = "lick")
        ax5.set_ylabel("# licks", fontsize = 20)
        ax5.set_ylim (min_cumlick, max_cumlick)
        ax5.set_xlim (-pre_min, post_min)
        
        ax4   = ax3.twinx()
        ax4.set_ylim (min_cumlick, max_cumlick)
        ax4.axes.get_yaxis().set_visible(False)
        
        
        lns    = lns4 + lns6 + lns7
        lns     = lns4 + lns5 + lns6 + lns7
        labels  = [l.get_label() for l in lns]        
          # ax3.legend(lns, labels, loc="upper right", bbox_to_anchor = (1, 1), fontsize = 18)
        ax3.set_title("aligned by time 0 = first bout", fontsize = 20)        
        
        suptitle = data.at[i,"mouseID"] + "_" + data.at[i,"expdate"] + "_" + stimulus + "_" + fluid + "_" + str(int(data.at[i,"age_mo"])) + "mo" + "_zs_lickfit"        
        plt.suptitle("%s" %suptitle, fontsize = 18) 
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])   
        
        save_file_path = analysis_path + "/" + suptitle + ".png"
        plt.savefig(save_file_path)
        plt.show()


#%#######################################################################  
    print ("files analyzed are: %s" %file_list)       
 
       
   
    #%%###################################################################
    #  get the smallest number of bouts  
    #####################################################################
    bout_number_list = [1000]
    for i in range(len(data)):
        bout_number_list.append(len(data.at[i, "bout_on"]))
             
    min_num_bouts    = np.nanmin(bout_number_list)
    # print (min_num_bouts)
    # bout_number      = np.arange(1, min_num_bouts+1, 1)
    # interbout_number = bout_number[:-1]
    
    for i in np.arange(len(data)):
        data.at[i, "bout_AUC"]              = np.array(data.at[i, "bout_AUC"][:min_num_bouts])
        data.at[i, "bout_size"]             = np.array(data.at[i, "bout_size"][:min_num_bouts])
        data.at[i, "bout_AUCperlick"]       = np.array(data.at[i, "bout_AUCperlick"][:min_num_bouts])

        data.at[i, "bout_meanzs"]           = np.array(data.at[i, "bout_meanzs"][:min_num_bouts])
        data.at[i, "bout_on"]               = np.array(data.at[i, "bout_on"][:min_num_bouts])
        data.at[i, "bout_off"]              = np.array(data.at[i, "bout_off"][:min_num_bouts])
        
        try:
            data.at[i, "interbout_AUC"]         = np.array(data.at[i, "interbout_AUC"][:min_num_bouts-1])
            # data.at[i, "interbout_size"]        = np.array(data.at[i, "interbout_size"][:min_num_bouts-1])
            data.at[i, "interbout_AUCperlick"]  = np.array(data.at[i, "interbout_AUCperlick"][:min_num_bouts-1])
            data.at[i, "interbout_meanzs"]      = np.array(data.at[i, "interbout_meanzs"][:min_num_bouts-1])
        except:
            data.at[i, "interbout_AUC"]         = np.array([np.nan])
            data.at[i, "interbout_AUCperlick"]  = np.array([np.nan])
            data.at[i, "interbout_meanzs"]      = np.array([np.nan])
        
    #%%#############################################################################################
    #   PROCESS DATA: "Groupby" by "mouseID", save the mean zs into the new dataframe
    ################################################################################################
    len_zscore_firstlick = []
    len_zscore_firstbout = []
    for i in np.arange(0, len(data), 1):
        len_zscore_firstbout.append(len(data.at[i, "zscore_firstbout"]))
    maxlen_zscore_firstbout = np.max(np.array(len_zscore_firstbout))
    
    for i in np.arange(0, len(data), 1):
        data.at[i, "zscore_firstbout"] = np.append(data.at[i, "zscore_firstbout"], np.zeros(maxlen_zscore_firstbout-len_zscore_firstbout[i]) + np.nan)
        
    meandata   = pd.DataFrame({ "BW":data.groupby("mouseID")["BW"].apply(np.mean),
                                "zscore_ds":data.groupby("mouseID")["zscore_ds"].apply(np.mean),
                                "zscore_heatmap":data.groupby("mouseID")["zscore_heatmap"].apply(np.mean),
                                "zscore_firstbout": data.groupby("mouseID")["zscore_firstbout"].apply(np.mean), 
                                
                                "zscore_heatmap_sm":data.groupby("mouseID")["zscore_heatmap_sm"].apply(np.mean),
                                "zscore_firstbout_sm": data.groupby("mouseID")["zscore_firstbout_sm"].apply(np.mean), 
    
                                "lickrate_sec_heatmap": data.groupby("mouseID")["lickrate_sec_heatmap"].apply(np.mean), 
                                "lickrate_sec_firstbout": data.groupby("mouseID")["lickrate_sec_firstbout"].apply(np.mean), 
                                "cumlick_sec_firstbout": data.groupby("mouseID")["cumlick_sec_firstbout"].apply(np.mean), 
                                "cumlever_sec_firstbout": data.groupby("mouseID")["cumlever_sec_firstbout"].apply(np.mean),
                                "no_leverpress": data.groupby("mouseID")["no_leverpress"].apply(np.mean),

                                # "cumwater_sec_firstlick_BW": data.groupby("mouseID")["cumwater_sec_firstlick_BW"].apply(np.mean), 
                                # "cumwater_sec_firstbout_BW": data.groupby("mouseID")["cumwater_sec_firstbout_BW"].apply(np.mean), 
                                
                                "bout_start_zscore": data.groupby("mouseID")["bout_start_zscore"].apply(np.mean),                                    
                                "zscore_firstbout_AUC":data.groupby("mouseID")["zscore_firstbout_AUC"].apply(np.mean),
                                
                                "lickrate_sec_firstlick_AUC": data.groupby("mouseID")["lickrate_sec_firstlick_AUC"].apply(np.mean),
                                "cumlick_sec_firstlick_AUC": data.groupby("mouseID")["cumlick_sec_firstlick_AUC"].apply(np.mean),

                                # "cumwater_firstlick_BW_AUC": data.groupby("mouseID")["cumwater_firstlick_BW_AUC"].apply(np.mean),
                                # "cumwater_firstbout_BW_AUC": data.groupby("mouseID")["cumwater_firstbout_BW_AUC"].apply(np.mean),
                            
                                "cumlick_firstlick": data.groupby("mouseID")["cumlick_firstlick"].apply(np.mean),
                                "cumlick_firstbout": data.groupby("mouseID")["cumlick_firstbout"].apply(np.mean),

                                # "cumwaterBW_firstlick": data.groupby("mouseID")["cumwaterBW_firstlick"].apply(np.mean),    
                                # "cumwaterBW_firstbout": data.groupby("mouseID")["cumwaterBW_firstbout"].apply(np.mean),
                                  
                                "bout_AUC": data.groupby("mouseID")["bout_AUC"].apply(np.mean),
                                "bout_size": data.groupby("mouseID")["bout_size"].apply(np.mean),
                                "bout_AUCperlick": data.groupby("mouseID")["bout_AUCperlick"].apply(np.mean),
                                "bout_meanzs": data.groupby("mouseID")["bout_meanzs"].apply(np.mean),
    
                                "interbout_AUC": data.groupby("mouseID")["interbout_AUC"].apply(np.mean),
                                # "interbout_size": data.groupby("mouseID")["interbout_size"].apply(np.mean),
                                "interbout_AUCperlick": data.groupby("mouseID")["interbout_AUCperlick"].apply(np.mean),
                                "interbout_meanzs": data.groupby("mouseID")["interbout_meanzs"].apply(np.mean),
                               
                                "all_lick_meanzs": data.groupby("mouseID")["all_lick_meanzs"].apply(np.mean),
                                "all_bout_meanzs": data.groupby("mouseID")["all_bout_meanzs"].apply(np.mean),
                                "all_interbout_meanzs": data.groupby("mouseID")["all_interbout_meanzs"].apply(np.mean),                                  
                               })#.reset_index()   

    #%###################################################
    #   Save data and meandata as .csv files
    #####################################################        
    
    meandata_to_save     = meandata[["BW", "cumlever_sec_firstbout",  "no_leverpress", 
                                    "all_lick_meanzs", "all_lever_meanzs", "all_bout_meanzs", "all_interbout_meanzs"                                                         
                                    ]]    

    csv_path         = analysis_path + "/meandata_" + stimulus + "_" + fluid + "_" + add_to_csv + ".csv"
    meandata_to_save.to_csv(csv_path, index = True)

    df1 = pd.DataFrame()
    for i in range(len(data)):
        df1[data.at[i, "expID"]] = pd.Series(data.at[i, "bout_meanzs"])
        
    csv_path              = analysis_path + "/meanBOUTzs" + "_" + stimulus + "_" + fluid + "_" + add_to_csv + ".csv"
    df1.to_csv(csv_path, index = True)
        
    df2 = pd.DataFrame()
    for i in range(len(data)):
        df2[data.at[i, "expID"]] = pd.Series(data.at[i, "interbout_meanzs"])
        
    csv_path              = analysis_path + "/meanINTERBOUTzs" + "_" + stimulus + "_" + fluid + "_" + add_to_csv + ".csv"
    df2.to_csv(csv_path, index = True)        


    #%%#############################################################
    ###     PSTH for first lick of individual bouts 
    ##############################################################                   
    psth_width      = math.floor(sampling_rate*psth_width_sec/ds_factor)
    psth_num_points = 2*psth_width+1
    psth_plot_times = np.linspace(-psth_width_sec, psth_width_sec, psth_num_points, endpoint = True, dtype = None)       
                   
    for i, file_path in enumerate (all_files):                                          
        psth            = np.empty([len(data.at[i,"bout_size"]), psth_num_points])
        psth[:]         = np.nan
        psth_y_min_list = []
        psth_y_max_list = []
        
        for k, bout_on in enumerate (data.at[i, "bout_on"]):
            bout_on     = bout_on.reshape(1, -1)
            ind         = tree.query(bout_on)[1][0][0] 
    
            ## Process without error if PSTH post-width exceeds the recording, fill with nan.
            if ind + psth_width + 1 >= len(data.at[i, "zscore_ds"]):  
                last_post_psth                          = len(data.at[i, "zscore_ds"]) - ind
                try: 
                    psth[k, : psth_width + last_post_psth]  = data.at[i, "zscore_ds"][ind - psth_width : ] - np.nanmedian(data.at[i, "zscore_ds"][ind - psth_width : ind]) 
                except:
                    continue
            else:    
                psth[k, :]  = data.at[i, "zscore_ds"][ind - psth_width : ind + psth_width +1] - np.nanmedian(data.at[i, "zscore_ds"][ind - psth_width : ind])
    
            psth_y_min_list.append(np.nanmin(psth[k, :]))
            psth_y_max_list.append(np.nanmax(psth[k, :]))
        data.at[i, "bout_on_PSTH"] = pd.DataFrame(psth)
        

    if "GRAB-DA" in path:
    
        len_zscore_firstlick = []
        len_zscore_firstbout = []
        for i in np.arange(0, len(data), 1):
            len_zscore_firstlick.append(len(data.at[i, "zscore_firstlick"]))
            len_zscore_firstbout.append(len(data.at[i, "zscore_firstbout"]))
        maxlen_zscore_firstlick = np.max(np.array(len_zscore_firstlick))
        maxlen_zscore_firstbout = np.max(np.array(len_zscore_firstbout))
        
        for i in np.arange(0, len(data), 1):
            data.at[i, "zscore_firstlick"] = np.append(data.at[i, "zscore_firstlick"], np.zeros(maxlen_zscore_firstlick-len_zscore_firstlick[i]) + np.nan)
            data.at[i, "zscore_firstbout"] = np.append(data.at[i, "zscore_firstbout"], np.zeros(maxlen_zscore_firstbout-len_zscore_firstbout[i]) + np.nan)
            # print ("%d %s %d" %(i, data.at[i, "expID"], len(data.at[i, "zscore_firstbout"])))
            
        meandata   = pd.DataFrame({ "BW":data.groupby("mouseID")["BW"].apply(np.mean),
                                    "zscore_ds":data.groupby("mouseID")["zscore_ds"].apply(np.mean),
                                    "zscore_heatmap":data.groupby("mouseID")["zscore_heatmap"].apply(np.mean),
                                    "zscore_firstlick": data.groupby("mouseID")["zscore_firstlick"].apply(np.mean), 
                                    "zscore_firstbout": data.groupby("mouseID")["zscore_firstbout"].apply(np.mean), 

                                    "zscore_heatmap_sm":data.groupby("mouseID")["zscore_heatmap_sm"].apply(np.mean),
                                    "zscore_firstbout_sm": data.groupby("mouseID")["zscore_firstbout_sm"].apply(np.mean), 
        
                                    "lickrate_sec_heatmap": data.groupby("mouseID")["lickrate_sec_heatmap"].apply(np.mean), 
                                    "lickrate_sec_firstlick": data.groupby("mouseID")["lickrate_sec_firstlick"].apply(np.mean), 
                                    "cumlick_sec_firstlick": data.groupby("mouseID")["cumlick_sec_firstlick"].apply(np.mean), 
                                    "lickrate_sec_firstbout": data.groupby("mouseID")["lickrate_sec_firstbout"].apply(np.mean), 
                                    "cumlick_sec_firstbout": data.groupby("mouseID")["cumlick_sec_firstbout"].apply(np.mean), 

                                    "cumwater_sec_firstlick_BW": data.groupby("mouseID")["cumwater_sec_firstlick_BW"].apply(np.mean), 
                                    "cumwater_sec_firstbout_BW": data.groupby("mouseID")["cumwater_sec_firstbout_BW"].apply(np.mean), 
                                    
                                    "T_halfdrop_sm_firstbout":data.groupby("mouseID")["T_halfdrop_sm_firstbout"].apply(np.mean), 
                                    "T_halfdrop_sm":data.groupby("mouseID")["T_halfdrop_sm"].apply(np.mean), 
                                    
                                    
                                    "bout_start_zscore": data.groupby("mouseID")["bout_start_zscore"].apply(np.mean),
                                    
                                    "zscore_rise_AUC": data.groupby("mouseID")["zscore_rise_AUC"].apply(np.mean),
                                    "zscore_rise_meanz": data.groupby("mouseID")["zscore_rise_meanz"].apply(np.mean),
                                    "zscore_firstlick_AUC":data.groupby("mouseID")["zscore_firstlick_AUC"].apply(np.mean),
                                    "zscore_firstbout_AUC":data.groupby("mouseID")["zscore_firstbout_AUC"].apply(np.mean),
                                    "zscore_firstbout_meanz":data.groupby("mouseID")["zscore_firstbout_meanz"].apply(np.mean),
                                    "zscore_peak":data.groupby("mouseID")["zscore_peak"].apply(np.mean),
                                    
                                    "lickrate_sec_firstlick_AUC": data.groupby("mouseID")["lickrate_sec_firstlick_AUC"].apply(np.mean),
                                    "cumlick_sec_firstlick_AUC": data.groupby("mouseID")["cumlick_sec_firstlick_AUC"].apply(np.mean),

                                    "cumwater_firstlick_BW_AUC": data.groupby("mouseID")["cumwater_firstlick_BW_AUC"].apply(np.mean),
                                    "cumwater_firstbout_BW_AUC": data.groupby("mouseID")["cumwater_firstbout_BW_AUC"].apply(np.mean),
                                
                                    "cumlick_firstlick": data.groupby("mouseID")["cumlick_firstlick"].apply(np.mean),
                                    "cumlick_firstbout": data.groupby("mouseID")["cumlick_firstbout"].apply(np.mean),
    
                                    "cumwaterBW_firstlick": data.groupby("mouseID")["cumwaterBW_firstlick"].apply(np.mean),    
                                    "cumwaterBW_firstbout": data.groupby("mouseID")["cumwaterBW_firstbout"].apply(np.mean),
                                      
                                    "bout_AUC": data.groupby("mouseID")["bout_AUC"].apply(np.mean),
                                    "bout_size": data.groupby("mouseID")["bout_size"].apply(np.mean),
                                    "bout_AUCperlick": data.groupby("mouseID")["bout_AUCperlick"].apply(np.mean),
                                    "bout_meanzs": data.groupby("mouseID")["bout_meanzs"].apply(np.mean),
        
                                    "interbout_AUC": data.groupby("mouseID")["interbout_AUC"].apply(np.mean),
                                    "interbout_size": data.groupby("mouseID")["interbout_size"].apply(np.mean),
                                    "interbout_AUCperlick": data.groupby("mouseID")["interbout_AUCperlick"].apply(np.mean),
                                    "interbout_meanzs": data.groupby("mouseID")["interbout_meanzs"].apply(np.mean),
                                   
                                    "all_lick_meanzs": data.groupby("mouseID")["all_lick_meanzs"].apply(np.mean),
                                    "all_bout_meanzs": data.groupby("mouseID")["all_bout_meanzs"].apply(np.mean),
                                    "all_interbout_meanzs": data.groupby("mouseID")["all_interbout_meanzs"].apply(np.mean),                                  
                                   })#.reset_index()   

        #%###################################################
        #   Save data and meandata as .csv files
        #####################################################        
        meandata_to_save     = meandata[["BW", "cumlick_firstlick", "cumwater_firstlick_BW_AUC",  
                                        "cumlick_firstbout", "cumwater_firstbout_BW_AUC", 
                                        "all_lick_meanzs", "all_bout_meanzs", "all_interbout_meanzs", "zscore_firstbout_meanz", "zscore_peak", "T_halfdrop_sm_firstbout"                                                      
                                        ]]    

        csv_path         = analysis_path + "/meandata_" + stimulus + "_" + fluid + "_" + add_to_csv + ".csv"
        meandata_to_save.to_csv(csv_path, index = True)
    
        # df = data["bout_meanzs"] 
        # df.columns = data["expID"]
        df1 = pd.DataFrame()
        for i in range(len(data)):
            df1[data.at[i, "expID"]] = pd.Series(data.at[i, "bout_meanzs"])
            
        csv_path              = analysis_path + "/meanBOUTzs" + "_" + stimulus + "_" + fluid + "_" + add_to_csv + ".csv"
        df1.to_csv(csv_path, index = True)
            
        df2 = pd.DataFrame()
        for i in range(len(data)):
            df2[data.at[i, "expID"]] = pd.Series(data.at[i, "interbout_meanzs"])
            
        csv_path              = analysis_path + "/meanINTERBOUTzs" + "_" + stimulus + "_" + fluid + "_" + add_to_csv + ".csv"
        df2.to_csv(csv_path, index = True)            
        
    else:
        len_zscore_firstlick = []
        len_zscore_firstbout = []
        for i in np.arange(0, len(data), 1):
            len_zscore_firstlick.append(len(data.at[i, "zscore_firstlick"]))
            len_zscore_firstbout.append(len(data.at[i, "zscore_firstbout"]))
        maxlen_zscore_firstlick = np.max(np.array(len_zscore_firstlick))
        maxlen_zscore_firstbout = np.max(np.array(len_zscore_firstbout))
        
        for i in np.arange(0, len(data), 1):
            data.at[i, "zscore_firstlick"] = np.append(data.at[i, "zscore_firstlick"], np.zeros(maxlen_zscore_firstlick-len_zscore_firstlick[i]) + np.nan)
            data.at[i, "zscore_firstbout"] = np.append(data.at[i, "zscore_firstbout"], np.zeros(maxlen_zscore_firstbout-len_zscore_firstbout[i]) + np.nan)
            # print ("%d %s %d" %(i, data.at[i, "expID"], len(data.at[i, "zscore_firstbout"])))
            
        meandata   = pd.DataFrame({ "BW":data.groupby("mouseID")["BW"].apply(np.mean),
                                    "tau_half_firstlick":data.groupby("mouseID")["tau_half_firstlick"].apply(np.mean),
                                    "tau_half_firstbout":data.groupby("mouseID")["tau_half_firstbout"].apply(np.mean),
                                    "zscore_ds":data.groupby("mouseID")["zscore_ds"].apply(np.mean),
                                    "zscore_heatmap":data.groupby("mouseID")["zscore_heatmap"].apply(np.mean),
                                    "zscore_firstlick": data.groupby("mouseID")["zscore_firstlick"].apply(np.mean), 
                                    "zscore_firstbout": data.groupby("mouseID")["zscore_firstbout"].apply(np.mean), 

                                    "zscore_heatmap_sm":data.groupby("mouseID")["zscore_heatmap_sm"].apply(np.mean),
                                    "zscore_firstbout_sm": data.groupby("mouseID")["zscore_firstbout_sm"].apply(np.mean), 
        
                                    "lickrate_sec_heatmap": data.groupby("mouseID")["lickrate_sec_heatmap"].apply(np.mean), 
                                   "lickrate_sec_firstlick": data.groupby("mouseID")["lickrate_sec_firstlick"].apply(np.mean), 
                                   "cumlick_sec_firstlick": data.groupby("mouseID")["cumlick_sec_firstlick"].apply(np.mean), 
                                   # "cumbout_sec_firstlick": data.groupby("mouseID")["cumbout_sec_firstlick"].apply(np.mean), 
        
                                   # "lickrate_min_firstbout": data.groupby("mouseID")["lickrate_min_firstbout"].apply(np.mean), 
                                    "lickrate_sec_firstbout": data.groupby("mouseID")["lickrate_sec_firstbout"].apply(np.mean), 
                                    "cumlick_sec_firstbout": data.groupby("mouseID")["cumlick_sec_firstbout"].apply(np.mean), 
                                   # "cumbout_sec_firstbout": data.groupby("mouseID")["cumbout_sec_firstbout"].apply(np.mean), 
        
                                   # "lickrate_min_firstlick_BW": data.groupby("mouseID")["lickrate_min_firstlick_BW"].apply(np.mean), 
                                   # "lickrate_sec_firstlick_BW": data.groupby("mouseID")["lickrate_sec_firstlick_BW"].apply(np.mean), 
                                    "cumwater_sec_firstlick_BW": data.groupby("mouseID")["cumwater_sec_firstlick_BW"].apply(np.mean), 
                                    "cumwater_sec_firstbout_BW": data.groupby("mouseID")["cumwater_sec_firstbout_BW"].apply(np.mean), 
            
                                    "taulick_half_firstlick": data.groupby("mouseID")["taulick_half_firstlick"].apply(np.mean),
                                    "taulick_half_firstbout": data.groupby("mouseID")["taulick_half_firstbout"].apply(np.mean),
                                    
                                    "T_halfdrop_sm_firstbout": data.groupby("mouseID")["T_halfdrop_sm_firstbout"].apply(np.mean),
                                    "T_halfdrop_sm": data.groupby("mouseID")["T_halfdrop_sm"].apply(np.mean),
                                    "T_halfrise_sm": data.groupby("mouseID")["T_halfrise_sm"].apply(np.mean),
                                    "T_halfdrop_cumlick_firstbout": data.groupby("mouseID")["T_halfdrop_cumlick_firstbout"].apply(np.mean),
                                    "bout_start_zscore": data.groupby("mouseID")["bout_start_zscore"].apply(np.mean),
                                    "T_halfdrop_zscoreDiff_firstbout": data.groupby("mouseID")["T_halfdrop_zscoreDiff_firstbout"].apply(np.mean),
                                    "T_halfdrop_zscoreDiff_over_licks": data.groupby("mouseID")["T_halfdrop_zscoreDiff_over_licks"].apply(np.mean),
                                    "tau_half_cumlick_firstbout": data.groupby("mouseID")["tau_half_cumlick_firstbout"].apply(np.mean),
                                    "tau_half_zscore_over_lick": data.groupby("mouseID")["tau_half_zscore_over_lick"].apply(np.mean),
                                    
                                    "zscore_rise_AUC": data.groupby("mouseID")["zscore_rise_AUC"].apply(np.mean),
                                    "zscore_rise_meanz": data.groupby("mouseID")["zscore_rise_meanz"].apply(np.mean),
                                                                       
                                    "zscore_firstlick_AUC":data.groupby("mouseID")["zscore_firstlick_AUC"].apply(np.mean),
                                    "zscore_firstbout_AUC":data.groupby("mouseID")["zscore_firstbout_AUC"].apply(np.mean),
    
                                    "zscore_firstlick_tauAUC":data.groupby("mouseID")["zscore_firstlick_tauAUC"].apply(np.mean),
                                    "zscore_firstbout_tauAUC":data.groupby("mouseID")["zscore_firstbout_tauAUC"].apply(np.mean),
    
                                    "zscore_firstlick_tauDiff":data.groupby("mouseID")["zscore_firstlick_tauDiff"].apply(np.mean),
                                    "zscore_firstbout_tauDiff":data.groupby("mouseID")["zscore_firstbout_tauDiff"].apply(np.mean),
                                   
                                    # "lickrate_min_firstlick_AUC": data.groupby("mouseID")["lickrate_min_firstlick_AUC"].apply(np.mean),
                                   "lickrate_sec_firstlick_AUC": data.groupby("mouseID")["lickrate_sec_firstlick_AUC"].apply(np.mean),
                                   "cumlick_sec_firstlick_AUC": data.groupby("mouseID")["cumlick_sec_firstlick_AUC"].apply(np.mean),
                                   # "lickrate_min_firstlick_BW_AUC": data.groupby("mouseID")["lickrate_min_firstlick_BW_AUC"].apply(np.mean),
                                   # "lickrate_sec_firstlick_BW_AUC": data.groupby("mouseID")["lickrate_sec_firstlick_BW_AUC"].apply(np.mean),
                                    "cumwater_firstlick_BW_AUC": data.groupby("mouseID")["cumwater_firstlick_BW_AUC"].apply(np.mean),
                                    "cumwater_firstbout_BW_AUC": data.groupby("mouseID")["cumwater_firstbout_BW_AUC"].apply(np.mean),
                                    
                                    "lickrate_sec_firstbout_AUC": data.groupby("mouseID")["lickrate_sec_firstbout_AUC"].apply(np.mean),
                                    "cumlick_sec_firstbout_AUC": data.groupby("mouseID")["cumlick_sec_firstbout_AUC"].apply(np.mean),
     
                                    "lickrate_sec_firstlick_tauAUC": data.groupby("mouseID")["lickrate_sec_firstlick_tauAUC"].apply(np.mean),
                                    "cumlick_sec_firstlick_tauAUC": data.groupby("mouseID")["cumlick_sec_firstlick_tauAUC"].apply(np.mean),  
                                    "cumwaterBW_firstlick_tauAUC": data.groupby("mouseID")["cumwaterBW_firstlick_tauAUC"].apply(np.mean),  
                                    "cumwaterBW_firstbout_tauAUC": data.groupby("mouseID")["cumwaterBW_firstbout_tauAUC"].apply(np.mean),  
    
                                    "lickrate_sec_firstbout_tauAUC": data.groupby("mouseID")["lickrate_sec_firstbout_tauAUC"].apply(np.mean),
                                    "cumlick_sec_firstbout_tauAUC": data.groupby("mouseID")["cumlick_sec_firstbout_tauAUC"].apply(np.mean),
                                                                    
                                    "cumlick_firstlick": data.groupby("mouseID")["cumlick_firstlick"].apply(np.mean),
                                    "cumlick_firstbout": data.groupby("mouseID")["cumlick_firstbout"].apply(np.mean),
        
                                    "tau_cumlick_firstlick":data.groupby("mouseID")["tau_cumlick_firstlick"].apply(np.mean),
                                    "tau_cumlick_firstbout":data.groupby("mouseID")["tau_cumlick_firstbout"].apply(np.mean),
    
                                    "nontau_cumlick_firstlick":data.groupby("mouseID")["nontau_cumlick_firstlick"].apply(np.mean),
                                    "nontau_cumlick_firstbout":data.groupby("mouseID")["nontau_cumlick_firstbout"].apply(np.mean),
    
                                    "cumwaterBW_firstlick": data.groupby("mouseID")["cumwaterBW_firstlick"].apply(np.mean),    
                                    "cumwaterBW_firstbout": data.groupby("mouseID")["cumwaterBW_firstbout"].apply(np.mean),
        
                                    "zscore_firstlick_setnumlickAUC":data.groupby("mouseID")["zscore_firstlick_setnumlickAUC"].apply(np.mean),
                                    "zscore_firstbout_setnumlickAUC":data.groupby("mouseID")["zscore_firstbout_setnumlickAUC"].apply(np.mean),                               
                
                                    "zscore_firstlick_setnumlickDiff":data.groupby("mouseID")["zscore_firstlick_setnumlickDiff"].apply(np.mean),
                                    "zscore_firstbout_setnumlickDiff":data.groupby("mouseID")["zscore_firstbout_setnumlickDiff"].apply(np.mean),                               
                                              
                                    "bout_AUC": data.groupby("mouseID")["bout_AUC"].apply(np.mean),
                                    "bout_size": data.groupby("mouseID")["bout_size"].apply(np.mean),
                                    "bout_AUCperlick": data.groupby("mouseID")["bout_AUCperlick"].apply(np.mean),
                                    "bout_meanzs": data.groupby("mouseID")["bout_meanzs"].apply(np.mean),
        
                                    "interbout_AUC": data.groupby("mouseID")["interbout_AUC"].apply(np.mean),
                                    # "interbout_size": data.groupby("mouseID")["interbout_size"].apply(np.mean),
                                    "interbout_AUCperlick": data.groupby("mouseID")["interbout_AUCperlick"].apply(np.mean),
                                    "interbout_meanzs": data.groupby("mouseID")["interbout_meanzs"].apply(np.mean),
                                   
                                    "all_lick_meanzs": data.groupby("mouseID")["all_lick_meanzs"].apply(np.mean),
                                    "all_bout_meanzs": data.groupby("mouseID")["all_bout_meanzs"].apply(np.mean),
                                    "all_interbout_meanzs": data.groupby("mouseID")["all_interbout_meanzs"].apply(np.mean),
                                   
                                   })#.reset_index()   
        
        #%###################################################
        #   Save data and meandata as .csv files
        #####################################################
        data_to_save     = data[["expID", "BW", "tau_half_firstlick", "tau_half_firstbout", "T_halfdrop_sm_firstbout", "T_halfdrop_sm", "T_halfrise_sm", "T_halfdrop_cumlick_firstbout",
                                 "bout_start_zscore", "T_halfdrop_zscoreDiff_firstbout", "T_halfdrop_zscoreDiff_over_licks",
                                 #"zscore_firstlick_AUC", 
                                 # "lickrate_min_firstlick_AUC", "lickrate_sec_firstlick_AUC","cumlick_sec_firstlick_AUC", 
                                  "cumlick_firstlick", "tau_cumlick_firstlick", "nontau_cumlick_firstlick",
                                  "cumwater_firstlick_BW_AUC", "cumwaterBW_firstlick_tauAUC",  
                                  "zscore_rise_AUC",     "zscore_rise_meanz",
                                  "zscore_firstlick_tauAUC", "zscore_firstlick_tauDiff",
                                  "zscore_firstlick_setnumlickDiff",  "taulick_half_firstlick", 
    
                                   # "lickrate_min_firstlick_BW_AUC","lickrate_sec_firstlick_BW_AUC",
    
                                    # "zscore_firstbout_AUC", #"lickrate_sec_firstbout_AUC","cumlick_sec_firstbout_AUC", 
                                    "cumlick_firstbout", "tau_cumlick_firstbout", "nontau_cumlick_firstbout",
                                    "cumwater_firstbout_BW_AUC", "cumwaterBW_firstbout_tauAUC",
                                    "zscore_firstbout_tauAUC", "zscore_firstbout_tauDiff", 
                                    "tau_half_cumlick_firstbout","tau_half_zscore_over_lick",
                                    # "zscore_firstbout_setnumlickAUC", 
                                    "zscore_firstbout_setnumlickDiff", "taulick_half_firstbout"                                                                  
                                  ]]
        # stimulus = "24hrWD"
        # fluid = "water"
        csv_path         = analysis_path + "/data_" + stimulus + "_" + fluid + "_" + add_to_csv + ".csv"
        data_to_save.to_csv(csv_path, index = True)
        
        meandata_to_save     = meandata[["BW", 
                                         "tau_half_firstlick", "tau_half_firstbout", 
                                         "T_halfdrop_sm_firstbout", "T_halfdrop_sm", "T_halfrise_sm", "T_halfdrop_cumlick_firstbout",
                                         "bout_start_zscore", "T_halfdrop_zscoreDiff_firstbout", "T_halfdrop_zscoreDiff_over_licks",
                                         #"zscore_firstlick_AUC", 
                                         # "lickrate_min_firstlick_AUC", "lickrate_sec_firstlick_AUC","cumlick_sec_firstlick_AUC", 
                                         "cumlick_firstlick", "tau_cumlick_firstlick", "nontau_cumlick_firstlick",
                                         "cumwater_firstlick_BW_AUC", "cumwaterBW_firstlick_tauAUC",              
                                         "zscore_rise_AUC", "zscore_rise_meanz",
                                         "zscore_firstlick_tauAUC", "zscore_firstlick_tauDiff",
                                         "zscore_firstlick_setnumlickDiff", "taulick_half_firstlick",
                                         # "lickrate_min_firstlick_BW_AUC","lickrate_sec_firstlick_BW_AUC",
    
                                        # "zscore_firstbout_AUC", #"lickrate_sec_firstbout_AUC","cumlick_sec_firstbout_AUC", 
                                        "cumlick_firstbout", "tau_cumlick_firstbout", "nontau_cumlick_firstbout",
                                        "cumwater_firstbout_BW_AUC", "cumwaterBW_firstbout_tauAUC",
                                        "zscore_firstbout_tauAUC", "zscore_firstbout_tauDiff", 
                                        "tau_half_cumlick_firstbout","tau_half_zscore_over_lick",
                                        # "zscore_firstbout_setnumlickAUC", 
                                        "zscore_firstbout_setnumlickDiff", "taulick_half_firstbout",
                                        "all_lick_meanzs", "all_bout_meanzs", "all_interbout_meanzs"                                                         
                                  ]]    
    
        csv_path         = analysis_path + "/meandata_" + stimulus + "_" + fluid + "_" + add_to_csv + ".csv"
        meandata_to_save.to_csv(csv_path, index = True)
    
        # df = data["bout_meanzs"] 
        # df.columns = data["expID"]
        df1 = pd.DataFrame()
        for i in range(len(data)):
            df1[data.at[i, "expID"]] = pd.Series(data.at[i, "bout_meanzs"])
            
        csv_path              = analysis_path + "/meanBOUTzs" + "_" + stimulus + "_" + fluid + "_" + add_to_csv + ".csv"
        df1.to_csv(csv_path, index = True)
            
        df2 = pd.DataFrame()
        for i in range(len(data)):
            df2[data.at[i, "expID"]] = pd.Series(data.at[i, "interbout_meanzs"])
            
        csv_path              = analysis_path + "/meanINTERBOUTzs" + "_" + stimulus + "_" + fluid + "_" + add_to_csv + ".csv"
        df2.to_csv(csv_path, index = True)        
        
    #%%#################################################################################
    #   Curve fitting of mean z score trace and cumlick plot
    ####################################################################################      

    zs_values       = meandata["zscore_firstbout_sm"].values
    # zs_mean_sm      = smooth(zs_values.mean(), window_len = 101, window="hanning")
    zs_mean_sm      = zs_values.mean()
    
    x_firstbout = time_sec_firstbout[int(pre_min*60*sampling_rate/ds_factor):int((pre_min+post_min)*60*sampling_rate/ds_factor)]

    y_firstbout = zs_mean_sm[int(pre_min*60*sampling_rate/ds_factor):int((pre_min+post_min)*60*sampling_rate/ds_factor)]
    # valid       = ~np.isnan(y) | np.isnan(y2)
    
    if "GRAB-DA" in path:
        pw0                 = (.05, -2, 0.1, 3) # a guess for slope, intercept, slope, intercept       
    else:
        pw0                 = (-0.005, -1, 0.0001, -2)
    pw_firstbout, cov   = curve_fit(two_lines_max, x_firstbout, y_firstbout, pw0, check_finite = False)
    crossover_firstbout = (pw_firstbout[3] - pw_firstbout[1]) / (pw_firstbout[0] - pw_firstbout[2]) 

    cumlick_values      = meandata["cumlick_sec_firstbout"].values
    cumlick_mean        = cumlick_values.mean()
    
    xlick_firstbout = time_sec_bin[600:]
    ylick_firstbout = cumlick_mean[600:]
        
    pw0 = (.67, 0, 0.1, 150) # a guess for slope, intercept, slope, intercept
    # pwlick_firstbout, cov = curve_fit(two_lines_min, xlick_firstbout, ylick_firstbout, pw0)
    # crossoverlick_firstbout = (pwlick_firstbout[3] - pwlick_firstbout[1]) / (pwlick_firstbout[0] - pwlick_firstbout[2]) 


    bout_start_zscore   = zs_mean_sm[int(pre_min*60*sampling_rate/ds_factor)]
    diff                = np.absolute(zs_mean_sm[int(pre_min*60*sampling_rate/ds_factor):int(pre_min*3*60*sampling_rate/ds_factor)] - bout_start_zscore/3*2)        
    ind_T_halfdrop_sm   = diff.argmin() + int(pre_min*60*sampling_rate/ds_factor)
    T_halfdrop_sm_firstbout               = time_sec_firstbout[ind_T_halfdrop_sm] 

    # T_halfdrop_values            = meandata["T_halfdrop_sm_firstbout"].values
    # T_halfdrop_sm_firstbout      = T_halfdrop_values.mean()
    

    
    #%%###############################################################################
    #   FIGURE. zscore plot for entire timecourse with raw lick rates 
    ################################################################################

    fig, ax     = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 8))
    ax3         = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)          
    
    zs_values       = meandata["zscore_heatmap"].values
    size            = len(meandata)
    title           = stimulus + "_" + " (n=%d)" %(size)        
         
    if size > 1:
        lns3 = ax3.plot(time_sec_raw/60, smooth(zs_values.mean(), window_len = 11),  linewidth = 2, color = palette_mean_GCaMP[0], label = "GCaMP")
        yerr = zs_values.std()/math.sqrt(size-1)
        ax3.fill_between(time_sec_raw/60, zs_values.mean() - yerr, zs_values.mean() + yerr, color = palette_ste_GCaMP[0], alpha = 0.5) ,
    else:
        lns3 = ax3.plot(time_sec_raw/60, smooth(zs_values.mean(), window_len = 11), linewidth = 2, color = palette_mean_GCaMP[0], label = "GCaMP")   
    
    ax3.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
    ax3.axvline (x = water_delay_min, linewidth = 2, linestyle = "--", color = "gray")
    ax3.set_ylabel("z score", fontsize = 32)
    ax3.set_xlabel("time (min)", fontsize = 32)
    plt.xticks(size = 24)
    plt.yticks(size = 24)
    ax3.set_ylim (min_zs, max_zs)
    ax3.set_xlim (-pre_min_heatmap, post_min_heatmap)
    
    ax4                  = ax3.twinx()       
    lickrate_values      = meandata["lickrate_sec_heatmap"].values
    
    if size > 1:
        lns4 = ax4.plot(time_sec_bin_heatmap/60, lickrate_values.mean()[:len(time_sec_bin_heatmap)], linewidth = 1, color ="blue", label = "lick rate") #
        yerr = lickrate_values.std()/math.sqrt(size-1)
        yerr = yerr[:len(time_sec_bin_heatmap)]
        ax4.fill_between(time_sec_bin_heatmap/60, lickrate_values.mean()[:len(time_sec_bin_heatmap)] - yerr, lickrate_values.mean()[:len(time_sec_bin_heatmap)] + yerr, color = "lightgrey", alpha = 0.5) #
    else:
        lns4 = ax4.plot(time_sec_bin_heatmap/60, lickrate_values.mean()[:len(time_sec_bin_heatmap)], linewidth = 1, color ="blue", label = "lick rate") #color = "black",       
    ax4.set_ylabel("lick rate (Hz)", fontsize = 32)
    ax4.set_ylim (min_lickrate, max_lickrate + 8)
    ax4.set_xlim (-pre_min_heatmap, post_min_heatmap)
    
    
    lns     = lns3 + lns4
    labels  = [l.get_label() for l in lns]  
    ax3.legend(lns, labels, loc="upper right", fontsize = 24)

    ax3.set_title(title, fontsize = 32)  
    suptitle = stimulus + "_" + fluid 
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
    if fixed_ylim_mean == 1:
        plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim.png")
        # plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim_noID.pdf", dpi = 1200, format = 'pdf', bbox_inches = 'tight')   
    
    else:
        plt.savefig(analysis_path + "/" + suptitle + "_mean.png")
        # plt.savefig(analysis_path + "/" + suptitle + "_mean_noID.pdf", dpi = 1200, format = 'pdf', bbox_inches = 'tight')
 
    #%%#############################################################
    ###     PSTH for first lick of individual bouts 
    ##############################################################                   
    
    psth_width      = math.floor(sampling_rate*psth_width_sec/ds_factor)
    psth_num_points = 2*psth_width+1
    psth_plot_times = np.linspace(-psth_width_sec, psth_width_sec, psth_num_points, endpoint = True, dtype = None)
    
    for i, file_path in enumerate (all_files):                                   
        psth            = np.empty([len(data.at[i,"bout_size"]), psth_num_points])
        psth[:]         = np.nan
        psth_y_min_list = []
        psth_y_max_list = []
        
        for k, bout_on in enumerate (data.at[i, "bout_on"]):
            bout_on     = bout_on.reshape(1, -1)
            ind         = tree.query(bout_on)[1][0][0] 
            # print (ind + psth_width+1)
            # print(len(data.at[i, "zscore_ds"]))
            # print ("\n")
            
            ## Process without error if PSTH post-width exceeds the recording, fill with nan.
            if ind + psth_width + 1 >= len(data.at[i, "zscore_heatmap"]):  
                last_post_psth                          = len(data.at[i, "zscore_heatmap"]) - ind
                try: 
                    psth[k, : psth_width + last_post_psth]  = data.at[i, "zscore_heatmap"][ind - psth_width : ] - np.nanmedian(data.at[i, "zscore_heatmap"][ind - psth_width : ind]) 
                except:
                    continue
            else:    
                psth[k, :]  = data.at[i, "zscore_heatmap"][ind - psth_width : ind + psth_width +1] - np.nanmedian(data.at[i, "zscore_heatmap"][ind - psth_width : ind])

    
            psth_y_min_list.append(np.nanmin(psth[k, :]))
            psth_y_max_list.append(np.nanmax(psth[k, :]))
        data.at[i, "bout_on_PSTH"] = pd.DataFrame(psth)



# ##################################################
# #   POP open the analysis folder
# ##################################################
os.startfile(path)
