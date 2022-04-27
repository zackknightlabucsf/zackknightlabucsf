# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 17:47:40 2020

This script is for group analysis of non-TTL data.

**** skips photometry_mouseplot.py **** 

1) reads ALL .pkl files under a specified directory generated from photometry_rawplot.py
2) generates mean dFF for the same "mouseID" to generate a mean data for a mouse
3) groups mean data by age and calculates mean for the age group

"""
#%%
from scipy.integrate import simps
from scipy import stats
from scipy import signal
import pandas as pd
import os
import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle as pk
#import seaborn as sns
from datetime import date
#import csv
import sys
import math
import seaborn as sns
# import scipy.stats as ss 
matplotlib.rcParams['font.size'] = 18
#
print("Enter the full path for directory with .pkl files from rawplot: ")
path      = input()

chow  = 0 # 1 if stimulus is chow, 0 if not chow
                        
##  Import chow_g and BW data file into a dataframe
if chow == 1:
   df_chow_g = pd.read_csv(r"C:\Users\heeun\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Mouse_metadata\photometrymiceTTL_wchow.csv")

age_breakpoints = {}
grouped         = {}
group_key_list  = {}
meandFF         = {}
stddFF          = {}


pre_length            = 10  # minutes of pre treatment to plot 
length_cutoff         = 55  # minutes after pre treatment length 

## edit here to set the starting minute for AUC. Default is pre-length
min_start_for_AUC     = pre_length         
min_end_for_AUC       = 4

frame_rate            = 1017.2527
ds_factor1            = 5
ds_factor2            = 5
ds_factor3            = 2
ds_factor             = ds_factor1*ds_factor2*ds_factor3
## ds_factor = 50 reduces 1000 Hz to 50 Hz (50 frames per second)

sm_factor             = 61  ## smoothing over 60 frames (about 1 second, after downsampling)

if (ds_factor == 50) and (sm_factor == 61):
    analysis_path = path + "/analysis"
else:
    analysis_path = path + "/analysis_ds" + str(ds_factor) + "_sm" + str(sm_factor)
    
if not os.path.exists(analysis_path):
    os.makedirs(analysis_path)            


frame_start_for_AUC  = int(min_start_for_AUC*60*frame_rate/ds_factor) 
frame_end_for_AUC    = int((min_start_for_AUC+min_end_for_AUC)*60*frame_rate/ds_factor)

## y range for individual traces and mean
fixed_ylim_ind   = 0    # 1 to use pre-defined y axis, 0 if no 
if fixed_ylim_ind  == 1:
    fixed_ylim_ind_range = [-2, 10]  ## Edit here

## y range for mean and ste 
fixed_ylim_mean    = 0  # 1 to use pre-defined y axis, 0 if no 
if fixed_ylim_mean  == 1:
    fixed_ylim_mean_range = [-10, 20]  ## Edit here    

all_files             = glob.glob(path + "/*.pkl")
file_list             = []
expID_list            = []
tempdata_list         = []
min_length_list       = []


## smoothing 
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


###################################################################################################
#   Input: entire pd.DataFrame from all experiments of a single mouse saved as .pkl file
################################################################################################### 


if chow == 1:
    data = pd.DataFrame(columns = ["expID", "mouseID", "expdate", "BW", "stimulus",
                                   "dFF_ds_sm", "dFF_ds_sm_AUC", "dFF_ds_sm_peak", "time_sec_ds",
                                   "chow_g", "dFF_ds_sm_chow_g", "dFF_ds_sm_chow_BW", 
                                   "dFF_ds_sm_AUC_chow_g", "dFF_ds_sm_AUC_chow_BW"])

    data["dFF_ds_sm_chow_g"]    = data["dFF_ds_sm_chow_g"].astype(object)
    data["dFF_ds_sm_chow_BW"]   = data["dFF_ds_sm_chow_BW"].astype(object)


else:    
    data = pd.DataFrame(columns = ["expID", "mouseID", "expdate", "BW", "stimulus", 
                               "dFF_ds_sm", "dFF_ds_sm_AUC", "dFF_ds_sm_peak", "time_sec_ds"])

data["dFF_ds_sm"]           = data["dFF_ds_sm"].astype(object) 
data["time_sec_ds"]         = data["time_sec_ds"].astype(object)


## Go through each file in the path to save the data into a single dataframe
for i, file_path in enumerate (all_files):
    expID            = file_path[(len(path)+1):-4] ## previous called as "file_name_cut
    file_list.append(expID)  
    inputfile        = open(file_path, 'rb')
    tempdata         = pk.load(inputfile)
    inputfile.close()

    if i == 0:
        stim  = tempdata["stimulus"]
    data.at[i, "expID"]         = expID
    expID_list.append(expID)
    
    mouseID                     = tempdata["mouseID"]
    data.at[i, "mouseID"]       = mouseID
    data.at[i, "expdate"]       = tempdata["date"]
 
    ## downsample all time series by a factor of 50 (1000Hz >> 20Hz. 1ms >> 50ms exposure time )
    dFF                       = tempdata["dFF"]
    temp_trace                = signal.decimate(dFF, ds_factor1)
    temp_trace                = signal.decimate(temp_trace, ds_factor2)
    temp_trace                = signal.decimate(temp_trace, ds_factor3)
    data.at[i, "dFF_ds_sm"]   = smooth(temp_trace, window_len = sm_factor)    
    data.at[i, "time_sec_ds"] = tempdata["time_sec"][::ds_factor] ## take every 50th element from the list

    data.at[i, "dFF_ds_sm_AUC"]     = simps(data.at[i, "dFF_ds_sm"][frame_start_for_AUC:frame_end_for_AUC])   
    data.at[i, "dFF_ds_sm_peak"]    = np.max(data.at[i, "dFF_ds_sm"][frame_start_for_AUC:frame_end_for_AUC])
# 

    min_length_temp   = len(data.at[i, "time_sec_ds"])
    min_length_list.append(min_length_temp)
    
    
    expdate_temp = tempdata["date"]     
    expdate      = date(int(expdate_temp[0:4]), int(expdate_temp[4:6]), int(expdate_temp[6:8])) 

    if chow == 1:
    # find chow_g for the mouse ID and experiment date    
        for j in range(len(df_chow_g)):
            if (df_chow_g.iloc[j]["mouseID"] == mouseID and df_chow_g.iloc[j]["ExpDate"] == int(data.at[i, "expdate"])):
                data.at[i, "chow_g"]  = df_chow_g.iloc[j]["chow_g"]
                data.at[i, "BW"]      = df_chow_g.iloc[j]["BW"]
                break
        try:
            data.at[i, "chow_g"]
        except:
            print ("%s does not have chow(g)" %expID)

        data.at[i, "dFF_ds_sm_chow_g"]   = data.at[i, "dFF_ds_sm"] /data.at[i, "chow_g"]
        data.at[i, "dFF_ds_sm_chow_BW"]  = data.at[i, "dFF_ds_sm_chow_g"]/data.at[i, "BW"]    

        data.at[i, "dFF_ds_sm_AUC_chow_g"]  = simps(data.at[i, "dFF_ds_sm_chow_g"][frame_start_for_AUC:frame_end_for_AUC].clip(min=0), dx=1)
        data.at[i, "dFF_ds_sm_AUC_chow_BW"] = simps(data.at[i, "dFF_ds_sm_chow_BW"][frame_start_for_AUC:frame_end_for_AUC].clip(min=0), dx=1)

print ("files analyzed are: %s" %file_list)        


## If trace length is slightly off, trim the long one
min_length  = np.min(np.array(min_length_list))
for i in np.arange(len(data)):
    data.at[i, "time_sec_ds"]      = data.at[i, "time_sec_ds"][:min_length]
    data.at[i, "dFF_ds_sm"]        = data.at[i, "dFF_ds_sm"][:min_length]


    # data.at[i, "dFF_sm_scaled"] = data.at[i, "dFF_sm_scaled"][:minlendFF]
    if chow == 1:
        data.at[i, "dFF_ds_sm_chow_g"]  = data.at[i, "dFF_ds_sm_chow_g"][:min_length]
        data.at[i, "dFF_ds_sm_chow_BW"] = data.at[i, "dFF_ds_sm_chow_BW"][:min_length]
        
#  Import metadata file and get the start feed time 
if chow == 1:    
    for i, expID in enumerate (file_list):   
      data.at[i, "startfeed_sec"] = np.nan
      for j in range(len(df_chow_g)):   ## look for the identical expdate and mouseID to find BW
          if (str(df_chow_g.at[j, "ExpDate"]) == str(data.at[i, "expdate"])) & (df_chow_g.at[j, "mouseID"] == data.at[i, "mouseID"]):
              timezero_sec   = df_chow_g.at[j, "TimeZero_mm"]*60 + df_chow_g.at[j, "TimeZero_ss"]
              startfeed_temp = df_chow_g.at[j, "StartFeed_mm"]*60 + df_chow_g.at[j, "StartFeed_ss"]
              startfeed_sec  = startfeed_temp - timezero_sec 
              data.at[i, "startfeed_sec"]   = startfeed_sec
              break
          
      if data.at[i, "startfeed_sec"] == np.nan:   ## if no BW data, then find +/- 1 week from the expdate
          print("%s has no chow_g data" %data.at[i, "expID"])


#%%#############################################################################################
#   PROCESS DATA   
#   1) "Groupby" by "mouseID", save the mean dFF, mean AUC and age into the new dataframe
#   2) Then "Groupby" by "age break points", save the mean and std into new variables 
################################################################################################

#   1) "Groupby" by "mouseID", save the mean dFF, mean AUC and age into the new dataframe

if chow == 1:  
    meandata   = pd.DataFrame({"dFF_ds_sm_mean": data.groupby("mouseID")["dFF_ds_sm"].apply(np.mean), 
                           # "dFF_ds_sm_scaled_mean": data.groupby("mouseID")["dFF_ds_sm_scaled"].apply(np.mean),
                           "dFF_ds_sm_mean_chow_g": data.groupby("mouseID")["dFF_ds_sm_chow_g"].apply(np.mean),
                           "dFF_ds_sm_mean_chow_BW": data.groupby("mouseID")["dFF_ds_sm_chow_BW"].apply(np.mean),
                           "dFF_ds_sm_AUC": data.groupby("mouseID")["dFF_ds_sm_AUC"].apply(np.mean),
                           "dFF_ds_sm_AUC_chow_g": data.groupby("mouseID")["dFF_ds_sm_AUC_chow_g"].apply(np.mean),
                           "dFF_ds_sm_AUC_chow_BW": data.groupby("mouseID")["dFF_ds_sm_AUC_chow_BW"].apply(np.mean),
                           "dFF_ds_sm_peak": data.groupby("mouseID")["dFF_ds_sm_peak"].apply(np.mean),
                           "num_trials": data.groupby("mouseID").size()
                           }).reset_index()
    meandFF_chow       = {}
    meandFF_chow_BW    = {}
    
else:
    meandata   = pd.DataFrame({"dFF_ds_sm_mean": data.groupby("mouseID")["dFF_ds_sm"].apply(np.mean), 
                           # "dFF_ds_sm_scaled_mean": data.groupby("mouseID")["dFF_ds_sm_scaled"].apply(np.mean),
                           "dFF_ds_sm_AUC": data.groupby("mouseID")["dFF_ds_sm_AUC"].apply(np.mean),
                           "dFF_ds_sm_peak": data.groupby("mouseID")["dFF_ds_sm_peak"].apply(np.mean),
                           "num_trials": data.groupby("mouseID").size()
                           }).reset_index()   
    
time_sec     = data.at[0, "time_sec_ds"]

#   2) Then "Groupby" by "age break points", save the mean and std into new variables 
# time_sec_raw
## Group the dataframe by age groups
for i in range(len(age_breakpoints)):
    grouped[i]        = meandata.groupby(pd.cut(meandata["age_mo_min"], np.array(age_breakpoints[i])))
    group_key_list[i] = list(grouped[i].groups.keys())
#
    meandFF[i]        = grouped[i]["dFF_ds_sm_mean"].apply(np.mean)
    if chow == 1:
        meandFF_chow[i]     = grouped[i]["dFF_ds_sm_mean_chow_g"].apply(np.mean)
        meandFF_chow_BW[i]  = grouped[i]["dFF_ds_sm_mean_chow_BW"].apply(np.mean)




#%%
################################################################################
#   FIGURE 1-1. dFF plot with individual + mean by 1) three age groups 2) two age groups
################################################################################

################################################################################
#   Set max and mix for y-axis
################################################################################
maxdFF_list  = [np.max(meandata["dFF_ds_sm_mean"][x]) for x in range(len(meandata))]    
mindFF_list  = [np.min(meandata["dFF_ds_sm_mean"][x]) for x in range(len(meandata))]
maxAUC_list  = [np.max(meandata["dFF_ds_sm_AUC"][x]) for x in range(len(meandata))]
maxpeak_list = [np.max(meandata["dFF_ds_sm_peak"][x]) for x in range(len(meandata))]

maxdFF       = np.max(maxdFF_list)
mindFF       = np.min(mindFF_list)
maxAUC       = np.max(maxAUC_list)
maxpeak      = np.max(maxpeak_list)

if fixed_ylim_ind == 1:
    mindFF = fixed_ylim_ind_range[0]
    maxdFF = fixed_ylim_ind_range[1]


################################################################################
rowperpage = len(age_breakpoints)
colperpage = np.max([len(age_breakpoints[x]) for x in range(len(age_breakpoints))]) - 1 

    
fig, ax    = plt.subplots(squeeze = False, figsize = (6*colperpage, 5*rowperpage))       

for l in range(len(age_breakpoints)):
    ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
    ind_list   = [-i-1 for i in range(len(grouped[l]))]
    ind_list.sort(reverse = False)
    
    for i, k in enumerate(ind_list):   
        ax          = plt.subplot2grid((rowperpage, colperpage), (l, i), rowspan=1, colspan=1)       
        maxage      = grouped[l]["age_mo_max"].max()[k]
        minage      = grouped[l]["age_mo_min"].min()[k]   
        size        = grouped[l].size()[k]
#        if size > 1:
#            ste   = stddFF[l][k]/math.sqrt(size-1)
#        else:
#            ste   = stddFF[l][i]*0

        color       = age_palette[i]
        legend      = stim #+ " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
    #    markers, caps, bars = plt.errorbar(time_sec/60, meandFF[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
    #    markers.set_color(color)
    #    [bar.set_alpha(0.5) for bar in bars]
        plt.plot(time_sec/60, meandFF[l][k], linewidth = 2, color = color)

        for j in np.arange(size):
            plt.plot(time_sec/60, grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_mean"].values[j], linewidth = 2, color = "black", alpha = 0.5)

        mouseID_list = grouped[l].get_group(group_key_list[l][i])["mouseID"].values
        for j, mouseID in enumerate(mouseID_list):
            # plt.annotate("%s\n" %mouseID, xy = (time_sec[-1]/60*0.75, maxdFF-(maxdFF-mindFF)/12*(j+2)), fontsize = 12)
            plt.annotate("%s\n" %mouseID, xy = (time_sec[0]/60, maxdFF-(maxdFF-mindFF)/14*(j+1)), fontsize = 10)
        plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
        plt.title(legend, fontsize = 20)
        plt.tick_params(axis = "both", which = "major", labelsize = 18)
        plt.xlabel("time (min)", fontsize = 18)
        plt.ylabel("dFF (%)", fontsize = 18)
        plt.ylim(mindFF*1.05, maxdFF*1.5)
        plt.xlim(-pre_length, length_cutoff)

suptitle = stim + "_dFF" 
plt.suptitle("%s" %suptitle, fontsize = 20)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
if fixed_ylim_ind == 1:
    plt.savefig(analysis_path + "/" + suptitle + "_indv_ylim.png")
    plt.savefig(analysis_path + "/" + suptitle + "_indv_ylim.pdf")
else:
    plt.savefig(analysis_path + "/" + suptitle + "_indv.png")
    plt.savefig(analysis_path + "/" + suptitle + "_indv.pdf")    
plt.show()   


#%%###############################################################################
#   FIGURE 1-2. dFF plot with individual + mean, with no mouseID 
##################################################################################
    
fig, ax    = plt.subplots(squeeze = False, figsize = (6*colperpage, 5*rowperpage))       

for l in range(len(age_breakpoints)):
    ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
    ind_list   = [-i-1 for i in range(len(grouped[l]))]
    ind_list.sort(reverse = False)
    
    for i, k in enumerate(ind_list):   
        ax          = plt.subplot2grid((rowperpage, colperpage), (l, i), rowspan=1, colspan=1)       
        maxage      = grouped[l]["age_mo_max"].max()[k]
        minage      = grouped[l]["age_mo_min"].min()[k]   
        size        = grouped[l].size()[k]
#        if size > 1:
#            ste   = stddFF[l][k]/math.sqrt(size-1)
#        else:
#            ste   = stddFF[l][i]*0

        color       = age_palette[i]
        legend      = stim #+ " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
    #    markers, caps, bars = plt.errorbar(time_sec/60, meandFF[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
    #    markers.set_color(color)
    #    [bar.set_alpha(0.5) for bar in bars]
        plt.plot(time_sec/60, meandFF[l][k], linewidth = 2, color = color)

        for j in np.arange(size):
            plt.plot(time_sec/60, grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_mean"].values[j], linewidth = 2, color = "black", alpha = 0.5)

        mouseID_list = grouped[l].get_group(group_key_list[l][i])["mouseID"].values
#        for j, mouseID in enumerate(mouseID_list):
            # plt.annotate("%s\n" %mouseID, xy = (time_sec[-1]/60*0.75, maxdFF-(maxdFF-mindFF)/12*(j+2)), fontsize = 12)
#            plt.annotate("%s\n" %mouseID, xy = (time_sec[0]/60, maxdFF-(maxdFF-mindFF)/14*(j+1)), fontsize = 10)
        plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
        plt.title(legend, fontsize = 20)
        plt.tick_params(axis = "both", which = "major", labelsize = 18)
        plt.xlabel("time (min)", fontsize = 18)
        plt.ylabel("dFF (%)", fontsize = 18)
        # plt.ylim(mindFF*1.05, maxdFF*1.05) 
        plt.xlim(-pre_length, length_cutoff)

suptitle = stim + "_dFF" 
plt.suptitle("%s" %suptitle, fontsize = 20)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
if fixed_ylim_ind == 1:
    plt.savefig(analysis_path + "/" + suptitle + "_indv_ylim_noID.png")
    plt.savefig(analysis_path + "/" + suptitle + "_indv_ylim_noID.pdf")
else:
    plt.savefig(analysis_path + "/" + suptitle + "_indv_noID.png")
    plt.savefig(analysis_path + "/" + suptitle + "_indv_noID.pdf")   
plt.show()   


#%%###############################################################################
#   FIGURE 2-1. dFF plot with mean + stderr by 1) three age groups 2) two age groups
################################################################################

################################################################################
#   Set max and mix for y-axis
################################################################################
maxdFF_list  = []

for l in range(len(age_breakpoints)):
    ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
    ind_list   = [-i-1 for i in range(len(grouped[l]))]
    ind_list.sort(reverse = False)
    
    for i, k in enumerate(ind_list):   
        size        = grouped[l].size()[k]        
        dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_mean"].values    
        
        if size > 1:
            ste         = dFF_values.std()/math.sqrt(size-1)
        else:
            ste         = [np.nan]*len(dFF_values.mean())

        maxdFF_list.append(np.max(dFF_values.mean()) + np.max(ste))
    
maxdFF       = np.nanmax(maxdFF_list)

if fixed_ylim_mean == 1:
    mindFF = fixed_ylim_mean_range[0]
    maxdFF = fixed_ylim_mean_range[1]


#%%##############################################################################
rowperpage = len(age_breakpoints)
colperpage = np.max([len(age_breakpoints[x]) for x in range(len(age_breakpoints))]) - 1 

fig, ax    = plt.subplots(squeeze = False, figsize = (6*colperpage, 5*rowperpage))       

for l in range(len(age_breakpoints)):
    ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
    ind_list   = [-i-1 for i in range(len(grouped[l]))]
    ind_list.sort(reverse = False)
    
    for i, k in enumerate(ind_list):   
        ax          = plt.subplot2grid((rowperpage, colperpage), (l, i), rowspan=1, colspan=1)       
        maxage      = grouped[l]["age_mo_max"].max()[k]
        minage      = grouped[l]["age_mo_min"].min()[k]   
        size        = grouped[l].size()[k]
           
        dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_mean"].values    

        if size > 1:
            ste     = dFF_values.std()/math.sqrt(size-1)
        else:
            ste     = [np.nan]*len(dFF_values.mean())

        color   = age_palette[i] 
        legend      = stim + " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
    #    markers, caps, bars = plt.errorbar(time_sec/60, meandFF[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
    #    markers.set_color(color)
    #    [bar.set_alpha(0.5) for bar in bars]
        plt.plot(time_sec/60, dFF_values.mean(), linewidth = 2, color = color)
        plt.fill_between(time_sec/60, dFF_values.mean() - ste, dFF_values.mean() + ste, color = age_palette_ste[i], alpha = 0.5)
            
        mouseID_list = grouped[l].get_group(group_key_list[l][i])["mouseID"].values
        for j, mouseID in enumerate(mouseID_list):
            # plt.annotate("%s\n" %mouseID, xy = (time_sec[-1]/60*0.75, maxdFF-(maxdFF-mindFF)/14*(j+1)), fontsize = 10)
            plt.annotate("%s\n" %mouseID, xy = (time_sec[0]/60, maxdFF*1.15-(maxdFF*1.15-mindFF)/14*(j+1)), fontsize = 10)
        plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
        plt.title(legend, fontsize = 20)
        plt.tick_params(axis = "both", which = "major", labelsize = 18)
        plt.xlabel("time (min)", fontsize = 18)
        plt.ylabel("dFF (%)", fontsize = 18)
        plt.ylim(mindFF*1.05, maxdFF*1.2) 
        plt.xlim(-pre_length, length_cutoff)

suptitle = stim + "_dFF" 
plt.suptitle("%s" %suptitle, fontsize = 20)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
if fixed_ylim_mean == 1:
    plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim.png")
    plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim.pdf")
    # plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim.tiff")
    # plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim.jpeg")
else:
    plt.savefig(analysis_path + "/" + suptitle + "_mean.png")
    plt.savefig(analysis_path + "/" + suptitle + "_mean.pdf")
    # plt.savefig(analysis_path + "/" + suptitle + "_mean.tiff")  
    # plt.savefig(analysis_path + "/" + suptitle + "_mean.jpeg") 
plt.show()     


#%%###############################################################################
#   FIGURE 2-2. dFF plot with mean + stderr by 1) three age groups 2) two age groups, 
#                without mouseID
################################################################################

################################################################################
rowperpage = len(age_breakpoints)
colperpage = np.max([len(age_breakpoints[x]) for x in range(len(age_breakpoints))]) - 1 

fig, ax    = plt.subplots(squeeze = False, figsize = (6*colperpage, 5*rowperpage))       

for l in range(len(age_breakpoints)):
    ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
    ind_list   = [-i-1 for i in range(len(grouped[l]))]
    ind_list.sort(reverse = False)
    
    for i, k in enumerate(ind_list):   
        ax          = plt.subplot2grid((rowperpage, colperpage), (l, i), rowspan=1, colspan=1)       
        maxage      = grouped[l]["age_mo_max"].max()[k]
        minage      = grouped[l]["age_mo_min"].min()[k]   
        size        = grouped[l].size()[k]
           
        dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_mean"].values    

        if size > 1:
            ste     = dFF_values.std()/math.sqrt(size-1)
        else:
            ste     = [np.nan]*len(dFF_values.mean())

        color   = age_palette[i] 
        legend      = stim + " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
    #    markers, caps, bars = plt.errorbar(time_sec/60, meandFF[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
    #    markers.set_color(color)
    #    [bar.set_alpha(0.5) for bar in bars]
        plt.plot(time_sec/60, dFF_values.mean(), linewidth = 2, color = color)

        plt.fill_between(time_sec/60, dFF_values.mean() - ste, dFF_values.mean() + ste, color = age_palette_ste[i], alpha = 0.5)
        mouseID_list = grouped[l].get_group(group_key_list[l][i])["mouseID"].values
        # for j, mouseID in enumerate(mouseID_list):
            # plt.annotate("%s\n" %mouseID, xy = (time_sec[-1]/60*0.75, maxdFF-(maxdFF-mindFF)/14*(j+1)), fontsize = 10)
            # plt.annotate("%s\n" %mouseID, xy = (time_sec[0]/60, maxdFF*1.15-(maxdFF*1.15-mindFF)/14*(j+1)), fontsize = 10)
        plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
        plt.title(legend, fontsize = 20)
        plt.tick_params(axis = "both", which = "major", labelsize = 18)
        plt.xlabel("time (min)", fontsize = 18)
        plt.ylabel("dFF (%)", fontsize = 18)
        plt.ylim(mindFF*1.05, maxdFF*1.5)
        plt.ylim(-15, 75)
        # plt.ylim(-3, 8)
        # plt.xlim(-pre_length, length_cutoff)
        # plt.xlim(-5,30)
# 
suptitle = stim + "_dFF" 
plt.suptitle("%s" %suptitle, fontsize = 20)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
if fixed_ylim_mean == 1:
    plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim_noID.png")
    plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim_noID.pdf", dpi = 1200, format = 'pdf', bbox_inches = 'tight')   

else:
    plt.savefig(analysis_path + "/" + suptitle + "_mean_noID.png")
    plt.savefig(analysis_path + "/" + suptitle + "_mean_noID.pdf", dpi = 1200, format = 'pdf', bbox_inches = 'tight')


plt.show()    



# # #%%
# ################################################################################
# #   FIGURE 3-1. SCALED dFF plot with individual + mean by 1) three age groups 2) two age groups
# ################################################################################

# ################################################################################
# #   Set max and mix for y-axis
# ################################################################################
# maxdFF_list  = [np.max(meandata["dFF_ds_sm_scaled_mean"][x]) for x in range(len(meandata))]    
# mindFF_list  = [np.min(meandata["dFF_ds_sm_scaled_mean"][x]) for x in range(len(meandata))]
# # maxAUC_list  = [np.max(meandata["dFF_ds_sm_AUC"][x]) for x in range(len(meandata))]
# # maxpeak_list = [np.max(meandata["dFF_ds_sm_peak"][x]) for x in range(len(meandata))]

# maxdFF       = np.max(maxdFF_list)
# mindFF       = np.min(mindFF_list)
# # maxAUC       = np.max(maxAUC_list)
# # maxpeak      = np.max(maxpeak_list)

# if fixed_ylim_ind == 1:
#     mindFF = fixed_ylim_ind_range[0]
#     maxdFF = fixed_ylim_ind_range[1]


# for l in range(len(age_breakpoints)):
#     ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
#     ind_list   = [-i-1 for i in range(len(grouped[l]))]
#     ind_list.sort(reverse = False)
    
#     for i, k in enumerate(ind_list):         
#         dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_scaled_mean"].values    
        
# ################################################################################
# rowperpage = len(age_breakpoints)
# colperpage = np.max([len(age_breakpoints[x]) for x in range(len(age_breakpoints))]) - 1 

    
# fig, ax    = plt.subplots(squeeze = False, figsize = (6*colperpage, 5*rowperpage))       

# for l in range(len(age_breakpoints)):
#     ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
#     ind_list   = [-i-1 for i in range(len(grouped[l]))]
#     ind_list.sort(reverse = False)
    
#     for i, k in enumerate(ind_list):   
#         ax          = plt.subplot2grid((rowperpage, colperpage), (l, i), rowspan=1, colspan=1)       
#         maxage      = grouped[l]["age_mo_max"].max()[k]
#         minage      = grouped[l]["age_mo_min"].min()[k]   
#         size        = grouped[l].size()[k]
        
#         dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_scaled_mean"].values    

#         if stim_color == 1:
#             color       = palette[palette["stimulus"] == stim]["color"].values[0]
#         else: 
#             color       = age_palette[i]
#         legend      = stim #+ " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
#     #    markers, caps, bars = plt.errorbar(time_sec/60, meandFF[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
#     #    markers.set_color(color)
#     #    [bar.set_alpha(0.5) for bar in bars]
#         plt.plot(time_sec/60, dFF_values.mean(), linewidth = 4, color = color)

#         for j in np.arange(size):
#             plt.plot(time_sec/60, grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_scaled_mean"].values[j], linewidth = 2, color = "black", alpha = 0.5)

#         mouseID_list = grouped[l].get_group(group_key_list[l][i])["mouseID"].values
#         for j, mouseID in enumerate(mouseID_list):
#             # plt.annotate("%s\n" %mouseID, xy = (time_sec[-1]/60*0.75, maxdFF-(maxdFF-mindFF)/12*(j+2)), fontsize = 12)
#             plt.annotate("%s\n" %mouseID, xy = (time_sec[0]/60, maxdFF-(maxdFF-mindFF)/14*(j+1)), fontsize = 10)
#         plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
#         plt.title(legend, fontsize = 20)
#         plt.tick_params(axis = "both", which = "major", labelsize = 18)
#         plt.xlabel("time (min)", fontsize = 18)
#         plt.ylabel("scaled dFF", fontsize = 18)
#         plt.ylim(mindFF*1.05, maxdFF*1.5) 

# suptitle = stim + "_dFF" 
# plt.suptitle("%s" %suptitle, fontsize = 20)
# plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
# if fixed_ylim_ind == 1:
#     plt.savefig(analysis_path + "/" + suptitle + "_indv_ylim_scaled.png")
#     plt.savefig(analysis_path + "/" + suptitle + "_indv_ylim_scaled.pdf")
# else:
#     plt.savefig(analysis_path + "/" + suptitle + "_indv_scaled.png")
#     plt.savefig(analysis_path + "/" + suptitle + "_indv_scaled.pdf")    
# plt.show()   


# #%%###############################################################################
# #   FIGURE 3-2. SCALED dFF plot with individual + mean, with no mouseID 
# ##################################################################################
    
# fig, ax    = plt.subplots(squeeze = False, figsize = (6*colperpage, 5*rowperpage))       

# for l in range(len(age_breakpoints)):
#     ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
#     ind_list   = [-i-1 for i in range(len(grouped[l]))]
#     ind_list.sort(reverse = False)
    
#     for i, k in enumerate(ind_list):   
#         ax          = plt.subplot2grid((rowperpage, colperpage), (l, i), rowspan=1, colspan=1)       
#         maxage      = grouped[l]["age_mo_max"].max()[k]
#         minage      = grouped[l]["age_mo_min"].min()[k]   
#         size        = grouped[l].size()[k]
        
#         dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_scaled_mean"].values    

#         if stim_color == 1:
#             color       = palette[palette["stimulus"] == stim]["color"].values[0]
#         else: 
#             color       = age_palette[i]
#         legend      = stim #+ " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
#     #    markers, caps, bars = plt.errorbar(time_sec/60, meandFF[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
#     #    markers.set_color(color)
#     #    [bar.set_alpha(0.5) for bar in bars]
#         plt.plot(time_sec/60, dFF_values.mean(), linewidth = 4, color = color)

#         for j in np.arange(size):
#             plt.plot(time_sec/60, grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_scaled_mean"].values[j], linewidth = 2, color = "black", alpha = 0.5)

#         mouseID_list = grouped[l].get_group(group_key_list[l][i])["mouseID"].values
# #        for j, mouseID in enumerate(mouseID_list):
#             # plt.annotate("%s\n" %mouseID, xy = (time_sec[-1]/60*0.75, maxdFF-(maxdFF-mindFF)/12*(j+2)), fontsize = 12)
# #            plt.annotate("%s\n" %mouseID, xy = (time_sec[0]/60, maxdFF-(maxdFF-mindFF)/14*(j+1)), fontsize = 10)
#         plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
#         plt.title(legend, fontsize = 20)
#         plt.tick_params(axis = "both", which = "major", labelsize = 18)
#         plt.xlabel("time (min)", fontsize = 18)
#         plt.ylabel("scaled dFF", fontsize = 18)
#         plt.ylim(mindFF*1.05, maxdFF*1.05) 

# suptitle = stim + "_dFF" 
# plt.suptitle("%s" %suptitle, fontsize = 20)
# plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
# if fixed_ylim_ind == 1:
#     plt.savefig(analysis_path + "/" + suptitle + "_indv_ylim_scaled_noID.png")
#     plt.savefig(analysis_path + "/" + suptitle + "_indv_ylim_scaled_noID.pdf")
# else:
#     plt.savefig(analysis_path + "/" + suptitle + "_indv_scaled_noID.png")
#     plt.savefig(analysis_path + "/" + suptitle + "_indv_scaled_noID.pdf")   
# plt.show()   


# #%%###############################################################################
# #   FIGURE 4-1. SCALED dFF plot with mean + stderr by 1) three age groups 2) two age groups
# ################################################################################
# ################################################################################
# #   Set max and mix for y-axis
# ################################################################################
# maxdFF_list  = []

# for l in range(len(age_breakpoints)):
#     ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
#     ind_list   = [-i-1 for i in range(len(grouped[l]))]
#     ind_list.sort(reverse = False)
    
#     for i, k in enumerate(ind_list):   
#         size        = grouped[l].size()[k]        
#         dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_scaled_mean"].values    
        
#         if size > 1:
#             ste         = dFF_values.std()/math.sqrt(size-1)
#         else:
#             ste         = [np.nan]*len(dFF_values.mean())

#         maxdFF_list.append(np.max(dFF_values.mean()) + np.max(ste))
    
# maxdFF       = np.max(maxdFF_list)

# if fixed_ylim_mean == 1:
#     mindFF = fixed_ylim_mean_range[0]
#     maxdFF = fixed_ylim_mean_range[1]


# ################################################################################
# rowperpage = len(age_breakpoints)
# colperpage = np.max([len(age_breakpoints[x]) for x in range(len(age_breakpoints))]) - 1 

# fig, ax    = plt.subplots(squeeze = False, figsize = (6*colperpage, 5*rowperpage))       

# for l in range(len(age_breakpoints)):
#     ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
#     ind_list   = [-i-1 for i in range(len(grouped[l]))]
#     ind_list.sort(reverse = False)
    
#     for i, k in enumerate(ind_list):   
#         ax          = plt.subplot2grid((rowperpage, colperpage), (l, i), rowspan=1, colspan=1)       
#         maxage      = grouped[l]["age_mo_max"].max()[k]
#         minage      = grouped[l]["age_mo_min"].min()[k]   
#         size        = grouped[l].size()[k]
           
#         dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_scaled_mean"].values    

#         if size > 1:
#             ste     = dFF_values.std()/math.sqrt(size-1)
#         else:
#             ste     = [np.nan]*len(dFF_values.mean())
#         if stim_color == 1:
#             color   = palette[palette["stimulus"] == stim]["color"].values[0]
#         else:
#             color   = age_palette[i] 
#         legend      = stim + " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
#     #    markers, caps, bars = plt.errorbar(time_sec/60, meandFF[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
#     #    markers.set_color(color)
#     #    [bar.set_alpha(0.5) for bar in bars]
#         plt.plot(time_sec/60, dFF_values.mean(), linewidth = 4, color = color)
#         if stim_color == 1:
#             plt.fill_between(time_sec/60, dFF_values.mean() - ste, dFF_values.mean() + ste, color = "grey", alpha = 0.5)
#         else:
#             plt.fill_between(time_sec/60, dFF_values.mean() - ste, dFF_values.mean() + ste, color = age_palette_ste[i], alpha = 0.5)
            
#         mouseID_list = grouped[l].get_group(group_key_list[l][i])["mouseID"].values
#         for j, mouseID in enumerate(mouseID_list):
#             # plt.annotate("%s\n" %mouseID, xy = (time_sec[-1]/60*0.75, maxdFF-(maxdFF-mindFF)/14*(j+1)), fontsize = 10)
#             plt.annotate("%s\n" %mouseID, xy = (time_sec[0]/60, maxdFF*1.15-(maxdFF*1.15-mindFF)/14*(j+1)), fontsize = 10)
#         plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
#         plt.title(legend, fontsize = 20)
#         plt.tick_params(axis = "both", which = "major", labelsize = 18)
#         plt.xlabel("time (min)", fontsize = 18)
#         plt.ylabel("dFF (%)", fontsize = 18)
#         plt.ylim(mindFF*1.05, maxdFF*1.2) 

# suptitle = stim + "_dFF" 
# plt.suptitle("%s" %suptitle, fontsize = 20)
# plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
# if fixed_ylim_mean == 1:
#     plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim_scaled.png")
#     plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim_scaled.pdf")
# else:
#     plt.savefig(analysis_path + "/" + suptitle + "_mean_scaled.png")
#     plt.savefig(analysis_path + "/" + suptitle + "_mean_scaled.pdf")
# plt.show()     


# #%%###############################################################################
# #   FIGURE 4-2. SCALED dFF plot with mean + stderr by 1) three age groups 2) two age groups, 
# #                without mouseID
# ################################################################################

# ################################################################################
# rowperpage = len(age_breakpoints)
# colperpage = np.max([len(age_breakpoints[x]) for x in range(len(age_breakpoints))]) - 1 

# fig, ax    = plt.subplots(squeeze = False, figsize = (6*colperpage, 5*rowperpage))       

# for l in range(len(age_breakpoints)):
#     ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
#     ind_list   = [-i-1 for i in range(len(grouped[l]))]
#     ind_list.sort(reverse = False)
    
#     for i, k in enumerate(ind_list):   
#         ax          = plt.subplot2grid((rowperpage, colperpage), (l, i), rowspan=1, colspan=1)       
#         maxage      = grouped[l]["age_mo_max"].max()[k]
#         minage      = grouped[l]["age_mo_min"].min()[k]   
#         size        = grouped[l].size()[k]
           
#         dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_scaled_mean"].values    

#         if size > 1:
#             ste     = dFF_values.std()/math.sqrt(size-1)
#         else:
#             ste     = [np.nan]*len(dFF_values.mean())
#         if stim_color == 1:
#             color   = palette[palette["stimulus"] == stim]["color"].values[0]
#         else:
#             color   = age_palette[i] 
#         legend      = stim + " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
#     #    markers, caps, bars = plt.errorbar(time_sec/60, meandFF[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
#     #    markers.set_color(color)
#     #    [bar.set_alpha(0.5) for bar in bars]
#         plt.plot(time_sec/60, dFF_values.mean(), linewidth = 4, color = color)
#         if stim_color == 1:
#             plt.fill_between(time_sec/60, dFF_values.mean() - ste, dFF_values.mean() + ste, color = "grey", alpha = 0.5)
#         else:
#             plt.fill_between(time_sec/60, dFF_values.mean() - ste, dFF_values.mean() + ste, color = age_palette_ste[i], alpha = 0.5)
            
#         mouseID_list = grouped[l].get_group(group_key_list[l][i])["mouseID"].values
#         # for j, mouseID in enumerate(mouseID_list):
#             # plt.annotate("%s\n" %mouseID, xy = (time_sec[-1]/60*0.75, maxdFF-(maxdFF-mindFF)/14*(j+1)), fontsize = 10)
#             # plt.annotate("%s\n" %mouseID, xy = (time_sec[0]/60, maxdFF*1.15-(maxdFF*1.15-mindFF)/14*(j+1)), fontsize = 10)
#         plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
#         plt.title(legend, fontsize = 20)
#         plt.tick_params(axis = "both", which = "major", labelsize = 18)
#         plt.xlabel("time (min)", fontsize = 18)
#         plt.ylabel("dFF (%)", fontsize = 18)
#         plt.ylim(mindFF*1.05, maxdFF*1.2) 

# suptitle = stim + "_dFF" 
# plt.suptitle("%s" %suptitle, fontsize = 20)
# plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
# if fixed_ylim_mean == 1:
#     plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim_scaled_noID.png")
#     plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim_scaled_noID.pdf")
# else:
#     plt.savefig(analysis_path + "/" + suptitle + "_mean_scaled_noID.png")
#     plt.savefig(analysis_path + "/" + suptitle + "_mean_scaled_noID.pdf")
# plt.show()    


#%%################################################################################
# for chow, plot i) dFF/chow intake (g)
###################################################################################
if chow == 1:
    
    ###############################################################################
    #   Set max and mix for y-axis
    ################################################################################
    maxdFF_list  = []
    mindFF_list  = []
    
    for l in range(len(age_breakpoints)):
        ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
        ind_list   = [-i-1 for i in range(len(grouped[l]))]
        ind_list.sort(reverse = False)
        
        for i, k in enumerate(ind_list):   
            size        = grouped[l].size()[k]        
            dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_mean_chow_g"].values    
            
            if size > 1:
                ste         = dFF_values.std()/math.sqrt(size-1)
            else:
                ste   = [np.nan]*len(dFF_values.mean())
        
            maxdFF_list.append(np.max(dFF_values.mean()) + np.max(ste))
            mindFF_list.append(np.min(dFF_values.mean()) - np.max(ste))
        
    maxdFF       = np.nanmax(maxdFF_list)
    mindFF       = np.nanmin(mindFF_list)
    
    if fixed_ylim_ind == 1:
        mindFF = fixed_ylim_ind_range[0]
        maxdFF = fixed_ylim_ind_range[1]
    
    
    ################################################################################
    rowperpage = len(age_breakpoints)
    colperpage = np.max([len(age_breakpoints[x]) for x in range(len(age_breakpoints))]) - 1 
    
    fig, ax    = plt.subplots(squeeze = False, figsize = (6*colperpage, 5*rowperpage))       
    
    for l in range(len(age_breakpoints)):
        ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
        ind_list   = [-i-1 for i in range(len(grouped[l]))]
        ind_list.sort(reverse = False)
        
        for i, k in enumerate(ind_list):   
            ax          = plt.subplot2grid((rowperpage, colperpage), (l, i), rowspan=1, colspan=1)       
            maxage      = grouped[l]["age_mo_max"].max()[k]
            minage      = grouped[l]["age_mo_min"].min()[k]   
            size        = grouped[l].size()[k]
               
            dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_mean_chow_g"].values    
    
            if size > 1:
                ste     = dFF_values.std()/math.sqrt(size-1)
            else:
                ste     = [np.nan]*len(dFF_values.mean())

            color   = age_palette[i] 
            legend      = stim + " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
        #    markers, caps, bars = plt.errorbar(time_sec/60, meandFF[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
        #    markers.set_color(color)
        #    [bar.set_alpha(0.5) for bar in bars]
            plt.plot(time_sec/60, dFF_values.mean(), linewidth = 2, color = color)

            plt.fill_between(time_sec/60, dFF_values.mean() - ste, dFF_values.mean() + ste, color = age_palette_ste[i], alpha = 0.5)
                
            # mouseID_list = grouped[l].get_group(group_key_list[l][i])["mouseID"].values
            # for j, mouseID in enumerate(mouseID_list):
            #     # plt.annotate("%s\n" %mouseID, xy = (time_sec[-1]/60*0.75, maxdFF-(maxdFF-mindFF)/14*(j+1)), fontsize = 10)
            #     plt.annotate("%s\n" %mouseID, xy = (time_sec[0]/60, maxdFF*1.15-(maxdFF*1.15-mindFF)/14*(j+1)), fontsize = 10)
            plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
            plt.title(legend, fontsize = 20)
            plt.tick_params(axis = "both", which = "major", labelsize = 18)
            plt.xlabel("time (min)", fontsize = 18)
            plt.ylabel("dFF(%)/chow intake(g)", fontsize = 18)
            plt.ylim(mindFF*1.05, maxdFF*1.2) 
    
    suptitle = stim + "_dFF" 
    plt.suptitle("%s" %suptitle, fontsize = 20)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
    if fixed_ylim_mean == 1:
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g_ylim.png")
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g_ylim.pdf")
        # plt.savefig(analysis_path + "/" + suptitle + "_mean_g_ylim.tiff")
        # plt.savefig(analysis_path + "/" + suptitle + "_mean_g_ylim.jpeg")
    else:
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g.png")
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g.pdf")
        # plt.savefig(analysis_path + "/" + suptitle + "_mean_g.tiff")
        # plt.savefig(analysis_path + "/" + suptitle + "_mean_g.jpeg")
    # plt.show()     


    #%%################################################################################
    # for chow, plot ii) dFF/chow intake (g)/BW (g))
    ###################################################################################
        
    ###############################################################################
    #   Set max and mix for y-axis
    ################################################################################
    maxdFF_list  = []
    mindFF_list  = []    
    
    for l in range(len(age_breakpoints)):
        ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
        ind_list   = [-i-1 for i in range(len(grouped[l]))]
        ind_list.sort(reverse = False)
        
        for i, k in enumerate(ind_list):   
            size        = grouped[l].size()[k]        
            dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_mean_chow_BW"].values    
            
            if size > 1:
                ste         = dFF_values.std()/math.sqrt(size-1)
            else:
                ste   = [np.nan]*len(dFF_values.mean())
        
            maxdFF_list.append(np.max(dFF_values.mean()) + np.max(ste))
            mindFF_list.append(np.min(dFF_values.mean()) - np.max(ste))
            
    maxdFF       = np.nanmax(maxdFF_list)
    mindFF       = np.nanmin(mindFF_list)
    
    # if fixed_ylim_mean == 1:
    #     mindFF = fixed_ylim_mean_range[0]
    #     maxdFF = fixed_ylim_mean_range[1]
    
    
    ################################################################################
    rowperpage = len(age_breakpoints)
    colperpage = np.max([len(age_breakpoints[x]) for x in range(len(age_breakpoints))]) - 1 
    
    fig, ax    = plt.subplots(squeeze = False, figsize = (6*colperpage, 5*rowperpage))       
    
    for l in range(len(age_breakpoints)):
        ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
        ind_list   = [-i-1 for i in range(len(grouped[l]))]
        ind_list.sort(reverse = False)
        
        for i, k in enumerate(ind_list):   
            ax          = plt.subplot2grid((rowperpage, colperpage), (l, i), rowspan=1, colspan=1)       
            maxage      = grouped[l]["age_mo_max"].max()[k]
            minage      = grouped[l]["age_mo_min"].min()[k]   
            size        = grouped[l].size()[k]
               
            dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_mean_chow_BW"].values    
    
            if size > 1:
                ste     = dFF_values.std()/math.sqrt(size-1)
            else:
                ste     = [np.nan]*len(dFF_values.mean())

            color   = age_palette[i] 
            legend      = stim + " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
        #    markers, caps, bars = plt.errorbar(time_sec/60, meandFF[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
        #    markers.set_color(color)
        #    [bar.set_alpha(0.5) for bar in bars]
            plt.plot(time_sec/60, dFF_values.mean(), linewidth = 2, color = color)

            plt.fill_between(time_sec/60, dFF_values.mean() - ste, dFF_values.mean() + ste, color = age_palette_ste[i], alpha = 0.5)
                
            mouseID_list = grouped[l].get_group(group_key_list[l][i])["mouseID"].values
            for j, mouseID in enumerate(mouseID_list):
                # plt.annotate("%s\n" %mouseID, xy = (time_sec[-1]/60*0.75, maxdFF-(maxdFF-mindFF)/14*(j+1)), fontsize = 10)
                plt.annotate("%s\n" %mouseID, xy = (time_sec[0]/60, maxdFF*1.15-(maxdFF*1.15-mindFF)/14*(j+1)), fontsize = 10)
            plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
            plt.title(legend, fontsize = 20)
            plt.tick_params(axis = "both", which = "major", labelsize = 18)
            plt.xlabel("time (min)", fontsize = 18)
            plt.ylabel("dFF(%)/chow intake(g) per BW(g)", fontsize = 18)
            plt.ylim(mindFF*1.05, maxdFF*1.2) 
    
    suptitle = stim + "_dFF" 
    plt.suptitle("%s" %suptitle, fontsize = 20)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
    if fixed_ylim_ind == 1:
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g_BW_ylim.png")
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g_BW_ylim.pdf")
    else:
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g_BW.png")
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g_BW.pdf")
    plt.show()     


    ################################################################################
    rowperpage = len(age_breakpoints)
    colperpage = np.max([len(age_breakpoints[x]) for x in range(len(age_breakpoints))]) - 1 
    
    fig, ax    = plt.subplots(squeeze = False, figsize = (6*colperpage, 5*rowperpage))       
    
    for l in range(len(age_breakpoints)):
        ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
        ind_list   = [-i-1 for i in range(len(grouped[l]))]
        ind_list.sort(reverse = False)
        
        for i, k in enumerate(ind_list):   
            ax          = plt.subplot2grid((rowperpage, colperpage), (l, i), rowspan=1, colspan=1)       
            maxage      = grouped[l]["age_mo_max"].max()[k]
            minage      = grouped[l]["age_mo_min"].min()[k]   
            size        = grouped[l].size()[k]
               
            dFF_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_mean_chow_BW"].values    
    
            if size > 1:
                ste     = dFF_values.std()/math.sqrt(size-1)
            else:
                ste     = [np.nan]*len(dFF_values.mean())

            color       = age_palette[i] 
            legend      = stim + " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
        #    markers, caps, bars = plt.errorbar(time_sec/60, meandFF[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
        #    markers.set_color(color)
        #    [bar.set_alpha(0.5) for bar in bars]
            plt.plot(time_sec/60, dFF_values.mean(), linewidth = 2, color = color)

            plt.fill_between(time_sec/60, dFF_values.mean() - ste, dFF_values.mean() + ste, color = age_palette_ste[i], alpha = 0.5)
            plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
            plt.title(legend, fontsize = 20)
            plt.tick_params(axis = "both", which = "major", labelsize = 18)
            plt.xlabel("time (min)", fontsize = 18)
            plt.ylabel("dFF(%)/chow intake(g) per BW(g)", fontsize = 18)
            plt.ylim(mindFF*1.05, maxdFF*1.2) 

    
    suptitle = stim + "_dFF" 
    plt.suptitle("%s" %suptitle, fontsize = 20)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
    if fixed_ylim_ind == 1:
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g_BW_ylim_noID.png")
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g_BW_ylim_noID.pdf")
    else:
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g_BW_noID.png")
        plt.savefig(analysis_path + "/" + suptitle + "_mean_g_BW_noID.pdf")
    plt.show()     


 
#%%##############################################################################
#   FIGURE 3. AUC and Peak 
################################################################################
rowperpage = len(age_breakpoints)
colperpage = 2
bar_width  = 0.75     

x_labels   = {}
ind_list   = {}
fig, ax    = plt.subplots(rowperpage, colperpage, figsize = (4*colperpage, 5*rowperpage))       
for l in range(len(age_breakpoints)):
    ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
    ind_list[l]     = [-i-1 for i in range(len(grouped[l]))]
    ind_list[l].sort(reverse = False)

    x_positions     = np.arange(1, len(ind_list[l])+1, 1)
    x_labels[l]     = []
    ax1  = plt.subplot2grid((rowperpage, colperpage), (l, 0), rowspan=1, colspan=1)       
    for i, k in enumerate(ind_list[l]):   
        AUC_values  = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_AUC"].values
        maxage      = grouped[l]["age_mo_max"].max()[k]
        minage      = grouped[l]["age_mo_min"].min()[k]   
        size        = grouped[l].size()[k]
        legend      = "%d-%d" %(minage, maxage)
        x_labels[l].append(legend)
        ax1.bar(x_positions[i], AUC_values.mean(), yerr = stats.sem(AUC_values), width = bar_width, color = color, alpha = 0.6)
        ax1.scatter([x_positions[i]]*size, AUC_values, c="black")

    ax1.set_xticks(x_positions)
    ax1.set_xticklabels(x_labels[l], fontsize = 16)
    ax1.set_xlabel("age (mo)")
    ax1.tick_params(axis = "both", which = "major", labelsize = 16)
    ax1.set_ylim(0, maxAUC*1.3)
    ax1.set_ylabel("AUC (a.u.)", fontsize = 16)

    ax2  = plt.subplot2grid((rowperpage, colperpage), (l, 1), rowspan=1, colspan=1)       
    for i, k in enumerate(ind_list[l]):   
        peak_values = grouped[l].get_group(group_key_list[l][i])["dFF_ds_sm_peak"].values
        size        = grouped[l].size()[k]
        ax2.bar(x_positions[i], peak_values.mean(), yerr = stats.sem(peak_values), width = bar_width, color = color, alpha = 0.6)
        ax2.scatter([x_positions[i]]*size, peak_values, c="black")

    ax2.set_xticks(x_positions)
    ax2.set_xticklabels(x_labels[l], fontsize = 16)
    ax2.set_xlabel("age (mo)")
    ax2.tick_params(axis = "both", which = "major", labelsize = 16)
    ax2.set_ylim(0, maxpeak*1.3)
    ax2.set_ylabel("Peak fluorescence (a.u.)", fontsize = 16)

suptitle = stim + "_AUC_peak"
plt.suptitle("%s" %suptitle, fontsize = 18)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
plt.savefig(analysis_path + "/" + suptitle + ".png")
plt.show()    



#####################################################
#   Save data and meandata as .csv files
#####################################################
if chow == 1:
    data_to_save     = data[["expID", "dFF_ds_sm_AUC", "dFF_ds_sm_AUC_chow_g", "dFF_ds_sm_AUC_chow_BW", "dFF_ds_sm_peak", "age_mo"]]
else:
    data_to_save     = data[["expID", "dFF_ds_sm_AUC", "dFF_ds_sm_peak", "age_mo"]]    
csv_path         = analysis_path + "/data_" + stim + ".csv"
data_to_save.to_csv(csv_path, index = True)

if chow == 1:
    meandata_to_save = meandata[["mouseID", "dFF_ds_sm_AUC", "dFF_ds_sm_AUC_chow_g", "dFF_ds_sm_AUC_chow_BW", "dFF_ds_sm_peak", "age_mo_min", "age_mo_max", "num_trials"]]
    
else:
    meandata_to_save = meandata[["mouseID", "dFF_ds_sm_AUC", "dFF_ds_sm_peak", "age_mo_min", "age_mo_max", "num_trials"]]    
       
csv_path         = analysis_path + "/meandata_" + stim + ".csv"
meandata_to_save.to_csv(csv_path, index = True)

