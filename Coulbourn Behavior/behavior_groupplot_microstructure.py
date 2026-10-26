# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:23:30 2020

Group-level analysis of mouse drinking behavior and licking microstructure.

Loads processed behavioral data (.pkl), groups mice by age and sex,
and averages repeated recordings from individual mice.

Calculates cumulative licking, estimated water intake, lick rates,
first-lick and first-bout latencies, and drinking-bout microstructure,
including bout number, size, duration, and interbout intervals.

Generates group comparison plots and exports summary measurements
as CSV, pickle, PNG, and PDF files.


"""
#%%

import os
import glob
import pickle as pk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats as ss 

import itertools
import subprocess
#import csv


### To input the path from the command line
print("Enter the full path for the directory with .pkl files from behavior_rawplot.py: ")
path             = input()


# stim      = "24hrWD"
stim      = "2MNaCl"
# stim      = "PEG"
# stim      = "BL"
# stim      = "FD"


# feeder     = 1   # 1 if feeder was used, 0 if feeder was not used
feeder     = 0

lickometer = 1  # 1 if lickometer was used, 0 if lickometer was not used

age_cutoff     = 32  # maximum age (months)

lick_start_AUC  = 0 # mins
# lick_end_AUC    = 1440/2
lick_end_AUC    = 30
# lick_end_AUC      = 4*60

feed_start_AUC  = 0 # mins
feed_end_AUC    = 1440/2
# feed_end_AUC    = 30
# feed_end_AUC      = 4*60

# length_cutoff  = 1440 # mins = 24hrs
# length_cutoff   = 720
length_cutoff = 30
# length_cutoff = 4*60  # 4hrs (for FD then feeder+lickometer experiment)

##  Import meatadata file into a dataframe (mouseID, expdate, BW)
BWpath = r"D:\Data_Dropbox_Sync\Mouse_metadata\allmiceBW.csv"
BWdata = pd.read_csv(BWpath)    

#if stim == "FD":
if feeder == 1:    
    # foodspill_g_path =  r"C:\Users\heeun\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Mouse_metadata\foodspill_g.csv"
    foodspill_g_path = r"D:\Data_Dropbox_Sync\Mouse_metadata\foodspill_g.csv"
    foodspill_g_data = pd.read_csv(foodspill_g_path)
    
analysis_path = path + "/analysis_/" + stim + "_" + str(length_cutoff) + "min_microstructure"
if not os.path.exists(analysis_path):
    os.makedirs(analysis_path)   
#
palette        = ["dimgray", "maroon"]
palette2       = ["silver", "pink"]


male_list   = ["M3", "M26"]
female_list = ["F3", "F26"]  

sexes       = ["M", "F"]


df_lickvol = {}
for i, sex_age in enumerate (male_list):
#    df_lickvol[sex_age] = 1.205697403329695/1000 
    df_lickvol[sex_age] = 1.21/1000

for i, sex_age in enumerate (female_list):
    df_lickvol[sex_age] = 1.21/1000
    
pellet_g = 0.018 
lick_g   = 1.21/1000
    


###########################################################################
###     Define functions      
###########################################################################

def adjust_age (age):
#    if 29 <= age <= 31:
#        adjusted_age = 30
    if 25 <= age <= 32: # 25, 26, 27 months old
        adjusted_age = 26                         
#    if 21 <= age <= 23: # 21, 22, 23 months old
#        adjusted_age = 22               
#    elif 17 <= age <= 19: # 17, 18, 19 months old
#        adjusted_age = 18                           
    elif 2  <= age <= 7: # 2, 3, 4 months old
        adjusted_age = 3                
    else: 
        adjusted_age = 0
    return (adjusted_age)
        
def lick_rate_by_sec (lick_on, time_sec):

    lick_on_ceil = np.ceil(lick_on)
    lick_on_ceil = lick_on_ceil.astype(int)
        
    lick_rate_sec_ = np.histogram(lick_on_ceil, bins = time_sec)   
    lick_rate_sec  = np.concatenate([[0], lick_rate_sec_[0]])
    lick_rate_cum  = np.cumsum(lick_rate_sec)

    return (lick_rate_sec, lick_rate_cum)

## Bout analysis
def bout_analysis (ontimes, bout_interval_cutoff, bout_size_cutoff):

    # Define time (s) that separates bouts and minimum number of licks for a bout
    # bout_interval_cutoff = amount of time (sec) between two licks that terminates a bout (I use 1 sec)
    # bout_size_cutoff = minimum # of licks to be considered a bout (I use 10)

    # while loop finds bout ons, offs, and licks per bout
    bout_ons        = np.asarray([])
    bout_offs       = np.asarray([])
    bout_size       = np.asarray([])
    bout_lengths    = np.asarray([])

    k = 0 # index for while loop
    while k < (len(ontimes) - bout_size_cutoff):
        licks     = 1
        while (ontimes[k+licks] - ontimes[k+licks-1]) < bout_interval_cutoff:
            licks += 1
            if (k+licks) >= len(ontimes):
                break        
        if licks >= bout_size_cutoff:
            bout_ons     = np.append(bout_ons, ontimes[k])
            bout_offs    = np.append(bout_offs, ontimes[k+licks - 1])
            bout_size    = np.append(bout_size, licks)
            bout_lengths = np.append(bout_lengths, ontimes[k+licks-1] - ontimes[k])
        k = k + licks      
    return (bout_ons, bout_offs, bout_size, bout_lengths)
         

all_files             = glob.glob(path + "/*.pkl")    
    
file_list             = []
data                  = {}

###########################################################################
###     Process data     
###########################################################################

# Go through each file in the directory to get .pkl data
for i, file_path in enumerate (all_files):
    expID        = file_path[(len(path)+1):-4] ## previous called as "file_name_cut
    stimulus     = expID.split("_")[4]
    length       = int((expID.split("_")[-1]).split("min")[0])
    age          = int((expID.split("_")[2]).split("mo")[0])
    if stimulus  == stim: 
        if length >= length_cutoff:  ## select files only for given stimulus and experiment duration longer than cutoff
            if (age-1 <= age_cutoff): 
                adjusted_age = adjust_age(age)
                if (adjusted_age != 0):
                    file_list.append(expID)  
                    inputfile    = open(file_path, 'rb')
                    data[expID]  = pk.load(inputfile)
                    inputfile.close()                      

time_sec    = np.arange(1, length_cutoff*60+1)  
time_min    = np.arange(1, length_cutoff*60+1, 60)

df_BW                       = {}

if lickometer == 1:
    df_firstlick_temp            = {}
    df_firstbout_temp            = {}
    df_lick_cum_raw_temp         = {}
    df_lick_cum_firstlick_temp   = {}
    df_lick_cum_firstbout_temp   = {}
    df_lick_rate_raw_temp        = {}
    df_lick_rate_firstlick_temp  = {}
    df_lick_rate_firstbout_temp  = {}
    df_lick_rate_raw_min_temp    = {}
    df_lick_rate_firstlick_min_temp = {}
    df_lick_rate_firstbout_min_temp = {}
    
    df_lick_on_firstlick         = {}
    
    df_firstlick                 = {}
    df_firstbout                 = {}
    df_lick_cum_raw              = {}
    df_lick_cum_firstlick        = {}
    df_lick_cum_firstbout        = {}
    df_lick_rate_raw             = {}
    df_lick_rate_firstlick       = {}
    df_lick_rate_firstbout       = {}
    df_lick_rate_raw_min         = {}
    df_lick_rate_firstlick_min   = {}
    df_lick_rate_firstbout_min   = {}
    
    max_firstlick                = np.array([0])
    max_firstbout                = np.array([0])
    max_cumlick                  = np.array([0])
    max_cumwater_g               = np.array([0])
    max_cumlick_BW               = np.array([0])
    max_cumwater_food_g          = np.array([0])
        
    df_lick_cum_raw_temp_BW     = {}
    df_lick_cum_raw_BW          = {}
    df_lick_rate_raw_temp_BW    = {}
    df_lick_rate_raw_BW         = {}
    df_lick_rate_raw_min_BW     = {}
    
    df_lick_cum_firstlick_temp_BW  = {}
    df_lick_cum_firstlick_BW       = {}
    df_lick_rate_firstlick_temp_BW = {}
    df_lick_rate_firstlick_BW      = {}
    df_lick_rate_firstlick_min_BW  = {}
    
if feeder == 1:
    df_feed_cum_raw_temp     = {}
    df_feed_cum_raw          = {}
    df_feed_rate_raw_temp    = {}
    df_feed_rate_raw         = {}
    df_feed_rate_raw_min_temp= {}
    df_feed_rate_raw_min     = {}
    df_feed_on_raw           = {}
    
    max_cumfeed              = np.array([0])
    max_cumfood_g            = np.array([0])
    max_cumfood_g_BW         = np.array([0])
    
    df_feed_cum_raw_temp_BW     = {}
    df_feed_cum_raw_BW          = {}
    df_feed_rate_raw_temp_BW    = {}
    df_feed_rate_raw_BW         = {}
    df_feed_rate_raw_min_temp_BW = {}
    df_feed_rate_raw_min_BW     = {}

#    if stim == "FD":
        
    df_food_g                       = {}
    df_food_cum_corrected           = {}
    df_food_cum_corrected_temp      = {}
    df_food_cum_corrected_BW        = {}
    df_food_cum_corrected_temp_BW   = {}

    if lickometer == 1:
        df_lick_cum_raw_temp_BW     = {}
        df_lick_cum_raw_BW          = {}
        df_lick_rate_raw_temp_BW    = {}
        df_lick_rate_raw_BW         = {}
        df_lick_rate_raw_min_BW     = {}
    
        df_lick_per_feed         = {}
        df_lick_on_raw           = {}
        max_cumlick_perfeed      = np.array([0])
    
        df_water_to_food_g              = {}
        df_water_to_food_g_temp         = {}
        df_water_to_food_g_BW           = {}
        df_water_to_food_g_BW_temp      = {}

# Initialize dataseries and dataframes
for sex_list in [male_list, female_list]: 
    for i, sex_age in enumerate (sex_list):  
        df_BW[sex_age]                   = pd.Series(0, index = ["dummy"])
        
        if lickometer == 1:
            df_firstlick[sex_age]           = pd.Series(0, index = ["dummy"])
            df_firstbout[sex_age]           = pd.Series(0, index = ["dummy"])
    
            df_lick_cum_raw_temp[sex_age]           = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_lick_cum_firstlick_temp[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_lick_cum_firstbout_temp[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
    
            df_lick_cum_raw[sex_age]         = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_lick_cum_firstlick[sex_age]   = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_lick_cum_firstbout[sex_age]   = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
    
            df_lick_rate_raw_temp[sex_age]          = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_lick_rate_firstlick_temp[sex_age]    = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])           
            df_lick_rate_firstbout_temp[sex_age]    = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
    
            df_lick_rate_raw[sex_age]          = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_lick_rate_firstlick[sex_age]    = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])           
            df_lick_rate_firstbout[sex_age]    = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            
            df_lick_rate_raw_min_temp[sex_age]           = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])
            df_lick_rate_firstlick_min_temp[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])
            df_lick_rate_firstbout_min_temp[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])       
            
            df_lick_rate_raw_min[sex_age]           = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])
            df_lick_rate_firstlick_min[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])
            df_lick_rate_firstbout_min[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])
           
            df_lick_on_firstlick[sex_age]    = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_lick_cum_firstlick_temp_BW[sex_age] = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_lick_cum_firstlick_BW[sex_age]      = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            
            df_lick_rate_firstlick_temp_BW[sex_age]= pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_lick_rate_firstlick_BW[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_lick_rate_firstlick_min_BW[sex_age] = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])

        if feeder == 1:
            df_feed_cum_raw_temp[sex_age]    = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_feed_cum_raw[sex_age]         = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])

            df_feed_rate_raw_temp[sex_age]      = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_feed_rate_raw[sex_age]           = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_feed_rate_raw_min_temp[sex_age]  = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])
            df_feed_rate_raw_min[sex_age]       = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])

            df_feed_on_raw[sex_age]          = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_feed_cum_raw_temp_BW[sex_age] = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])            
            df_feed_cum_raw_BW[sex_age]      = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])

            df_feed_rate_raw_temp_BW[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy']) 
            df_feed_rate_raw_BW[sex_age]          = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_feed_rate_raw_min_temp_BW[sex_age] = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])
            df_feed_rate_raw_min_BW[sex_age]      = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])
            
#            if stim == "FD":
            df_food_g[sex_age]                      = pd.Series(0, index = ["dummy"])  
            df_food_cum_corrected[sex_age]          = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_food_cum_corrected_temp[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_food_cum_corrected_BW[sex_age]       = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
            df_food_cum_corrected_temp_BW[sex_age]  = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
  

            if lickometer == 1:
                df_lick_per_feed[sex_age]        = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
    
                df_lick_on_raw[sex_age]          = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
    #            df_BW[sex_age]                   = pd.Series(0, index = ["dummy"])
                
                df_lick_cum_raw_temp_BW[sex_age] = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
                df_lick_cum_raw_BW[sex_age]      = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
                
                df_lick_rate_raw_temp_BW[sex_age]= pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
                df_lick_rate_raw_BW[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])
                df_lick_rate_raw_min_BW[sex_age] = pd.DataFrame(np.nan, index = np.arange(length_cutoff), columns = ['dummy'])
                  
                df_water_to_food_g[sex_age]          = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])                 
                df_water_to_food_g_temp[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy']) 
                df_water_to_food_g_BW[sex_age]       = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])                   
                df_water_to_food_g_BW_temp[sex_age]  = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ['dummy'])                   
                                    
              

#%% Read data and save it to the dataseries and dataframes by age group and by gender
for i, expID in enumerate (file_list):   
    adjusted_age = adjust_age(data[expID]["age"])
#    print ("processing %s" %expID)
    
    data[expID]["BW"] = np.nan
    for j in range(len(BWdata)):   ## look for the identical expdate and mouseID to find BW
        if (str(BWdata.at[j, "expdate"]) == str(data[expID]["date"])) & (BWdata.at[j, "mouseID"] == data[expID]["mouseID"]):
            data[expID]["BW"]  = BWdata.at[j, "BW"]
            break
        
    if data[expID]["BW"] == np.nan:   
        print("%s has no BW data" %expID)
        
#    if stim == "FD":
    if feeder == 1:
        for j in range(len(foodspill_g_data)):   ## look for the identical expdate and mouseID to find BW
#            expdate = foodspill_g_data.at[j, "expID"].split("_")[1]
#            mouseID = foodspill_g_data.at[j, "expID"].split("_")[0]
            expdate = foodspill_g_data.at[j, "expdate"]
            mouseID = foodspill_g_data.at[j, "mouseID"]
            
            if (str(expdate) == str(data[expID]["date"])) & (mouseID == data[expID]["mouseID"]):
                data[expID]["foodspill_g"]  = foodspill_g_data.at[j, "foodspill_g"]
                break       

        if data[expID]["foodspill_g"] == np.nan:   
            print("%s has no foodspill_g data" %expID)
        
            
    for sex in sexes:
        if data[expID]["sex"]  == sex:  
            # get cumulative licks  
            sex_age                                     = sex + str(adjusted_age)              
            BW_temp                                     = pd.Series(data[expID]["BW"], index = [expID])
            df_BW[sex_age]                              = df_BW[sex_age].append(BW_temp)
            time_sec_from_pkl                              = np.arange(1, data[expID]["totalmins"]*60+1)   
            time_min_from_pkl                              = np.arange(1, (data[expID]["totalmins"]+1)*60, 60)   


            if lickometer == 1: 
                firstlick_temp                              = pd.Series(data[expID]["first_lick_sec"], index = [expID])
                df_firstlick[sex_age]                       = df_firstlick[sex_age].append(firstlick_temp) 
                firstbout_temp                              = pd.Series(data[expID]["first_bout_sec"], index = [expID])
                df_firstbout[sex_age]                       = df_firstbout[sex_age].append(firstbout_temp)    
     
                lickrate_sec_raw, lickrate_cum_raw             = lick_rate_by_sec(data[expID]["lickontimes_raw"], time_sec_from_pkl)   
                lickrate_sec_firstlick, lickrate_cum_firstlick = lick_rate_by_sec(data[expID]["lickontimes_firstlick"], time_sec_from_pkl)
                lickrate_sec_firstbout, lickrate_cum_firstbout = lick_rate_by_sec(data[expID]["lickontimes_firstbout"], time_sec_from_pkl)              
                
                lickrate_min_raw, dump                         = lick_rate_by_sec(data[expID]["lickontimes_raw"], time_min) 
                lickrate_min_firstlick, dump                   = lick_rate_by_sec(data[expID]["lickontimes_firstlick"], time_min) 
                lickrate_min_firstbout, dump                   = lick_rate_by_sec(data[expID]["lickontimes_firstbout"], time_min) 
                
                df_lick_cum_raw_temp[sex_age][expID]           = pd.DataFrame(lickrate_cum_raw+1) 
                df_lick_cum_firstlick_temp[sex_age][expID]     = pd.DataFrame(lickrate_cum_firstlick+1) 
                df_lick_cum_firstbout_temp[sex_age][expID]     = pd.DataFrame(lickrate_cum_firstbout+1) 
           
                df_lick_rate_raw_temp[sex_age][expID]          = pd.DataFrame(lickrate_sec_raw+1)
                df_lick_rate_firstlick_temp[sex_age][expID]    = pd.DataFrame(lickrate_sec_firstlick+1)            
                df_lick_rate_firstbout_temp[sex_age][expID]    = pd.DataFrame(lickrate_sec_firstbout+1)
                
                df_lick_rate_raw_min_temp[sex_age][expID]      = pd.DataFrame(lickrate_min_raw)
                df_lick_rate_firstlick_min_temp[sex_age][expID]= pd.DataFrame(lickrate_min_firstlick)
                df_lick_rate_firstbout_min_temp[sex_age][expID]= pd.DataFrame(lickrate_min_firstbout)
                 
                df_lick_on_firstlick[sex_age][expID]           = pd.DataFrame(data[expID]["lickontimes_firstlick"])   
               
                df_lick_cum_firstlick_temp_BW[sex_age][expID]  = df_lick_cum_firstlick_temp[sex_age][expID]*df_lickvol[sex_age]/(df_BW[sex_age][expID]**0.667) 
                df_lick_rate_firstlick_temp_BW[sex_age][expID] = df_lick_rate_firstlick_temp[sex_age][expID]*df_lickvol[sex_age]/(df_BW[sex_age][expID]**0.667)

            if feeder == 1:  
                feedrate_sec_raw, feedrate_cum_raw         = lick_rate_by_sec(data[expID]["feedontimes_raw"], time_sec_from_pkl)  
                df_feed_cum_raw_temp[sex_age][expID]       = pd.DataFrame(feedrate_cum_raw+1) 
                df_feed_rate_raw_temp[sex_age][expID]      = pd.DataFrame(feedrate_sec_raw+1)
                feedrate_min_raw, dump                     = lick_rate_by_sec(data[expID]["feedontimes_raw"], time_min_from_pkl)
                df_feed_rate_raw_min_temp[sex_age][expID]  = pd.DataFrame(feedrate_min_raw+1)
                
#                if stim == "FD":
                raw_food_g_temp                            = feedrate_cum_raw[-1]*pellet_g
                corrected_food_g_temp                      = pd.Series(raw_food_g_temp - data[expID]["foodspill_g"], index = [expID])
                df_food_g[sex_age]                         = df_food_g[sex_age].append(corrected_food_g_temp)

                df_feed_cum_raw_temp_BW[sex_age][expID]      = df_feed_cum_raw_temp[sex_age][expID]/df_BW[sex_age][expID]
                df_food_cum_corrected_temp[sex_age][expID]    = df_feed_cum_raw_temp[sex_age][expID]*pellet_g/raw_food_g_temp*corrected_food_g_temp[expID]
                df_food_cum_corrected_temp_BW[sex_age][expID] = df_food_cum_corrected_temp[sex_age][expID]/df_BW[sex_age][expID]

                if lickometer == 1:              
                    df_water_to_food_g_temp[sex_age][expID]    = df_lick_cum_raw_temp[sex_age][expID]*lick_g/df_food_g[sex_age][expID] 
                    df_lick_cum_raw_temp_BW[sex_age][expID]    = df_lick_cum_raw_temp[sex_age][expID]*lick_g/(df_BW[sex_age][expID]**0.667) 
    #                df_lick_rate_raw_temp_BW[sex_age][expID]   = df_lick_rate_raw_temp[sex_age][expID]/df_BW[sex_age][expID]
    
    #                df_feed_rate_raw_temp_BW[sex_age][expID]   = df_feed_rate_raw_temp[sex_age][expID]/df_BW[sex_age][expID]
    #                df_feed_rate_raw_min_temp_BW[sex_age][expID]   = df_feed_rate_raw_min_temp[sex_age][expID]/df_BW[sex_age][expID]
                    df_water_to_food_g_BW_temp[sex_age][expID] = df_lick_cum_raw_temp_BW[sex_age][expID]/df_food_cum_corrected_temp_BW[sex_age][expID]  
                    
    if lickometer == 1:        
        max_firstlick = max(max_firstlick, max(df_firstlick[sex_age]))
        max_firstbout = max(max_firstbout, max(df_firstbout[sex_age])) 
    
#for sex_age in male_list + female_list:
#    print ("keys in df_lick_cum_firstlick_temp[%s] are:" %sex_age)
#    print (df_lick_cum_firstlick_temp[sex_age].keys())

#%%#######################################################################################
###     TAKE THE MEAN if there is any other technical replicates from the same mouse     
##########################################################################################
mouseID_sex_age_list  = []
num_replicates        = 0

for i, expID in enumerate (file_list):  
#    print ("processing %s" %expID)
    sex_age                 = data[expID]["sex"] + str(adjust_age(data[expID]["age"]))  
    mouseID_sex_age         = data[expID]["mouseID"] + "_" + sex_age  
    if mouseID_sex_age in mouseID_sex_age_list:
        continue   # skip to the next loop if this mouse has been updated with replicates already
    mouseID_sex_age_list.append(mouseID_sex_age) 

    if lickometer == 1:
    #    firstlick_temp                  = pd.DataFrame(df_firstlick_temp[sex_age][expID], index = [expID])    
        lick_cum_raw_temp               = pd.DataFrame(df_lick_cum_raw_temp[sex_age][expID], columns = [expID])      
        lick_cum_firstlick_temp         = pd.DataFrame(df_lick_cum_firstlick_temp[sex_age][expID], columns = [expID]) 
        lick_cum_firstbout_temp         = pd.DataFrame(df_lick_cum_firstbout_temp[sex_age][expID], columns = [expID]) 
        lick_rate_raw_temp              = pd.DataFrame(df_lick_rate_raw_temp[sex_age][expID], columns = [expID])      
        lick_rate_firstlick_temp        = pd.DataFrame(df_lick_rate_firstlick_temp[sex_age][expID], columns = [expID]) 
        lick_rate_firstbout_temp        = pd.DataFrame(df_lick_rate_firstbout_temp[sex_age][expID], columns = [expID])
     
        lick_rate_raw_min_temp          = pd.DataFrame(df_lick_rate_raw_min_temp[sex_age][expID], columns = [expID])
        lick_rate_firstlick_min_temp    = pd.DataFrame(df_lick_rate_firstlick_min_temp[sex_age][expID], columns = [expID])
        lick_rate_firstbout_min_temp    = pd.DataFrame(df_lick_rate_firstbout_min_temp[sex_age][expID], columns = [expID])
    #
        lick_cum_firstlick_temp_BW        = pd.DataFrame(df_lick_cum_firstlick_temp_BW[sex_age][expID], columns = [expID])
        lick_rate_firstlick_temp_BW       = pd.DataFrame(df_lick_rate_firstlick_temp_BW[sex_age][expID], columns = [expID])
        lick_rate_firstlick_min_temp_BW   = pd.DataFrame(df_lick_rate_firstlick_temp_BW[sex_age][expID], columns = [expID])

    if feeder == 1:
        feed_cum_raw_temp           = pd.DataFrame(df_feed_cum_raw_temp[sex_age][expID], columns = [expID])  
        feed_rate_raw_temp          = pd.DataFrame(df_feed_rate_raw_temp[sex_age][expID], columns = [expID])   
        feed_rate_raw_min_temp      = pd.DataFrame(df_feed_rate_raw_min_temp[sex_age][expID], columns = [expID])   
        
#        lick_rate_raw_temp_BW       = pd.DataFrame(df_lick_rate_raw_temp_BW[sex_age][expID], columns = [expID])
#        lick_rate_raw_min_temp_BW   = pd.DataFrame(df_lick_rate_raw_temp_BW[sex_age][expID], columns = [expID])
        feed_cum_raw_temp_BW        = pd.DataFrame(df_feed_cum_raw_temp_BW[sex_age][expID], columns = [expID])
#        feed_rate_raw_temp_BW       = pd.DataFrame(df_feed_rate_raw_temp_BW[sex_age][expID], columns = [expID])
#        feed_rate_raw_min_temp_BW   = pd.DataFrame(df_feed_rate_raw_temp_BW[sex_age][expID], columns = [expID])
        food_cum_corrected_temp    = pd.DataFrame(df_food_cum_corrected_temp[sex_age][expID], columns = [expID])
        food_cum_corrected_temp_BW = pd.DataFrame(df_food_cum_corrected_temp_BW[sex_age][expID], columns = [expID])  

        if lickometer == 1:
            lick_cum_raw_temp_BW       = pd.DataFrame(df_lick_cum_raw_temp_BW[sex_age][expID], columns = [expID])
            water_to_food_g_temp       = pd.DataFrame(df_water_to_food_g_temp[sex_age][expID], columns = [expID])
            water_to_food_g_BW_temp    = pd.DataFrame(df_water_to_food_g_BW_temp[sex_age][expID], columns = [expID])
                
    for j, expID2 in enumerate (file_list):
        mouseID_sex_age2   = data[expID2]["mouseID"] + "_" + data[expID2]["sex"] + str(adjust_age(data[expID2]["age"]))
        if ((mouseID_sex_age == mouseID_sex_age2) and (expID != expID2)):            
#            firstlick_temp[expID2]              = pd.Series(df_firstlick_temp[sex_age][expID2]) 
            if lickometer == 1:
                lick_cum_raw_temp[expID2]           = pd.Series(df_lick_cum_raw_temp[sex_age][expID2])
                lick_cum_firstlick_temp[expID2]     = pd.Series(df_lick_cum_firstlick_temp[sex_age][expID2]) 
                lick_cum_firstbout_temp[expID2]     = pd.Series(df_lick_cum_firstbout_temp[sex_age][expID2])

#            lick_rate_raw_temp[expID2]           = pd.Series(df_lick_rate_raw_temp[sex_age][expID2])
#            lick_rate_firstlick_temp[expID2]     = pd.Series(df_lick_rate_firstlick_temp[sex_age][expID2]) 
#            lick_rate_firstbout_temp[expID2]     = pd.Series(df_lick_rate_firstbout_temp[sex_age][expID2])                       
#
#            lick_rate_raw_min_temp[expID2]           = pd.Series(df_lick_rate_raw_min_temp[sex_age][expID2])
#            lick_rate_firstlick_min_temp[expID2]     = pd.Series(df_lick_rate_firstlick_min_temp[sex_age][expID2]) 
#            lick_rate_firstbout_min_temp[expID2]     = pd.Series(df_lick_rate_firstbout_min_temp[sex_age][expID2])                       

                print ("%s has technical replicates" %mouseID_sex_age) 
        
                lick_cum_firstlick_temp_BW[expID2]        = pd.Series(df_lick_cum_firstlick_temp_BW[sex_age][expID2])
#            lick_rate_firstlick_temp_BW[expID2]       = pd.Series(df_lick_rate_raw_temp_BW[sex_age][expID2])
            
            if feeder == 1:
                feed_cum_raw_temp[expID2]       = pd.Series(df_feed_cum_raw_temp[sex_age][expID2])
                feed_rate_raw_temp[expID2]      = pd.Series(df_feed_rate_raw_temp[sex_age][expID2])
                feed_rate_raw_min_temp[expID2]  = pd.Series(df_feed_rate_raw_min_temp[sex_age][expID2])

                feed_cum_raw_temp_BW[expID2]      = pd.Series(df_feed_cum_raw_temp_BW[sex_age][expID2])
#                feed_rate_raw_temp_BW[expID2]     = pd.Series(df_feed_rate_raw_temp[sex_age][expID2])
#                feed_rate_raw_min_temp_BW[expID2] = pd.Series(df_feed_rate_raw_temp[sex_age][expID2])
                food_cum_corrected_temp[expID2]    = pd.Series(df_food_cum_corrected_temp[sex_age][expID2])
                food_cum_corrected_temp_BW[expID2] = pd.Series(df_food_cum_corrected_temp_BW[sex_age][expID2]) 

                if lickometer == 1:
                    water_to_food_g_temp[expID2]    = pd.Series(df_water_to_food_g_temp[sex_age][expID2]) 
                    lick_cum_raw_temp_BW[expID2]    = pd.Series(df_lick_cum_raw_temp_BW[sex_age][expID2])

#                lick_rate_raw_temp_BW[expID2]   = pd.Series(df_lick_rate_raw_temp[sex_age][expID2])
#                lick_rate_raw_min_temp_BW[expID2]   = pd.Series(df_lick_rate_raw_temp[sex_age][expID2])
                    water_to_food_g_BW_temp[expID2]    = pd.Series(df_water_to_food_g_BW_temp[sex_age][expID2])

    if lickometer == 1:
    #    df_firstlick[sex_age][mouseID_sex_age]                = firstlick_temp.mean(axis=1)
        df_lick_cum_raw[sex_age][mouseID_sex_age]             = lick_cum_raw_temp.mean(axis=1) 
        df_lick_cum_firstlick[sex_age][mouseID_sex_age]       = lick_cum_firstlick_temp.mean(axis=1)          
        df_lick_cum_firstbout[sex_age][mouseID_sex_age]       = lick_cum_firstbout_temp.mean(axis=1)
        
        df_lick_rate_raw[sex_age][mouseID_sex_age]             = lick_rate_raw_temp.mean(axis=1) 
        df_lick_rate_firstlick[sex_age][mouseID_sex_age]       = lick_rate_firstlick_temp.mean(axis=1)          
        df_lick_rate_firstbout[sex_age][mouseID_sex_age]       = lick_rate_firstbout_temp.mean(axis=1)
     
        df_lick_rate_raw_min[sex_age][mouseID_sex_age]         = lick_rate_raw_min_temp.mean(axis=1) 
        df_lick_rate_firstlick_min[sex_age][mouseID_sex_age]   = lick_rate_firstlick_min_temp.mean(axis=1)          
        df_lick_rate_firstbout_min[sex_age][mouseID_sex_age]   = lick_rate_firstbout_min_temp.mean(axis=1)
        
    #    #% process mean from BW adjusted data     
        df_lick_cum_firstlick_BW[sex_age][mouseID_sex_age]      = lick_cum_firstlick_temp_BW.mean(axis=1)
    #    df_lick_rate_firstlick_BW[sex_age][mouseID_sex_age]     = lick_rate_firstlick_temp_BW.mean(axis=1)
    #    df_lick_rate_firstlick_min_BW[sex_age][mouseID_sex_age] = lick_rate_firstlick_min_temp_BW.mean(axis=1)
        if len(df_lick_cum_raw) > 1:
            num_replicates = num_replicates +1 

    #% process mean from feed data 
    if feeder == 1:
        df_feed_cum_raw[sex_age][mouseID_sex_age]         = feed_cum_raw_temp.mean(axis=1)
        df_feed_rate_raw[sex_age][mouseID_sex_age]        = feed_rate_raw_temp.mean(axis=1)
        df_feed_rate_raw_min[sex_age][mouseID_sex_age]    = feed_rate_raw_min_temp.mean(axis=1)
       
        df_feed_cum_raw_BW[sex_age][mouseID_sex_age]      = feed_cum_raw_temp_BW.mean(axis=1)
#        print ("df_feed_cum_raw_BW[%s][%s]" %(sex_age, mouseID_sex_age))
#        df_feed_rate_raw_BW[sex_age][mouseID_sex_age]     = feed_rate_raw_temp_BW.mean(axis=1)
#        df_feed_rate_raw_min_BW[sex_age][mouseID_sex_age] = feed_rate_raw_min_temp_BW.mean(axis=1)
        df_food_cum_corrected[sex_age][mouseID_sex_age]    = food_cum_corrected_temp.mean(axis=1)
        df_food_cum_corrected_BW[sex_age][mouseID_sex_age] = food_cum_corrected_temp_BW.mean(axis=1)      

        if len(df_feed_cum_raw) > 1:
            num_replicates = num_replicates +1 

        if lickometer == 1:
            df_lick_cum_raw_BW[sex_age][mouseID_sex_age]      = lick_cum_raw_temp_BW.mean(axis=1)

    #        print ("df_lick_cum_raw_BW[%s][%s]" %(sex_age, mouseID_sex_age))
    #        df_lick_rate_raw_BW[sex_age][mouseID_sex_age]     = lick_rate_raw_temp_BW.mean(axis=1)
    #        df_lick_rate_raw_min_BW[sex_age][mouseID_sex_age] = lick_rate_raw_min_temp_BW.mean(axis=1)
    
            ## caclulate licks per pellet consumption along the time series        
            try:
                df_lick_per_feed[sex_age][mouseID_sex_age]    = df_lick_cum_raw[sex_age][mouseID_sex_age]/df_feed_cum_raw[sex_age][mouseID_sex_age]
            except:
                df_lick_per_feed[sex_age][mouseID_sex_age]    = np.nan
            
            df_water_to_food_g[sex_age][mouseID_sex_age]       = water_to_food_g_temp.mean(axis=1)
            df_water_to_food_g_BW[sex_age][mouseID_sex_age]    = water_to_food_g_BW_temp.mean(axis=1)
        

    
#print ("mice included in this analysis are: ")
#print (mouseID_sex_age_list)    
print ("number of mice with technical replicates are: %d" %num_replicates)

#%%################################################################################
#   feeder all cumulative lick (mean value if same mouse repeated))
###################################################################################

if feeder == 1:
    for sex_age in male_list + female_list:
        try:    df_lick_cum_raw[sex_age]                = df_lick_cum_raw[sex_age].drop(columns = ['dummy'])    
        except: pass #print ("dummy column already removed from df_lick_cum_raw[%s]" %sex_age)        
        try:    df_feed_cum_raw[sex_age]                = df_feed_cum_raw[sex_age].drop(columns = ['dummy'])    
        except: pass #print ("dummy column already removed from df_feed_cum_raw[%s]" %sex_age)    
        try:    df_lick_cum_raw_temp[sex_age]           = df_lick_cum_raw_temp[sex_age].drop(columns = ['dummy'])
        except: pass
        try:    df_feed_cum_raw_temp[sex_age]           = df_feed_cum_raw_temp[sex_age].drop(columns = ['dummy'])
        except: pass
        try:    df_water_to_food_g[sex_age]            = df_water_to_food_g[sex_age].drop(columns = ['dummy'])
        except: pass

    
#        try:    df_lick_per_feed[sex_age]               = df_lick_per_feed[sex_age].drop(columns = ['dummy'])
#        except: pass #print("dummy column already removed from df_lick_per_feed[%s]" %sex_age)
#        try:    df_lick_rate_raw[sex_age]               = df_lick_cum_raw[sex_age].drop(columns = ['dummy'])    
#        except: pass #print ("dummy column already removed from df_lick_rate_raw[%s]" %sex_age)        
#        try:    df_feed_rate_raw[sex_age]               = df_feed_cum_raw[sex_age].drop(columns = ['dummy'])    
#        except: pass #print ("dummy column already removed from df_feed_rate_raw[%s]" %sex_age)    
#    
#        try:    df_lick_rate_raw_min[sex_age]           = df_lick_rate_raw_min[sex_age].drop(columns = ['dummy'])    
#        except: pass #print ("dummy column already removed from df_lick_rate_raw_min[%s]" %sex_age)        
#        try:    df_feed_rate_raw_min[sex_age]           = df_feed_rate_raw_min[sex_age].drop(columns = ['dummy'])    
#        except: pass #print ("dummy column already removed from df_feed_rate_raw_min[%s]" %sex_age)    
      
        try:    df_BW[sex_age]                          = df_BW[sex_age].drop(labels = ["dummy"])
        except: pass #print ("dummy column already removed from df_BW[%s]" %sex_age) 
        try:    df_lick_cum_raw_BW[sex_age]             = df_lick_cum_raw_BW[sex_age].drop(columns = ['dummy'])    
        except: pass #print ("dummy column already removed from df_lick_cum_raw[%s]" %sex_age)        
        try:    df_feed_cum_raw_BW[sex_age]             = df_feed_cum_raw_BW[sex_age].drop(columns = ['dummy'])    
        except: pass #print ("dummy column already removed from df_feed_cum_raw[%s]" %sex_age)    

        try:    df_water_to_food_g[sex_age]             = df_water_to_food_g[sex_age].drop(columns = ["dummy"])
        except: pass
        try:    df_food_cum_corrected[sex_age]          = df_food_cum_corrected[sex_age].drop(columns = ["dummy"])
        except: pass
        try:    df_food_cum_corrected_BW[sex_age]       = df_food_cum_corrected_BW[sex_age].drop(columns = ["dummy"])
        except: pass
        try:    df_water_to_food_g_BW[sex_age]          = df_water_to_food_g_BW[sex_age].drop(columns = ["dummy"])
        except: pass
    
#        try:    df_lick_rate_raw_BW[sex_age]            = df_lick_rate_raw_BW[sex_age].drop(columns = ['dummy'])    
#        except: pass #print ("dummy column already removed from df_lick_rate_raw[%s]" %sex_age)        
#        try:    df_feed_rate_raw_BW[sex_age]            = df_feed_rate_raw_BW[sex_age].drop(columns = ['dummy'])    
#        except: pass #print ("dummy column already removed from df_feed_rate_raw[%s]" %sex_age)    
#        try:    df_lick_rate_raw_min_BW[sex_age]            = df_lick_rate_raw_min_BW[sex_age].drop(columns = ['dummy'])    
#        except: pass #print ("dummy column already removed from df_lick_rate_raw[%s]" %sex_age)        
#        try:    df_feed_rate_raw_min_BW[sex_age]            = df_feed_rate_raw_min_BW[sex_age].drop(columns = ['dummy'])    
#        except: pass #print ("dummy column already removed from df_feed_rate_raw[%s]" %sex_age) 
   
else:
    for sex_age in male_list + female_list:
        try:    df_lick_cum_raw[sex_age]                = df_lick_cum_raw[sex_age].drop(columns = ['dummy'])    
        except: pass # print ("dummy column already removed from df_lick_cum_raw[%s]" %sex_age)        
        try:    df_lick_cum_firstlick[sex_age]          = df_lick_cum_firstlick[sex_age].drop(columns = ['dummy'])
        except: pass #print ("dummy column already removed from df_lick_cum_firstlick[%s]" %sex_age)    
        try:    df_lick_cum_firstbout[sex_age]          = df_lick_cum_firstbout[sex_age].drop(columns = ['dummy'])
        except: pass #print ("dummy column already removed from df_lick_cum_firstbout[%s]" %sex_age)     


#%%########################################################################################
#   FEEDER different types of cumulative lick plots 
#########################################################################################
#   1) cumulative lick (mean value if same mouse repeated))

# subplot(nrows, ncols, index, **kwargs)
# plt.subplot2grid(shape, location, rowspan, colspan)
if feeder == 1: 
    ##################################
    #   FEEDER PLOT
    ##################################    
    # set maximum value for y-axis 
    for k, sex_list in enumerate ([male_list, female_list]):  
        for i, sex_age in enumerate (sex_list):
            feed_cum_mean       = df_feed_cum_raw[sex_age].mean(axis=1)
            feed_cum_stderr     = ss.sem(df_feed_cum_raw[sex_age], axis=1)
            max_cumfeed         = max(max_cumfeed, float(feed_cum_mean[-1:] + feed_cum_stderr[-1:])) 
        
    for k, sex_list in enumerate ([male_list, female_list]):  
        plt.figure(figsize = (8, 6))
        for i, sex_age in enumerate (sex_list):
            feed_cum_mean       = df_feed_cum_raw[sex_age].mean(axis=1)
            feed_cum_stderr     = ss.sem(df_feed_cum_raw[sex_age], axis=1)
    #        # Print and save mean figures 
            legend              = sex_age + " " + "n=%d" %len(df_feed_cum_raw[sex_age].columns)
    #            markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5, label = legend)
            markers, caps, bars = plt.errorbar(time_sec/60, feed_cum_mean, feed_cum_stderr, fmt='-', linewidth = 4, ecolor = palette2[i], label = legend)
            markers.set_color(palette[i])
            [bar.set_alpha(0.05) for bar in bars]
            [cap.set_alpha(0.05) for cap in caps]   
        plt.legend(loc="upper left", bbox_to_anchor = (1,1), fontsize = 18)
        plt.tick_params(axis = "both", which = "major", labelsize = 18)
        plt.xlabel("time (min)", fontsize = 18)
        plt.ylabel("cum. # pellets", fontsize = 18)
        plt.ylim(0, np.max(max_cumfeed)*1.05)

        suptitle = "cumulative_feed_" + stim + "_" + sex_list[0][0] + "_" + str(length_cutoff) + "min"
        plt.title("%s" %suptitle, fontsize = 18)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
        plt.savefig(analysis_path + "/" + suptitle + ".png")
        plt.savefig(analysis_path + "/" + suptitle + ".pdf")
        # plt.show()
        plt.close()    
    
        
    ##################################
    #   FOOD INTAKE PLOT
    ##################################    
    # set maximum value for y-axis 
    for k, sex_list in enumerate ([male_list, female_list]):  
        for i, sex_age in enumerate (sex_list):
            feed_cum_mean       = df_feed_cum_raw[sex_age].mean(axis=1)*pellet_g
            feed_cum_stderr     = ss.sem(df_feed_cum_raw[sex_age]*pellet_g, axis=1)
            max_cumfood_g         = max(max_cumfood_g, float(feed_cum_mean[-1:] + feed_cum_stderr[-1:])) 
        
    for k, sex_list in enumerate ([male_list, female_list]):  
        plt.figure(figsize = (8, 6))
        for i, sex_age in enumerate (sex_list):
            feed_cum_mean       = df_feed_cum_raw[sex_age].mean(axis=1)*pellet_g
            feed_cum_stderr     = ss.sem(df_feed_cum_raw[sex_age]*pellet_g, axis=1)
    #        # Print and save mean figures 
            legend              = sex_age + " " + "n=%d" %len(df_feed_cum_raw[sex_age].columns)
    #            markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5, label = legend)
            markers, caps, bars = plt.errorbar(time_sec/60, feed_cum_mean, feed_cum_stderr, fmt='-', linewidth = 4, ecolor = palette2[i], label = legend)
            markers.set_color(palette[i])
            [bar.set_alpha(0.05) for bar in bars]
            [cap.set_alpha(0.05) for cap in caps]   
        plt.legend(loc="upper left", bbox_to_anchor = (1,1), fontsize = 18)
        plt.tick_params(axis = "both", which = "major", labelsize = 18)
        plt.xlabel("time (min)", fontsize = 18)
        plt.ylabel("food intake (g)", fontsize = 18)
        plt.ylim(0, np.max(max_cumfood_g)*1.05)

        suptitle = "cumulative_food_g_" + stim + "_" + sex_list[0][0] + "_" + str(length_cutoff) + "min"
        plt.title("%s" %suptitle, fontsize = 18)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
        plt.savefig(analysis_path + "/" + suptitle + ".png")
        plt.savefig(analysis_path + "/" + suptitle + ".pdf")
    #    plt.show()
        plt.close()    
        
    ###########################################################################################################################
    #   cumulative raw food intake/BW 
    ###########################################################################################################################
    
    # set maximum value for y-axis 
    for k, sex_list in enumerate ([male_list, female_list]):  
        for i, sex_age in enumerate (sex_list):
            feed_cum_BW_mean    = df_feed_cum_raw_BW[sex_age].mean(axis=1)
            feed_cum_BW_stderr  = ss.sem(df_feed_cum_raw_BW[sex_age], axis=1)
    #        lick_cum_BW_stderr  = df_lick_cum_firstlick_BW[sex_age].std(axis=1)/math.sqrt(len(df_lick_cum_firstlick_BW[sex_age].columns)-1)
            max_cumfood_g_BW       = max(max_cumfood_g_BW, float(feed_cum_BW_mean[-1:] + feed_cum_BW_stderr[-1:]) )
            
    nonnan_BW   = {}
    ind_nonnan  = {}
    for k, sex_list in enumerate ([male_list, female_list]):  
        plt.figure(figsize = (8, 6))
        for i, sex_age in enumerate (sex_list):
            feed_cum_mean       = df_feed_cum_raw_BW[sex_age].mean(axis=1)
            feed_cum_stderr     = df_feed_cum_raw_BW[sex_age].std(axis=1)
            nonnan_BW[sex_age]  = [x for x in df_BW[sex_age].values if x >0]
            try:
                feed_cum_stderr = feed_cum_stderr/math.sqrt(len(nonnan_BW[sex_age])-1)
            except:
                feed_cum_stderr = pd.Series([np.nan]*len(time_sec))        # Print and save mean figures 
            BW_columns          = list(df_feed_cum_raw_BW[sex_age].columns)
            BW_columns_range    = np.arange(0, len(BW_columns), 1)
            ind_nonnan[sex_age] = [x for x in BW_columns_range if df_feed_cum_raw_BW[sex_age][BW_columns[x]][0] >=0] # count the # of columns with nonnan values
            legend              = sex_age + " " + "n=%d" %(len(ind_nonnan[sex_age]))
    #            markers, caps, bars = plt.errorbar(time_sec/60, feed_cum_mean, feed_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5,  label = legend)
            markers, caps, bars = plt.errorbar(time_sec/60, feed_cum_mean, feed_cum_stderr, fmt='-', linewidth = 4, ecolor = palette2[i], label = legend)
            markers.set_color(palette[i])
            [bar.set_alpha(0.05) for bar in bars]
            [cap.set_alpha(0.05) for cap in caps]
        #ax9.set_xticklabels(m, fontsize = 16)
        plt.legend(loc="upper left", bbox_to_anchor = (1, 1), fontsize = 18)
        plt.tick_params(axis = "both", which = "major", labelsize = 18)
        plt.xlabel("time (min)", fontsize = 18)
        plt.ylabel("food intake(g)/BW(g)", fontsize = 18)
        plt.ylim(0, np.max(max_cumfood_g_BW)*1.1)

        suptitle = "cum_raw_food_g_BW_" + stim + "_" + sex_list[0][0] + "_" + str(length_cutoff) + "min"
        plt.title ("%s" %suptitle, fontsize = 18)
        plt.suptitle("body weight adjusted", fontsize = 18)     
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
        plt.savefig(analysis_path + "/" + suptitle + ".png")
        plt.savefig(analysis_path + "/" + suptitle + ".pdf")
        # plt.show()
        plt.close()

    #%%##########################################################################################################################
    #   cumulative corrected food intake 
    ###########################################################################################################################
    
    # set maximum value for y-axis 
    max_cumfood_g   =  np.array([0])
    for k, sex_list in enumerate ([male_list, female_list]):  
        for i, sex_age in enumerate (sex_list):
            feed_cum_mean    = df_food_cum_corrected[sex_age].mean(axis=1)
            feed_cum_stderr  = ss.sem(df_food_cum_corrected[sex_age], axis=1)
    #        lick_cum_BW_stderr  = df_lick_cum_firstlick_BW[sex_age].std(axis=1)/math.sqrt(len(df_lick_cum_firstlick_BW[sex_age].columns)-1)
            max_cumfood_g       = max(max_cumfood_g, float(feed_cum_mean[-1:] + feed_cum_stderr[-1:]) )
            
    nonnan_BW   = {}
    ind_nonnan  = {}
    for k, sex_list in enumerate ([male_list, female_list]):  
        plt.figure(figsize = (8, 6))
        for i, sex_age in enumerate (sex_list):
            feed_cum_mean       = df_food_cum_corrected[sex_age].mean(axis=1)
            feed_cum_stderr     = df_food_cum_corrected[sex_age].std(axis=1)
            nonnan_BW[sex_age]  = [x for x in df_BW[sex_age].values if x >0]
            try:
                feed_cum_stderr = feed_cum_stderr/math.sqrt(len(nonnan_BW[sex_age])-1)
            except:
                feed_cum_stderr = pd.Series([np.nan]*len(time_sec))        # Print and save mean figures 
            BW_columns          = list(df_food_cum_corrected[sex_age].columns)
            BW_columns_range    = np.arange(0, len(BW_columns), 1)
            ind_nonnan[sex_age] = [x for x in BW_columns_range if df_food_cum_corrected[sex_age][BW_columns[x]][0] >=0] # count the # of columns with nonnan values
            legend              = sex_age + " " + "n=%d" %(len(ind_nonnan[sex_age]))
    #            markers, caps, bars = plt.errorbar(time_sec/60, feed_cum_mean, feed_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5,  label = legend)
            markers, caps, bars = plt.errorbar(time_sec/60, feed_cum_mean, feed_cum_stderr, fmt='-', linewidth = 4, ecolor = palette2[i], label = legend)
            markers.set_color(palette[i])
            [bar.set_alpha(0.05) for bar in bars]
            [cap.set_alpha(0.05) for cap in caps]
        #ax9.set_xticklabels(m, fontsize = 16)
        plt.legend(loc="upper left", bbox_to_anchor = (1, 1), fontsize = 18)
        plt.tick_params(axis = "both", which = "major", labelsize = 18)
        plt.xlabel("time (min)", fontsize = 18)
        plt.ylabel("food intake(g)", fontsize = 18)
        plt.ylim(0, np.max(max_cumfood_g)*1.1)

        suptitle = "cum_corr_food_g_" + stim + "_" + sex_list[0][0] + "_" + str(length_cutoff) + "min"
        plt.title ("%s" %suptitle, fontsize = 18)
        plt.suptitle("food spill corrected", fontsize = 18)     
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
        plt.savefig(analysis_path + "/" + suptitle + ".png")
        plt.savefig(analysis_path + "/" + suptitle + ".pdf")
        # plt.show()
        plt.close()

    ###########################################################################################################################
    #   cumulative corrected food intake/BW
    ###########################################################################################################################
    
    # set maximum value for y-axis 
    max_cumfood_g_BW   =  np.array([0])
    for k, sex_list in enumerate ([male_list, female_list]):  
        for i, sex_age in enumerate (sex_list):
            feed_cum_mean    = df_food_cum_corrected_BW[sex_age].mean(axis=1)
            feed_cum_stderr  = ss.sem(df_food_cum_corrected_BW[sex_age], axis=1)
    #        lick_cum_BW_stderr  = df_lick_cum_firstlick_BW[sex_age].std(axis=1)/math.sqrt(len(df_lick_cum_firstlick_BW[sex_age].columns)-1)
            max_cumfood_g_BW      = max(max_cumfood_g_BW, float(feed_cum_mean[-1:] + feed_cum_stderr[-1:]) )
            
    nonnan_BW   = {}
    ind_nonnan  = {}
    for k, sex_list in enumerate ([male_list, female_list]):  
        plt.figure(figsize = (8, 6))
        for i, sex_age in enumerate (sex_list):
            feed_cum_mean       = df_food_cum_corrected_BW[sex_age].mean(axis=1)
            feed_cum_stderr     = df_food_cum_corrected_BW[sex_age].std(axis=1)
            nonnan_BW[sex_age]  = [x for x in df_BW[sex_age].values if x >0]
            try:
                feed_cum_stderr = feed_cum_stderr/math.sqrt(len(nonnan_BW[sex_age])-1)
            except:
                feed_cum_stderr = pd.Series([np.nan]*len(time_sec))        # Print and save mean figures 
            BW_columns          = list(df_food_cum_corrected_BW[sex_age].columns)
            BW_columns_range    = np.arange(0, len(BW_columns), 1)
            ind_nonnan[sex_age] = [x for x in BW_columns_range if df_food_cum_corrected_BW[sex_age][BW_columns[x]][0] >=0] # count the # of columns with nonnan values
            legend              = sex_age + " " + "n=%d" %(len(ind_nonnan[sex_age]))
    #            markers, caps, bars = plt.errorbar(time_sec/60, feed_cum_mean, feed_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5,  label = legend)
            markers, caps, bars = plt.errorbar(time_sec/60, feed_cum_mean, feed_cum_stderr, fmt='-', linewidth = 4, ecolor = palette2[i], label = legend)
            markers.set_color(palette[i])
            [bar.set_alpha(0.05) for bar in bars]
            [cap.set_alpha(0.05) for cap in caps]
        #ax9.set_xticklabels(m, fontsize = 16)
        plt.legend(loc="upper left", bbox_to_anchor = (1, 1), fontsize = 18)
        plt.tick_params(axis = "both", which = "major", labelsize = 18)
        plt.xlabel("time (min)", fontsize = 18)
        plt.ylabel("food intake(g)/BW(g)", fontsize = 18)
        plt.ylim(0, np.max(max_cumfood_g_BW)*1.1)


        suptitle = "cum_corr_food_g_BW_" + stim + "_" + sex_list[0][0] + "_" + str(length_cutoff) + "min"
        plt.title ("%s" %suptitle, fontsize = 18)
        plt.suptitle("food spill corrected, BW adjusted", fontsize = 18)     
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
        plt.savefig(analysis_path + "/" + suptitle + ".png")
        plt.savefig(analysis_path + "/" + suptitle + ".pdf")
        plt.show()
        # plt.close()
        
    if lickometer == 1:    
        
        ###########################################################################################################################
        #  Divide water intake by the total pellet consumption (g) and plot by different time lengths
        ###########################################################################################################################
        # set maximum value for y-axis 
        for k, sex_list in enumerate ([male_list, female_list]):  
            for i, sex_age in enumerate (sex_list):
                lick_cum_mean       = df_water_to_food_g[sex_age].mean(axis=1)
                lick_cum_stderr     = ss.sem(df_water_to_food_g[sex_age], axis=1)                
    #                if len(df_water_to_food_g[sex_age].columns) > 2:
    #                    lick_cum_stderr     = df_water_to_food_g[sex_age].std(axis=1)/math.sqrt(len(df_water_to_food_g[sex_age].columns)-2)  # subtracting 1 and also 1 for "dummy" column
    #                else:
    #                    continue
                max_cumwater_food_g      = max(max_cumwater_food_g, float(lick_cum_mean[-1:] + lick_cum_stderr[-1:]))
            
        for k, sex_list in enumerate ([male_list, female_list]):  
            plt.figure(figsize = (8, 6))
            for i, sex_age in enumerate (sex_list):
                lick_cum_mean       = df_water_to_food_g[sex_age].mean(axis=1)
                lick_cum_stderr     = ss.sem(df_water_to_food_g[sex_age], axis=1)
    #                lick_cum_stderr     = df_water_to_food_g[sex_age].std(axis=1)/math.sqrt(len(df_water_to_food_g[sex_age].columns)-2)  # subtracting 1 and also 1 for "dummy" column
                
        #        # Print and save mean figures 
                legend              = sex_age + " " + "n=%d" %(len(df_water_to_food_g[sex_age].columns))
        #            markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5, label = legend)
                markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='-', linewidth = 4, ecolor = palette2[i], label = legend)
                markers.set_color(palette[i])
                [bar.set_alpha(0.05) for bar in bars]
                [cap.set_alpha(0.05) for cap in caps]   
            plt.legend(loc="upper left", bbox_to_anchor = (1,1), fontsize = 18)
            plt.tick_params(axis = "both", which = "major", labelsize = 18)
            plt.xlabel("time (min)", fontsize = 18)
            plt.ylabel("water intake (g)/total food intake (g)", fontsize = 18)
            plt.ylim(0, np.max(max_cumwater_food_g)*1.05)

            suptitle = "water_to_food_g_" + stim + "_" + sex_list[0][0] + "_" + str(length_cutoff) + "min"
            plt.title("%s" %suptitle, fontsize = 18)
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
            plt.savefig(analysis_path + "/" + suptitle + ".png")
            plt.savefig(analysis_path + "/" + suptitle + ".pdf")
        #    plt.show()
            plt.close()    
    

        ###########################################################################################################################
        #   cumulative body surface area corrected water intake / BW corrected, food-spill corrected food intake
        ###########################################################################################################################
        
        # set maximum value for y-axis 
        max_cumfood_g_BW   =  np.array([0])
        for k, sex_list in enumerate ([male_list, female_list]):  
            for i, sex_age in enumerate (sex_list):
                feed_cum_mean    = df_water_to_food_g_BW[sex_age].mean(axis=1)
                feed_cum_stderr  = ss.sem(df_water_to_food_g_BW[sex_age], axis=1)
        #        lick_cum_BW_stderr  = df_lick_cum_firstlick_BW[sex_age].std(axis=1)/math.sqrt(len(df_lick_cum_firstlick_BW[sex_age].columns)-1)
                max_cumfood_g_BW      = max(max_cumfood_g_BW, float(feed_cum_mean[-1:] + feed_cum_stderr[-1:]) )
                
        nonnan_BW   = {}
        ind_nonnan  = {}
        for k, sex_list in enumerate ([male_list, female_list]):  
            plt.figure(figsize = (8, 6))
            for i, sex_age in enumerate (sex_list):
                feed_cum_mean       = df_water_to_food_g_BW[sex_age].mean(axis=1)
                feed_cum_stderr     = df_water_to_food_g_BW[sex_age].std(axis=1)
                nonnan_BW[sex_age]  = [x for x in df_BW[sex_age].values if x >0]
                try:
                    feed_cum_stderr = feed_cum_stderr/math.sqrt(len(nonnan_BW[sex_age])-1)
                except:
                    feed_cum_stderr = pd.Series([np.nan]*len(time_sec))        # Print and save mean figures 
                BW_columns          = list(df_water_to_food_g_BW[sex_age].columns)
                BW_columns_range    = np.arange(0, len(BW_columns), 1)
                ind_nonnan[sex_age] = [x for x in BW_columns_range if df_water_to_food_g_BW[sex_age][BW_columns[x]][0] >=0] # count the # of columns with nonnan values
                legend              = sex_age + " " + "n=%d" %(len(ind_nonnan[sex_age]))
        #            markers, caps, bars = plt.errorbar(time_sec/60, feed_cum_mean, feed_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5,  label = legend)
                markers, caps, bars = plt.errorbar(time_sec/60, feed_cum_mean, feed_cum_stderr, fmt='-', linewidth = 4, ecolor = palette2[i], label = legend)
                markers.set_color(palette[i])
                [bar.set_alpha(0.05) for bar in bars]
                [cap.set_alpha(0.05) for cap in caps]
            #ax9.set_xticklabels(m, fontsize = 16)
            plt.legend(loc="upper left", bbox_to_anchor = (1, 1), fontsize = 18)
            plt.tick_params(axis = "both", which = "major", labelsize = 18)
            plt.xlabel("time (min)", fontsize = 18)
            plt.ylabel("(water intake(g)/BW(g)^0.667)/ \n (food intake(g)/BW(g))", fontsize = 18)
       
            suptitle = "cum_corr_water_food_g_BW_" + stim + "_" + sex_list[0][0] + "_" + str(length_cutoff) + "min"
            plt.title ("%s" %suptitle, fontsize = 18)
            plt.suptitle("food spill corrected, BW adjusted", fontsize = 18)     
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
            plt.savefig(analysis_path + "/" + suptitle + ".png")
            plt.savefig(analysis_path + "/" + suptitle + ".pdf")
            plt.show()
            # plt.close()

if lickometer == 1:
     #################################
     #   LICK PLOT
     #################################
     # set maximum value for y-axis 
     for k, sex_list in enumerate ([male_list, female_list]):  
         for i, sex_age in enumerate (sex_list):
             lick_cum_mean       = df_lick_cum_raw[sex_age].mean(axis=1)
             lick_cum_stderr     = ss.sem(df_lick_cum_raw[sex_age], axis=1)
             max_cumlick         = max(max_cumlick, float(lick_cum_mean[-1:] + lick_cum_stderr[-1:])) 
         
     for k, sex_list in enumerate ([male_list, female_list]):  
         plt.figure(figsize = (8, 6))
         for i, sex_age in enumerate (sex_list):
             lick_cum_mean       = df_lick_cum_raw[sex_age].mean(axis=1)
             lick_cum_stderr     = ss.sem(df_lick_cum_raw[sex_age], axis=1)
     #        # Print and save mean figures 
             legend              = sex_age + " " + "n=%d" %len(df_lick_cum_raw[sex_age].columns)
     #            markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5, label = legend)
             markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='-', linewidth = 4, ecolor = palette2[i], label = legend)
             markers.set_color(palette[i])
             [bar.set_alpha(0.05) for bar in bars]
             [cap.set_alpha(0.05) for cap in caps]   
         plt.legend(loc="upper left", bbox_to_anchor = (1,1), fontsize = 18)
         plt.tick_params(axis = "both", which = "major", labelsize = 18)
         plt.xlabel("time (min)", fontsize = 18)
         plt.ylabel("cumulative # licks", fontsize = 18)
         plt.ylim(0, np.max(max_cumlick)*1.05)

         suptitle = "cumulative_lick_" + stim + "_" + sex_list[0][0] + "_" + str(length_cutoff) + "min"
         plt.title("%s" %suptitle, fontsize = 18)
         plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
         plt.savefig(analysis_path + "/" + suptitle + ".png")
         plt.savefig(analysis_path + "/" + suptitle + ".pdf")
     #    plt.show()
         plt.close()        
         
     #################################
     #   WATER INTAKE PLOT
     #################################
     # set maximum value for y-axis 
     for k, sex_list in enumerate ([male_list, female_list]):  
         for i, sex_age in enumerate (sex_list):
             lick_cum_mean       = df_lick_cum_raw[sex_age].mean(axis=1)*lick_g
             lick_cum_stderr     = ss.sem(df_lick_cum_raw[sex_age]*lick_g, axis=1)
             max_cumwater_g         = max(max_cumwater_g, float(lick_cum_mean[-1:] + lick_cum_stderr[-1:])) 
         
     for k, sex_list in enumerate ([male_list, female_list]):  
         plt.figure(figsize = (8, 6))
         for i, sex_age in enumerate (sex_list):
             lick_cum_mean       = df_lick_cum_raw[sex_age].mean(axis=1)*lick_g
             lick_cum_stderr     = ss.sem(df_lick_cum_raw[sex_age]*lick_g, axis=1)
     #        # Print and save mean figures 
             legend              = sex_age + " " + "n=%d" %len(df_lick_cum_raw[sex_age].columns)
     #            markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5, label = legend)
             markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='-', linewidth = 4, ecolor = palette2[i], label = legend)
             markers.set_color(palette[i])
             [bar.set_alpha(0.05) for bar in bars]
             [cap.set_alpha(0.05) for cap in caps]   
         plt.legend(loc="upper left", bbox_to_anchor = (1,1), fontsize = 18)
         plt.tick_params(axis = "both", which = "major", labelsize = 18)
         plt.xlabel("time (min)", fontsize = 18)
         plt.ylabel("water intake (g)", fontsize = 18)
         plt.ylim(0, np.max(max_cumwater_g)*1.05)

         suptitle = "cumulative_water_g_" + stim + "_" + sex_list[0][0] + "_" + str(length_cutoff) + "min"
         plt.title("%s" %suptitle, fontsize = 18)
         plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
         plt.savefig(analysis_path + "/" + suptitle + ".png")
         plt.savefig(analysis_path + "/" + suptitle + ".pdf")
     #    plt.show()
         plt.close()    



     # ###########################################################################################################################
     # #   cumulative water intake/BW^0.667 
     # ###########################################################################################################################
     
     # # set maximum value for y-axis 
     # for k, sex_list in enumerate ([male_list, female_list]):  
     #     for i, sex_age in enumerate (sex_list):
     #         lick_cum_BW_mean    = df_lick_cum_raw_BW[sex_age].mean(axis=1)
     #         lick_cum_BW_stderr  = ss.sem(df_lick_cum_raw_BW[sex_age], axis=1)
     # #        lick_cum_BW_stderr  = df_lick_cum_firstlick_BW[sex_age].std(axis=1)/math.sqrt(len(df_lick_cum_firstlick_BW[sex_age].columns)-1)
     #         max_cumlick_BW      = max(max_cumlick_BW, float(lick_cum_BW_mean[-1:] + lick_cum_BW_stderr[-1:]) )
             
     # nonnan_BW   = {}
     # ind_nonnan  = {}
     # for k, sex_list in enumerate ([male_list, female_list]):  
     #     plt.figure(figsize = (8, 6))
     #     for i, sex_age in enumerate (sex_list):
     #         lick_cum_BW_mean       = df_lick_cum_raw_BW[sex_age].mean(axis=1)
     #         lick_cum_BW_stderr     = df_lick_cum_raw_BW[sex_age].std(axis=1)
     #         nonnan_BW[sex_age]  = [x for x in df_BW[sex_age].values if x >0]
     #         try:
     #             lick_cum_BW_stderr = lick_cum_BW_stderr/math.sqrt(len(nonnan_BW[sex_age])-1)
     #         except:
     #             lick_cum_BW_stderr = pd.Series([np.nan]*len(time_sec))        # Print and save mean figures 
     #         BW_columns          = list(df_lick_cum_raw_BW[sex_age].columns)
     #         BW_columns_range    = np.arange(0, len(BW_columns), 1)
     #         ind_nonnan[sex_age] = [x for x in BW_columns_range if df_lick_cum_raw_BW[sex_age][BW_columns[x]][0] >=0] # count the # of columns with nonnan values
     #         legend              = sex_age + " " + "n=%d" %(len(ind_nonnan[sex_age]))
     # #            markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5,  label = legend)
     #         markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_BW_mean, lick_cum_BW_stderr, fmt='-', linewidth = 4, ecolor = palette2[i], label = legend)
     #         markers.set_color(palette[i])
     #         [bar.set_alpha(0.05) for bar in bars]
     #         [cap.set_alpha(0.05) for cap in caps]
     #     #ax9.set_xticklabels(m, fontsize = 16)sim
     #     plt.legend(loc="upper left", bbox_to_anchor = (1, 1), fontsize = 18)
     #     plt.tick_params(axis = "both", which = "major", labelsize = 18)
     #     plt.xlabel("time (min)", fontsize = 18)
     #     plt.ylabel("water intake(g)/BW(g)^0.667", fontsize = 18)
     #     # plt.ylim(0, np.max(max_cumlick_BW)*1.1)

     #     if opto == 1:       
     #         plt.axvspan(laser_on_mins[0], laser_on_mins[1], color = "red", alpha=0.1)
     # #    if stim == "2MNaCl" or length_cutoff >= 80:
     # #        ax.set_ylim(0, 60)
     # #        
     # #    # elif length_cutoff >= 80:
     # #        # ax.set_ylim(0, 60)
     # #
     # #    else:
     # #        ax.set_ylim(0, 30)
     # #    plt.ylim(0, 0.07)
     #     suptitle = "cum_water_BW_" + stim + "_" + sex_list[0][0] + "_" + str(length_cutoff) + "min"
     #     plt.title ("%s" %suptitle, fontsize = 18)
     #     plt.suptitle("body surface area adjusted", fontsize = 18)     
     #     plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
     #     plt.savefig(analysis_path + "/" + suptitle + ".png")
     #     plt.savefig(analysis_path + "/" + suptitle + ".pdf")
     #     plt.show()
     #     # plt.close()


##########################################################################################
#   endpoint values for water and food intake, save into dictionary and export to csv file 
##########################################################################################          
         
if feeder == 1:
    dic_cumfood         = {}
    dic_cumfood_BW      = {}
    dic_cumwater        = {}
    dic_cumwater_BSA    = {}
    dic_cumwater_BSA_to_cumfood_BW = {}
    
    if lickometer == 1:
        dic_cumwater_BSA_to_cumfood_BW  = {}
    
    for k, sex_age in enumerate (male_list + female_list): 
        for j, mouseID_sex_age in enumerate(df_feed_cum_raw[sex_age].columns):      
            dic_cumfood[mouseID_sex_age]                    = df_food_cum_corrected[sex_age][mouseID_sex_age][length_cutoff*60-1]   
            dic_cumfood_BW[mouseID_sex_age]                 = df_food_cum_corrected_BW[sex_age][mouseID_sex_age][length_cutoff*60-1]               

            if lickometer == 1:      
                dic_cumwater[mouseID_sex_age]                      = df_lick_cum_raw[sex_age][mouseID_sex_age][length_cutoff*60-1]*lick_g
                dic_cumwater_BSA[mouseID_sex_age]                  = df_lick_cum_raw_BW[sex_age][mouseID_sex_age][length_cutoff*60-1]
                dic_cumwater_BSA_to_cumfood_BW[mouseID_sex_age]    = df_water_to_food_g_BW[sex_age][mouseID_sex_age][length_cutoff*60-1]   
    

    csv_path    = analysis_path + "/cumfood_" + stim + "_" + str(length_cutoff) + "min" + ".csv"      
    with open(csv_path, 'w') as f:
        for key in dic_cumfood.keys():
            f.write("%s, %s\n" % (key, dic_cumfood[key]))
                        
    csv_path    = analysis_path + "/cumfood_BW_" + stim + "_" + str(length_cutoff) + "min" + ".csv"      
    with open(csv_path, 'w') as f:
        for key in dic_cumfood_BW.keys():
            f.write("%s, %s\n" % (key, dic_cumfood_BW[key]))

    if lickometer == 1:
        csv_path    = analysis_path + "/cumwater_" + stim + "_" + str(length_cutoff) + "min" + ".csv"      
        with open(csv_path, 'w') as f:
            for key in dic_cumwater.keys():
                f.write("%s, %s\n" % (key, dic_cumwater[key]))   
              
        csv_path    = analysis_path + "/cumwater_BW_" + stim + "_" + str(length_cutoff) + "min" + ".csv"      
        with open(csv_path, 'w') as f:
            for key in dic_cumwater_BSA.keys():
                f.write("%s, %s\n" % (key, dic_cumwater_BSA[key]))  

        csv_path    = analysis_path + "/cumwater_BSA_to_cumfood_BW_" + stim + "_" + str(length_cutoff) + "min" + ".csv"      
        with open(csv_path, 'w') as f:
            for key in dic_cumwater_BSA_to_cumfood_BW.keys():
                f.write("%s, %s\n" % (key, dic_cumwater_BSA_to_cumfood_BW[key]))   


elif lickometer == 1:
    dic_cumwater        = {}
    dic_cumwater_BSA    = {}    
    for k, sex_age in enumerate (male_list + female_list): 
        for j, mouseID_sex_age in enumerate(df_lick_cum_raw[sex_age].columns): 
            dic_cumwater[mouseID_sex_age]           = df_lick_cum_firstlick[sex_age][mouseID_sex_age][length_cutoff*60-1]*lick_g
            dic_cumwater_BSA[mouseID_sex_age]       = df_lick_cum_firstlick_BW[sex_age][mouseID_sex_age][length_cutoff*60-1]
            
    csv_path    = analysis_path + "/cumwater_" + stim + "_" + str(length_cutoff) + "min" + ".csv"      
    with open(csv_path, 'w') as f:
        for key in dic_cumwater.keys():
            f.write("%s, %s\n" % (key, dic_cumwater[key]))   
        
    csv_path    = analysis_path + "/cumwater_BW_" + stim + "_" + str(length_cutoff) + "min" + ".csv"      
    with open(csv_path, 'w') as f:
        for key in dic_cumwater_BSA.keys():
            f.write("%s, %s\n" % (key, dic_cumwater_BSA[key]))       
            
                 
 


# #%%################################################################
# #   Area under the curve and peak values
# ###################################################################
     
# AUC         = {}
# AUC_        = {}
# AUC_BW      = {}

       
# if feeder == 1:
#     AUC         = {}
#     AUC_BW      = {}
#     for k, sex_age in enumerate (male_list + female_list):    
#         for j, mouseID_sex_age in enumerate(df_food_cum_corrected[sex_age].columns):
#             AUC[mouseID_sex_age]    = simps(df_food_cum_corrected[sex_age][mouseID_sex_age][feed_start_AUC*60:feed_end_AUC*60])
#             AUC_BW[mouseID_sex_age] = simps(df_food_cum_corrected_BW[sex_age][mouseID_sex_age][feed_start_AUC*60:feed_end_AUC*60])

#     csv_path    = analysis_path + "/AUC_cumfood_" + stim + "_" + str(length_cutoff) + "min" + ".csv"       
#     with open(csv_path, 'w') as f:
#         for key in AUC.keys():
#             f.write("%s, %s\n" % (key, AUC[key]))    
    
#     csv_path    = analysis_path + "/AUC_cumfoodBW_" + stim + "_" + str(length_cutoff) + "min" + ".csv"       
#     with open(csv_path, 'w') as f:
#         for key in AUC.keys():
#             f.write("%s, %s\n" % (key, AUC_BW[key]))    
            
#     if lickometer == 1:
#         AUC         = {}
#         AUC_BW      = {}
#         for k, sex_age in enumerate (male_list + female_list):    
#             for j, mouseID_sex_age in enumerate(df_lick_cum_raw[sex_age].columns):
#                 df_temp = df_lick_cum_raw[sex_age][mouseID_sex_age]*lick_g
#                 AUC[mouseID_sex_age]    = simps(df_temp[lick_start_AUC*60:lick_end_AUC*60])
#                 AUC_BW[mouseID_sex_age] = simps(df_lick_cum_raw_BW[sex_age][mouseID_sex_age][lick_start_AUC*60:lick_end_AUC*60])
        
#         csv_path    = analysis_path + "/AUC_raw_cumwater_" + stim + "_" + str(length_cutoff) + "min" + ".csv"       
#         with open(csv_path, 'w') as f:
#             for key in AUC.keys():
#                 f.write("%s, %s\n" % (key, AUC[key]))    
        
#         csv_path    = analysis_path + "/AUC_raw_cumwaterBSA_" + stim + "_" + str(length_cutoff) + "min" + ".csv"       
#         with open(csv_path, 'w') as f:
#             for key in AUC_BW.keys():
#                 f.write("%s, %s\n" % (key, AUC_BW[key]))      
            
# elif lickometer == 1:
#     AUC      = {}
#     AUC_     = {}
#     AUC_BW   = {}
#     for k, sex_age in enumerate (male_list + female_list):    
#         for j, mouseID_sex_age in enumerate(df_lick_cum_firstlick[sex_age].columns):
#             AUC[mouseID_sex_age]     = simps(df_lick_cum_firstlick[sex_age][mouseID_sex_age][lick_start_AUC*60:lick_end_AUC*60])
#             temp                     = df_lick_cum_firstlick[sex_age][mouseID_sex_age]*lick_g
#             AUC_[mouseID_sex_age]    = simps(temp[lick_start_AUC*60:lick_end_AUC*60])
#             AUC_BW[mouseID_sex_age]  = simps(df_lick_cum_firstlick_BW[sex_age][mouseID_sex_age][lick_start_AUC*60:lick_end_AUC*60])
    
#     csv_path    = analysis_path + "/AUC_firstlick_cumlick_" + stim + "_" + str(length_cutoff) + "min" + ".csv"       
#     with open(csv_path, 'w') as f:
#         for key in AUC.keys():
#             f.write("%s, %s\n" % (key, AUC[key]))    

#     csv_path    = analysis_path + "/AUC_firstlick_cumwater_" + stim + "_" + str(length_cutoff) + "min" + ".csv"       
#     with open(csv_path, 'w') as f:
#         for key in AUC_.keys():
#             f.write("%s, %s\n" % (key, AUC_[key]))    
    
#     csv_path    = analysis_path + "/AUC_firstlick_cumwaterBSA_" + stim + "_" + str(length_cutoff) + "min" + ".csv"       
#     with open(csv_path, 'w') as f:
#         for key in AUC_BW.keys():
#             f.write("%s, %s\n" % (key, AUC_BW[key]))

# sys.exit()

        
   

#%%##########################################################
#   BOUT STRUCTURE ANALYSIS
#############################################################

# data processing for bout structure histograms
df_bout_on      = {}
df_bout_off     = {}
df_bout_size    = {}
df_bout_length  = {}
df_bout_numbers = {}
df_interbout_interval = {}
df_bout_numbers = {}

for sex_age in male_list + female_list:
    df_bout_on[sex_age]            = {}
    df_bout_off[sex_age]           = {}
    df_bout_size[sex_age]          = {}
    df_bout_length[sex_age]        = {}
    df_interbout_interval[sex_age] = {}
    df_bout_numbers[sex_age]       = {}

for i, expID in enumerate (file_list):
    adjusted_age = adjust_age(data[expID]["age"])
    sex_age = data[expID]["sex"] + str(adjusted_age)
    bout_on, bout_off, bout_size, bout_length = bout_analysis(data[expID]["lickontimes_firstlick"], 1, 10)

    df_bout_on[sex_age][expID]      = bout_on
    df_bout_off[sex_age][expID]     = bout_off
    df_bout_length[sex_age][expID]  = bout_length
    df_bout_size[sex_age][expID]    = bout_size

    interbout_interval              = np.diff(bout_on) 
    df_interbout_interval[sex_age][expID]   = interbout_interval
    df_bout_numbers[sex_age][expID]         = [len(df_bout_on[sex_age][expID])]

df_bout_length_merged        = {}
df_bout_size_merged          = {}
df_interbout_interval_merged = {}
df_bout_numbers_merged       = {}

for sex, sex_list in zip (sexes, [male_list, female_list]):
    for sex_age in sex_list:
        df_bout_length_merged[sex_age]    = list(itertools.chain(*df_bout_length[sex_age].values()))
        df_bout_size_merged[sex_age]      = list(itertools.chain(*df_bout_size[sex_age].values())) 
        df_interbout_interval_merged[sex_age]      = list(itertools.chain(*df_interbout_interval[sex_age].values()))
        df_bout_numbers_merged[sex_age]            = list(itertools.chain(*df_bout_numbers[sex_age].values()))
    

#%%###############################################
#   POP open the analysis folder
##################################################
subprocess.Popen(f'explorer {os.path.realpath(analysis_path)}')    
    
