# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:23:30 2020

Group-level analysis of optogenetically modulated drinking behavior.

Loads processed mouse licking data (.pkl), groups recordings by
age and sex, and averages repeated recordings from individual mice.

Calculates cumulative licks, lick rates, estimated water intake,
and body-size-adjusted drinking measurements during optogenetic
stimulation experiments.

Generates group comparison plots with laser stimulation periods,
licking raster plots, and exports summary measurements as CSV,
PNG, and PDF files.

"""
#%%

import os
import glob
import pickle as pk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss 
import subprocess





print("Enter the full path for the directory with individual .pkl files: ")
path             = input()

stim       = "hyd"

total_number_laser = 1
laser_duration     = 15  # min
laser_freq         = 20  # Hz


protocol        = "15h15o15f"    
length_min      = 45 # minutes 
laser_on_times  = [15]
laser_off_times = [30]
water_on_time   = 15

BWpath = r"D:\Data_Dropbox_Sync\Mouse_metadata\allmiceBW.csv"
BWdata = pd.read_csv(BWpath)    
    
analysis_path = path + "/analysis/" + stim + "_" + str(laser_freq) + "Hz_" + protocol
if not os.path.exists(analysis_path):
    os.makedirs(analysis_path)   
#
palette        = ["dimgray", "maroon"]
palette2       = ["silver", "pink"]


male_list   = ["M3", "M26"]
female_list = ["F3", "F26"]  
sexes       = ["M", "F"]



lick_g      = 1.21/1000  # g or mL 
age_cutoff = 30

###########################################################################
###     Define functions      
###########################################################################

def adjust_age (age):
    if 25 <= age <= 29: # 25, 26, 27 months old
        adjusted_age = 26                         
                         
    elif 2  <= age <= 6: # 2, 3, 4 months old
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
tempdata_list         = []
data                  = {}

###########################################################################
###     Process data     
###########################################################################

# Go through each file in the directory to get .pkl data
for i, file_path in enumerate (all_files):
    expID             = file_path[(len(path)+1):-4] ## previous called as "file_name_cut
    stimulus          = expID.split("_")[4]
    # protocol_         = expID.split("_")[6].replace("n", "o")
    protocol_         = expID.split("_")[5].replace("n", "o")
    length            = int((expID.split("_")[-1]).split("min")[0])
    try:
        laser_frequency   = str((expID.split("_")[6]).split("Hz")[0])
    except:
        continue
    age               = int((expID.split("_")[2]).split("mo")[0])  
    # if ((stimulus  == stim) and (laser_frequency == str(laser_freq)) and (length >= length_min)): 
        # ((laser_freq == 0) or (protocol_ == protocol)) and (length >= length_min)): 
    if (age-1 <= age_cutoff): 
        adjusted_age = adjust_age(age)
        if (adjusted_age != 0):
            file_list.append(expID)  
            inputfile               = open(file_path, 'rb')
            data[expID]             = pk.load(inputfile)
            data[expID]["protocol"] = protocol_
            inputfile.close()
           

time_sec    = np.arange(1, length_min*60+1)  
time_min    = np.arange(1, length_min*60+1, 60)

df_firstlick_temp            = {}
df_firstbout_temp            = {}
df_lick_cum_raw_temp         = {}
df_lick_rate_raw_temp        = {}
df_lick_rate_raw_min_temp    = {}

df_lick_on_raw               = {}

df_firstlick                 = {}
df_firstbout                 = {}
df_lick_cum_raw              = {}
df_lick_rate_raw             = {}
df_lick_rate_raw_min         = {}

max_firstlick                = np.array([0])
max_firstbout                = np.array([0])
max_cumlick                  = np.array([0])
max_cumwater_g               = np.array([0])
max_cumlick_BW               = np.array([0])
    
df_BW                       = {}
df_lick_cum_raw_temp_BW     = {}
df_lick_cum_raw_BW          = {}
df_lick_rate_raw_temp_BW    = {}
df_lick_rate_raw_BW         = {}
df_lick_rate_raw_min_BW     = {}

          
# Initialize dataseries and dataframes
for sex_list in [male_list, female_list]: 
    for i, sex_age in enumerate (sex_list):  
        df_firstlick[sex_age]           = pd.Series(0, index = ["dummy"])
        df_firstbout[sex_age]           = pd.Series(0, index = ["dummy"])

        df_lick_cum_raw_temp[sex_age]    = pd.DataFrame(np.nan, index = np.arange(length_min*60), columns = ['dummy'])
        df_lick_cum_raw[sex_age]         = pd.DataFrame(np.nan, index = np.arange(length_min*60), columns = ['dummy'])

        df_lick_rate_raw_temp[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_min*60), columns = ['dummy'])
        df_lick_rate_raw[sex_age]          = pd.DataFrame(np.nan, index = np.arange(length_min*60), columns = ['dummy'])
       
        df_lick_rate_raw_min_temp[sex_age]      = pd.DataFrame(np.nan, index = np.arange(length_min), columns = ['dummy'])        
        df_lick_rate_raw_min[sex_age]           = pd.DataFrame(np.nan, index = np.arange(length_min), columns = ['dummy'])
       
        df_lick_on_raw[sex_age]          = pd.DataFrame(np.nan, index = np.arange(length_min*60), columns = ['dummy'])
        
        df_BW[sex_age]                   = pd.Series(0, index = ["dummy"])
                
        df_lick_cum_raw_temp_BW[sex_age] = pd.DataFrame(np.nan, index = np.arange(length_min*60), columns = ['dummy'])
        df_lick_cum_raw_BW[sex_age]      = pd.DataFrame(np.nan, index = np.arange(length_min*60), columns = ['dummy'])
        
        df_lick_rate_raw_temp_BW[sex_age]= pd.DataFrame(np.nan, index = np.arange(length_min*60), columns = ['dummy'])
        df_lick_rate_raw_BW[sex_age]     = pd.DataFrame(np.nan, index = np.arange(length_min*60), columns = ['dummy'])
        df_lick_rate_raw_min_BW[sex_age] = pd.DataFrame(np.nan, index = np.arange(length_min), columns = ['dummy'])

               
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
                   
    for sex in sexes:
        if data[expID]["sex"]  == sex:  
            # get cumulative licks  
            sex_age                                     = sex + str(adjusted_age) 
            firstlick_temp                              = pd.Series(data[expID]["first_lick_sec"], index = [expID])
            df_firstlick[sex_age]                       = df_firstlick[sex_age].append(firstlick_temp) 
            firstbout_temp                              = pd.Series(data[expID]["first_bout_sec"], index = [expID])
            df_firstbout[sex_age]                       = df_firstbout[sex_age].append(firstbout_temp)    
 
            time_sec_from_pkl                              = np.arange(1, data[expID]["totalmins"]*60+1)   
            lickrate_sec_raw, lickrate_cum_raw             = lick_rate_by_sec(data[expID]["lickontimes_hab_corr"], time_sec_from_pkl)   
           
            time_min_from_pkl                              = np.arange(1, (data[expID]["totalmins"]+1)*60, 60)   
            lickrate_min_raw, dump                         = lick_rate_by_sec(data[expID]["lickontimes_hab_corr"], time_min) 
            
            df_lick_cum_raw_temp[sex_age][expID]           = pd.DataFrame(lickrate_cum_raw+1)       
            df_lick_rate_raw_temp[sex_age][expID]          = pd.DataFrame(lickrate_sec_raw+1)           
            df_lick_rate_raw_min_temp[sex_age][expID]      = pd.DataFrame(lickrate_min_raw)           
            df_lick_on_raw[sex_age][expID]                 = pd.DataFrame(data[expID]["lickontimes_hab_corr"])   
              
            BW_temp                                        = pd.Series(data[expID]["BW"], index = [expID])
            df_BW[sex_age]                                 = df_BW[sex_age].append(BW_temp)
            
            df_lick_cum_raw_temp_BW[sex_age][expID]  = df_lick_cum_raw_temp[sex_age][expID]*lick_g/(df_BW[sex_age][expID]**0.667) 
            df_lick_rate_raw_temp_BW[sex_age][expID] = df_lick_rate_raw_temp[sex_age][expID]*lick_g/(df_BW[sex_age][expID]**0.667)
               
        
    max_firstlick = max(max_firstlick, max(df_firstlick[sex_age]))
    max_firstbout = max(max_firstbout, max(df_firstbout[sex_age])) 
    


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

#    firstlick_temp                  = pd.DataFrame(df_firstlick_temp[sex_age][expID], index = [expID])    
    lick_cum_raw_temp           = pd.DataFrame(df_lick_cum_raw_temp[sex_age][expID], columns = [expID])      
    lick_rate_raw_temp          = pd.DataFrame(df_lick_rate_raw_temp[sex_age][expID], columns = [expID])      
 
    lick_rate_raw_min_temp      = pd.DataFrame(df_lick_rate_raw_min_temp[sex_age][expID], columns = [expID])

    lick_cum_raw_temp_BW        = pd.DataFrame(df_lick_cum_raw_temp_BW[sex_age][expID], columns = [expID])
    lick_rate_raw_temp_BW       = pd.DataFrame(df_lick_rate_raw_temp_BW[sex_age][expID], columns = [expID])
    lick_rate_raw_min_temp_BW   = pd.DataFrame(df_lick_rate_raw_temp_BW[sex_age][expID], columns = [expID])
          
    for j, expID2 in enumerate (file_list):
        mouseID_sex_age2   = data[expID2]["mouseID"] + "_" + data[expID2]["sex"] + str(adjust_age(data[expID2]["age"]))
        if ((mouseID_sex_age == mouseID_sex_age2) and (expID != expID2) and ((laser_freq == 0) or (data[expID]["protocol"] == data[expID2]["protocol"]))):
            num_replicates = num_replicates +1 
            
#            firstlick_temp[mouseID_sex_age2]              = pd.Series(df_firstlick_temp[sex_age][expID2]) 
            lick_cum_raw_temp[mouseID_sex_age2]     = pd.Series(df_lick_cum_raw_temp[sex_age][expID2])
            lick_cum_raw_temp[mouseID_sex_age2]     = pd.Series(df_lick_cum_raw_temp[sex_age][expID2]) 
            lick_cum_raw_temp[mouseID_sex_age2]     = pd.Series(df_lick_cum_raw_temp[sex_age][expID2])

            lick_rate_raw_temp[mouseID_sex_age2]       = pd.Series(df_lick_rate_raw_temp[sex_age][expID2])
            lick_rate_raw_min_temp[mouseID_sex_age2]   = pd.Series(df_lick_rate_raw_min_temp[sex_age][expID2])

            print ("%s has technical replicates" %mouseID_sex_age) 
        
            lick_cum_raw_temp_BW[mouseID_sex_age2]        = pd.Series(df_lick_cum_raw_temp_BW[sex_age][expID2])
            lick_rate_raw_temp_BW[mouseID_sex_age2]       = pd.Series(df_lick_rate_raw_temp_BW[sex_age][expID2])
    

         
#    df_firstlick[sex_age][mouseID_sex_age]                = firstlick_temp.mean(axis=1)
    df_lick_cum_raw[sex_age][mouseID_sex_age]             = lick_cum_raw_temp.mean(axis=1)     
    df_lick_rate_raw[sex_age][mouseID_sex_age]            = lick_rate_raw_temp.mean(axis=1)  
    df_lick_rate_raw_min[sex_age][mouseID_sex_age]        = lick_rate_raw_min_temp.mean(axis=1) 
    
#    #% process mean from BW adjusted data     
    df_lick_cum_raw_BW[sex_age][mouseID_sex_age]      = lick_cum_raw_temp_BW.mean(axis=1)

#print ("mice included in this analysis are: ")
#print (mouseID_sex_age_list)    
print ("number of mice with technical replicates are: %d" %num_replicates)


#%%########################################################################################
#   Cumulative lick (mean value if same mouse repeated))
#########################################################################################

# subplot(nrows, ncols, index, **kwargs)
# plt.subplot2grid(shape, location, rowspan, colspan)

for sex_age in male_list + female_list:
    try:    df_lick_cum_raw[sex_age]          = df_lick_cum_raw[sex_age].drop(columns = ['dummy'])
    except: pass #print ("dummy column already removed from df_lick_cum_raw[%s]" %sex_age)    
    try:    df_lick_cum_raw_BW[sex_age]       = df_lick_cum_raw_BW[sex_age].drop(columns = ['dummy'])
    except: pass #print ("dummy column already removed from df_lick_cum_firstbout[%s]" %sex_age)  
 
# set maximum value for y-axis 
for k, sex_list in enumerate ([male_list, female_list]):  
    for i, sex_age in enumerate (sex_list):
        lick_cum_mean       = df_lick_cum_raw[sex_age].mean(axis=1)
        lick_cum_stderr     = ss.sem(df_lick_cum_raw[sex_age], axis=1)
        max_cumlick         = max(max_cumlick, float(lick_cum_mean[-1:] + lick_cum_stderr[-1:])) 
        
for k, sex_list in enumerate ([male_list, female_list]):  
    plt.figure(figsize = (10, 8))
    for i, sex_age in enumerate (sex_list):
        lick_cum_mean       = df_lick_cum_raw[sex_age].mean(axis=1)
        lick_cum_stderr     = ss.sem(df_lick_cum_raw[sex_age], axis=1)
#        # Print and save mean figures 
        legend              = sex_age + " " + "n=%d" %len(df_lick_cum_raw[sex_age].columns)
#            markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5, label = legend)
        markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='-', linewidth = 2, ecolor = palette2[i], label = legend)
        markers.set_color(palette[i])
        [bar.set_alpha(0.02) for bar in bars]
        [cap.set_alpha(0.02) for cap in caps]   
    plt.legend(loc="upper left", bbox_to_anchor = (1,1), fontsize = 18)
    plt.tick_params(axis = "both", which = "major", labelsize = 18)
    plt.xlabel("time (min)", fontsize = 18)
    plt.ylabel("cumulative # licks", fontsize = 18)
#     plt.ylim(0, np.max(max_cumlick)*1.05)
    plt.ylim(0, 2000)
    for j in np.arange(total_number_laser):
        plt.axvspan(laser_on_times[j], laser_on_times[j]+laser_duration, color = "red", alpha=0.2)
    plt.axvline (x = water_on_time, linewidth = 1, linestyle = "--", color = "gray")
    plt.xlim(0, 50)
#    ax.set_ylim(0, 800)
#    plt.title("raw", fontsize = 18) 
    suptitle = "cum_licks_" + stim + "_" + sex_list[0][0] + "_" + str(laser_freq) + "Hz_" + protocol
    plt.title("%s" %suptitle, fontsize = 18)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
    plt.savefig(analysis_path + "/" + suptitle + ".png")
    plt.savefig(analysis_path + "/" + suptitle + ".pdf")
    plt.show()
#    plt.close()    
    

   
#%%################################
#   Lick rate
###################################   

for sex_age in male_list + female_list:
    try:    df_lick_rate_raw[sex_age]          = df_lick_rate_raw[sex_age].drop(columns = ['dummy'])
    except: pass #print ("dummy column already removed from df_lick_cum_raw[%s]" %sex_age)    



for k, sex_list in enumerate ([male_list, female_list]):  
    plt.figure(figsize = (16, 6))
    for i, sex_age in enumerate (sex_list):
        lick_rate_mean       = df_lick_rate_raw[sex_age].mean(axis=1)
        lick_rate_stderr     = ss.sem(df_lick_rate_raw[sex_age], axis=1)
#        # Print and save mean figures 
        legend              = sex_age + " " + "n=%d" %len(df_lick_rate_raw[sex_age].columns)
#            markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5, label = legend)
        markers, caps, bars = plt.errorbar(time_sec/60, lick_rate_mean, lick_rate_stderr, fmt='-', linewidth = 2, ecolor = palette2[i], label = legend)
        markers.set_color(palette[i])
        [bar.set_alpha(0.02) for bar in bars]
        [cap.set_alpha(0.02) for cap in caps]   
    plt.legend(loc="upper left", bbox_to_anchor = (1,1), fontsize = 18)
    plt.tick_params(axis = "both", which = "major", labelsize = 18)
    plt.xlabel("time (min)", fontsize = 18)
    plt.ylabel("lick rate (Hz)", fontsize = 18)
    plt.ylim(-1, 23)
    
    for j in np.arange(total_number_laser):
        plt.axvspan(laser_on_times[j], laser_on_times[j]+laser_duration, color = "red", alpha=0.2)
    plt.axvline (x = water_on_time, linewidth = 1, linestyle = "--", color = "gray")

    plt.xlim(0, 50)
    plt.ylim(0, 12)
#    plt.title("raw", fontsize = 18) 
    suptitle = "lickrate_" + stim + "_" + sex_list[0][0] + "_" + str(laser_freq) + "Hz_" + protocol
    plt.title("%s" %suptitle, fontsize = 18)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
    plt.savefig(analysis_path + "/" + suptitle + ".png")
    plt.savefig(analysis_path + "/" + suptitle + ".pdf")
    plt.show()
#    plt.close()    

#%%################################
#   Cumulative water intake (g)
###################################   

for sex_age in male_list + female_list:
    try:    df_lick_cum_raw[sex_age]          = df_lick_cum_raw[sex_age].drop(columns = ['dummy'])
    except: pass #print ("dummy column already removed from df_lick_cum_raw[%s]" %sex_age)    
    try:    df_lick_cum_raw_BW[sex_age]       = df_lick_cum_raw_BW[sex_age].drop(columns = ['dummy'])
    except: pass #print ("dummy column already removed from df_lick_cum_firstbout[%s]" %sex_age)  
 
# set maximum value for y-axis 
for k, sex_list in enumerate ([male_list, female_list]):  
    for i, sex_age in enumerate (sex_list):
        lick_cum_mean       = df_lick_cum_raw[sex_age].mean(axis=1)
        lick_cum_stderr     = ss.sem(df_lick_cum_raw[sex_age], axis=1)
        max_cumlick         = max(max_cumlick, float(lick_cum_mean[-1:] + lick_cum_stderr[-1:]))
        
for k, sex_list in enumerate ([male_list, female_list]):  
    plt.figure(figsize = (10, 8))
    for i, sex_age in enumerate (sex_list):
        lick_cum_mean       = df_lick_cum_raw[sex_age].mean(axis=1)*lick_g
        lick_cum_stderr     = ss.sem(df_lick_cum_raw[sex_age], axis=1)*lick_g
#        # Print and save mean figures 
        legend              = sex_age + " " + "n=%d" %len(df_lick_cum_raw[sex_age].columns)
#            markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='.', ms=8, mfc = palette[i], mec = palette[i], ecolor = palette2[i], alpha = 0.5, label = legend)
        markers, caps, bars = plt.errorbar(time_sec/60, lick_cum_mean, lick_cum_stderr, fmt='-', linewidth = 2, ecolor = palette2[i], label = legend)
        markers.set_color(palette[i])
        [bar.set_alpha(0.02) for bar in bars]
        [cap.set_alpha(0.02) for cap in caps]   
    plt.legend(loc="upper left", bbox_to_anchor = (1,1), fontsize = 18)
    plt.tick_params(axis = "both", which = "major", labelsize = 18)
    plt.xlabel("time (min)", fontsize = 18)
    plt.ylabel("water intake (mg)", fontsize = 18)
    plt.ylim(0, np.max(max_cumlick)*lick_g*1.05)
    # plt.ylim(0, 2000)
    for j in np.arange(total_number_laser):
        plt.axvspan(laser_on_times[j], laser_on_times[j]+laser_duration, color = "red", alpha=0.2)
    plt.axvline (x = water_on_time, linewidth = 1, linestyle = "--", color = "gray")
    plt.xlim(0, 50)
#    ax.set_ylim(0, 800)
#    plt.title("raw", fontsize = 18) 
    suptitle = "cumwater_" + stim + "_" + sex_list[0][0] + "_" + str(laser_freq) + "Hz_" + protocol
    plt.title("%s" %suptitle, fontsize = 18)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
    plt.savefig(analysis_path + "/" + suptitle + ".png")
    plt.savefig(analysis_path + "/" + suptitle + ".pdf")
    plt.show()
#    plt.close()     
    


#%%################################################################
#   total cumulative number of licks at the end of the experiment
###################################################################
for k, sex_age in enumerate (male_list + female_list):  
    try:    df_lick_cum_raw[sex_age]                = df_lick_cum_raw[sex_age].drop(columns = ['dummy'])    
    except:
        pass

    try:    df_lick_cum_raw_BW[sex_age]             = df_lick_cum_raw_BW[sex_age].drop(columns = ['dummy'])    
    except:
        pass
    
    try:    df_lick_cum_raw_temp[sex_age]           = df_lick_cum_raw_temp[sex_age].drop(columns = ['dummy'])
    except:
        pass
    
    try:    df_lick_cum_raw_temp_BW[sex_age]        = df_lick_cum_raw_temp_BW[sex_age].drop(columns = ['dummy'])
    except:
        pass
    
dic_total_cumlick    = {}
dic_total_cumlick_BW = {}
for k, sex_age in enumerate (male_list + female_list):    
    for j, expID in enumerate(df_lick_cum_raw_temp[sex_age].columns):
        dic_total_cumlick[expID] = df_lick_cum_raw_temp[sex_age][expID][length_min*60-1]
    
    for j, expID in enumerate(df_lick_cum_raw_temp_BW[sex_age].columns):
        dic_total_cumlick_BW[expID] = df_lick_cum_raw_temp_BW[sex_age][expID][length_min*60-1]

csv_path    = analysis_path + "/total_cumlick_" + stim + "_" + str(laser_freq) + "Hz_" + protocol + ".csv"      
with open(csv_path, 'w') as f:
    for key in dic_total_cumlick.keys():
        f.write("%s, %s\n" % (key, dic_total_cumlick[key]))  

csv_path    = analysis_path + "/total_cumlick_BS_" + stim + "_" + str(laser_freq) + "Hz_" + protocol + ".csv"      
with open(csv_path, 'w') as f:
    for key in dic_total_cumlick_BW.keys():
        f.write("%s, %s\n" % (key, dic_total_cumlick_BW[key]))



#%%###############################################
#   POP open the analysis folder
##################################################
subprocess.Popen(f'explorer {os.path.realpath(analysis_path)}')    
    
