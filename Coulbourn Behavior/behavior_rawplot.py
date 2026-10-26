# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 16:31:57 2020

Mouse drinking behavior analysis and visualization.

Reads raw licking behavior data exported from Coulbourn Graphic State
and integrates mouse metadata, including age, sex, and body weight.

Identifies individual licks and drinking bouts, calculates lick rates
and cumulative licking, and determines first-lick and first-bout
latencies.

Generates cumulative licking plots aligned to recording onset,
first lick, and first bout, and saves processed behavioral data
as .pkl files and figures as .png files.

"""

import numpy as np
import pandas as pd
#import scipy as sc
import matplotlib.pyplot as plt
#%matplotlib inline
import os
import glob
import csv
from datetime import date
import math
import pickle as pk
 
# path = r"C:\Users\heeun\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Mouse_Behavior_to_analyze\20210115_opto1"         


### To input the path from the command line
print("Enter the full path for the directory with raw graphic state data: ")
path             = input()

baseline = 0
# baseline = 1  # set as 0 if NOT baseline


## Define directories 
analysis_path = path + "/analysis" 
if not os.path.exists(analysis_path):
    os.makedirs(analysis_path)


## lick rate by second and cumulative licks
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


#%% 0. Read the metadata file first. 
#df_metadata       = pd.read_csv("@@@@")
# with open(r'C:\Users\heeun\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Mouse_metadata\allmiceDOB.csv', mode='r') as infile:
# with open(r"F:\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Mouse_metadata\allmiceDOB.csv", mode = "r") as infile:
with open(r"D:\Data_Dropbox_Sync\Mouse_metadata\allmiceDOB.csv", mode='r') as infile:    
    reader = csv.reader(infile)
    metadata = {rows[0]:rows[1] for rows in reader}

#%% 0-1. Read another metadata for BW for each experiment date.
# df_BW = pd.read_csv(r"C:\Users\heeun\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Mouse_metadata\allmiceBW.csv")
# df_BW = pd.read_csv(r"F:\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Mouse_metadata\allmiceBW.csv")
df_BW = pd.read_csv(r"D:\Data_Dropbox_Sync\Mouse_metadata\allmiceBW.csv")

#%% 1. Read all files under the directory 
# if baseline == 1:
# #    all_files       = glob.iglob(path + "**/*.csv", recursive = True)
all_files       = glob.glob(path + "/*.csv")
# else:
# all_files       = glob.iglob(path + "**/*/*.csv", recursive = True)
file_list       = [] 
mouseID_list    = []
rawfile         = {} # Dictionary of keys = file name, values = dataframe of raw csv.
data            = {}


## Go through each file in the directory to get data    
for i, file_name in enumerate(all_files):
    if baseline == 1:
        file_name_cut            = file_name[(len(path)+1):-4]
    else:
        file_name_cut            = (file_name[(len(path)+1):-4]).split("\\")[0]  
    file_list.append(file_name_cut)
    rawfile[file_name_cut]   = pd.read_csv(file_name, index_col = None, header = 0)
    data["mouseID"]            = file_name_cut.split("_")[7][:-3].upper() + file_name_cut.split("_")[7][-3:]
#    data["mouseID"]            = file_name_cut.split("_")[7][:-3] + file_name_cut.split("_")[7][-3:]

    if "W1378" in data["mouseID"]:
        data["mouseID"] = data["mouseID"].replace("W1378", "HW34") ## unify mouseID with multiple names
    elif "W1355" in data["mouseID"]:
        data["mouseID"] = data["mouseID"].replace("W1355", "HW111")  
    data["date"]             = file_name_cut.split("_")[0] + file_name_cut.split("_")[1] + file_name_cut.split("_")[2]
    data["time"]             = file_name_cut.split("_")[4] + file_name_cut.split("_")[5] + file_name_cut.split("_")[6]   
    if "m" in data["mouseID"]:
        data["sex"]       = "M"
    else: 
        data["sex"]       = "F"

    # baseline
    if baseline == 1:
        data["stimulus"]       = "baseline"
        data["totalmins"]      = 60*24
    # stimulus 
    else:
        data["stimulus"]         = file_name_cut.split("_")[-2] 
        if "WD" in data["stimulus"] or "wd" in data["stimulus"] or "Wd" in data["stimulus"] or "wD" in data["stimulus"]:
            data["stimulus"] = "24hrWD"
        elif "Saline" in data["stimulus"]: 
            data["stimulus"] = "saline"
            
        protocol = file_name_cut.split("_")[-1]
        if "hr" in protocol: ## if protocal name is e.g. 3hrlick
            data["totalmins"] = int(protocol.split("hr")[0])*60
        else:  ## if protocol name is e.g. 45minlick 
            data["totalmins"] = int(protocol.split("min")[0])

    time_sec = np.arange(1, data["totalmins"]*60+1)    
    # find DOB from the metadata file    
    for key in metadata.keys():
        if key == data["mouseID"]:
            data["DOB"] = str(metadata[key])
            break  
    # Print no DOB info!!
   
    # find mouse age based on DOB and experiment date        
    DOB     = date(int(data["DOB"][0:4]), int(data["DOB"][4:6]), int(data["DOB"][6:8])) 
    expdate = date(int(data["date"][0:4]), int(data["date"][4:6]), int(data["date"][6:8])) 
    age_days    = (expdate - DOB).days  
    age_months  = math.floor(age_days*0.032876676)
    data["age"] = age_months
       
    # find body weight info for the mouse ID and experiment date    
    for i in range(len(df_BW)):
        if (df_BW.iloc[i]["mouseID"] == data["mouseID"] and df_BW.iloc[i]["expdate"] == int(data["date"])):
            data["BW"]  = df_BW.iloc[i]["BW"]
            break
    try:
        data["BW"]
    except:
        print ("%s does not have BW" %file_name_cut)
        
    timestamp       = np.array(rawfile[file_name_cut]["Time"])
    if baseline == 1:
        feed_on_mask       = [i == 1 for i in rawfile[file_name_cut]["On1A3"]] 
        feedontimes_raw    = timestamp[feed_on_mask]
    
    try:
        lick_on_mask    = [i == 1 for i in rawfile[file_name_cut]["optical"]]
        lick_off_mask   = [i == 1 for i in rawfile[file_name_cut]["opticaloff"]] 

        lickontimes_raw    = timestamp[lick_on_mask]
        lickofftimes_raw   = timestamp[lick_off_mask]
    
    except:  

#    if len(lickontimes_raw) <1:
#        print("%s has no lick data by \"optical\"" %file_name_cut)
        try:
            lick_on_mask    = [i == 1 for i in rawfile[file_name_cut]["On1A1"]]
            lick_off_mask   = [i == 1 for i in rawfile[file_name_cut]["Off1A1"]] 
        
            lickontimes_raw    = timestamp[lick_on_mask]
            lickofftimes_raw   = timestamp[lick_off_mask]
        except: 
            print("%s has no lick data by \"On1A1\"" %file_name_cut)
            
#    if len(lickontimes_raw) <1:
            try:
                lick_on_mask    = [i == 1 for i in rawfile[file_name_cut]["contact"]]
                lick_off_mask   = [i == 1 for i in rawfile[file_name_cut]["contactoff"]] 
        
                lickontimes_raw    = timestamp[lick_on_mask]
                lickofftimes_raw   = timestamp[lick_off_mask]
            except:
                print("%s has no lick data by \"contact\"" %file_name_cut)
            
    if len(lickontimes_raw) <1:
        continue
    
    data["lickontimes_raw"]         = lickontimes_raw
    data["lickofftimes_raw"]        = lickofftimes_raw
    data["interlick_interval"]      = np.diff(lickontimes_raw)
    
    data["first_lick_sec"]          = lickontimes_raw[0]
    data["lickontimes_firstlick"]   = lickontimes_raw  - data["first_lick_sec"]
    data["lickofftimes_firstlick"]  = lickofftimes_raw - data["first_lick_sec"]
 
    # bout analysis
    boutontimes_raw, boutofftimes_raw, boutsize, boutlengths = bout_analysis(lickontimes_raw, 1, 10)
    nobout     = 0
    try: 
        data["first_bout_sec"]      = boutontimes_raw[0]   
    except:
        print("%s has licks but no bouts" %file_name_cut)
        data["first_bout_sec"]      = lickontimes_raw[0]
        nobout                      = 1
        
    data["lickontimes_firstbout"]   = lickontimes_raw  - data["first_bout_sec"]
    data["interbout_interval"]      = np.diff(boutontimes_raw)   

    # if baseline, process feed data
    if data["stimulus"]  == "baseline":
        data["feedontimes_raw"]         = feedontimes_raw - 10

    # calculate lick rate and cumulative lick to make a figure
    lickrate_sec_raw, lickrate_cum_raw             = lick_rate_by_sec(lickontimes_raw, time_sec)   
    lickrate_sec_firstlick, lickrate_cum_firstlick = lick_rate_by_sec(data["lickontimes_firstlick"], time_sec)
    lickrate_sec_firstbout, lickrate_cum_firstbout = lick_rate_by_sec(data["lickontimes_firstbout"], time_sec)  

    # if baseline, process feed data
    if data["stimulus"]  == "baseline":
        feedrate_sec_raw, feedrate_cum_raw             = lick_rate_by_sec(data["feedontimes_raw"], time_sec) 

    suptitle = data["mouseID"] +"_"+ data["date"] +"_"+ str(data["age"]) + "mo_" + data["sex"] + "_" + data["stimulus"] + "_" + str(data["totalmins"]) + "min"
    
    ## save pkl data for each experiment 
    save_file_path  = analysis_path + "/" + suptitle + ".pkl"
    output = open(save_file_path, 'wb')
    pk.dump(data, output)
    output.close()   
    
    #%% Figure
    if data["stimulus"]  == "baseline":
        ## print raw cumulative lick
        fig, ax = plt.subplots(2, 1, figsize = (12,12)) 
        #    Plt.subplot2grid(shape, location, rowspan, colspan)
        ax[0] = plt.subplot2grid((2, 1), (0, 0), rowspan = 2)
        ax[0].plot(time_sec/60, lickrate_cum_raw, color = "blue", linewidth = 2, label = "lick")
        ax[0].set_xlabel("time (min)", fontsize = 16)
        ax[0].set_ylabel("cumulative licks", fontsize = 16)
        ax[0].set_ylim(-lickrate_cum_raw[-1]/25, 1.05*lickrate_cum_raw[-1])
        ax[0].tick_params(labelsize = 16) 
        ax[0].annotate("total cumulative licks = %d\ntotal cumulative feeds = %d\nBW = %.1f (g)" %(lickrate_cum_raw[-1], feedrate_cum_raw[-1], data["BW"]), xy = (data["totalmins"]*0.4, 0.5), fontsize = 18)           
    
        ax[1] = ax[0].twinx()
        ax[1].plot(time_sec/60, feedrate_cum_raw, color = "black", linewidth = 2, label = "feed")
        ax[1].set_ylabel("cumulative feeds", fontsize = 16)
        ax[1].set_ylim(-feedrate_cum_raw[-1]/25, 1.05*feedrate_cum_raw[-1])
        ax[1].tick_params(labelsize = 16)
        ax[1].legend(loc='upper left', fontsize = 16)         

#        ax[1] = plt.subplot2grid((5, 1), (0, 1), rowspan = 1)        
#        ax[1].eventplot(feedontimes_raw/60, color="black", lineoffsets = 2, linelengths = 2)
#        ax[1].set_xlim(-10, time_sec[-1]/60)
#        ax[1].set_xlabel("feeding times (sec)")
#        ax[1].set_title("feed rastor")
#        
#        ax[2] = plt.subplot2grid((5, 1), (0, 2), rowspan = 1)     
#        ax[2].eventplot(lickontimes_raw/60, color="blue", lineoffsets = 2, linelengths = 2)
#        ax[2].set_xlim(-10, time_sec[-1]/60)
#        ax[2].set_xlabel("lick times (sec)")
#        ax[2].set_title("lick rastor")
          
        plt.suptitle("%s" %suptitle, fontsize = 18) 
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])        
        save_file_path  = analysis_path + "/" + suptitle + ".png"
        plt.savefig(save_file_path)
        plt.close(plt.gcf())
        
    else: 
        fig, ax = plt.subplots(3, 1, figsize = (8,12))
    #    ax[0] = plt.subplot2grid((3, 1), (0, 0), colspan=8)
        ax[0].plot(time_sec/60, lickrate_cum_raw, color = "black", linewidth = 2, label = "raw")
        ax[0].set_xlabel("time (min)", fontsize = 16)
        ax[0].set_ylabel("cumulative licks", fontsize = 16)
        ax[0].set_ylim(-lickrate_cum_raw[-1]/25, 1.05*lickrate_cum_raw[-1])
        ax[0].tick_params(labelsize = 16)
        ax[0].legend(loc='upper left', fontsize = 16)   
        if nobout == 1:
            ax[0].annotate("first lick: %.2f sec\nno bouts\ntotal cumulative licks = %d" %(data["first_lick_sec"], lickrate_cum_raw[-1]), xy = (data["totalmins"]*0.4, 0.5), fontsize = 18)           
        else:
            ax[0].annotate("first lick: %.2f sec\nfirst bout %.2f sec\ntotal cumulative licks = %d" %(data["first_lick_sec"], data["first_bout_sec"], lickrate_cum_raw[-1]), xy = (data["totalmins"]*0.4, 0.5), fontsize = 18)    
        
        ## print cumulative lick aligned to first lick 
    #    ax[1] = plt.subplot2grid((3, 1), (1, 0), colspan=8)
        ax[1].plot(time_sec/60, lickrate_cum_firstlick, color = "green", linewidth = 2, label = "aligned to first lick" )
        ax[1].set_xlabel("time (min)", fontsize = 16)
        ax[1].set_ylabel("cumulative licks", fontsize = 16)
        ax[1].set_ylim(-lickrate_cum_firstlick[-1]/25, 1.05*lickrate_cum_firstlick[-1])
        ax[1].tick_params(labelsize = 16)
        ax[1].legend(loc='lower right', fontsize = 16)
    
        ## print cumulative lick aligned to first bout      
    #    ax[2] = plt.subplot2grid((3, 1), (2, 0), colspan=8)
        ax[2].plot(time_sec/60, lickrate_cum_firstbout, color = "blue", linewidth = 2, label = "aligned to first bout")
        ax[2].set_xlabel("time (min)", fontsize = 16)
        ax[2].set_ylabel("cumulative licks", fontsize = 16)
        ax[2].set_ylim(-lickrate_cum_firstlick[-1]/25, 1.05*lickrate_cum_firstbout[-1])
        ax[2].tick_params(labelsize = 16)
        ax[2].legend(loc='lower right', fontsize = 16)
          
        plt.suptitle("%s" %suptitle, fontsize = 18) 
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])        
        save_file_path  = analysis_path + "/" + suptitle + ".png"
        plt.savefig(save_file_path)
        plt.close(plt.gcf())

        ## lickrate checking
    #    ax[0] = plt.subplot2grid((3, 1), (0, 0), colspan=8)
        plt.figure()    
        plt.plot(time_sec/60, lickrate_sec_firstlick, color = "black", linewidth = 1)
        plt.xlabel("time (min)", fontsize = 16)
        plt.ylabel("lick rate(/sec)", fontsize = 16)
        plt.tick_params(labelsize = 16)
        plt.annotate("max lickrate = %d" %np.max(lickrate_sec_firstlick), xy = (data["totalmins"]*0.4, np.max(lickrate_sec_firstlick)*0.9), fontsize = 18)           
    
        suptitle = "lickrate_" + suptitle
        plt.suptitle("%s" %suptitle, fontsize = 18) 
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])   
              
        save_file_path  = analysis_path + "/" + suptitle + ".png"
        plt.savefig(save_file_path)
        plt.close(plt.gcf())  
    
    

        