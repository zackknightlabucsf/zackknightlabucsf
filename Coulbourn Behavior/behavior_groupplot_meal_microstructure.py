# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:23:30 2020

Mouse feeding and drinking microstructure analysis.

Loads processed behavioral data (.pkl) and analyzes feeding and
licking events during selected dark or light phases.

Identifies feeding meals and drinking bouts, calculates meal size,
duration, frequency, and intermeal intervals, and measures the
latency from pellet consumption or meal completion to subsequent
drinking bouts.

Quantifies temporal differences between normalized cumulative
food and water intake and exports behavioral measurements
as .csv files.

"""
#%%
import os
import glob
import pickle as pk
import pandas as pd
import numpy as np
import math
from scipy.integrate import simps
from datetime import datetime
import subprocess


### To input the path from the command line
print("Enter the full path for the directory with .pkl files from behavior_rawplot.py: ")
path             = input()

# stim      = "FD"
stim       = "BL"

feeder    = 1

dark      = 1
light     = 0


# length_cutoff     = 4*60
length_cutoff   = 24*60 # do not change this for dark or light cycle


if feeder == 1:    

    foodspill_g_path = r"D:\Data_Dropbox_Sync\Mouse_metadata\foodspill_g.csv"
    foodspill_g_data = pd.read_csv(foodspill_g_path)
    
    

##############################################
###     time window for individual raster plot
##############################################
window_hr           = 1
window_sec          = int(window_hr*60*60)
window_sec_range    = np.arange(0, int(length_cutoff/60/window_hr), 1) 
# window_sec_range    = [5, 17]

today = datetime.today().strftime('%Y%m%d')


analysis_path = path + "/analysis_" + today + "_raster"

if dark == 1:
    analysis_path = analysis_path + "_dark"

elif light == 1:
    analysis_path = analysis_path + "_light"
    
if not os.path.exists(analysis_path):
    os.makedirs(analysis_path)   

time_min  = np.arange(1, length_cutoff+1)  
if dark == 1:
    time_min = time_min[:int(len(time_min)/2)]
    
elif light == 1:
    time_min = time_min[int(len(time_min)/2):]
    
# ##############################################
# #   pellet in grams 
# #   water volume per lick in mL 
# ##############################################
pellet_g    = 0.018
lick_g      = 1.21/1000


# ##############################################
# #   definition of a meal in mice = 50 mg  (ref PMID: 1523228, PMID: 33539467)
# ##############################################

meal_definition_g  = 0.05 


# ##############################################
# #   definition of minimum intermeal interval to be considered as a separate meal (ref PMID: 1523228, PMID: 33539467)
# ##############################################

if light  == 1:
    intermeal_interval_sec = 8*60   # 8 minutes

else: # if 24hrs or dark cycle
    intermeal_interval_sec = 5*60   # 5 minutes

###########################################################################
###     Define functions      
###########################################################################
        
def lick_rate (lick_on, time_bins):

    lick_on_ceil = np.ceil(lick_on)
    lick_on_ceil = lick_on_ceil.astype(int)
        
    lick_rate_ = np.histogram(lick_on_ceil, bins = time_bins)   
    lick_rate  = np.concatenate([[0], lick_rate_[0]])
    lick_cum   = np.cumsum(lick_rate)

    return (lick_rate, lick_cum)

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

#########################################################################
###     Process data     
###########################################################################

data = pd.DataFrame(columns = ["expID", "mouseID", "expdate", "BW","sex", "age", "foodspill_g", "total_pellets",
                               "actual_food_g_per_pellet", "num_pellet_for_meal", "roundup_num_pellet_for_meal",
                               "lick_on_raw", "feed_on_raw", "bout_on", "fbout_on", "fbout_off", 
                               "meal_size", "meal_length", "meal_numbers", "intermeal_interval",
                               "abs_norm_food_minus_water_cum",
                               "bout_on_after_pellet","pellet_bout_interval",
                               "bout_on_after_fbout","fbout_bout_interval"])

data["lick_on_raw"]                   = data["lick_on_raw"].astype(object)
data["feed_on_raw"]                   = data["feed_on_raw"].astype(object)
data["bout_on"]                       = data["bout_on"].astype(object)
data["fbout_on"]                      = data["fbout_on"].astype(object)
data["fbout_off"]                     = data["fbout_off"].astype(object)
data["abs_norm_food_minus_water_cum"] = data["abs_norm_food_minus_water_cum"].astype(object)
data["bout_on_after_pellet"]          = data["bout_on_after_pellet"].astype(object)   
data["pellet_bout_interval"]          = data["pellet_bout_interval"].astype(object)
data["bout_on_after_fbout"]           = data["bout_on_after_fbout"].astype(object)   
data["fbout_bout_interval"]           = data["fbout_bout_interval"].astype(object)


data["meal_size"]                     = data["meal_size"].astype(object)
data["meal_length"]                   = data["meal_length"].astype(object)
data["meal_numbers"]                  = data["meal_numbers"].astype(object)
data["intermeal_interval"]            = data["intermeal_interval"].astype(object)


# Go through each file in the directory to get .pkl data
for i, file_path in enumerate (all_files):
    expID        = file_path[(len(path)+1):-4] ## previous called as "file_name_cut
    print ("processing %s" %expID)
    file_list.append(expID) 
    inputfile        = open(file_path, 'rb')
    tempdata         = pk.load(inputfile)
    inputfile.close()
    
    data.at[i, "expID"]         = expID
    data.at[i, "mouseID"]       = tempdata["mouseID"]
    data.at[i, "age"]           = tempdata["age"]
    data.at[i, "sex"]           = expID.split("_")[3]
    data.at[i, "expdate"]       = expID.split("_")[1]
    data.at[i, "BW"]            = tempdata["BW"]        
    data.at[i, "lick_on_raw"]   = tempdata["lickontimes_raw"]
    data.at[i, "feed_on_raw"]   = tempdata["feedontimes_raw"]
    total_pellets_24hrs         = len(data.at[i, "feed_on_raw"])

      
    if dark == 1:
            data.at[i, "lick_on_raw"] = data.at[i, "lick_on_raw"][data.at[i, "lick_on_raw"] < 24*60*60/2]
            data.at[i, "feed_on_raw"] = data.at[i, "feed_on_raw"][data.at[i, "feed_on_raw"] < 24*60*60/2]
            
    elif light == 1:
            data.at[i, "lick_on_raw"] = data.at[i, "lick_on_raw"][data.at[i, "lick_on_raw"] >= 24*60*60/2]
            data.at[i, "feed_on_raw"] = data.at[i, "feed_on_raw"][data.at[i, "feed_on_raw"] >= 24*60*60/2]

      
    if feeder == 1:
        for j in range(len(foodspill_g_data)):   ## look for the identical expdate and mouseID to find BW
            expdate = str(foodspill_g_data.at[j, "expdate"])
            mouseID = str(foodspill_g_data.at[j, "mouseID"])
            if (str(expdate) == str(data.at[i, "expdate"])) & (mouseID == data.at[i, "mouseID"]):
                data.at[i, "foodspill_g"]  = foodspill_g_data.at[j, "foodspill_g"]
                print (data.at[i, "foodspill_g"])
                break       

            try: 
                data.at[i, "foodspill_g"]
            except:
                print("%s has no foodspill_g data" %expID)

    data.at[i, "total_pellets"]               = len(data.at[i, "feed_on_raw"])          
    data.at[i, "actual_food_g_per_pellet"]    = (data.at[i, "total_pellets"]*pellet_g - data.at[i, "foodspill_g"]*data.at[i, "total_pellets"]/total_pellets_24hrs)/(data.at[i, "total_pellets"]*pellet_g)*pellet_g
    data.at[i, "num_pellet_for_meal"]         = meal_definition_g /data.at[i, "actual_food_g_per_pellet"] 
    data.at[i, "roundup_num_pellet_for_meal"] = math.ceil(data.at[i, "num_pellet_for_meal"])    
              
    time_sec = np.arange(1, length_cutoff*60+1)  

    df_raster         = pd.DataFrame(np.nan, index = np.arange(length_cutoff*60), columns = ["dummy"])
    df_raster["feed"] = pd.DataFrame(data.at[i, "feed_on_raw"])
    df_raster["lick"] = pd.DataFrame(data.at[i, "lick_on_raw"])
    
    try:    
        df_raster     = df_raster.drop(columns = ["dummy"])
    
    except: 
        pass
    
    array_raster_temp = df_raster.values
    array_raster      = array_raster_temp.transpose()

    
    lickrate_min_raw, lickcum_min_raw   = lick_rate(data.at[i, "lick_on_raw"]/60, time_min)  
    feedrate_min_raw, feedcum_min_raw   = lick_rate(data.at[i, "feed_on_raw"]/60, time_min) 
    
    foodrate_min_raw        = feedrate_min_raw*pellet_g
    waterrate_min_raw       = lickrate_min_raw*lick_g
    food_minus_water_min    = foodrate_min_raw - waterrate_min_raw
    
    feedrate_min_norm           = feedrate_min_raw/feedcum_min_raw[-1]
    lickrate_min_norm           = lickrate_min_raw/lickcum_min_raw[-1]
    feed_minus_lick_min_norm    = feedrate_min_norm  - lickrate_min_norm
    
    feedcum_min_norm             = feedcum_min_raw/feedcum_min_raw[-1]
    lickcum_min_norm             = lickcum_min_raw/lickcum_min_raw[-1]
    feed_minus_lick_min_cum_norm = feedcum_min_norm - lickcum_min_norm
    
    data.at[i, "abs_norm_food_minus_water_cum"] = abs(feed_minus_lick_min_cum_norm)    
    
    # lickrate_sec_raw, lickcum_sec_raw   = lick_rate(tempdata["lickontimes_raw"], time_sec) 
    # feedrate_sec_raw, feedcum_sec_raw   = lick_rate(tempdata["feedontimes_raw"], time_sec) 
   
    # foodrate_sec_raw        = feedrate_sec_raw*pellet_g
    # waterrate_sec_raw       = lickrate_sec_raw*lick_g
    # food_minus_water_sec    = foodrate_sec_raw - waterrate_sec_raw    



    # ###################################################################
    # Analysis: time interval from each pellet to the first lick of first bout      
    # ####################################################################  

    bout_on, bout_off, bout_size, bout_length = bout_analysis(data.at[i, "lick_on_raw"], 1, 10)
    data.at[i, "bout_on"]                     = bout_on

    list_bout_on_after_pellet = []
    list_pellet_bout_interval = []
    
    for k, pellet_sec in enumerate(data.at[i, "feed_on_raw"]):
        pellet_sec_mask         = [j > pellet_sec for j in data.at[i, "bout_on"]]
        bout_on_masked          = data.at[i, "bout_on"][pellet_sec_mask]  
        for j, bout_on in enumerate(bout_on_masked):
            if bout_on > 0:
                bout_on_after_pellet    = bout_on_masked[j]
                break
        pellet_bout_interval    = bout_on_after_pellet - pellet_sec
        list_bout_on_after_pellet.append(bout_on_after_pellet)
        list_pellet_bout_interval.append(pellet_bout_interval)
        
    data.at[i, "bout_on_after_pellet"] = np.array(list_bout_on_after_pellet)    
    data.at[i, "pellet_bout_interval"] = np.array(list_pellet_bout_interval)
    data.at[i, "pellet_bout_interval"] = data.at[i, "pellet_bout_interval"][data.at[i, "pellet_bout_interval"] > 0]    
    
 
    
    ###################################################################
    # Analysis: time interval from last pellet of meal (defined by an individual animal) to the first lick of first bout          
    ####################################################################  

    fbout_on, fbout_off, fbout_size, fbout_length = bout_analysis(data.at[i, "feed_on_raw"], intermeal_interval_sec, data.at[i, "roundup_num_pellet_for_meal"])    

    data.at[i, "fbout_on"]                        = fbout_on
    data.at[i, "fbout_off"]                       = fbout_off    
    
    data.at[i, "meal_length"]                     = fbout_length
    data.at[i, "meal_size"]                       = fbout_size*data.at[i, "actual_food_g_per_pellet"] # meal size in g

    data.at[i, "intermeal_interval"]              = np.diff(data.at[i, "fbout_on"]) 
    data.at[i, "meal_numbers"]                    = [len(data.at[i, "fbout_on"])]






    list_bout_on_after_fbout = []
    list_fbout_bout_interval = []

    for k, pellet_sec in enumerate(data.at[i, "fbout_off"]):
        pellet_sec_mask         = [j > pellet_sec for j in data.at[i, "bout_on"]]
        bout_on_masked          = data.at[i, "bout_on"][pellet_sec_mask]  
        for j, bout_on in enumerate(bout_on_masked):
            if bout_on > 0:
                bout_on_after_fbout    = bout_on_masked[j]
                break
        fbout_bout_interval    = bout_on_after_fbout - pellet_sec
        # print (bout_on_after_fbout)
        list_bout_on_after_fbout.append(bout_on_after_fbout)
        list_fbout_bout_interval.append(fbout_bout_interval)
        
    data.at[i, "bout_on_after_fbout"] = np.array(list_bout_on_after_fbout)    
    data.at[i, "fbout_bout_interval"] = np.array(list_fbout_bout_interval)
    data.at[i, "fbout_bout_interval"] = data.at[i, "fbout_bout_interval"][data.at[i, "fbout_bout_interval"] > 0]      
    
    ###################################################################
#    Analysis: time interval from last pellet of meal (defined by an individual animal) to the first lick of first bout      
    ####################################################################  

    data.at[i, "fbout_on"]                        = fbout_on
    data.at[i, "fbout_off"]                       = fbout_off    
    
    list_bout_on_after_fbout = []
    list_fbout_bout_interval = []

    for k, pellet_sec in enumerate(data.at[i, "fbout_off"]):
        pellet_sec_mask         = [j > pellet_sec for j in data.at[i, "bout_on"]]
        bout_on_masked          = data.at[i, "bout_on"][pellet_sec_mask]  
        for j, bout_on in enumerate(bout_on_masked):
            if bout_on > 0:
                bout_on_after_fbout    = bout_on_masked[j]
                break
        fbout_bout_interval    = bout_on_after_fbout - pellet_sec
        # print (bout_on_after_fbout)
        list_bout_on_after_fbout.append(bout_on_after_fbout)
        list_fbout_bout_interval.append(fbout_bout_interval)
        
    data.at[i, "bout_on_after_fbout"] = np.array(list_bout_on_after_fbout)    
    data.at[i, "fbout_bout_interval"] = np.array(list_fbout_bout_interval)
    data.at[i, "fbout_bout_interval"] = data.at[i, "fbout_bout_interval"][data.at[i, "fbout_bout_interval"] > 0]        
    


#################################################################################################
#    SAVE DATA: filter selected columns of "data" dataframe into a new dataframe and save into a .csv file
#################################################################################################
abs_AUC              = pd.DataFrame()
interval_median      = pd.DataFrame()

abs_AUC              = {}
interval_median      = {}
interval_mean        = {}
finterval_median     = {}
finterval_mean       = {}
numpellet            = {}

for i in np.arange(len(data)):
    expID                   = data["expID"][i]
    abs_AUC[expID]          = simps(data.at[i, "abs_norm_food_minus_water_cum"])
    interval_median[expID]  = np.median(data.at[i, "pellet_bout_interval"])
    interval_mean[expID]    = np.mean(data.at[i, "pellet_bout_interval"])
    finterval_median[expID] = np.median(data.at[i, "fbout_bout_interval"])
    finterval_mean[expID]   = np.mean(data.at[i, "fbout_bout_interval"])
    numpellet[expID]        = data.at[i, "roundup_num_pellet_for_meal"]

suptitle       = stim + "_abs_AUC"
save_file_path = analysis_path + "/" + suptitle  + ".csv"
with open(save_file_path, 'w') as f:
    for key in abs_AUC.keys():
        f.write("%s, %s\n" % (key, abs_AUC[key]))

suptitle       = stim + "_pellet_bout_interval_med"
save_file_path = analysis_path + "/" + suptitle  + ".csv"
with open(save_file_path, 'w') as f:
    for key in interval_median.keys():
        f.write("%s, %s\n" % (key, interval_median[key]))

suptitle       = stim + "_pellet_bout_interval_mean"
save_file_path = analysis_path + "/" + suptitle  + ".csv"
with open(save_file_path, 'w') as f:
    for key in interval_mean.keys():
        f.write("%s, %s\n" % (key, interval_mean[key]))

suptitle       = stim + "_meal_waterbout_interval_med"
save_file_path = analysis_path + "/" + suptitle  + ".csv"
with open(save_file_path, 'w') as f:
    for key in finterval_median.keys():
        f.write("%s, %s\n" % (key, finterval_median[key]))

suptitle       = stim + "_meal_waterbout_interval_mean"
save_file_path = analysis_path + "/" + suptitle  + ".csv"
with open(save_file_path, 'w') as f:
    for key in finterval_mean.keys():
        f.write("%s, %s\n" % (key, finterval_mean[key]))


suptitle       = stim + "_roundup_num_pellet_for_meal"
save_file_path = analysis_path + "/" + suptitle  + ".csv"
with open(save_file_path, 'w') as f:
    for key in finterval_mean.keys():
        f.write("%s, %s\n" % (key, numpellet[key]))




df_bout_length_forcsv        = pd.DataFrame(np.nan, index = np.arange(1000), columns = ["dummy"])
df_bout_size_forcsv          = pd.DataFrame(np.nan, index = np.arange(1000), columns = ["dummy"])
df_interbout_interval_forcsv = pd.DataFrame(np.nan, index = np.arange(1000), columns = ["dummy"])
df_bout_numbers_forcsv       = pd.DataFrame(np.nan, index = np.arange(1000), columns = ["dummy"])
df_bite_numbers_forcsv       = pd.DataFrame(np.nan, index = np.arange(1000), columns = ["dummy"])
df_T_halfrise_forcsv         = pd.DataFrame(np.nan, index = np.arange(1000), columns = ["dummy"])

for i, expID in enumerate (file_list):
    df_bout_length_forcsv[expID]        = pd.Series(data.at[i, "meal_length"])
    df_bout_size_forcsv[expID]          = pd.Series(data.at[i, "meal_size"])       
    df_interbout_interval_forcsv[expID] = pd.Series(data.at[i, "intermeal_interval"])
    df_bout_numbers_forcsv[expID]       = pd.Series(data.at[i, "meal_numbers"])
        
        
try:    
    df_bout_length_forcsv           = df_bout_length_forcsv.drop(columns = ["dummy"])
    df_bout_size_forcsv             = df_bout_size_forcsv.drop(columns = ["dummy"])
    df_interbout_interval_forcsv    = df_interbout_interval_forcsv.drop(columns = ["dummy"])    
    df_bout_numbers_forcsv          = df_bout_numbers_forcsv.drop(columns = ["dummy"])
  
  
    
except: 
    pass        
        
csv_path    = analysis_path + "/meal_length.csv"
df_bout_length_forcsv.T.to_csv(csv_path, index = True)
        
csv_path    = analysis_path + "/meal_size_mg.csv"    
df_bout_size_forcsv.T.to_csv(csv_path, index = True)

csv_path    = analysis_path + "/intermeal_interval.csv"
df_interbout_interval_forcsv.T.to_csv(csv_path, index = True)

csv_path    = analysis_path + "/meal_numbers.csv"
df_bout_numbers_forcsv.T.to_csv(csv_path, index = True)


##################################################
#   POP open the analysis folder
##################################################
subprocess.Popen(f'explorer {os.path.realpath(analysis_path)}')  






