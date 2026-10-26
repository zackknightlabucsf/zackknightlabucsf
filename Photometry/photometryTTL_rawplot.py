# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 17:04:00 2020

Raw fiber photometry and TTL licking behavior analysis.

Reads two-chamber TDT recordings and mouse metadata, aligns GCaMP signals and TTL lick events to stimulus onset, and calculates dF/F0 and z-scores.

Identifies licking bouts, generates individual photometry and behavior plots, and calculates bout-aligned PSTHs and within-bout changes in neural activity.

Exports processed traces, licking and bout measurements, and PSTH results as .pkl, .csv, and figure files.

""" 

# Jupyter magic
#%matplotlib inline

import os
import pickle as pk
import numpy as np
import pandas as pd
import csv
# from tkinter.filedialog import askdirectory 
from datetime import date
import math
from scipy import signal
from sklearn.neighbors import KDTree
import matplotlib.pyplot as plt  # standard Python plotting library
import matplotlib
matplotlib.rcParams['font.size'] = 18 # set font size for all figures

# import the read_block function from the tdt package
from tdt import read_block #, StructType
# import subprocess


#%% Smoothing 
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



#%% 

def main (BLOCKPATH):
# if 1== 1:
    ##  Import DOB data file into a dictionary
    with open (r"D:\Data_Dropbox_Sync\Mouse_metadata\allmiceDOB.csv", mode='r') as infile:
        reader = csv.reader(infile)
        DOBdata = {rows[0]:rows[1] for rows in reader}
        
    TTLpath = r"D:\Data_Dropbox_Sync\Mouse_metadata\photometrymiceTTL.csv" # __.csv is for entire timecourse    
    TTLdata = pd.read_csv(TTLpath)    
    
    #%% Importing the Data
    print("Enter the full path for directory with raw photometry files: ")
    BLOCKPATH      = input()

    setup = 2
    if "setup3" in BLOCKPATH:
        setup = 3
 

    # fluid_open_min   = 0
    
    Z_pre_post_min   = [10, 30] ## put the stimuluss application/drug injection time as [pre_mm, post_mm]
    K_pre_post_min   = [10, 30] ## 10, 55 for injection + water exposure (inj, 10min, water open, 45min)
                              ## 10, 45 for WD
    Z_F0_window_min  = [5, 10]  ## setting F0 window as a median value of 5-10 mins (equals -5 t0 0 pre-stimulus)
    K_F0_window_min  = [5, 10]
    F0_window_min    = [5, 10]
    
    ds_factor1            = 5
    ds_factor2            = 5
    ds_factor3            = 5   
    ds_factor             = int(ds_factor1*ds_factor2*ds_factor3)
    sm_factor             = 61    ## smoothing over 300 frames / 20 Hz (fps) = 15 secs
            
    psth_width_sec  = 120      
    
    ds_factor              = 4
    # baseline              = 10  # in minutes 
    # baseline_index        = np.arange(0, baseline*60*ds_factor, 1)   ## make an array of indices for baseline 
        
    data  = read_block(BLOCKPATH)
    trace = {}   
    tree  = {} 
    df_bout_zs = {} 
    blockname = data.info.blockname
    
    ## for mouse in Z
    mouseZ           = blockname.split('_')[0]
    mouseK           = (blockname.split('_')[1]).split('-')[0]
    mouse_list       = [mouseZ, mouseK]
    # mouse_list       = [mouseZ]
    
    # mouseZ = "none"
    
    for i, mouseID in enumerate (mouse_list):
        trace[i]                   = {}
        # if i == 0:
        #     continue
        if mouseID == "none":
    #    if mouse == "W1764f22":
            if i == 0:  
                print("No mouse in Z chamber")
            else:
                print("No mouse in K chamber")
            continue
            
        else:
            trace[i]["date"]           = "20" + (blockname.split('_')[1]).split('-')[1]    ## saves as a string 20100408
            trace[i]["time"]           = data.info.utc_start_time  ## saves hh:mm:ss "10:04:10" as a string
            trace[i]["mouseID"]        = mouseID
    
            # find DOB from the metadata file    
            for key in DOBdata.keys():
                if key == mouseID:
                    trace[i]["DOB"] = str(DOBdata[key])
                    break  
            try: 
                trace[i]["DOB"]
            except:
                print ("**************%s has no DOB" %mouseID)
                continue
                    
            # find mouse age based on DOB and experiment date        
            DOB                   = date(int(trace[i]["DOB"][0:4]), int(trace[i]["DOB"][4:6]), int(trace[i]["DOB"][6:8])) 
            expdate               = date(int(trace[i]["date"][0:4]), int(trace[i]["date"][4:6]), int(trace[i]["date"][6:8])) 
            age_days              = (expdate - DOB).days  
            trace[i]["age_wk"]      = math.floor(age_days/7)
            trace[i]["age_mo"]      = math.floor(age_days*0.032876676)
            if trace[i]["age_mo"]   > 19:
                age_color = "maroon"
            else:
                age_color = "dimgrey"
            # print ("%s has age_wk %d and age_mo %d" %(mouseID, trace[i]["age_wk"], trace[i]["age_mo"]))      
    
            for j in range(len(TTLdata)):
                if (str(int(TTLdata.at[j, "ExpDate"])) == str(trace[i]["date"])) & (TTLdata.at[j, "mouseID"] == trace[i]["mouseID"]):
                    trace[i]["BW"]          = TTLdata.at[j, "BW"]
                    trace[i]["stimulus"]    = TTLdata.at[j, "State"]
                    trace[i]["LICK_fluid"]  = TTLdata.at[j, "Access"] 
                    stim_sec                = int(TTLdata.at[j, "Timezero_mm"])*60 + int(TTLdata.at[j, "Timezero_ss"])
                    print ("mouseID is %s, Expdate is %s, stimulus is %s, access is %s" %(trace[i]["mouseID"], trace[i]["date"], trace[i]["stimulus"], trace[i]["LICK_fluid"]))
                    break
            try: 
                trace[i]["stimulus"]
            except:
                print ("%s_%s does not have a matching TTL metadata" %(mouseID, expdate))
    
            if i == 0:
                sampling_rate       = data.streams["G__Z"].fs ## in Hz
                raw_GCaMP           = data.streams["G__Z"].data
                raw_UV              = data.streams["U__Z"].data
                pre_sec             = int(Z_pre_post_min[0])*60
                post_sec            = int(Z_pre_post_min[1])*60
    
            else:
                sampling_rate       = data.streams["G__K"].fs ## in Hz
                raw_GCaMP           = data.streams["G__K"].data
                raw_UV              = data.streams["U__K"].data
                pre_sec             = int(K_pre_post_min[0])*60
                post_sec            = int(K_pre_post_min[1])*60        
                
            ## Crop the trace[i] around the stimulus time
            crop_start_sec   = stim_sec - pre_sec
            crop_end_sec     = stim_sec + post_sec
            
            crop_start_frame = int(np.ceil(crop_start_sec*sampling_rate))
            crop_end_frame   = int(np.floor(crop_end_sec*sampling_rate))
            
            trace[i]["GCaMP"]   = raw_GCaMP[crop_start_frame:crop_end_frame+1]
            trace[i]["UV"]      = raw_UV[crop_start_frame:crop_end_frame+1]        
            # trace[i]["Fratio"]   = trace[i]["GCaMP"]/trace[i]["UV"]
            trace[i]["Fratio"]   = trace[i]["GCaMP"]
            num_frames           = len(trace[i]["GCaMP"])
            trace[i]["time_sec"] = np.linspace(1, num_frames, num_frames)/sampling_rate
            trace[i]["time_sec"] = trace[i]["time_sec"] - pre_sec
            
            
            ## F0 is the median value from -5 min to 0 pre-injection 
            if i == 0:           
                F0_window_start_frame = int(np.ceil(int(Z_F0_window_min[0])*60*sampling_rate))
                F0_window_end_frame   = int(np.floor(int(Z_F0_window_min[1])*60*sampling_rate))
            else:
                F0_window_start_frame = int(np.ceil(int(K_F0_window_min[0])*60*sampling_rate))
                F0_window_end_frame   = int(np.floor(int(K_F0_window_min[1])*60*sampling_rate))
    
            # F0                        = np.median(trace[i]["Fratio"][F0_window_start_frame:F0_window_end_frame])
            F0                        = np.median(trace[i]["Fratio"][F0_window_start_frame:F0_window_end_frame])
          
            trace[i]["dFF"]           = (trace[i]["Fratio"] - F0)/F0*100 
            trace[i]["frame"]         = np.arange(0, len(trace[i]["GCaMP"]), 1)
            
    
            ##  downsample all time series by a factor of 50 (1000Hz >> 20Hz. 1ms >> 50ms exposure time )       		
            temp_trace                = signal.decimate(trace[i]["dFF"], ds_factor1)
            temp_trace                = signal.decimate(temp_trace, ds_factor2)
            trace[i]["dFF_ds"]        = signal.decimate(temp_trace, ds_factor3)
            trace[i]["dFF_ds_sm"]     = smooth(trace[i]["dFF_ds"], window_len = sm_factor)
            trace[i]["time_sec_ds"]   = trace[i]["time_sec"][::ds_factor] ## take every 50th element from the list
    	    
            trace[i]["dFF_ds_sm_peak"] = np.max(trace[i]["dFF_ds_sm"])
            trace[i]["dFF_ds_sm_scaled"]  = trace[i]["dFF_ds_sm"]/trace[i]["dFF_ds_sm_peak"]
    
            ####### plot dFF with all dictionary keys      
            fig, ax = plt.subplots(5,1, figsize = (8,12))
            ax[0].plot(trace[i]["time_sec"]/60, trace[i]["GCaMP"], color = "green", linewidth = 2, label = "GCaMP")
            ax[0].set_ylim(0.95*min(trace[i]["GCaMP"]), 1.05*max(trace[i]["GCaMP"]))
            ax[0].set_ylabel("raw (a.u.)")
            ax[0].legend(loc='upper right', fontsize = 16)
        #    ax[1] = ax[0].twinx()      
        #    ax[1].plot(trace[i]["time_sec"]/60, trace[i]["UV"], color = "red", linewidth = 2, label = "UV")
        #    ax[1].set_ylim(0.95*min(trace[i]["UV"]), 1.05*max(trace[i]["UV"]))
        #    ax[1].set_ylabel("raw fluorescence (a.u.)")
        #         
            ax[1].plot(trace[i]["time_sec"]/60, trace[i]["UV"], color = "red", linewidth = 2, label = "UV")
            ax[1].set_ylim(0.95*min(trace[i]["UV"]), 1.05*max(trace[i]["UV"]))
            ax[1].set_ylabel("raw (a.u.)")
            ax[1].legend(loc='upper right', fontsize = 16)
            
        #    ax[2] = plt.subplot2grid((1, 1), (1, 0), colspan=8)
            ax[2].plot(trace[i]["time_sec"]/60, trace[i]["Fratio"], color = "black", linewidth = 2, label = "GCaMP/UV")
            ax[2].set_ylim(0.95*min(trace[i]["Fratio"]), 1.05*max(trace[i]["Fratio"]))
            ax[2].set_ylabel("ratio")
            ax[2].legend(loc='upper right', fontsize = 16)
            
            ax[3].plot(trace[i]["time_sec"]/60, trace[i]["dFF"], color = "blue", linewidth = 2, label = "dFF")
            if min(trace[i]["dFF"]) > 0:
                ax[3].set_ylim(0.95*min(trace[i]["dFF"]), 1.05*max(trace[i]["dFF"]))
            else: 
                ax[3].set_ylim(1.05*min(trace[i]["dFF"]), 1.05*max(trace[i]["dFF"]))
            ax[3].set_ylabel("dF/F0 (%)")  
            ax[3].legend(loc='upper right', fontsize = 16)
              
            ax[4].plot(trace[i]["time_sec_ds"]/60, trace[i]["dFF_ds_sm"], color = "navy", linewidth = 2, label = "dFF_ds_sm")
            if min(trace[i]["dFF_ds_sm"]) > 0:
                ax[4].set_ylim(0.95*min(trace[i]["dFF_ds_sm"]), 1.05*max(trace[i]["dFF_ds_sm"]))
            else: 
                ax[4].set_ylim(1.05*min(trace[i]["dFF_ds_sm"]), 1.05*max(trace[i]["dFF_ds_sm"]))
            ax[4].set_ylabel("dF/F0 (%)") 
            ax[4].set_xlabel("time (min)")
            ax[4].legend(loc='upper right', fontsize = 16)
            
            suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_" + str(trace[i]["age_mo"]) + "mo" + "_noTTL"          
            plt.suptitle("%s" %suptitle) 
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])    
            
            analysis_path  = BLOCKPATH + "/analysis_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] #+ "_post5min"       
            if not os.path.exists(analysis_path):
                os.makedirs(analysis_path)    
    
            save_file_path = analysis_path + "/" + suptitle + ".png"
            plt.savefig(save_file_path)
            # save_file_path = analysis_path + "/" + suptitle + ".pdf"
            # plt.savefig(save_file_path)
            plt.close()
            #%%#############################################################################
            #   z score and figure
            ##################################################################################

            ## z score based on simple G/U 
            F0_window_start_frame = int(np.ceil(int(F0_window_min[0])*60*sampling_rate))
            F0_window_end_frame   = int(np.floor(int(F0_window_min[1])*60*sampling_rate))
                    
            trace[i]["zscore"]       = (trace[i]["Fratio"] - np.mean(trace[i]["Fratio"][F0_window_start_frame:F0_window_end_frame]))/np.std(trace[i]["Fratio"])
        
            temp_trace                = signal.decimate(trace[i]["zscore"], ds_factor1)
            temp_trace                = signal.decimate(temp_trace, ds_factor2)
            trace[i]["zscore_ds"]        = signal.decimate(temp_trace, ds_factor3)
            trace[i]["time_sec_ds_zs"]   = trace[i]["time_sec"][::ds_factor] ## take every 50th element from the list    
                
                
            ##  plot z score, not fitted, simple G/U       
            fig, ax = plt.subplots(5, 1, figsize = (12, 20))
            
            ## subplot 0. raw GCaMP
            ax[0].plot(trace[i]["time_sec"]/60, trace[i]["GCaMP"], color = "green", linewidth = 2, label = "GCaMP")
            ax[0].set_ylim(0.95*min(trace[i]["GCaMP"]), 1.05*max(trace[i]["GCaMP"]))
            ax[0].set_ylabel("raw (a.u.)")
            ax[0].legend(loc='upper right', fontsize = 14)
        
            ## subplot 1. raw UV
            ax[1].plot(trace[i]["time_sec"]/60, trace[i]["UV"], color = "red", linewidth = 2, label = "UV")
            ax[1].set_ylim(0.95*min(trace[i]["UV"]), 1.05*max(trace[i]["UV"]))
            ax[1].set_ylabel("raw (a.u.)")
            ax[1].legend(loc='upper right', fontsize = 14)
            
            ## subplot 2. F = raw GCaMP/raw UV
            ax[2].plot(trace[i]["time_sec"]/60, trace[i]["Fratio"], color = "black", linewidth = 2, label = "GCaMP/UV")
            ax[2].set_ylim(0.95*min(trace[i]["Fratio"]), 1.05*max(trace[i]["Fratio"]))
            ax[2].set_ylabel("ratio")
            ax[2].legend(loc='upper right', fontsize = 14)
            
            ## subplot 3. z score 
            ax[3].plot(trace[i]["time_sec"]/60, trace[i]["zscore"], color = "blue", linewidth = 2, label = "zscore")
            if min(trace[i]["zscore"]) > 0:
                ax[3].set_ylim(0.95*min(trace[i]["zscore"]), 1.05*max(trace[i]["zscore"]))
            else: 
                ax[3].set_ylim(1.05*min(trace[i]["zscore"]), 1.05*max(trace[i]["zscore"]))
            ax[3].set_ylabel("z score")  
            ax[3].legend(loc='upper right', fontsize = 14)
            
            ## subplot 4. zscore_ds = downsampled dFF  
            ax[4].plot(trace[i]["time_sec_ds_zs"]/60, trace[i]["zscore_ds"], color = "purple", linewidth = 2, label = "zscore ds")
            if min(trace[i]["zscore_ds"]) > 0:
                ax[4].set_ylim(0.95*min(trace[i]["zscore_ds"]), 1.05*max(trace[i]["zscore_ds"]))
            else: 				
                ax[4].set_ylim(1.05*min(trace[i]["zscore_ds"]), 1.05*max(trace[i]["zscore_ds"]))
            ax[4].set_ylabel("z score") 
            ax[4].set_xlabel("time (min)")
            ax[4].legend(loc='upper right', fontsize = 14)
            
            
            suptitle = trace[i]["mouseID"] +"_"+ trace[i]["date"] +"_"+ trace[i]["stimulus"] + "_zs"
            plt.suptitle("%s" %suptitle) 
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])    
          
            
            save_file_path = analysis_path + "/" + suptitle + ".png"
            plt.savefig(save_file_path)
            # save_file_path = analysis_path + "/" + suptitle + ".pdf"
            # plt.savefig(save_file_path)
            # plt.show()
            plt.close()   
        
            if "lickpellet" in BLOCKPATH:                
                if i == 0:
                    activelever_on_temp   = data.epocs["Epo1"]['onset']
                    activelever_off_temp  = data.epocs["Epo1"]['onset']
                    LICK_on_temp          = data.epocs["Epo2"]['onset']
                    LICK_off_temp         = data.epocs["Epo2"]['onset']
    
                else:
                    activelever_on_temp   = data.epocs["Epo3"]['onset']
                    activelever_off_temp  = data.epocs["Epo3"]['onset']
                    LICK_on_temp          = data.epocs["Epo4"]['onset']
                    LICK_off_temp         = data.epocs["Epo4"]['onset']        

                        
                activelever_on_temp  = activelever_on_temp[activelever_on_temp <1E308]
                activelever_off_temp = activelever_off_temp[activelever_off_temp <1E308]
               
                activelever_on_temp  = activelever_on_temp - stim_sec
                activelever_off_temp = activelever_off_temp - stim_sec
        
                activelever_on_temp  = activelever_on_temp[(activelever_on_temp > 0)]
                activelever_off_temp = activelever_off_temp[(activelever_off_temp > 0)]
        
                activelever_on_temp     = activelever_on_temp[(activelever_on_temp < trace[i]["time_sec"][-1])]
                activelever_off_temp    = activelever_off_temp[(activelever_off_temp < trace[i]["time_sec"][-1])]
        
           
            else:     
                if i == 0:
                    if setup == 3:
                        LICK_on_temp   = data.epocs["Ep2_"]['onset']
                        LICK_off_temp  = data.epocs["Ep2_"]['offset']
                    
                    else:
                        # LICK_on_temp   = data.epocs["PC1_"]['onset']
                        # LICK_off_temp  = data.epocs["PC1_"]['offset']
                        LICK_on_temp   = data.epocs["Ep1_"]['onset']
                        LICK_off_temp  = data.epocs["Ep1_"]['offset']
                        # LICK_on_temp   = data.epocs["Epo1"]['onset']
                        # LICK_off_temp  = data.epocs["Epo1"]['offset']
        
                else:
        
                    if setup == 3:
                        LICK_on_temp   = data.epocs["Ep3_"]['onset']
                        LICK_off_temp  = data.epocs["Ep3_"]['offset']
        
                    else:
                        # LICK_on_temp   = data.epocs["PC3_"]['onset']
                        # LICK_off_temp  = data.epocs["PC3_"]['offset']
                        LICK_on_temp   = data.epocs["Ep2_"]['onset']
                        LICK_off_temp  = data.epocs["Ep2_"]['offset']        
                        # LICK_on_temp   = data.epocs["Epo3"]['onset']
                        # LICK_off_temp  = data.epocs["Epo3"]['offset'] 
                   
            
            LICK_on_temp  = LICK_on_temp[LICK_on_temp <1E308]
            LICK_off_temp = LICK_off_temp[LICK_off_temp <1E308]
           
            LICK_on_temp  = LICK_on_temp - stim_sec
            LICK_off_temp = LICK_off_temp - stim_sec
    
            LICK_on_temp  = LICK_on_temp[(LICK_on_temp > 0)]
            LICK_off_temp = LICK_off_temp[(LICK_off_temp > 0)]
    
            LICK_on_temp     = LICK_on_temp[(LICK_on_temp < trace[i]["time_sec"][-1])]
            LICK_off_temp    = LICK_off_temp[(LICK_off_temp < trace[i]["time_sec"][-1])]
    
            trace[i]["LICK_on"]  = LICK_on_temp
            trace[i]["LICK_off"] = LICK_off_temp[:len(trace[i]["LICK_on"])]
    #        trace[i]["LICK_off"] = LICK_on_temp
    
            # Add the first and last time stamps to make tails on the TTL stream
            try:
                trace[i]["LICK_x"]   = np.append(np.append(trace[i]["time_sec"][0], np.reshape(np.kron([trace[i]["LICK_on"], trace[i]["LICK_off"]],
                                                       np.array([[1], [1]])).T, [1,-1])[0]), trace[i]["time_sec"][-1])
                licks_exist = 1                
                LICK_size            = len(trace[i]["LICK_on"])
                d                    = np.ones(LICK_size)
            
            # Add zeros to beginning and end of 0,1 value array to match len of LICK_x
                trace[i]["LICK_y"]   = np.append(np.append(0,np.reshape(np.vstack([np.zeros(LICK_size),
                                   d, d, np.zeros(LICK_size)]).T, [1, -1])[0]),0)
                
            except: 
                print(mouseID, "has no licks")
                licks_exist = 0
            
            
            if "lickpellet" in BLOCKPATH:
                
                trace[i]["lever_on"]  = activelever_on_temp
                trace[i]["lever_off"] = activelever_off_temp
        
                trace[i]["lever_x"]   = np.append(np.append(trace[i]["time_sec"][0], np.reshape(np.kron([trace[i]["lever_on"], trace[i]["lever_off"]],
                                                           np.array([[1], [1]])).T, [1,-1])[0]), trace[i]["time_sec"][-1])
        
                lever_size            = len(trace[i]["lever_on"])
                d                    = np.ones(lever_size)
                
                # Add zeros to beginning and end of 0,1 value array to match len of LICK_x
                trace[i]["lever_y"]   = np.append(np.append(0,np.reshape(np.vstack([np.zeros(lever_size),
                                       d, d, np.zeros(lever_size)]).T, [1, -1])[0]),0)
        
                    
                
    #%%########################################################
    #       BOUT ANALYSIS 
    ###########################################################      
            if licks_exist == 0:
                continue
                    
            
            else:
            # Define time (s) that separates bouts and minimum number of licks for a bout
            # bout_interval_cutoff = amount of time (sec) between two licks that terminates a bout (I use 1 sec)
            # bout_size_cutoff = minimum # of licks to be considered a bout (I use 10)
                bout_interval_cutoff_sec = 1
                bout_size_cutoff         = 10
            
                bout_ons        = np.asarray([])
                bout_offs       = np.asarray([])
                bout_size       = np.asarray([])
                bout_length     = np.asarray([])
                bout_data       = np.asarray([])   
        
                lickon          = trace[i]["LICK_on"]
                # lickoff         = trace[i]["LICK_off"]
                # lickoff          = trace[i]["LICK_on"]
        #        bout_licktimes  = []
                k = 0 # index for while loop
                # while loop finds bout ons, offs, and licks per bout
                while k < (len(lickon) - bout_size_cutoff):
                    # bout_licktimes_temp  = []
                    licks     = 1
                    while (lickon[k+licks] - lickon[k+licks-1]) < bout_interval_cutoff_sec:
                        licks += 1
                        
                        if (k+licks) >= len(lickon):
                            break
                
                    if licks >= bout_size_cutoff:
                        bout_ons     = np.append(bout_ons, lickon[k])
                        bout_offs    = np.append(bout_offs, lickon[k+licks - 1]) ## bout off is the last on time
                        bout_size    = np.append(bout_size, licks)
                        bout_length  = np.append(bout_length, lickon[k+licks-1] - lickon[k])   
                        bout_data    = np.append(bout_data, 1)
                    
                    k = k + licks
                
                # Make a continuous time series for lick BOUTS for plotting
                try: 
                    LICK_EVENT_x    = np.append(trace[i]["time_sec"][0], np.append(
                                    np.reshape(np.kron([bout_ons , bout_offs], np.array([[1], [1]])).T, [1,-1])[0], trace[i]["time_sec"][-1]))          
                    LICK_EVENT_sz   = len(bout_ons)
                    LICK_EVENT_y = np.append(np.append(
                        0, np.reshape(np.vstack([np.zeros(LICK_EVENT_sz), bout_data, bout_data, np.zeros(LICK_EVENT_sz)]).T, [1 ,-1])[0]), 0)
                    
                    
                    trace[i]["BOUT_on"]     = bout_ons
                    trace[i]["BOUT_off"]    = bout_offs
                    trace[i]["BOUT_x"]      = LICK_EVENT_x
                    trace[i]["BOUT_y"]      = LICK_EVENT_y
                    trace[i]["BOUT_size"]   = bout_size
                    trace[i]["BOUT_length"] = bout_length 
                    
                    bouts_exist    = 1 
                except:     
                    bouts_exist    = 0 
                    print (mouseID, " has no bouts")


            
            
        #%%######################################################################
        #       Entire timecourse: dFF and lick +  bout
        #########################################################################
                #################################################################
                #       FIGURE 1 dFF + licks and bouts 
                #################################################################
               
                if bouts_exist == 1:
                
                    y_scale = 1
                    y_shift = 0
                    
                    fig, ax = plt.subplots(3, 1, figsize = (10, 8))
                    ax1     = plt.subplot2grid((6, 1), (0, 0), rowspan=4, colspan=1)   
                    ax1.plot(trace[i]["time_sec_ds"]/60, trace[i]["zscore_ds"], linewidth=2, color='green', label='GCaMP/UV deltaF/F')
                    ax1.set_ylabel("zscore", fontsize = 20)
                    ax1.set_xlabel("time (min)")
                    # ax1.axes.get_xaxis().set_visible(False)        
                    
                    ax2     = plt.subplot2grid((6, 1), (4, 0), rowspan=1, colspan=1)   
                    try: 
                        ax2.plot(trace[i]["LICK_x"], y_scale*trace[i]["LICK_y"] + y_shift, linewidth=0.8, color='dodgerblue', label='Lick Event')
                    except:
                        
                        ax2.plot(trace[i]["LICK_x"], y_scale*trace[i]["LICK_y"] + y_shift, linewidth=0.8, color='dodgerblue', label='Lick Event')
            
                    ax2.axes.get_xaxis().set_visible(False)    
                    ax2.axis("off")
                    
                    ax3     = plt.subplot2grid((6, 1), (5, 0), rowspan=1, colspan=1)   
                    ax3.plot(LICK_EVENT_x, y_scale*LICK_EVENT_y+y_shift, linewidth=1, color='blue', label='Lick Bout')
                    ax3.axis("off")
                    
                    suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_" + str(trace[i]["age_mo"]) + "mo" + "2"       
                    plt.suptitle("%s" %suptitle, fontsize = 20) 
                    plt.tight_layout(rect=[0, 0.03, 1, 0.95])           
                    save_file_path = analysis_path + "/" + suptitle + ".png"
                    plt.savefig(save_file_path)
                    # save_file_path = analysis_path + "/" + suptitle + ".pdf"
                    # plt.savefig(save_file_path)
                    plt.show()
             
            
            
                    #################################################################
                    #       FIGURE 2  dFF + licks and bouts     
                    #################################################################
                    # Cosmetics for plotting
                    if "WD" in trace[i]["stimulus"]: 
                        y_scale    = -np.min(trace[i]["dFF"])/3
                        y_shift    = np.min(trace[i]["dFF"])/3*4
                    else:
                        y_scale    = np.max(trace[i]["dFF"])/3
                        y_shift    = -np.max(trace[i]["dFF"])  
            
                    fig, ax = plt.subplots(2, 1, figsize = (10, 8))
                    ax1     = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)   
                    lns1    = ax1.plot(trace[i]["time_sec_ds"]/60, trace[i]["zscore_ds"], linewidth=1.5, color='green', label='GCaMP')
                    ax1.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
                    ax1.set_ylabel("zscore", fontsize = 20)
                    ax1.set_xlabel("time (min)")
                
                    ax2     = ax1.twinx()       
                    lns2    = ax2.plot(trace[i]["LICK_x"]/60, y_scale*trace[i]["LICK_y"] + y_shift, linewidth=1, color='dodgerblue', alpha = 0.5, label='Lick')
                    ax2.axes.get_yaxis().set_visible(False)
            
                    lns     = lns1 + lns2
                    labels  = [l.get_label() for l in lns]        
                    ax1.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16)
                    
            
                    ax3     = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)   
                    lns3    = ax3.plot(trace[i]["time_sec_ds"]/60, trace[i]["zscore_ds"], linewidth=1.5, color='green', label='GCaMP')
                    ax3.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
                    # ax3.axvline (x = stim_sec/60, linewidth = 2, linestyle = "--", color = "magenta")
                    ax3.set_ylabel("zscore", fontsize = 18)
                    ax3.set_xlabel("time (min)")
                    
                    ax4     = ax3.twinx()
                    lns4    = ax4.plot(LICK_EVENT_x/60, y_scale*LICK_EVENT_y+y_shift, linewidth=2, color='blue', alpha = 0.5, label='Bout')
                    ax4.axis("off")
            
                    lns     = lns3 + lns4
                    labels  = [l.get_label() for l in lns]        
                    ax3.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16)
            
                    suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_" + str(trace[i]["age_mo"]) + "mo"          
                    plt.suptitle("%s" %suptitle, fontsize = 18) 
                    plt.tight_layout(rect=[0, 0.03, 1, 0.95])           
                    save_file_path = analysis_path + "/" + suptitle + ".png"
                    plt.savefig(save_file_path)
                    # save_file_path = analysis_path + "/" + suptitle + ".pdf"
                    # plt.savefig(save_file_path)
                    plt.show()
            
                    #%%#################################################################
                    #       FIGURE 3  dFF + licks and bouts     
                    #################################################################
                    # Cosmetics for plotting
                    if "WD" in trace[i]["stimulus"]: 
                        y_scale    = -np.min(trace[i]["dFF"])/3
                        y_shift    = np.min(trace[i]["dFF"])/3*4
                    else:
                        y_scale    = np.max(trace[i]["dFF"])/3
                        y_shift    = -np.max(trace[i]["dFF"])  
            
                    fig, ax = plt.subplots(2, 1, figsize = (10, 8))
                    ax1     = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)   
                    lns1    = ax1.plot(trace[i]["time_sec_ds"]/60, trace[i]["zscore_ds"], linewidth=1.5, color=age_color, label='GCaMP')
                    ax1.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
                    ax1.set_ylabel("zscore", fontsize = 20)
                    ax1.set_xlabel("time (min)")
                
                    ax2     = ax1.twinx()       
                    lns2    = ax2.plot(trace[i]["LICK_x"]/60, y_scale*trace[i]["LICK_y"] + y_shift, linewidth=1, color='dodgerblue', alpha = 0.5, label='Lick')
                    ax2.axes.get_yaxis().set_visible(False)
            
                    lns     = lns1 + lns2
                    labels  = [l.get_label() for l in lns]        
                    ax1.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16)
                    
            
                    ax3     = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)   
                    lns3    = ax3.plot(trace[i]["time_sec_ds"]/60, trace[i]["zscore_ds"], linewidth=1.5, color=age_color, label='GCaMP')
                    ax3.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
            #         ax3.axvline (x = stim_sec/60, linewidth = 2, linestyle = "--", color = "magenta")
                    ax3.set_ylabel("zscore", fontsize = 18)
                    ax3.set_xlabel("time (min)")
                    
                    ax4     = ax3.twinx()
                    lns4    = ax4.plot(LICK_EVENT_x/60, y_scale*LICK_EVENT_y+y_shift, linewidth=2, color='blue', alpha = 0.5, label='Bout')
                    ax4.axis("off")
            
                    lns     = lns3 + lns4
                    labels  = [l.get_label() for l in lns]        
                    ax3.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16)
            
                    suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_" + str(trace[i]["age_mo"]) + "mo"          
                    plt.suptitle("%s" %suptitle, fontsize = 18) 
                    plt.tight_layout(rect=[0, 0.03, 1, 0.95])           
                    save_file_path = analysis_path + "/" + suptitle + "3.png"
                    plt.savefig(save_file_path)
                    # save_file_path = analysis_path + "/" + suptitle + "3.pdf"
                    # plt.savefig(save_file_path)
                    plt.show()
            
            #%%##############################################################
            #       PSTH
            #################################################################
                    
                    # mask                    = trace[i]["time_sec"] > 0
            #        time_sec_dims               = np.expand_dims(trace[i]["time_sec"], axis = 1)
            #        tree[i]                     = KDTree(time_sec_dims, leaf_size = 2)  
                    # trace[i]["time_sec_KDtree"] = trace[i]["time_sec_KDtree"].astype(object)
            #        trace[i]["time_sec_KDtree"] = tree[i]     
            
            
                    time_sec_dims             = np.expand_dims(trace[i]["time_sec_ds"], axis = 1)
                    tree[i]                   = KDTree(time_sec_dims, leaf_size = 2)  
            
                    frame_rate      = 1017.2527                                       
                    psth_width      = math.floor(frame_rate*psth_width_sec/ds_factor)
                    psth_num_points = 2*psth_width+1
                    psth_plot_times = np.linspace(-psth_width_sec, psth_width_sec, psth_num_points, endpoint = True, dtype = None)
                            
                    psth            = np.empty([len(trace[i]["BOUT_size"]), psth_num_points])
                    psth[:]         = np.nan
                    psth_y_min_list = []
                    psth_y_max_list = []
                    
                    # print ("total number of bouts are %d" %len(trace[i]["BOUT_on"]))
                    for k, bout_on in enumerate (trace[i]["BOUT_on"]):
                        bout_on     = bout_on.reshape(1, -1)
                        ind         = tree[i].query(bout_on)[1][0][0] 
                        # print ("processing PSTH for bout # %d, bout_on time is %.2f" %(k, bout_on[0][0]))
                        ## Process without error if PSTH post-width exceeds the recording, fill with nan.
                        if ind + psth_width + 1 >= len(trace[i]["zscore_ds"]):  
                            last_post_psth                          = len(trace[i]["zscore_ds"]) - ind
                            psth[k, : psth_width + last_post_psth]  = trace[i]["zscore_ds"][ind - psth_width : ] - np.nanmedian(trace[i]["zscore_ds"][ind - psth_width : ind]) 
                        else:    
                            psth[k, :]  = trace[i]["zscore_ds"][ind - psth_width : ind + psth_width +1] - np.nanmedian(trace[i]["zscore_ds"][ind - psth_width : ind])
            
                        psth_y_min_list.append(np.nanmin(psth[k, :]))
                        psth_y_max_list.append(np.nanmax(psth[k, :]))
                    trace[i]["BOUT_on_PSTH_KDquery"] = psth 
            
            #%%####################################################################################        
            ##      FIGURE 1. PSTH of first lick of all individual bouts in the order 
            #####################################################################################
                    
                    psth_y_min  = np.nanmin(psth_y_min_list)
                    psth_y_max  = np.nanmax(psth_y_max_list)
                    
                    fig, ax     = plt.subplots(len(psth), 1, figsize = (6, 2*len(psth)))       
                    for k in range(len(psth)):  
                        ax      = plt.subplot2grid((len(psth), 1), (k, 0), rowspan=1, colspan=1)       
                        legend  = "bout # %d" %(k+1)
                        if k == 0:
                            plt.plot(psth_plot_times, psth[0], linewidth = 2, color = "red", label = legend)
                            plt.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18, color = "red")            
                        else:
                            plt.plot(psth_plot_times, psth[k], linewidth = 2, color = "black", label = legend)
                            plt.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18)           
                        plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
                        plt.xlim(-psth_width_sec, psth_width_sec)
                        plt.ylim(psth_y_min*1.05, psth_y_max*1.05)   
                        plt.yticks(size = 18)
                        plt.ylabel("zscore", size = 18)
                        if k == len(psth)-1 :
                            plt.xticks(size = 18)
                            plt.xlabel("Time (sec)", size = 18)
                        else: 
                            ax.axes.get_xaxis().set_visible(False)               
                    suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_" + str(trace[i]["age_mo"]) + "mo" + "_PSTH2"                   
                    plt.suptitle("%s\nfirst lick in each bout " %suptitle, size = 18)
                    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
                    save_file_path = analysis_path + "/" + suptitle + ".png"
                    plt.savefig(save_file_path)
                    plt.show()           
                     
            
            
            #####################################################################################        
            ##      FIGURE 2. PSTH of first lick in first bout VS first lick in other bouts
                    #####################################################################################
                    
                    psth_first      = psth[0,]
                    psth_mean       = np.nanmean(psth[1:,], axis =0)
                    psth_stderr     = np.nanstd(psth[1:,], axis =0)
                    if len(trace[i]["BOUT_on"]) >2:
                        psth_stderr     = psth_stderr/math.sqrt(len(trace[i]["BOUT_on"])-2)
                    else:
                        psth_stderr     = psth_stderr
                    
                    df_psth_Z       = pd.DataFrame(psth_plot_times, columns = ["time_sec"])
                    df_psth_Z["firstlick_firstbout"]        = pd.Series(psth_first)
                    df_psth_Z["firstlick_otherbout_mean"]   = pd.Series(psth_mean)
                    df_psth_Z["firstlick_otherbout_stderr"] = pd.Series(psth_stderr)
                    
                    plt.figure(figsize = (10,6))
                    plt.plot(psth_plot_times, psth_first, linewidth = 3, color = "red")
                    plt.errorbar(psth_plot_times, psth_mean, psth_stderr, linewidth = 3, color = "black", ecolor = "lightgray")
                    plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
                    plt.xlim(-psth_width_sec, psth_width_sec)
                    #plt.ylim(-2.1, 2.1)   
                    plt.xticks(size = 18)
                    plt.yticks(size = 18)
                    plt.xlabel("Time (sec)", size = 18)
                    plt.ylabel("zscore", size = 18)
                    suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_" + str(trace[i]["age_mo"]) + "mo" + "_PSTH1"                   
                    plt.title("%s\nfirst lick in all bouts " %suptitle, size = 18)
                    save_file_path = analysis_path + "/" + suptitle + ".png"
                    plt.savefig(save_file_path)
                    plt.show()
            
            
        
        
    
    #%%#################################################################
    #       ZOOM in by bouts 
    ###################################################################
            ## Zoomed-in bouts are saved into one dataframe.
            if bouts_exist == 1:
                df_bout_zs[i] = pd.DataFrame(columns = ["bout_on", "bout_off", "bout_size",
                                                       "bout_number", "bout_on_zs", "bout_off_zs",
                                                       "bout_delzs", "bout_delzs_perlick",
                                                       "bout_plot_times", "bout_plot_zs",
                                                       "bout_lick_x", "bout_lick_y"])
                
                for k, (bout_on, bout_off) in enumerate (zip(trace[i]["BOUT_on"], trace[i]["BOUT_off"])):
                     ## Use the KDtree created for PSTH to find bout on and bout off   
                     bout_on           = bout_on.reshape(1, -1)
                     ind_on            = tree[i].query(bout_on)[1][0][0] 
                     
                     bout_off          = bout_off.reshape(1, -1)
                     ind_off           = tree[i].query(bout_off)[1][0][0]
                
                     df_bout_zs[i].at[k, "bout_on_ind"]  = ind_on
                     df_bout_zs[i].at[k, "bout_off_ind"] = ind_off           
                     df_bout_zs[i].at[k, "bout_on_dFF"]  = trace[i]["zscore_ds"][ind_on]
                     df_bout_zs[i].at[k, "bout_off_dFF"] = trace[i]["zscore_ds"][ind_off]
                     
                     df_bout_zs[i].at[k, "bout_deldFF"]         = trace[i]["zscore_ds"][ind_off] - trace[i]["zscore_ds"][ind_on]
                     df_bout_zs[i].at[k, "bout_deldFF_perlick"] = df_bout_zs[i].at[k, "bout_deldFF"]/trace[i]["BOUT_size"][k]
                
                     df_bout_zs[i].at[k, "bout_number"]  = k+1
                     df_bout_zs[i].at[k, "bout_on"]      = trace[i]["BOUT_on"][k]
                     df_bout_zs[i].at[k, "bout_off"]     = trace[i]["BOUT_off"][k]
                     df_bout_zs[i].at[k, "bout_size"]    = trace[i]["BOUT_size"][k]
                     


                #%%###############################################################
                ## Figure: z score with licks, ZOOMED in for each bout + lever press 
                #################################################################
               
                if "lickpellet" in BLOCKPATH:
                    pre_post_sec    = 3.0   ## 3 sec pre- and post- bout
                    frame_rate      = 1017.2527
                    pre_post_width  = math.floor(pre_post_sec*frame_rate)
                      
                    fig, ax     = plt.subplots(len(df_bout_zs[i]), 1, figsize = (16, 4*len(df_bout_zs[i])))       
                    for k in range(len(df_bout_zs[i])):  
                         ax1     = plt.subplot2grid((len(df_bout_zs[i]), 1), (k, 0), rowspan=1, colspan=1)       
                         legend  = "bout # %d" %(k+1)
                    
                         ind_start             = int(df_bout_zs[i].at[k, "bout_on_ind"] - pre_post_width)
                         ind_end               = int(df_bout_zs[i].at[k, "bout_off_ind"] + pre_post_width+1)  
                         x                     = trace[i]["time_sec_ds_zs"][ind_start:ind_end]
                         y                     = trace[i]["zscore_ds"][ind_start:ind_end]
            
                         lick_x                = trace[i]["LICK_x"]
                         lick_x_masked         = lick_x[lick_x >= x[0]]
                         lick_x_masked_masked  = lick_x_masked[lick_x_masked <= x[-1]]
                         lick_y                = trace[i]["LICK_y"][1:1+len(lick_x_masked_masked)]
            
                         lever_x               = trace[i]["lever_x"]
                         lever_x_masked        = lever_x[lever_x >= x[0]]
                         lever_x_masked_masked = lever_x_masked[lever_x_masked <= x[-1]]
                         lever_y               = trace[i]["lever_y"][1:1+len(lever_x_masked_masked)]
                        
                         df_bout_zs[i].at[k,"bout_plot_times"] = x
                         df_bout_zs[i].at[k,"bout_plot_zs"]   = y
                         df_bout_zs[i].at[k,"bout_lick_x"] = lick_x_masked_masked
                         df_bout_zs[i].at[k,"bout_lick_y"] = lick_y
            
                         if k == 0:
                             lns1    = ax1.plot(x, y, linewidth=2, color = "black", label = legend)
            #                 ax1.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18, color = "red")            
                             ax2     = ax1.twinx() 
                             lns2    = ax2.plot(lick_x_masked_masked, lick_y, linewidth=1, color='dodgerblue', alpha = 0.5, label='lick')
                             lns3    = ax2.plot(lever_x_masked_masked, lever_y, linewidth=1.5, color='red', alpha = 0.8, label='lever')
                             ax2.axes.get_yaxis().set_visible(False)
                         else:
                             lns1    = ax1.plot(x, y, linewidth=2, color = "black", label = legend)
            #                 ax1.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18)            
                    
                             ax2     = ax1.twinx()       
                             lns2    = ax2.plot(lick_x_masked_masked, lick_y, linewidth=1, color='dodgerblue', alpha = 0.5, label='lick')
                             lns3    = ax2.plot(lever_x_masked_masked, lever_y, linewidth=1.5, color='red', alpha = 0.8, label='lever')
                             ax2.axes.get_yaxis().set_visible(False)
                              
                         lns     = lns1 + lns2 + lns3
                         labels  = [l.get_label() for l in lns]        
                         ax1.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16) 
                         ax1.set_ylabel("zscore", size = 18)
                         
                         if k == len(df_bout_zs[i]) -1 :
                             ax1.set_xlabel("Time (sec)", size = 18)
            
                         # plt.axvline(x = df_bout_zs[i]["bout_on"][k], linewidth = 2, linestyle = "--", color = "gray")
                         # plt.axvline(x = df_bout_zs[i]["bout_off"][k], linewidth = 2, linestyle = "--", color = "gray")
                         plt.yticks(size = 18)
                         plt.ylabel("zscore", size = 18)
                         plt.xticks(size = 18)
                         plt.xlabel("Time (sec)", size = 18)
                         
                    suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_" + str(trace[i]["age_mo"]) + "mo" + "_boutzoom"                   
                    plt.suptitle("%s\nfirst lick in each bout " %suptitle, size = 18)
                    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
                    save_file_path = analysis_path + "/" + suptitle + "_lever.png"
                    plt.savefig(save_file_path)
                    plt.show()   
                    
                
                    #%%###############################################################
                    ## Figure: z score with licks, FURTHER ZOOMED in for each bout + lever press 
                    #################################################################
                    pre_post_sec    = 0.1   ## 3 sec pre- and post- bout
                    frame_rate      = 1017.2527
                    pre_post_width  = math.floor(pre_post_sec*frame_rate)
                      
                    fig, ax     = plt.subplots(len(df_bout_zs[i]), 1, figsize = (16, 4*len(df_bout_zs[i])))       
                    for k in range(len(df_bout_zs[i])):  
                         ax1     = plt.subplot2grid((len(df_bout_zs[i]), 1), (k, 0), rowspan=1, colspan=1)       
                         legend  = "bout # %d" %(k+1)
                    
                         ind_start             = int(df_bout_zs[i].at[k, "bout_on_ind"] - pre_post_width)
                         ind_end               = int(df_bout_zs[i].at[k, "bout_off_ind"] + pre_post_width+1)  
                         x = trace[i]["time_sec_ds_zs"][ind_start:ind_end]
                         y = trace[i]["zscore_ds"][ind_start:ind_end]
            
                         lick_x               = trace[i]["LICK_x"]
                         lick_x_masked        = lick_x[lick_x >= x[0]]
                         lick_x_masked_masked = lick_x_masked[lick_x_masked <= x[-1]]
                         lick_y               = trace[i]["LICK_y"][1:1+len(lick_x_masked_masked)]
            
                         lever_x               = trace[i]["lever_x"]
                         lever_x_masked        = lever_x[lever_x >= x[0]]
                         lever_x_masked_masked = lever_x_masked[lever_x_masked <= x[-1]]
                         lever_y               = trace[i]["lever_y"][1:1+len(lever_x_masked_masked)]
                        
                         df_bout_zs[i].at[k,"bout_plot_times"] = x
                         df_bout_zs[i].at[k,"bout_plot_zs"]   = y
                         df_bout_zs[i].at[k,"bout_lick_x"] = lick_x_masked_masked
                         df_bout_zs[i].at[k,"bout_lick_y"] = lick_y
            
                         if k == 0:
                             lns1    = ax1.plot(x, y, linewidth=2, color = "black", label = legend)
            #                 ax1.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18, color = "red")            
                             ax2     = ax1.twinx() 
                             lns2    = ax2.plot(lick_x_masked_masked, lick_y, linewidth=1, color='dodgerblue', alpha = 0.5, label='lick')
                             lns3    = ax2.plot(lever_x_masked_masked, lever_y, linewidth=1.5, color='red', alpha = 0.8, label='lever')
                             ax2.axes.get_yaxis().set_visible(False)
                         else:
                             lns1    = ax1.plot(x, y, linewidth=2, color = "black", label = legend)
            #                 ax1.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18)            
                    
                             ax2     = ax1.twinx()       
                             lns2    = ax2.plot(lick_x_masked_masked, lick_y, linewidth=1, color='dodgerblue', alpha = 0.5, label='lick')
                             lns3    = ax2.plot(lever_x_masked_masked, lever_y, linewidth=1.5, color='red', alpha = 0.8, label='lever')
                             ax2.axes.get_yaxis().set_visible(False)
                              
                         lns     = lns1 + lns2 + lns3
                         labels  = [l.get_label() for l in lns]        
                         ax1.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16) 
                         ax1.set_ylabel("zscore", size = 18)
                         
                         if k == len(df_bout_zs[i]) -1 :
                             ax1.set_xlabel("Time (sec)", size = 18)
            
                         # plt.axvline(x = df_bout_zs[i]["bout_on"][k], linewidth = 2, linestyle = "--", color = "gray")
                         # plt.axvline(x = df_bout_zs[i]["bout_off"][k], linewidth = 2, linestyle = "--", color = "gray")
                         plt.yticks(size = 18)
                         plt.ylabel("zscore", size = 18)
                         plt.xticks(size = 18)
                         plt.xlabel("Time (sec)", size = 18)
                         
                    suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_" + str(trace[i]["age_mo"]) + "mo" + "_boutzoom2"                   
                    plt.suptitle("%s\nfirst lick in each bout " %suptitle, size = 18)
                    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
                    save_file_path = analysis_path + "/" + suptitle + "_leverzoom.png"
                    plt.savefig(save_file_path)
                    plt.show()  




                
                #%%#################################################################
                ## Figure: bar graph of z score change per lick for individual bout      
                ####################################################################
                else:
                    
                    fig, axes     = plt.subplots(nrows = 3, ncols = 1, figsize = (0.5*len(df_bout_zs[i]), 12))       
                    
                    df_bout_zs[i].plot.bar(x = "bout_number", y = "bout_delzs", rot = 0, color = "black", ax = axes[0])
                    axes[0].set_title("Changes in zscore") 
                    
                    #ax2      = plt.subplot2grid((3, 1), (1, 0), rowspan=1, colspan=1)           
                    df_bout_zs[i].plot.bar(x = "bout_number", y = "bout_size", rot = 0, color = "dimgray", ax = axes[1])
                    axes[1].set_title("number of licks in bout")        
                    
                    #ax3      = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1)   
                    df_bout_zs[i].plot.bar(x = "bout_number", y = "bout_delzs_perlick", rot = 0, color = "darkgray", ax = axes[2])
                    axes[2].set_title("zscore changes/lick")
                    
                    suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_" + str(trace[i]["age_mo"]) + "mo" + "_zslick"                   
                    plt.suptitle("%s\nzs change per lick in each bout" %suptitle, size = 18)
                    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
                    save_file_path = analysis_path + "/" + suptitle + ".png"
                    plt.savefig(save_file_path)
                    plt.show()           
                    
                    #%%###############################################################
                    ## Figure: z score with licks, ZOOMED in for each bout          
                    #################################################################
                    pre_post_sec    = 3.0   ## 3 sec pre- and post- bout
                    frame_rate      = 1017.2527
                    pre_post_width  = math.floor(pre_post_sec*frame_rate)
                      
                    fig, ax     = plt.subplots(len(df_bout_zs[i]), 1, figsize = (8, 2*len(df_bout_zs[i])))       
                    for k in range(len(df_bout_zs[i])):  
                         ax1     = plt.subplot2grid((len(df_bout_zs[i]), 1), (k, 0), rowspan=1, colspan=1)       
                         legend  = "bout # %d" %(k+1)
                    
                         ind_start             = int(df_bout_zs[i].at[k, "bout_on_ind"] - pre_post_width)
                         ind_end               = int(df_bout_zs[i].at[k, "bout_off_ind"] + pre_post_width+1)  
                         x = trace[i]["time_sec_ds_zs"][ind_start:ind_end]
                         y = trace[i]["zscore_ds"][ind_start:ind_end]
            
                         lick_x               = trace[i]["LICK_x"]
                         lick_x_masked        = lick_x[lick_x >= x[0]]
                         lick_x_masked_masked = lick_x_masked[lick_x_masked <= x[-1]]
                         lick_y               = trace[i]["LICK_y"][1:1+len(lick_x_masked_masked)]
                        
                         df_bout_zs[i].at[k,"bout_plot_times"] = x
                         df_bout_zs[i].at[k,"bout_plot_zs"]   = y
                         df_bout_zs[i].at[k,"bout_lick_x"] = lick_x_masked_masked
                         df_bout_zs[i].at[k,"bout_lick_y"] = lick_y
            
                         if k == 0:
                             lns1    = ax1.plot(x, y, linewidth=2, color = "red", label = legend)
            #                 ax1.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18, color = "red")            
                             ax2     = ax1.twinx() 
                             lns2    = ax2.plot(lick_x_masked_masked, lick_y, linewidth=1, color='dodgerblue', alpha = 0.5, label='lick')
                             ax2.axes.get_yaxis().set_visible(False)
                         else:
                             lns1    = ax1.plot(x, y, linewidth=2, color = "black", label = legend)
            #                 ax1.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18)            
                    
                             ax2     = ax1.twinx()       
                             lns2    = ax2.plot(lick_x_masked_masked, lick_y, linewidth=1, color='dodgerblue', alpha = 0.5, label='lick')
                             ax2.axes.get_yaxis().set_visible(False)
                              
                         lns     = lns1 + lns2
                         labels  = [l.get_label() for l in lns]        
                         ax1.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16) 
                         ax1.set_ylabel("zscore", size = 18)
                         
                         if k == len(df_bout_zs[i]) -1 :
                             ax1.set_xlabel("Time (sec)", size = 18)
            
                         plt.axvline(x = df_bout_zs[i]["bout_on"][k], linewidth = 2, linestyle = "--", color = "gray")
                         plt.axvline(x = df_bout_zs[i]["bout_off"][k], linewidth = 2, linestyle = "--", color = "gray")
                         plt.yticks(size = 18)
                         plt.ylabel("zscore", size = 18)
                         plt.xticks(size = 18)
                         plt.xlabel("Time (sec)", size = 18)
                         
                    suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_" + str(trace[i]["age_mo"]) + "mo" + "_boutzoom"                   
                    plt.suptitle("%s\nfirst lick in each bout " %suptitle, size = 18)
                    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
                    save_file_path = analysis_path + "/" + suptitle + ".png"
                    plt.savefig(save_file_path)
                    plt.show()    
        
    
    
    	#%%##############################################################
    	#       SAVE data files 
    	#################################################################
     
    		########################################################
    		## Save df_bout_zs dataframe as a .csv file
    		########################################################
            if bouts_exist == 1:

                suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_boutAUC"                   
                save_file_path = analysis_path + "/" + suptitle
                df_bout_zs[i].to_csv(save_file_path + ".csv")
                
        		########################################################
        		## Save df_bout_zs and df_lick_dFF dataframe as a part of trace dictionary 
        		########################################################
                trace[i]["df_bout_zs"]                   = df_bout_zs[i]
        # 		trace[i]["lick_meandFF_dFF_mean"]         = df_lick_dFF[i]["lick_meandFF"]       
                
                
        
    ##########################################################
    #       SAVE trace data as a pkl file
    #################################################################
            suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_" + str(trace[i]["age_mo"]) + "mo"          
            plt.suptitle("%s" %suptitle, fontsize = 18)        
            save_file_path  = analysis_path + "/" + suptitle + "_zs.pkl"
            output = open(save_file_path, 'wb')
            pk.dump(trace[i], output)
            output.close()  
     
            
    
#%%############################################################################################
if __name__ == '__main__':
    print("Enter the full path for the directory with .pkl files for group analysis: ")
    path        = input()
    # path = r"D:\Data_Dropbox_Sync\Photometry_DROPBOX\_SFO_GCaMP\TTL_WD_water\old_female"

    main (path)
    
    # ##################################################
    # #   POP open the analysis folder
    # ##################################################
    os.startfile(path)


