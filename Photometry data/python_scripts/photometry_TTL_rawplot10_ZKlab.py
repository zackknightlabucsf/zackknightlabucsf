# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 17:04:00 2020

Processing order for TTL lick fiber photometry data
1. photometry_TTL_raw_PSTH.py (or photometry_TTL_raw_PSTH_ds.py for downsampling dFF)
                      ********** this script
    1-1. photometry_TTL_groupplot_PSTH.py 
2. photometry_TTL_finddecay.py  
3. photometry_TTL_groupplot.py


Input: 
    1) individual directory with "raw" fiber photometry datafile from synapse

    *** NO NEED TO PROCESS with photometry_rawplot or photometry_mouseplot

    *** This script NEEDS at least 1 TTL lick input 

    2) DOB metadata file and TTL metadata file with stimulus and water exposure time, etc. 
            
Output: 
    1) .pkl file with: mouseID, age, dFF, date, stimulus, lick times and bout times
    2) plot with dFF trace and individual lick events and bout events
    3) PSTH plot for first lick in first bout vs first lick in other bouts
    
    
@author: heeun
""" 

# Jupyter magic
#%matplotlib inline

import os
import sys
import pickle as pk
import numpy as np
import pandas as pd
import csv
# from tkinter.filedialog import askdirectory 
from datetime import date
import math
from scipy import signal
import scipy.stats as ss
from scipy.integrate import simps
from sklearn.neighbors import KDTree
import matplotlib.pyplot as plt  # standard Python plotting library
import matplotlib
matplotlib.rcParams['font.size'] = 18 # set font size for all figures

# import the read_block function from the tdt package
from tdt import read_block #, StructType
# import subprocess


##  Import TTL data file into a dataframe

chow_to_TTLwater = 0 # To use chow as a pre-stimulus, water open as a stimulus and TTL water lick is post-stimulus  
WD_to_TTLwater   = 1

if chow_to_TTLwater == 1 or WD_to_TTLwater == 1:
    ##  1) If analyzing chow (pre-stim) - TTL water (post-stim) only: 
    TTLpath = r"C:\Users\heeun\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Mouse_metadata\photometrymiceTTL.csv"
    # TTLpath  = r"F:\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Mouse_metadata\photometrymiceTTL.csv"

else:    ## usually FD to chow to TTL water (entire time course)
    ##  2) If analyzing the entire timecourse from baseline - chow - TTL water:
    # TTLpath = r"C:\Users\heeun\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Mouse_metadata\photometrymiceTTL_wchow.csv"
    TTLpath  = r"F:\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Mouse_metadata\photometrymiceTTL_wchow.csv"

TTLdata = pd.read_csv(TTLpath)    

#%% Importing the Data
# BLOCKPATH = r"C:\Users\heeun\Dropbox\Garrison Lab, Buck Institute\Garrison Lab - Heeun\Data_Dropbox_Sync\Fiber_Photometry_new_ZK\Fiber_Photometry_Setup2_2020\20200128\W1764f22_W1764f01-200128-135727"
print("Enter the full path for directory with raw photometry files: ")
BLOCKPATH      = input()
# BLOCKPATH      = r"C:\Users\heeun\Dropbox\Garrison Lab\Heeun\Data_Dropbox_Sync\Photometry_DROPBOX\Photometry_TTL_water\PVH_AVP_GCaMP\Photometry_TTL_WD_water\AVP137f00_AVP137f01-210520-144811"

setup = 2
if "setup3" in BLOCKPATH:
    setup = 3

fluid_open_min   = 0

Z_pre_post_min   = [10, 45] ## put the stimuluss application/drug injection time as [pre_mm, post_mm]
K_pre_post_min   = [10, 45] ## 10, 55 for injection + water exposure (inj, 10min, water open, 45min)
                          ## 10, 45 for WD
Z_F0_window_min  = [5, 10]  ## setting F0 window as a median value of 5-10 mins (equals -5 t0 0 pre-stimulus)
K_F0_window_min  = [5, 10]

sm_factor        = 61
# sm_factor        = 3001     ## usually 1000 Hz (fps) so this is smoothing over 3000 frames = 3 sec window 
# sm_factor        = 30001      ## usually 1000 Hz (fps) so this is smoothing over 30000 frames = 30 sec window 
ds_factor        = 50
ds_factor1       = 5
ds_factor2       = 5
ds_factor3       = 2
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
data  = read_block(BLOCKPATH)
trace = {}   
tree  = {} 
df_bout_dFF = {} 
df_lick_dFF = {}
blockname = data.info.blockname

## for mouse in Z
mouseZ           = blockname.split('_')[0]
mouseK           = (blockname.split('_')[1]).split('-')[0]
mouse_list       = [mouseZ, mouseK]
# mouse_list       = [mouseK]

# mouseZ = "none"

for i, mouseID in enumerate (mouse_list):
    trace[i]                   = {}
    # if i ==1:
    #     continue
    if mouseID == "none":
        if i == 0:  
            print("No mouse in Z chamber")
        else:
            print("No mouse in K chamber")
        continue
        
    else:
        trace[i]["date"]           = "20" + (blockname.split('_')[1]).split('-')[1]    ## saves as a string 20100408
        trace[i]["time"]           = data.info.utc_start_time  ## saves hh:mm:ss "10:04:10" as a string
        trace[i]["mouseID"]        = mouseID
        
        for j in range(len(TTLdata)):
            if (str(int(TTLdata.at[j, "ExpDate"])) == str(trace[i]["date"])) & (TTLdata.at[j, "mouseID"] == trace[i]["mouseID"]):
                trace[i]["BW"]          = TTLdata.at[j, "BW"]
                trace[i]["stimulus"]    = TTLdata.at[j, "Stimulus"]
                trace[i]["LICK_fluid"]  = TTLdata.at[j, "LickAccess"] 
                stim_sec                = int(TTLdata.at[j, "TimeZero_mm"])*60 + int(TTLdata.at[j, "TimeZero_ss"])
                print ("mouseID is %s, Expdate is %s, stimulus is %s, access is %s" %(trace[i]["mouseID"], trace[i]["date"], trace[i]["stimulus"], trace[i]["LICK_fluid"]))
                break
        
        try:
            trace[i]["stimulus"]
        except:
            print("TTL data does not match")
        
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
        trace[i]["Fratio"]   = trace[i]["GCaMP"]/trace[i]["UV"]
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
        
        suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_noTTL"       

        plt.suptitle("%s" %suptitle) 
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])    

        analysis_path  = BLOCKPATH + "/analysis_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"]
       
        if not os.path.exists(analysis_path):
            os.makedirs(analysis_path)

        save_file_path = analysis_path + "/" + suptitle + ".png"
        plt.savefig(save_file_path)
        save_file_path = analysis_path + "/" + suptitle + ".pdf"
        plt.savefig(save_file_path)
        plt.show()
            
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
        # trace[i]["LICK_off"] = LICK_off_temp
        trace[i]["LICK_off"] = LICK_on_temp

    
            
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
            lickoff         = trace[i]["LICK_off"]
            # lickoff          = trace[i]["LICK_on"]
    #        bout_licktimes  = []
            k = 0 # index for while loop
            # while loop finds bout ons, offs, and licks per bout
            while k < (len(lickon) - bout_size_cutoff):
                bout_licktimes_temp  = []
                licks     = 1
                while (lickon[k+licks] - lickon[k+licks-1]) < bout_interval_cutoff_sec:
    #                bout_licktimes_temp.append(lickon[k+licks])
    #                print ("bout_licktimes of k %d, licks %d is %s" %(k, licks, bout_licktimes_temp))
                    licks += 1
                    
                    if (k+licks) >= len(lickon):
    #                    bout_licktimes.append(bout_licktimes_temp)
    #                    print ("end of bout, bout_licktimes is %s" %bout_licktimes)
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
            y_scale = 1
            y_shift = 0
            
            fig, ax = plt.subplots(3, 1, figsize = (10, 8))
            ax1     = plt.subplot2grid((6, 1), (0, 0), rowspan=4, colspan=1)   
            ax1.plot(trace[i]["time_sec"]/60, trace[i]["dFF"], linewidth=2, color='green', label='GCaMP/UV deltaF/F')
            ax1.set_ylabel("dFF (%)", fontsize = 20)
            ax1.set_xlabel("time (min)")
            # ax1.axes.get_xaxis().set_visible(False)        
            
            ax2     = plt.subplot2grid((6, 1), (4, 0), rowspan=1, colspan=1)   
            ax2.plot(trace[i]["LICK_x"], y_scale*trace[i]["LICK_y"] + y_shift, linewidth=0.8, color='dodgerblue', label='Lick Event')
            ax2.axes.get_xaxis().set_visible(False)    
            ax2.axis("off")
            
            ax3     = plt.subplot2grid((6, 1), (5, 0), rowspan=1, colspan=1)   
            if bouts_exist == 0:
                ax3.plot(trace[i]["LICK_x"], y_scale*trace[i]["LICK_y"] + y_shift, linewidth=0.8, color='dodgerblue', label='Lick Event')
    
            else:
                ax3.plot(LICK_EVENT_x, y_scale*LICK_EVENT_y+y_shift, linewidth=1, color='blue', label='Lick Bout')
            ax3.axis("off")
            
            suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_2"       
            plt.suptitle("%s" %suptitle, fontsize = 20) 
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])           
            save_file_path = analysis_path + "/" + suptitle + ".png"
            plt.savefig(save_file_path)
            save_file_path = analysis_path + "/" + suptitle + ".pdf"
            plt.savefig(save_file_path)
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
            lns1    = ax1.plot(trace[i]["time_sec_ds"]/60, trace[i]["dFF_ds_sm"], linewidth=2.5, color='black', alpha = 0.8, label='GCaMP')
            ax1.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
            ax1.set_ylabel("dFF (%)", fontsize = 20)
            # ax1.set_ylim(-60, 20)
            # ax1.set_xlim(-5, 30)
            ax1.set_xlabel("time (min)")
            
    
            ax2     = ax1.twinx()       
            lns2    = ax2.plot(trace[i]["LICK_x"]/60, y_scale*trace[i]["LICK_y"] + y_shift, linewidth=0.8, color='dodgerblue', alpha = 0.5, label='Lick')
            ax2.axes.get_yaxis().set_visible(False)
    
            lns     = lns1 + lns2
            labels  = [l.get_label() for l in lns]        
            ax1.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16)
            
    
            ax3     = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)   
            lns3    = ax3.plot(trace[i]["time_sec_ds"]/60, trace[i]["dFF_ds_sm"], linewidth=1.5, color='green', label='GCaMP')
            ax3.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
            # ax3.axvline (x = stim_sec/60, linewidth = 2, linestyle = "--", color = "magenta")
            ax3.set_ylabel("dFF (%)", fontsize = 18)
            ax3.set_xlabel("time (min)")
    
            if bouts_exist == 1:        
                ax4     = ax3.twinx()
                lns4    = ax4.plot(LICK_EVENT_x/60, y_scale*LICK_EVENT_y+y_shift, linewidth=2, color='blue', alpha = 0.5, label='Bout')

            else:
                lns4    = ax4.plot(trace[i]["LICK_x"]/60, y_scale*trace[i]["LICK_y"] + y_shift, linewidth=0.8, color='dodgerblue', alpha = 0.5, label='Lick')
        
            ax4.axis("off")
            lns     = lns3 + lns4 
            labels  = [l.get_label() for l in lns]        
            ax3.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16)
    
            suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"]         
            plt.suptitle("%s" %suptitle, fontsize = 18) 
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])           
            save_file_path = analysis_path + "/" + suptitle + ".png"
            plt.savefig(save_file_path)
            save_file_path = analysis_path + "/" + suptitle + ".pdf"
            plt.savefig(save_file_path)
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
            psth_width_sec  = 120                                         
            psth_width      = math.floor(frame_rate*psth_width_sec/ds_factor)
            psth_num_points = 2*psth_width+1
            psth_plot_times = np.linspace(-psth_width_sec, psth_width_sec, psth_num_points, endpoint = True, dtype = None)
                    
            psth            = np.empty([len(trace[i]["BOUT_size"]), psth_num_points])
            psth[:]         = np.nan
            psth_y_min_list = []
            psth_y_max_list = []
            
            print ("total number of bouts are %d" %len(trace[i]["BOUT_on"]))
            for k, bout_on in enumerate (trace[i]["BOUT_on"]):
                bout_on     = bout_on.reshape(1, -1)
                ind         = tree[i].query(bout_on)[1][0][0] 
                print ("processing PSTH for bout # %d, bout_on time is %.2f" %(k, bout_on[0][0]))
                ## Process without error if PSTH post-width exceeds the recording, fill with nan.
                if ind + psth_width + 1 >= len(trace[i]["dFF_ds_sm"]):  
                    last_post_psth                          = len(trace[i]["dFF_ds_sm"]) - ind
                    psth[k, : psth_width + last_post_psth]  = trace[i]["dFF_ds_sm"][ind - psth_width : ] - np.nanmedian(trace[i]["dFF_ds_sm"][ind - psth_width : ind]) 
                else:    
                    psth[k, :]  = trace[i]["dFF_ds_sm"][ind - psth_width : ind + psth_width +1] - np.nanmedian(trace[i]["dFF_ds_sm"][ind - psth_width : ind])
    
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
                plt.ylabel("dFF (%)", size = 18)
                if k == len(psth)-1 :
                    plt.xticks(size = 18)
                    plt.xlabel("Time (sec)", size = 18)
                else: 
                    ax.axes.get_xaxis().set_visible(False)               
            suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_PSTH2"                   
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
            plt.ylabel("dFF (%)", size = 18)
            suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_PSTH1"                   
            plt.title("%s\nfirst lick in all bouts " %suptitle, size = 18)
            save_file_path = analysis_path + "/" + suptitle + ".png"
            plt.savefig(save_file_path)
            plt.show()
    
    
    
    
    
    #%%###################################################################################################
    #       ZOOM in by bouts UPDATED for ZACK to INCLUDE AUC during individual BOUTs 5/24/2021 STARTS
    #####################################################################################################
    		## Zoomed-in bouts are saved into one dataframe
            if bouts_exist == 1:
                df_bout_dFF[i] = pd.DataFrame(columns = ["bout_on", "bout_off", "bout_size",
        											   "bout_number", "bout_on_dFF", "bout_off_dFF",
        											   "bout_deldFF", "bout_deldFF_perlick",
        											   "bout_plot_times", "bout_plot_dFF",
        												   "bout_lick_x", "bout_lick_y",
        												  
        												   "bout_AUC", "bout_AUCperlick", "bout_meandFF",
                                                           "interbout_number",
        												   "interbout_AUC", "interbout_licks", "interbout_AUCperlick",
        												   "interbout_meandFF"])
           
                df_lick_dFF[i] = pd.DataFrame(columns = ["lick_meandFF"])
        												
        			
                for k, (bout_on, bout_off) in enumerate (zip(trace[i]["BOUT_on"], trace[i]["BOUT_off"])):
        				 ## Use the KDtree created for PSTH to find bout on and bout off   
                    bout_on           = bout_on.reshape(1, -1)
                    ind_on            = tree[i].query(bout_on)[1][0][0] 
        				 
                    bout_off          = bout_off.reshape(1, -1)
                    ind_off           = tree[i].query(bout_off)[1][0][0]
        			
                    df_bout_dFF[i].at[k, "bout_on_ind"]  = ind_on
                    df_bout_dFF[i].at[k, "bout_off_ind"] = ind_off           
                    df_bout_dFF[i].at[k, "bout_on_dFF"]  = trace[i]["dFF_ds_sm"][ind_on]
                    df_bout_dFF[i].at[k, "bout_off_dFF"] = trace[i]["dFF_ds_sm"][ind_off]
        				 
                    df_bout_dFF[i].at[k, "bout_deldFF"]         = trace[i]["dFF_ds_sm"][ind_off] - trace[i]["dFF_ds_sm"][ind_on]
                    df_bout_dFF[i].at[k, "bout_deldFF_perlick"] = df_bout_dFF[i].at[k, "bout_deldFF"]/trace[i]["BOUT_size"][k]
        			
                    df_bout_dFF[i].at[k, "bout_number"]       = k+1
                    df_bout_dFF[i].at[k, "bout_on"]           = trace[i]["BOUT_on"][k]
                    df_bout_dFF[i].at[k, "bout_off"]          = trace[i]["BOUT_off"][k]
                    df_bout_dFF[i].at[k, "bout_size"]         = trace[i]["BOUT_size"][k]
        				 
                    df_bout_dFF[i].at[k, "bout_AUC"]        = simps(trace[i]["dFF_ds_sm"][ind_on:ind_off+1])
                    df_bout_dFF[i].at[k, "bout_AUCperlick"] = df_bout_dFF[i].at[k, "bout_AUC"]/df_bout_dFF[i].at[k, "bout_size"]
        				 
                    df_bout_dFF[i].at[k, "bout_meandFF"]     = np.mean(trace[i]["dFF_ds_sm"][ind_on:ind_off+1])
        	 
        																		
        			## Process data for interbout dFF
                    for k in range(len(df_bout_dFF[i])-1):       
                        df_bout_dFF[i].at[k, "interbout_number"]     = int(k+1)
                        df_bout_dFF[i].at[k, "interbout_size"]       = trace[i]["interbout_lick_size"][k]
        				 
                        interbout_on_ind                             = int(df_bout_dFF[i].at[k, "bout_off_ind"])
                        interbout_off_ind                            = int(df_bout_dFF[i].at[k+1, "bout_on_ind"])
                        df_bout_dFF[i].at[k, "interbout_AUC"]        = simps(trace[i]["dFF_ds_sm"][interbout_on_ind:interbout_off_ind])
                        df_bout_dFF[i].at[k, "interbout_AUCperlick"] = df_bout_dFF[i].at[k, "interbout_AUC"]/df_bout_dFF[i].at[k, "interbout_size"]
                        df_bout_dFF[i].at[k, "interbout_meandFF"]    = np.mean(trace[i]["dFF_ds_sm"][interbout_on_ind:interbout_off_ind+1])
        
        			# Get AUC and mean dFF for individual lick events 
                    for k, lick_on in enumerate (trace[i]["LICK_on"]):
                        lick_on                                = lick_on.reshape(1, -1)       				
                        ind_lick_on                            = tree[i].query(lick_on)[1][0][0]
                        df_lick_dFF[i].at[k, "lick_meandFF"]   = trace[i]["dFF_ds_sm"][ind_lick_on]
        				
                        
                    df_bout_dFF[i]["lick_meandFF"]  = df_lick_dFF[i]["lick_meandFF"]
        
    
    
    
        			####################################################################
        			## Additional Figure 1: bar graph of AUC and AUC per lick for individual bout      
        			####################################################################
                    fig, axes     = plt.subplots(nrows = 4, ncols = 1, figsize = (0.5*len(df_bout_dFF[i]), 4*4)) ## 4*of rows
        			
                    df_bout_dFF[i].plot.bar(x = "bout_number", y = "bout_AUC", legend = None, rot = 0, color = "black", ax = axes[0])
                    axes[0].set_title("dFF AUC") 
        			
        			#ax2      = plt.subplot2grid((3, 1), (1, 0), rowspan=1, colspan=1)           
                    df_bout_dFF[i].plot.bar(x = "bout_number", y = "bout_size", legend = None, rot = 0, color = "dimgray", ax = axes[1])
                    axes[1].set_title("number of licks in bout")        
        			
        			#ax3      = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1)   
                    df_bout_dFF[i].plot.bar(x = "bout_number", y = "bout_AUCperlick", legend = None, rot = 0, color = "darkgray", ax = axes[2])
                    axes[2].set_title("dFF AUC per lick")
        	   
        			#ax3      = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1)   
                    df_bout_dFF[i].plot.bar(x = "bout_number", y = "bout_meandFF", legend = None, rot = 0, color = "black", ax = axes[3])
                    axes[3].set_title("dFF mean")
        
        		 
                    suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_boutAUC"                   
                    plt.suptitle("%s\ndFF_AUC_ind_bout" %suptitle, size = 18)
                    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
                    save_file_path = analysis_path + "/" + suptitle + ".png"
                    plt.savefig(save_file_path)
                    save_file_path = analysis_path + "/" + suptitle + ".pdf"
                    plt.savefig(save_file_path)
                    plt.show()    
    		 
        			####################################################################
        			## Additional Figure 2: bar graph of AUC and AUC per lick for individual interbout interval
        			####################################################################
                    if len(trace[i]["BOUT_on"]) > 1:
                        fig, axes     = plt.subplots(nrows = 4, ncols = 1, figsize = (0.5*len(df_bout_dFF[i]), 4*4)) ## 4*of rows
        			
                        df_bout_dFF[i].plot.bar(x = "interbout_number", y = "interbout_AUC", legend = None, rot = 0, color = "black", ax = axes[0])
                        axes[0].set_title("dFF AUC") 
        				
        				#ax2      = plt.subplot2grid((3, 1), (1, 0), rowspan=1, colspan=1)           
                        df_bout_dFF[i].plot.bar(x = "interbout_number", y = "interbout_size", legend = None, rot = 0, color = "dimgray", ax = axes[1])
                        axes[1].set_title("number of licks during interbout")        
        				
        				#ax3      = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1)   
                        df_bout_dFF[i].plot.bar(x = "interbout_number", y = "interbout_AUCperlick", legend = None, rot = 0, color = "darkgray", ax = axes[2])
                        axes[2].set_title("dFF AUC per lick")
        	   
        				#ax3      = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1)   
                        df_bout_dFF[i].plot.bar(x = "interbout_number", y = "interbout_meandFF", legend = None, rot = 0, color = "black", ax = axes[3])
                        axes[3].set_title("dFF mean")
        
        		 
                        suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_interboutAUC"                   
                        plt.suptitle("%s\ndFF_AUC_ind_interbout" %suptitle, size = 18)
                        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
                        save_file_path = analysis_path + "/" + suptitle + ".png"
                        plt.savefig(save_file_path)
                        save_file_path = analysis_path + "/" + suptitle + ".pdf"
                        plt.savefig(save_file_path)
                        plt.show()    
        			
                    
    		   
            			# ###################################################################################
            			# ## Additional Figure 3: bar graph of (i) AUC dFF during lick vs during not lick    
                        #                                      (ii) mean dFF during lick vs during not lick   
            			# ###################################################################################
                        fig, ax       = plt.subplots(nrows = 1, ncols = 2, figsize = (8, 6)) ## 4*of rows
            # 			fig, axes     = plt.subplots(figsize = (3, 6)) ## 4*of rows
                        bar_width     = 0.75
            
                    
            
                        ## 1. mean dFF during all licks vs mean dFF during interbout interval 
                        values        = [df_lick_dFF[i]["lick_meandFF"].mean(), df_bout_dFF[i]["interbout_meandFF"].mean()]
                        labels        = ['all licks', 'no lick']
                        errorbars     = [ss.sem(df_lick_dFF[i]["lick_meandFF"]), ss.sem(df_bout_dFF[i]["interbout_meandFF"][:-1])]      
                        ax1           = plt.subplot2grid((1, 2), (0, 0), rowspan = 1, colspan = 1)
                        ax1.bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)
                        #    ax7   = plt.subplot2grid((4, 4), (0, k*2), rowspan=2, colspan=2)  
            
                        ax1.set_ylabel("dFF(%)")
                        ax1.set_title("mean dF/F0(%)")
            			
                        ## 2. mean dFF during all licks that are part of bouts vs mean dFF during interbout interval 
                        values        = [df_bout_dFF[i]["bout_meandFF"].mean(), df_bout_dFF[i]["interbout_meandFF"].mean()]
                        labels        = ['all bouts', 'interbout']
                        errorbars     = [ss.sem(df_bout_dFF[i]["bout_meandFF"]), ss.sem(df_bout_dFF[i]["interbout_meandFF"][:-1])]   
                        ax2           = plt.subplot2grid((1, 2), (0, 1), rowspan = 1, colspan = 1)            
                        ax2.bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)
            
                        ax2.set_ylabel("dFF(%)")
                        ax2.set_title("mean dF/F0(%)")
            
                        suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"]  
                        plt.suptitle("%s" %suptitle, size = 18)
                        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
                        save_file_path = analysis_path + "/" + suptitle + "_meandFF_lick" + ".png"
                        plt.savefig(save_file_path)
                        save_file_path = analysis_path + "/" + suptitle + "_meandFF_lick" + ".pdf"
                        plt.savefig(save_file_path)
                        plt.show()  
            			
                        ####################################################################
                        ## Figure 4: bar graph of dFF change per lick for individual bout      
                        ####################################################################
                        fig, axes     = plt.subplots(nrows = 3, ncols = 1, figsize = (0.5*len(df_bout_dFF[i]), 12))       
                        
                        df_bout_dFF[i].plot.bar(x = "bout_number", y = "bout_deldFF", rot = 0, color = "black", ax = axes[0])
                        axes[0].set_title("Changes in dFF(%)") 
                        
                        #ax2      = plt.subplot2grid((3, 1), (1, 0), rowspan=1, colspan=1)           
                        df_bout_dFF[i].plot.bar(x = "bout_number", y = "bout_size", rot = 0, color = "dimgray", ax = axes[1])
                        axes[1].set_title("number of licks in bout")        
                        
                        #ax3      = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1)   
                        df_bout_dFF[i].plot.bar(x = "bout_number", y = "bout_deldFF_perlick", rot = 0, color = "darkgray", ax = axes[2])
                        axes[2].set_title("dFF(%) changes/lick")
                        
                        suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_dFFlick"                   
                        plt.suptitle("%s\ndFF change per lick in each bout" %suptitle, size = 18)
                        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
                        save_file_path = analysis_path + "/" + suptitle + ".png"
                        plt.savefig(save_file_path)
                        plt.show()   
        
                         #%%###############################################################
                        ## Figure 5: dFF with licks, ZOOMED in for each bout          
                        #################################################################
                        pre_post_sec    = 3.0   ## 3 sec pre- and post- bout
                        frame_rate      = 1017.2527
                        pre_post_width  = math.floor(pre_post_sec*frame_rate)
                          
                        fig, ax     = plt.subplots(len(df_bout_dFF[i]), 1, figsize = (8, 2*len(df_bout_dFF[i])))       
                        for k in range(len(df_bout_dFF[i])):  
                             ax1     = plt.subplot2grid((len(df_bout_dFF[i]), 1), (k, 0), rowspan=1, colspan=1)       
                             legend  = "bout # %d" %(k+1)
                        
                             ind_start             = int(df_bout_dFF[i].at[k, "bout_on_ind"] - pre_post_width)
                             ind_end               = int(df_bout_dFF[i].at[k, "bout_off_ind"] + pre_post_width+1)  
                             x = trace[i]["time_sec_ds"][ind_start:ind_end]
                             y = trace[i]["dFF_ds_sm"][ind_start:ind_end]
                
                             lick_x               = trace[i]["LICK_x"]
                             lick_x_masked        = lick_x[lick_x >= x[0]]
                             lick_x_masked_masked = lick_x_masked[lick_x_masked <= x[-1]]
                             lick_y               = trace[i]["LICK_y"][1:1+len(lick_x_masked_masked)]
                            
                             df_bout_dFF[i].at[k,"bout_plot_times"] = x
                             df_bout_dFF[i].at[k,"bout_plot_dFF"]   = y
                             df_bout_dFF[i].at[k,"bout_lick_x"] = lick_x_masked_masked
                             df_bout_dFF[i].at[k,"bout_lick_y"] = lick_y
                
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
                             ax1.set_ylabel("dFF (%)", size = 18)
                             
                             if k == len(df_bout_dFF[i]) -1 :
                                 ax1.set_xlabel("Time (sec)", size = 18)
                
                             plt.axvline(x = df_bout_dFF[i]["bout_on"][k], linewidth = 2, linestyle = "--", color = "gray")
                             plt.axvline(x = df_bout_dFF[i]["bout_off"][k], linewidth = 2, linestyle = "--", color = "gray")
                             plt.yticks(size = 18)
                             plt.ylabel("dFF (%)", size = 18)
                             plt.xticks(size = 18)
                             plt.xlabel("Time (sec)", size = 18)
                             
                        suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_boutzoom"                   
                        plt.suptitle("%s\nfirst lick in each bout " %suptitle, size = 18)
                        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
                        save_file_path = analysis_path + "/" + suptitle + ".png"
                        plt.savefig(save_file_path)
                        plt.show()    
                    
     	#%%##############################################################
     	#       SAVE data files 
     	#################################################################
     
    		########################################################
    		## Save df_bout_dFF dataframe as a .csv file
    		########################################################
            suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"] + "_boutAUC"                   
            save_file_path = analysis_path + "/" + suptitle
            df_bout_dFF[i].to_csv(save_file_path + ".csv")
            
    		########################################################
    		## Save df_bout_dFF and df_lick_dFF dataframe as a part of trace dictionary 
    		########################################################
            trace[i]["df_bout_dFF"]                   = df_bout_dFF[i]
    # 		trace[i]["lick_meandFF_dFF_mean"]         = df_lick_dFF[i]["lick_meandFF"]       
    
        
    ##########################################################
    #       SAVE trace data as a pkl file
    #################################################################
            suptitle = trace[i]["mouseID"] + "_" + trace[i]["date"] + "_" + trace[i]["stimulus"] + "_" + trace[i]["LICK_fluid"]        
            plt.suptitle("%s" %suptitle, fontsize = 18)        
            save_file_path  = analysis_path + "/" + suptitle + ".pkl"
            output = open(save_file_path, 'wb')
            pk.dump(trace[i], output)
            output.close()  
     
            

   
# ##################################################
# #   POP open the analysis folder
# ##################################################
os.startfile(analysis_path)


