# -*- coding: utf-8 -*-
"""
Created on Mon Jan 6 20:20:16 2020

Raw fiber photometry processing and visualization.

Reads two-chamber (Z/K) TDT recordings, aligns GCaMP and UV signals
to stimulus onset, and calculates dF/F0 and z-scores using direct
GCaMP/UV ratios or UV-based linear fitting.

Downsamples and smooths traces, generates individual recording plots,
and exports processed data as .pkl files.

"""

#import scipy.io as sio
#from scipy.signal import decimate
#import pandas as pd
import os
#import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle as pk
#import tdt
from datetime import datetime
from tdt import read_block
from scipy import signal
import math
matplotlib.rcParams['font.size'] = 24
#
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
    

       
print("Enter the full path for the directory with raw photometry data: ")
path             = input()    


Z_stimulus       = "2MNaCl"
K_stimulus       = "2MNaCl"

Z_stim_time      = [14, 35]  ## put the stimulus application/drug injection time as [mm, ss]
K_stim_time      = [15, 10]

Z_pre_post_min   = [10, 45] ## put the stimulus application/drug injection time as [pre_mm, post_mm]
K_pre_post_min   = [10, 45]

 
Z_F0_window_min  = [5, 10]  ## setting F0 window as a median value of 5-10 mins (equals -5 t0 0 pre-stimulus)
K_F0_window_min  = [5, 10]
F0_window_min    = [5, 10]
   
ds_factor1            = 5
ds_factor2            = 5
ds_factor3            = 5   
ds_factor             = int(ds_factor1*ds_factor2*ds_factor3)
sm_factor             = 101
        
baseline              = 10  # in minutes 
baseline_index        = np.arange(0, baseline*60*ds_factor, 1)   ## make an array of indices for baseline 
## Read data
data  = read_block(path)
trace = {}
       
blockname = data.info.blockname

## for mouse in Z
mouseZ           = blockname.split('_')[0]
mouseK           = (blockname.split('_')[1]).split('-')[0]

mouse_list       = [mouseZ, mouseK]


today = datetime.today().strftime('%Y%m%d')
analysis_path = path + "/analysis_"  + today + "_zs"
if not os.path.exists(analysis_path):
    os.makedirs(analysis_path)


for i, mouseID in enumerate (mouse_list):
    # if i == 0:
    #     continue
    if mouseID == "none":
        if i == 0:  
            print("No mouse in Z chamber")
        else:
            print("No mouse in K chamber")
        continue
        
    else:
        trace["date"]           = "20" + (blockname.split('_')[1]).split('-')[1]    ## saves as a string 20100408
        trace["time"]           = data.info.utc_start_time  ## saves hh:mm:ss "10:04:10" as a string


        if i == 0:
            trace["mouseID"]    = mouseZ 
            trace["stimulus"]   = Z_stimulus
            # trace["BW"]         = ZBW
            # trace["age"]        = Zage
            sampling_rate       = data.streams["G__Z"].fs ## in Hz
            raw_GCaMP           = data.streams["G__Z"].data
            raw_UV              = data.streams["U__Z"].data
            pre_sec             = int(Z_pre_post_min[0])*60
            post_sec            = int(Z_pre_post_min[1])*60
            stim_sec            = int(Z_stim_time[0])*60 + int(Z_stim_time[1]) 

        else:
            trace["mouseID"]    = mouseK
            trace["stimulus"]   = K_stimulus
            # trace["BW"]         = KBW
            # trace["age"]        = Kage
            sampling_rate       = data.streams["G__K"].fs ## in Hz
            raw_GCaMP           = data.streams["G__K"].data
            raw_UV              = data.streams["U__K"].data
            pre_sec             = int(K_pre_post_min[0])*60
            post_sec            = int(K_pre_post_min[1])*60
            stim_sec            = int(K_stim_time[0])*60 + int(K_stim_time[1])             
            
    ## Crop the trace around the stimulus time
    crop_start_sec   = stim_sec - pre_sec
    crop_end_sec     = stim_sec + post_sec
    
    crop_start_frame = int(np.ceil(crop_start_sec*sampling_rate))
    crop_end_frame   = int(np.floor(crop_end_sec*sampling_rate))
    
    trace["GCaMP"]   = raw_GCaMP[crop_start_frame:crop_end_frame+1]
    trace["UV"]      = raw_UV[crop_start_frame:crop_end_frame+1]
    
    trace["Fratio"]  = trace["GCaMP"]/trace["UV"]
    num_frames       = len(trace["GCaMP"])
    trace["time_sec"] = np.linspace(1, num_frames, num_frames)/sampling_rate
    trace["time_sec"] = trace["time_sec"] - pre_sec
    
        
    ## F0 is the median value from -5 min to 0 pre-injection 
    if i == 0:           
        F0_window_start_frame = int(np.ceil(int(Z_F0_window_min[0])*60*sampling_rate))
        F0_window_end_frame   = int(np.floor(int(Z_F0_window_min[1])*60*sampling_rate))
    else:
        F0_window_start_frame = int(np.ceil(int(K_F0_window_min[0])*60*sampling_rate))
        F0_window_end_frame   = int(np.floor(int(K_F0_window_min[1])*60*sampling_rate))

    F0                   = np.median(trace["Fratio"][F0_window_start_frame:F0_window_end_frame])
   
    trace["dFF"]         = (trace["Fratio"] - F0)/F0*100 
    trace["frame"]       = np.arange(0, len(trace["GCaMP"]), 1)
		
    temp_trace           = signal.decimate(trace["dFF"], ds_factor1)
    temp_trace           = signal.decimate(temp_trace, ds_factor2)

    trace["dFF_ds"]        = signal.decimate(temp_trace, ds_factor3)
    trace["dFF_ds_sm"]     = smooth(trace["dFF_ds"], window_len = sm_factor)
    trace["time_sec_ds"]   = trace["time_sec"][::ds_factor] ## take every 50th element from the list
	    
    trace["dFF_sm_peak"]    = np.max(trace["dFF_ds_sm"])
    trace["dFF_sm_scaled"]  = trace["dFF_ds_sm"]/trace["dFF_sm_peak"]


    G             = trace["GCaMP"][::int(np.round(sampling_rate/ds_factor))]
    U             = trace["UV"][::int(np.round(sampling_rate/ds_factor))]
 
          
    ## G baseline is the function of U baseline, with the degree of 1 (first order equation)
    ## fitted Gfit_Z is the product of each U value from the fitted function   
    para_fit       = np.polyfit(U[0:len(baseline_index)], G[0:len(baseline_index)], 1)
    Gfit           = para_fit[0]*np.array(U) + np.array(para_fit[1])
    
    ##  y is the relative ratio of G_Z to the fitted Gfit_Z
    trace["GCaMP_zsfit"]    = G
    trace["UV_zsfit"]       = U
    trace["Ffit_ratio"]     = np.array(G/Gfit)   
    trace["zscore_fit"]     = (trace["Ffit_ratio"] - np.mean(trace["Ffit_ratio"][0:len(baseline_index)]))/np.std(trace["Ffit_ratio"])                                  
    trace["time_sec_zsfit"] = trace["time_sec"][::int(np.round(sampling_rate/ds_factor))]

    ## G baseline is the function of U baseline, with the degree of 1 (first order equation)
    ## fitted Gfit_Z is the product of each U value from the fitted function   
    para_fit2       = np.polyfit(U[int(math.ceil(len(baseline_index)/2)):len(baseline_index)], G[int(math.ceil(len(baseline_index)/2)):len(baseline_index)], 1)
    Gfit2           = para_fit2[0]*np.array(U) + np.array(para_fit2[1])
    
    ##  y is the relative ratio of G_Z to the fitted Gfit_Z
    trace["Ffit_ratio2"]     = np.array(G/Gfit2)   
    trace["zscore_fit2"]     = (trace["Ffit_ratio2"] - np.mean(trace["Ffit_ratio2"][0:len(baseline_index)]))/np.std(trace["Ffit_ratio2"])                                  

    trace["zscore_fit2"]     = (trace["Ffit_ratio2"] - np.mean(trace["Ffit_ratio2"][0:len(baseline_index)]))/np.std(trace["Ffit_ratio2"])                                 

               
    ## z score based on simple G/U 
    F0_window_start_frame = int(np.ceil(int(F0_window_min[0])*60*sampling_rate))
    F0_window_end_frame   = int(np.floor(int(F0_window_min[1])*60*sampling_rate))
            
    trace["zscore"]       = (trace["Fratio"] - np.mean(trace["Fratio"][F0_window_start_frame:F0_window_end_frame]))/np.std(trace["Fratio"])

    temp_trace                = signal.decimate(trace["zscore"], ds_factor1)
    temp_trace                = signal.decimate(temp_trace, ds_factor2)
    trace["zscore_ds"]        = signal.decimate(temp_trace, ds_factor3)
    trace["time_sec_ds_zs"]   = trace["time_sec"][::ds_factor] ## take every 50th element from the list    
        


    ##  plot dFF       
    fig, ax = plt.subplots(6, 1, figsize = (12, 24))
    
    ## subplot 0. raw GCaMP
    ax[0].plot(trace["time_sec"]/60, trace["GCaMP"], color = "green", linewidth = 2, label = "GCaMP")
    ax[0].set_ylim(0.95*min(trace["GCaMP"]), 1.05*max(trace["GCaMP"]))
    ax[0].set_ylabel("raw (a.u.)")
    ax[0].legend(loc='upper right', fontsize = 14)

    ## subplot 1. raw UV
    ax[1].plot(trace["time_sec"]/60, trace["UV"], color = "red", linewidth = 2, label = "UV")
    ax[1].set_ylim(0.95*min(trace["UV"]), 1.05*max(trace["UV"]))
    ax[1].set_ylabel("raw (a.u.)")
    ax[1].legend(loc='upper right', fontsize = 14)
    
    ## subplot 2. F = raw GCaMP/raw UV
    ax[2].plot(trace["time_sec"]/60, trace["Fratio"], color = "black", linewidth = 2, label = "GCaMP/UV")
    ax[2].set_ylim(0.95*min(trace["Fratio"]), 1.05*max(trace["Fratio"]))
    ax[2].set_ylabel("ratio")
    ax[2].legend(loc='upper right', fontsize = 14)
    
    ## subplot 3. dFF = (F-F0)/F0 (%)
    ax[3].plot(trace["time_sec"]/60, trace["dFF"], color = "blue", linewidth = 2, label = "dFF")
    if min(trace["dFF"]) > 0:
        ax[3].set_ylim(0.95*min(trace["dFF"]), 1.05*max(trace["dFF"]))
    else: 
        ax[3].set_ylim(1.05*min(trace["dFF"]), 1.05*max(trace["dFF"]))
    ax[3].set_ylabel("dF/F0 (%)")  
    ax[3].legend(loc='upper right', fontsize = 14)
    
    ## subplot 4. dFF_ds = downsampled dFF  
    ax[4].plot(trace["time_sec_ds"]/60, trace["dFF_ds"], color = "purple", linewidth = 2, label = "dFF_ds")
    if min(trace["dFF_ds"]) > 0:
        ax[4].set_ylim(0.95*min(trace["dFF_ds"]), 1.05*max(trace["dFF_ds"]))
    else: 				
        ax[4].set_ylim(1.05*min(trace["dFF_ds"]), 1.05*max(trace["dFF_ds"]))
    ax[4].set_ylabel("dF/F0 (%) ds") 
    ax[4].set_xlabel("time (min)")
    ax[4].legend(loc='upper right', fontsize = 14)
    
    ## subplot 5. zscore_ds = downsampled and smoothened dFF
    ax[5].plot(trace["time_sec_ds"]/60, trace["dFF_ds_sm"], color = "black", linewidth = 2, label = "dFF_ds_sm")
    if min(trace["dFF_ds_sm"]) > 0:
        ax[5].set_ylim(0.95*min(trace["dFF_ds_sm"]), 1.05*max(trace["dFF_ds_sm"]))
    else:
        ax[5].set_ylim(1.05*min(trace["dFF_ds_sm"]), 1.05*max(trace["dFF_ds_sm"]))
    ax[5].set_ylabel("dF/F0 (%) ds sm")
    ax[5].set_xlabel("time (min)")
    ax[5].legend(loc='upper right', fontsize = 14)
    
    suptitle = trace["mouseID"] +"_"+ trace["date"] +"_"+ trace["stimulus"] + "_dFF"
    plt.suptitle("%s" %suptitle) 
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])    
    
    save_file_path = analysis_path + "/" + suptitle + ".png"
    plt.savefig(save_file_path)
    # save_file_path = analysis_path + "/" + suptitle + ".pdf"
    # plt.savefig(save_file_path)
    # plt.show()
    plt.close()

    
        
    ##  plot z score, not fitted, simple G/U       
    fig, ax = plt.subplots(5, 1, figsize = (12, 20))
    
    ## subplot 0. raw GCaMP
    ax[0].plot(trace["time_sec"]/60, trace["GCaMP"], color = "green", linewidth = 2, label = "GCaMP")
    ax[0].set_ylim(0.95*min(trace["GCaMP"]), 1.05*max(trace["GCaMP"]))
    ax[0].set_ylabel("raw (a.u.)")
    ax[0].legend(loc='upper right', fontsize = 14)

    ## subplot 1. raw UV
    ax[1].plot(trace["time_sec"]/60, trace["UV"], color = "red", linewidth = 2, label = "UV")
    ax[1].set_ylim(0.95*min(trace["UV"]), 1.05*max(trace["UV"]))
    ax[1].set_ylabel("raw (a.u.)")
    ax[1].legend(loc='upper right', fontsize = 14)
    
    ## subplot 2. F = raw GCaMP/raw UV
    ax[2].plot(trace["time_sec"]/60, trace["Fratio"], color = "black", linewidth = 2, label = "GCaMP/UV")
    ax[2].set_ylim(0.95*min(trace["Fratio"]), 1.05*max(trace["Fratio"]))
    ax[2].set_ylabel("ratio")
    ax[2].legend(loc='upper right', fontsize = 14)
    
    ## subplot 3. z score 
    ax[3].plot(trace["time_sec"]/60, trace["zscore"], color = "blue", linewidth = 2, label = "zscore")
    if min(trace["zscore"]) > 0:
        ax[3].set_ylim(0.95*min(trace["zscore"]), 1.05*max(trace["zscore"]))
    else: 
        ax[3].set_ylim(1.05*min(trace["zscore"]), 1.05*max(trace["zscore"]))
    ax[3].set_ylabel("z score")  
    ax[3].legend(loc='upper right', fontsize = 14)
    
    ## subplot 4. zscore_ds = downsampled dFF  
    ax[4].plot(trace["time_sec_ds_zs"]/60, trace["zscore_ds"], color = "purple", linewidth = 2, label = "zscore ds")
    if min(trace["zscore_ds"]) > 0:
        ax[4].set_ylim(0.95*min(trace["zscore_ds"]), 1.05*max(trace["zscore_ds"]))
    else: 				
        ax[4].set_ylim(1.05*min(trace["zscore_ds"]), 1.05*max(trace["zscore_ds"]))
    ax[4].set_ylabel("z score") 
    ax[4].set_xlabel("time (min)")
    ax[4].legend(loc='upper right', fontsize = 14)
    
    
    suptitle = trace["mouseID"] +"_"+ trace["date"] +"_"+ trace["stimulus"] + "_zs"
    plt.suptitle("%s" %suptitle) 
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])    
  
    
    save_file_path = analysis_path + "/" + suptitle + ".png"
    plt.savefig(save_file_path)
    # save_file_path = analysis_path + "/" + suptitle + ".pdf"
    # plt.savefig(save_file_path)
    # plt.show()
    plt.close()   
    
    ##  plot z score, fitted 
    fig, ax = plt.subplots(5, 1, figsize = (12, 20))
    ## subplot 0. raw GCaMP
    ax[0].plot(trace["time_sec_zsfit"]/60, trace["GCaMP_zsfit"], color = "green", linewidth = 2, label = "GCaMP")
    ax[0].set_ylim(0.95*min(trace["GCaMP_zsfit"]), 1.05*max(trace["GCaMP_zsfit"]))
    ax[0].set_ylabel("raw (a.u.)")
    ax[0].legend(loc='upper right')

    ## subplot 1. raw UV
    ax[1].plot(trace["time_sec_zsfit"]/60, trace["UV_zsfit"], color = "red", linewidth = 2, label = "UV")
    ax[1].set_ylim(0.95*min(trace["UV_zsfit"]), 1.05*max(trace["UV_zsfit"]))
    ax[1].set_ylabel("raw (a.u.)")
    ax[1].legend(loc='upper right')
    
    ## subplot 2. F = raw GCaMP/fitted GCaMP
    ax[2].plot(trace["time_sec_zsfit"]/60, trace["Ffit_ratio"], color = "black", linewidth = 2, label = "GCaMP/UVfit")
    ax[2].set_ylim(0.95*min(trace["Ffit_ratio"]), 1.05*max(trace["Ffit_ratio"]))
    ax[2].set_ylabel("ratio")
    ax[2].legend(loc='upper right')
    
    ## subplot 3. z score  = 
    ax[3].plot(trace["time_sec_zsfit"]/60, trace["zscore_fit"], color = "blue", linewidth = 2, label = "zscore_fit")
    if min(trace["zscore_fit"]) > 0:
        ax[3].set_ylim(0.95*min(trace["zscore_fit"]), 1.05*max(trace["zscore_fit"]))
    else: 
        ax[3].set_ylim(1.05*min(trace["zscore_fit"]), 1.05*max(trace["zscore_fit"]))
    ax[3].set_ylabel("zscore_fit")  
    ax[3].legend(loc='upper right')
    
    ## subplot 3. z score  = 
    ax[4].plot(trace["time_sec_zsfit"]/60, trace["zscore_fit2"], color = "purple", linewidth = 2, label = "zscore_fit2")
    if min(trace["zscore_fit2"]) > 0:
        ax[4].set_ylim(0.95*min(trace["zscore_fit2"]), 1.05*max(trace["zscore_fit2"]))
    else: 
        ax[4].set_ylim(1.05*min(trace["zscore_fit2"]), 1.05*max(trace["zscore_fit2"]))
    ax[4].set_ylabel("zscore_fit2")  
    ax[4].legend(loc='upper right')
    
    
    suptitle = trace["mouseID"] +"_"+ trace["date"] +"_"+ trace["stimulus"] + "_zsfit"
    plt.suptitle("%s" %suptitle) 
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])    
         
    save_file_path = analysis_path + "/" + suptitle + ".png"
    plt.savefig(save_file_path)
    # save_file_path = analysis_path + "/" + suptitle + ".pdf"
    # plt.savefig(save_file_path)
    # plt.show()
    plt.close()
      
        #%%############################################################
    #      SAVE trace data as a pkl file
    ################################################################        
    suptitle        = trace["mouseID"] + "_" + trace["date"] + "_" + trace["stimulus"] + "_zs" #+ str(trace["age_mo"]) + "mo" + "_zs"          
    save_file_path  = analysis_path + "/" + suptitle + ".pkl"
    output          = open(save_file_path, 'wb')
    pk.dump(trace, output)
    output.close()  
   
