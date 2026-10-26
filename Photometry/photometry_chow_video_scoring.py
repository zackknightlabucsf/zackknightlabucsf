# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:37:24 2022

Fiber photometry analysis aligned to video-scored chow feeding.

Loads processed photometry z-score traces and video-scored chow-bite
timestamps to analyze feeding behavior and associated neural activity.

Identifies individual bites and feeding bouts, calculates bout-aligned
PSTHs, and quantifies z-score responses during bites, feeding bouts,
and interbout intervals.

Generates individual recording plots and exports processed data,
summary measurements, and figures as .pkl, .csv, PNG, and PDF files.


"""
import glob
import os
import pickle as pk
import numpy as np
import pandas as pd
from datetime import datetime
import math
import scipy.stats as ss
# from scipy.integrate import simps
from sklearn.neighbors import KDTree
import matplotlib.pyplot as plt  # standard Python plotting library
import matplotlib
matplotlib.rcParams['font.size'] = 18 # set font size for all figures


####################################################
#   import chow in time from the TTL file
####################################################
TTLpath  = r"D:\Data_Dropbox_Sync\Mouse_metadata\photometrymiceTTL_wchow.csv"
TTLdata = pd.read_csv(TTLpath)   

####################################################
#   import video timestamps for each chow bite
####################################################
videostamp_path       = r"D:\Data_Dropbox_Sync\Photometry_DROPBOX\HJ_videos_TTL\video_timestamp_files"
all_videostamp_files  = glob.glob(videostamp_path + "/*.csv")

#%% Importing the Data

state  = "FD"
access = "chow"

ds_factor1            = 5
ds_factor2            = 5
ds_factor3            = 5   
ds_factor             = int(ds_factor1*ds_factor2*ds_factor3)


print("Enter the full path for directory with .pkl files: ")
path      = input()


today = datetime.today().strftime('%Y%m%d')

analysis_path = path + "/analysis_bite_video_" + today
if not os.path.exists(analysis_path):
   os.makedirs(analysis_path)   


pre_length      = 5
post_length     = 30

fix_ylim = 1 # yes if 1, no if 0
if fix_ylim == 1:
    ylim = [-1, 4]

m               = 0
frame_rate      = 1017.2527


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



all_files                = glob.glob(path + "/*.pkl")
file_list                = []
expID_list               = []
df_bout_zs              = {}
df_lick_zs              = {}

data = pd.DataFrame(columns = ["expID", "mouseID", "date", "BW",
                               "zscore", "zscore_ds", "zscore_firstlick","zscore_firstbout", 
                               "state", "access", "bite_type", "bite_sec", 
                               
                               # "lickrate_sec", "lickrate_min", "lickrate_sec_sm60",
                               # "lickon_firstlick", "lickrate_sec_firstlick", 
                               # "lickrate_min_firstlick", "lickrate_sec_sm60_firstlick", 
                               # "lickrate_min_firstlick_BW", "lickrate_sec_sm60_firstlick_BW", 
                                 
                               # "lickon_firstBOUT", "lickrate_sec_firstBOUT", 
                               # "lickrate_min_firstBOUT", "lickrate_sec_sm60_firstBOUT",

                               "time_sec", "time_sec_ds_zs", 
                               # "time_sec_bin", "time_min_bin",
                               # "time_sec_ds_firstlick", "time_sec_ds_firstBOUT",
                                          
                               "LICK_on", "LICK_off", "LICK_x", "LICK_y",
                               "BOUT_on", "BOUT_off", "BOUT_x", "BOUT_y",
                               # "bout_on_ind", "bout_off_ind", "bout_on_zs", "bout_off_zs", 
                                "BOUT_size", "BOUT_length", "BOUT_number", 
                               # "bout_meanzs", "all_lick_meanzs",
                               # "all_bout_meanzs", "all_interbout_meanzs",
                               # "interbout_size", "interbout_meanzs",
                                
                               # "cumlick_sec", "cumbout_sec", 
                               # "cumlick_sec_firstlick", "cumbout_sec_firstlick", 
                               # "cumlick_sec_firstlick_BW",
                               # "cumlick_sec_firstbout", "cumbout_sec_firstbout", 
                               # "bout_on_PSTH", 
                                "BOUT_on_PSTH_KDquery"
                               
                               ])

# data["zscore"]             = data["zscore"].astype(object)
data["zscore_ds"]          = data["zscore_ds"].astype(object) 
data["time_sec"]           = data["time_sec"].astype(object) 
data["time_sec_ds_zs"]        = data["time_sec_ds_zs"].astype(object) 

data["LICK_on"]           = data["LICK_on"].astype(object) 
data["LICK_off"]          = data["LICK_off"].astype(object) 
data["LICK_x"]            = data["LICK_x"].astype(object) 
data["LICK_y"]            = data["LICK_y"].astype(object) 
data["BOUT_on"]           = data["BOUT_on"].astype(object) 
data["BOUT_off"]          = data["BOUT_off"].astype(object) 
data["BOUT_x"]            = data["BOUT_x"].astype(object) 
data["BOUT_y"]            = data["BOUT_y"].astype(object) 
data["BOUT_size"]         = data["BOUT_size"].astype(object)
data["BOUT_length"]       = data["BOUT_length"].astype(object)


#%% 
for i, file_path in enumerate (all_files):
    expID            = file_path[(len(path)+1):-4] ## previous called as "file_name_cut
    inputfile        = open(file_path, 'rb')
    tempdata         = pk.load(inputfile)
    inputfile.close()

    data.at[i, "expID"]             = expID    
    data.at[i, "mouseID"]           = tempdata["mouseID"]
    data.at[i, "date"]              = tempdata["date"]
    # data.at[i, "age"]               = tempdata["age_mo"]
    data.at[i, "state"]             = state
    data.at[i, "access"]            = access 
    data.at[i, "time_sec"]          = tempdata["time_sec"]
    # if "AVP" in path:
    #     data.at[i, "time_sec_ds_zs"]        = tempdata["time_sec"][::ds_factor]
    # else:
    try: 
        data.at[i, "time_sec_ds_zs"]    = tempdata["time_sec_ds_zs"]
    except:
        data.at[i, "time_sec_ds_zs"]    = tempdata["time_sec"][::ds_factor] ## take every 50th element from the list
    
    ####################################################
    #   import matching chow bite data from video
    ####################################################

    for j, file_path in enumerate(all_videostamp_files):
        expID_           = file_path[(len(videostamp_path)+1):-4]
        mouseID_         = expID_.split("_")[0]
        expdate_         = expID_.split("_")[2]
        if (expdate_ == str(data.at[i, "date"])) & (mouseID_ == str(data.at[i, "mouseID"])):
            video_timestamps            = pd.read_csv(file_path)
            data.at[i, "bite_type"]     = video_timestamps["Behavior"]
            data.at[i, "bite_sec"]      = video_timestamps["Start (s)"]
            video_chow_in               = float(data.at[i, "bite_sec"][0])
            print("Processing video chow bite data for %s" %expID)
            break
    
    try: 
        data.at[i, "bite_sec"][0]
    except:
        print("****Video chow bite data does not match for %s*****" %expID)
        continue

    file_list.append(expID)  
    expID_list.append(expID)
        
    
    bite_on_temp  = data.at[i, "bite_sec"] - video_chow_in
    bite_on_temp  = bite_on_temp[(bite_on_temp > 0)]
    bite_on_temp  = bite_on_temp[(bite_on_temp < data.at[i, "time_sec"][-1])]

    data.at[i, "LICK_on"]  = bite_on_temp
    # data.at[i, "LICK_off"] = LICK_off_temp
    data.at[i, "LICK_off"] = bite_on_temp
    
    bite_interval                    = np.diff(bite_on_temp)
    # data.at[i, "interlick_interval"] = bite_interval

    
    try:
        data.at[i, "BW"]    = tempdata["BW"]
    except:
        for j in range(len(TTLdata)):
            if (str(int(TTLdata.at[j, "ExpDate"])) == str(data.at[i, "date"])) & (TTLdata.at[j, "mouseID"] == data.at[i, "mouseID"]):
                data.at[i, "BW"]          = TTLdata.at[j, "BW"]
                break

    data.at[i, "zscore"]     = tempdata["zscore"]
    data.at[i, "zscore_ds"]  = tempdata["zscore_ds"]       
    # Add the first and last time stamps to make tails on the TTL stream
    try:
        data.at[i, "LICK_x"]   = np.append(np.append(data.at[i, "time_sec"][0], np.reshape(np.kron([data.at[i, "LICK_on"], data.at[i, "LICK_off"]],
                                               np.array([[1], [1]])).T, [1,-1])[0]), data.at[i, "time_sec"][-1])
        licks_exist = 1                
        LICK_size            = len(data.at[i, "LICK_on"])
        d                    = np.ones(LICK_size)
    
    # Add zeros to beginning and end of 0,1 value array to match len of LICK_x
        data.at[i, "LICK_y"]   = np.append(np.append(0,np.reshape(np.vstack([np.zeros(LICK_size),
                           d, d, np.zeros(LICK_size)]).T, [1, -1])[0]),0)
        
    except: 
        print(expID, "has no bites")
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
        bout_interval_cutoff_sec = 30
        bout_size_cutoff         = 2
    
        bout_ons        = np.asarray([])
        bout_offs       = np.asarray([])
        bout_size       = np.asarray([])
        bout_length     = np.asarray([])
        bout_data       = np.asarray([])   

        lickon          = data.at[i, "LICK_on"]
        lickoff         = data.at[i, "LICK_off"]
        # lickoff          = data.at[i, "LICK_on"]
#        bout_licktimes  = []
        k = 1 # index for while loop
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
            LICK_EVENT_x    = np.append(data.at[i, "time_sec"][0], np.append(
                            np.reshape(np.kron([bout_ons , bout_offs], np.array([[1], [1]])).T, [1,-1])[0], data.at[i, "time_sec"][-1]))          
            LICK_EVENT_sz   = len(bout_ons)
            LICK_EVENT_y = np.append(np.append(
                0, np.reshape(np.vstack([np.zeros(LICK_EVENT_sz), bout_data, bout_data, np.zeros(LICK_EVENT_sz)]).T, [1 ,-1])[0]), 0)
            
            
            data.at[i, "BOUT_on"]     = bout_ons
            data.at[i, "BOUT_off"]    = bout_offs
            data.at[i, "BOUT_x"]      = LICK_EVENT_x
            data.at[i, "BOUT_y"]      = LICK_EVENT_y
            data.at[i, "BOUT_size"]   = bout_size
            data.at[i, "BOUT_length"] = bout_length 
            
            # data.at[i,"interbout_lick_size"]    = np.asarray([])
            # data.at[i,"interbout_lick_length"]  = np.asarray([])
		
    		# ## save interbout lick on times 
      #       for j in range(len(bout_ons)-1):
      #           interbout_lickons = lickon[lickon > bout_offs[j]]
      #           interbout_lickons = interbout_lickons[interbout_lickons < bout_ons[j+1]]
      #           try: 
      #               lick_length       = interbout_lickons[interbout_lickons.index[-1]] - interbout_lickons[interbout_lickons.index[0]]
      #           except:
      #               lick_length       = 0
      #           data.at[i,"interbout_lick_size"]   = np.append(data.at[i,"interbout_lick_size"], len(interbout_lickons))
      #           data.at[i,"interbout_lick_length"] = np.append(data.at[i,"interbout_lick_length"], lick_length)
                              
            bouts_exist    = 1 
            
        except:     
            bouts_exist    = 0 
            print (expID, " has no bouts")
                
                
#%%######################################################################
#       Entire timecourse: zs and lick +  bout
#########################################################################
        #################################################################
        #       FIGURE 1 zs + bites
        #################################################################
        y_scale = 1
        y_shift = 0
        
        fig, ax = plt.subplots(3, 1, figsize = (10, 10))
        ax1     = plt.subplot2grid((6, 1), (0, 0), rowspan=4, colspan=1)   
        ax1.plot(data.at[i, "time_sec_ds_zs"]/60, data.at[i, "zscore_ds"], linewidth=2.5, color='black', label='GCaMP/UV deltaF/F')
        ax1.set_ylabel("zs (%)", fontsize = 20)
        ax1.set_xlim(-pre_length, post_length)
        ax1.set_xlabel("time (min)")
        if fix_ylim == 1:
            ax1.set_ylim(ylim[0], ylim[1])
        # ax1.axes.get_xaxis().set_visible(False)        
        
        ax2     = plt.subplot2grid((6, 1), (4, 0), rowspan=1, colspan=1)   
        ax2.plot(data.at[i, "LICK_x"], y_scale*data.at[i, "LICK_y"] + y_shift, linewidth=0.5, color='maroon', label='Bite Event')
        ax2.set_xlim(-pre_length, post_length)
        ax2.axes.get_xaxis().set_visible(False)    
        ax2.axis("off")
        
        ax3     = plt.subplot2grid((6, 1), (5, 0), rowspan=1, colspan=1)   
        if bouts_exist == 0:
            ax3.plot(data.at[i, "LICK_x"], y_scale*data.at[i, "LICK_y"] + y_shift, linewidth=0.5, color='maroon', label='Bite Event')

        else:
            ax3.plot(LICK_EVENT_x, y_scale*LICK_EVENT_y+y_shift, linewidth=1, color='maroon', label='Bite Bout')
        ax3.set_xlim(-pre_length, post_length)
        ax3.axis("off")
        
        suptitle = data.at[i, "mouseID"] + "_" + data.at[i, "date"] + "_" + data.at[i, "state"] + "_" + data.at[i, "access"]         
        plt.suptitle("%s" %suptitle, fontsize = 20) 
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])           
        save_file_path = analysis_path + "/" + suptitle + ".png"
        plt.savefig(save_file_path)
        save_file_path = analysis_path + "/" + suptitle + ".pdf"
        plt.savefig(save_file_path)
        plt.show()


        #################################################################
        #       FIGURE 2  zs + bites  
        #################################################################
        # Cosmetics for plotting
        if "WD" in data.at[i, "state"]: 
            y_scale    = -np.min(data.at[i, "zscore"])/3
            y_shift    = np.min(data.at[i, "zscore"])/3*4
        else:
            y_scale    = np.max(data.at[i, "zscore"])/3
            y_shift    = -np.max(data.at[i, "zscore"])  

        fig, ax = plt.subplots(2, 1, figsize = (10, 10))
        ax1     = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)   
        lns1    = ax1.plot(data.at[i, "time_sec_ds_zs"]/60, data.at[i, "zscore_ds"], linewidth=2.5, color='black', alpha = 0.8, label='GCaMP')
        ax1.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
        ax1.set_ylabel("zs (%)", fontsize = 20)
        ax1.set_xlim(-pre_length, post_length)
        ax1.set_xlabel("time (min)")
        if fix_ylim == 1:
            ax1.set_ylim(ylim[0], ylim[1])
        
        

        ax2     = ax1.twinx()       
        lns2    = ax2.plot(data.at[i, "LICK_x"]/60, y_scale*data.at[i, "LICK_y"] + y_shift, linewidth=0.5, color='maroon', alpha = 0.5, label='Bite')
        ax2.axes.get_yaxis().set_visible(False)
        ax2.set_xlim(-pre_length, post_length)

        lns     = lns1 + lns2
        labels  = [l.get_label() for l in lns]        
        ax1.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16)

        ax3     = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)   
        lns3    = ax3.plot(data.at[i, "time_sec_ds_zs"]/60, data.at[i, "zscore_ds"], linewidth=1.5, color='black', label='GCaMP')
        ax3.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
        # ax3.axvline (x = stim_sec/60, linewidth = 2, linestyle = "--", color = "magenta")
        ax3.set_ylabel("zs (%)", fontsize = 18)
        ax3.set_xlim(-pre_length, post_length)
        ax3.set_xlabel("time (min)")

        if bouts_exist == 1:        
            ax4     = ax3.twinx()
            lns4    = ax4.plot(LICK_EVENT_x/60, y_scale*LICK_EVENT_y+y_shift, linewidth=2, color='maroon', alpha = 0.5, label='Bite Bout')

        else:
            ax4     = ax3.twinx()
            lns4    = ax4.plot(data.at[i, "LICK_x"]/60, y_scale*data.at[i, "LICK_y"] + y_shift, linewidth=0.5, color='maroon', alpha = 0.5, label='Bite Bout')
        ax4.set_xlim(-pre_length, post_length) 
        ax4.axis("off")
        lns     = lns3 + lns4 
        labels  = [l.get_label() for l in lns]        
        ax3.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16)
       
        suptitle = data.at[i, "mouseID"] + "_" + data.at[i, "date"] + "_" + data.at[i, "state"] + "_" + data.at[i, "access"]         
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
        
        # mask                    = data.at[i, "time_sec"] > 0
#        time_sec_dims               = np.expand_dims(data.at[i, "time_sec"], axis = 1)
#        tree[i]                     = KDTree(time_sec_dims, leaf_size = 2)  
        # data.at[i, "time_sec_KDtree"] = data.at[i, "time_sec_KDtree"].astype(object)
#        data.at[i, "time_sec_KDtree"] = tree[i]     


        time_sec_dims             = np.expand_dims(data.at[i, "time_sec_ds_zs"], axis = 1)
        tree                      = KDTree(time_sec_dims, leaf_size = 2)  

        frame_rate      = 1017.2527
        psth_width_sec  = 120                                         
        psth_width      = math.floor(frame_rate*psth_width_sec/ds_factor)
        psth_num_points = 2*psth_width+1
        psth_plot_times = np.linspace(-psth_width_sec, psth_width_sec, psth_num_points, endpoint = True, dtype = None)
                
        psth            = np.empty([len(data.at[i, "BOUT_size"]), psth_num_points])
        psth[:]         = np.nan
        psth_y_min_list = []
        psth_y_max_list = []
        
        print ("total number of bouts are %d" %len(data.at[i, "BOUT_on"]))
        for k, bout_on in enumerate (data.at[i, "BOUT_on"]):
            bout_on     = bout_on.reshape(1, -1)
            ind         = tree.query(bout_on)[1][0][0] 
            # print ("processing PSTH for bout # %d, bout_on time is %.2f" %(k, bout_on[0][0]))
            ## Process without error if PSTH post-width exceeds the recording, fill with nan.
            if ind + psth_width + 1 >= len(data.at[i, "zscore_ds"]):  
                last_post_psth                          = len(data.at[i, "zscore_ds"]) - ind
                psth[k, : psth_width + last_post_psth]  = data.at[i, "zscore_ds"][ind - psth_width : ] - np.nanmedian(data.at[i, "zscore_ds"][ind - psth_width : ind]) 
            else:    
                psth[k, :]  = data.at[i, "zscore_ds"][ind - psth_width : ind + psth_width +1] - np.nanmedian(data.at[i, "zscore_ds"][ind - psth_width : ind])

            psth_y_min_list.append(np.nanmin(psth[k, :]))
            psth_y_max_list.append(np.nanmax(psth[k, :]))
        data.at[i, "BOUT_on_PSTH_KDquery"] = psth 

#%%#################################################################
#       ZOOM in by bouts 
###################################################################
        ## Zoomed-in bouts are saved into one dataframe.
        df_bout_zs[m] = pd.DataFrame(columns = ["bout_on", "bout_off", "bout_size",
                                               "bout_number", "bout_on_zs", "bout_off_zs",
                                               "bout_delzs", "bout_delzs_perlick",
                                               "bout_plot_times", "bout_plot_zs",
                                               "bout_lick_x", "bout_lick_y",
                                               "interbout_number",
                                                "interbout_licks", 
												"interbout_meanzs"])
        
        df_lick_zs[m] = pd.DataFrame(columns = ["lick_meanzs"])
        
             
        for k, (bout_on, bout_off) in enumerate (zip(data.at[i, "BOUT_on"], data.at[i, "BOUT_off"])):
             ## Use the KDtree created for PSTH to find bout on and bout off   
             bout_on           = bout_on.reshape(1, -1)
             ind_on            = tree.query(bout_on)[1][0][0] 
             
             bout_off          = bout_off.reshape(1, -1)
             ind_off           = tree.query(bout_off)[1][0][0]
        
             df_bout_zs[m].at[k, "bout_on_ind"]  = ind_on
             df_bout_zs[m].at[k, "bout_off_ind"] = ind_off           
             df_bout_zs[m].at[k, "bout_on_zs"]  = data.at[i, "zscore_ds"][ind_on]
             df_bout_zs[m].at[k, "bout_off_zs"] = data.at[i, "zscore_ds"][ind_off]
             
             df_bout_zs[m].at[k, "bout_delzs"]         = data.at[i, "zscore_ds"][ind_off] - data.at[i, "zscore_ds"][ind_on]
             df_bout_zs[m].at[k, "bout_delzs_perlick"] = df_bout_zs[m].at[k, "bout_delzs"]/data.at[i, "BOUT_size"][k]
        
             df_bout_zs[m].at[k, "bout_number"]  = k+1
             df_bout_zs[m].at[k, "bout_on"]      = data.at[i, "BOUT_on"][k]
             df_bout_zs[m].at[k, "bout_off"]     = data.at[i, "BOUT_off"][k]
             df_bout_zs[m].at[k, "bout_size"]    = data.at[i, "BOUT_size"][k]
             
             df_bout_zs[m].at[k, "bout_meanzs"]     = np.mean(data.at[i,"zscore_ds"][ind_on:ind_off+1])
	 
																		
			## Process data for interbout zs
        for k in range(len(df_bout_zs[m])-1):       
             df_bout_zs[m].at[k, "interbout_number"]     = int(k+1)
             # df_bout_zs[m].at[k, "interbout_size"]       = data.at[i,"interbout_lick_size"][k]
             
             interbout_on_ind                             = int(df_bout_zs[m].at[k, "bout_off_ind"])
             interbout_off_ind                            = int(df_bout_zs[m].at[k+1, "bout_on_ind"])
             df_bout_zs[m].at[k, "interbout_meanzs"]    = np.mean(data.at[i,"zscore_ds"][interbout_on_ind:interbout_off_ind+1])

			# Get mean zs for individual lick events
        for k, lick_on in enumerate (data.at[i,"LICK_on"]):
             lick_on                                = np.array(lick_on).reshape(1, -1)       			
             ind_lick_on                            = tree.query(lick_on)[1][0][0]
             df_lick_zs[m].at[k, "lick_meanzs"]     = data.at[i,"zscore_ds"][ind_lick_on]
             
             
        df_bout_zs[m]["lick_meanzs"]  = df_lick_zs[m]["lick_meanzs"]    
            
        #%%###############################################################
        ## Figure: zs with licks, ZOOMED in for each bout          
        #################################################################
        pre_post_sec    = 3.0   ## 3 sec pre- and post- bout
        frame_rate      = 1017.2527
        pre_post_width  = math.floor(pre_post_sec*frame_rate)
          
        fig, ax     = plt.subplots(len(df_bout_zs[m]), 1, figsize = (8, 2*len(df_bout_zs[m])))       
        for k in range(len(df_bout_zs[m])):  
             ax1     = plt.subplot2grid((len(df_bout_zs[m]), 1), (k, 0), rowspan=1, colspan=1)       
             legend  = "bout # %d" %(k+1)
        
             ind_start             = int(df_bout_zs[m].at[k, "bout_on_ind"] - pre_post_width)
             ind_end               = int(df_bout_zs[m].at[k, "bout_off_ind"] + pre_post_width+1)  
             x = data.at[i, "time_sec_ds_zs"][ind_start:ind_end]
             y = data.at[i, "zscore_ds"][ind_start:ind_end]

             lick_x               = data.at[i, "LICK_x"]
             lick_x_masked        = lick_x[lick_x >= x[0]]
             lick_x_masked_masked = lick_x_masked[lick_x_masked <= x[-1]]
             lick_y               = data.at[i, "LICK_y"][1:1+len(lick_x_masked_masked)]
            
             df_bout_zs[m].at[k,"bout_plot_times"] = x
             df_bout_zs[m].at[k,"bout_plot_zs"]   = y
             df_bout_zs[m].at[k,"bout_lick_x"] = lick_x_masked_masked
             df_bout_zs[m].at[k,"bout_lick_y"] = lick_y

             if k == 0:
                 lns1    = ax1.plot(x, y, linewidth=2, color = "red", label = legend)
#                 ax1.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18, color = "red")            
                 ax2     = ax1.twinx() 
                 lns2    = ax2.plot(lick_x_masked_masked, lick_y, linewidth=1, color='purple', alpha = 0.5, label='bite')
                 ax2.axes.get_yaxis().set_visible(False)
             else:
                 lns1    = ax1.plot(x, y, linewidth=2, color = "black", label = legend)
#                 ax1.annotate("%s" %legend, xy = (psth_plot_times[int(len(psth_plot_times)/10)], psth_y_min*0.9), fontsize = 18)            
        
                 ax2     = ax1.twinx()       
                 lns2    = ax2.plot(lick_x_masked_masked, lick_y, linewidth=1, color='purple', alpha = 0.5, label='bite')
                 ax2.axes.get_yaxis().set_visible(False)
                  
             lns     = lns1 + lns2
             labels  = [l.get_label() for l in lns]        
             ax1.legend(lns, labels, loc="upper left", bbox_to_anchor = (1, 1), fontsize = 16) 
             ax1.set_ylabel("zscore", size = 18)
             
             if k == len(df_bout_zs[m]) -1 :
                 ax1.set_xlabel("Time (sec)", size = 18)

             plt.axvline(x = df_bout_zs[m].at[k, "bout_on"], linewidth = 2, linestyle = "--", color = "gray")
             plt.axvline(x = df_bout_zs[m].at[k, "bout_off"], linewidth = 2, linestyle = "--", color = "gray")
             plt.yticks(size = 18)
             plt.ylabel("zscore", size = 18)
             plt.xticks(size = 18)
             plt.xlabel("Time (sec)", size = 18)
             
        suptitle = data.at[i, "mouseID"] + "_" + data.at[i, "date"] + "_" + data.at[i, "state"] + "_" + data.at[i, "access"] + "_boutzoom"                   
        plt.suptitle("%s\nfirst lick in each bout " %suptitle, size = 18)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
        save_file_path = analysis_path + "/" + suptitle + ".png"
        plt.savefig(save_file_path)
        plt.show()    
        
 
		   
		# ###################################################################################
		# ## Additional Figure: bar graph of  (i) mean zs all licks vs mean zs during interbout interval 
        #                                     (ii)  mean zs during all licks that are part of bouts vs mean zs during interbout interval 
		# ##################################################################################
        fig, ax       = plt.subplots(nrows = 1, ncols = 2, figsize = (8, 6)) ## 4*of rows
# 		fig, axes     = plt.subplots(figsize = (3, 6)) ## 4*of rows
        bar_width     = 0.75

        ## 1. mean zs during all licks vs mean zs during interbout interval
        values        = [df_lick_zs[m]["lick_meanzs"].mean(), df_bout_zs[m]["interbout_meanzs"].mean()]
        labels        = ['all bites', 'interbout']
        errorbars     = [ss.sem(df_lick_zs[m]["lick_meanzs"]), ss.sem(df_bout_zs[m]["interbout_meanzs"][:-1])]      
        ax1           = plt.subplot2grid((1, 2), (0, 0), rowspan = 1, colspan = 1)
        ax1.bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)
        #   ax7   = plt.subplot2grid((4, 4), (0, k*2), rowspan=2, colspan=2)  
        ax1.set_ylabel("zscore")
        ax1.set_title("mean zs")

        ## . mean zs during all licks that are part of bouts vs mean zs during interbout interval 
        values        = [df_bout_zs[m]["bout_meanzs"].mean(), df_bout_zs[m]["interbout_meanzs"].mean()]
        labels        = ['all bouts', 'interbout']
        errorbars     = [ss.sem(df_bout_zs[m]["bout_meanzs"]), ss.sem(df_bout_zs[m]["interbout_meanzs"][:-1])]  
        ax2           = plt.subplot2grid((1, 2), (0, 1), rowspan = 1, colspan = 1)            
        ax2.bar(labels, values, yerr = errorbars, width = bar_width, color = "black", ecolor='black', capsize = 10, alpha = 0.5)
        ax2.set_ylabel("zscore")
        ax2.set_title("mean zs")        
        suptitle = data.at[i, "mouseID"] + "_" + data.at[i, "date"] + "_" + data.at[i, "state"] + "_" + data.at[i, "access"]                            
        plt.suptitle("%s" %suptitle, size = 18)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
        save_file_path = analysis_path + "/" + suptitle + "_meanzs_lick" + ".png"
        plt.savefig(save_file_path)
        save_file_path = analysis_path + "/" + suptitle + "_meanzs_lick" + ".pdf"
        plt.savefig(save_file_path)
        plt.show()  
        
        m = m+1


    
#########################################################
#      SAVE data data as a pkl file
################################################################
        # print ("%s" %(data.at[i, "mouseID"] + "__" + str(data.at[i, "LICK_on"][1]) + "__" + str(data.at[i, "LICK_on"][10])))
        suptitle = data.at[i, "mouseID"] + "_" + data.at[i, "date"] + "_" + data.at[i, "state"] + "_" + data.at[i, "access"]        
        plt.suptitle("%s" %suptitle, fontsize = 18)        
        save_file_path  = analysis_path + "/" + suptitle + ".pkl"
        output = open(save_file_path, 'wb')
        pk.dump(data.iloc[i], output)
        output.close()  
     
#%%##############################################################
# #       SAVE data files 
# ################################################################# 
    
 
data_to_save = pd.DataFrame(columns = ["expID", # "age", 
                                       "lick_meanzs", "bout_meanzs", "interbout_meanzs", "bite_latency_sec"])

for i, expID in enumerate(expID_list):
    data_to_save.at[i, "expID"] 	          = expID
    # data_to_save.at[i, "age"]                 = data.at[i, "age"]
    data_to_save.at[i, "lick_meanzs"]         = df_lick_zs[i]["lick_meanzs"].mean()
    data_to_save.at[i, "bout_meanzs"]         = df_bout_zs[i]["bout_meanzs"].mean()
    data_to_save.at[i, "interbout_meanzs"]    = df_bout_zs[i]["interbout_meanzs"].mean()
    try: 
        data_to_save.at[i, "bite_latency_sec"]    = data.at[i, "LICK_on"][1]
    except:
        data_to_save.at[i, "bite_latency_sec"]   = np.nan
 
# csv_path         = analysis_path + "/meanlickboutzs" + ".csv"
csv_path         = analysis_path + "/meanbiteboutzs" + ".csv"
data_to_save.to_csv(csv_path, index = True)   



df0 = pd.DataFrame()
for i in range(len(expID_list)):
    df0[data.at[i, "expID"]] = pd.Series(df_bout_zs[i]["lick_meanzs"])

    
stimulus = "chowbite"
csv_path              = analysis_path + "/meanBITEzs" + "_" + stimulus + ".csv"
df0.to_csv(csv_path, index = True)
    

df1 = pd.DataFrame()
for i in range(len(expID_list)):
    df1[data.at[i, "expID"]] = pd.Series(df_bout_zs[i]["bout_meanzs"])
    
csv_path              = analysis_path + "/meanBOUTzs" + "_" + stimulus + ".csv"
df1.to_csv(csv_path, index = True)
    
df2 = pd.DataFrame()
for i in range(len(expID_list)):
    df2[data.at[i, "expID"]] = pd.Series(df_bout_zs[i]["interbout_meanzs"])
    
csv_path              = analysis_path + "/meanINTERBOUTzs" + "_" + stimulus + ".csv"
df2.to_csv(csv_path, index = True)
            


    
    
# ##################################################
# #   POP open the analysis folder
# ##################################################
os.startfile(analysis_path)
