# -*- coding: utf-8 -*-
"""
Created on Tue May 26 23:03:32 2020

Fiber photometry group analysis from processed z-score (.pkl) files.

Groups recordings by mouse ID and age, averaging repeated recordings
per mouse before calculating group-level neural activity.

Calculates z-score AUC, peak responses, and half-rise times.
For chow experiments, also normalizes signals and AUC to chow intake.

Generates age-group mean traces (± SEM), individual-recording heatmaps,
and CSV summaries of recording- and mouse-level measurements.

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
from datetime import date
from scipy.integrate import simps
import seaborn as sns
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


print("Enter the full path for directory with .pkl files from rawplot: ")
path      = input()

def main(path):
# if 1==1:

    
    # add_to_csv = path.split("\\")[-2]
    add_to_csv = path.split("\\")[-2] + "_" + path.split("\\")[-1]
    
    ##  Import DOB data file into a dictionary
    df_DOB    = pd.read_csv(r"D:\Data_Dropbox_Sync\Mouse_metadata\allmiceDOB.csv")
    
    ##  Import chow_g and BW data file into a dataframe
    if "chow" in path:
        chow = 1 # 1 if stimulus is chow, 0 if not chow
        df_chow_g = pd.read_csv(r"D:\Data_Dropbox_Sync\Mouse_metadata\photometrymiceTTL_wchow.csv")

    else:
        chow = 0
        
    age_breakpoints = {}
    grouped         = {}
    group_key_list  = {}
    meanzscore         = {}
     
    pre_min          = 10 # minutes
    post_min         = 45 # minutes

    age_breakpoints[0] = [0, 20, 32]

    
    min_start_for_AUC     = pre_min         ## edit here to set starting minute for AUC. Default is pre-length
    min_end_for_AUC       = post_min 
    
    sampling_rate         = 1017.2527

    ds_factor1            = 5
    ds_factor2            = 5
    ds_factor3            = 5
    ds_factor             = int(ds_factor1*ds_factor2*ds_factor3)   
    sm_factor             = 101
        
    today = datetime.today().strftime('%Y%m%d')
    analysis_path = path + "/groupAnalysis_"  + today  + "_nonTTLscript15_risetime"       
    if not os.path.exists(analysis_path):
        os.makedirs(analysis_path)            
    
    
    # this is to exclude small peak from injection itself  
    frame_start_for_AUC  = int(min_start_for_AUC*60*sampling_rate/ds_factor) 
    frame_end_for_AUC    = int((min_start_for_AUC+min_end_for_AUC)*60*sampling_rate/ds_factor)
    
    ## y range for individual traces and mean
    fixed_ylim_ind   = 0 # 1 if yes, 0 if no 
    if fixed_ylim_ind  == 1:
        fixed_ylim_ind_range = [-2, 10]  ## Edit here
    
    ## y range for mean and ste 
    fixed_ylim_mean    = 1
    if fixed_ylim_mean  == 1:
        fixed_ylim_mean_range = [-3, 6]  ## Edit here    
        

    age_palette           = ["black", "maroon"]
    age_palette_ste       = ["silver", "pink"]    
    
        
    all_files             = glob.glob(path + "/*.pkl")
    file_list             = []
    expID_list            = []
    # min_length_list       = []
    # 
    ###################################################################################################
    #   Input: entire pd.DataFrame from all experiments of a single mouse saved as .pkl file
    ###################################################################################################   
    
    if chow == 1:
        data = pd.DataFrame(columns = ["expID", "mouseID", "expdate", "BW",
                                       "sex", "stimulus", "age_wk", "age_mo", "zscore",
                                       "zscore_ds", "zscore_ds_AUC", "zscore_ds_peak", "time_sec_ds",
                                       "chow_g", "zscore_ds_chow_g", "zscore_ds_chow_BW", 
                                       "zscore_ds_AUC_chow_g", "zscore_ds_AUC_chow_BW", "T_halfrise_sm"])
    
        data["zscore_ds_chow_g"]    = data["zscore_ds_chow_g"].astype(object)
        data["zscore_ds_chow_BW"]   = data["zscore_ds_chow_BW"].astype(object)

    else:    
        data = pd.DataFrame(columns = ["expID", "mouseID", "expdate", "BW",
                                       "sex", "stimulus", "age_wk", "age_mo", "zscore",
                                   "zscore_ds", "zscore_ds_AUC", "zscore_ds_peak", "time_sec_ds", "T_halfrise_sm"])
    
    data["zscore"]            = data["zscore"].astype(object) 
    data["zscore_ds"]         = data["zscore_ds"].astype(object)    
    data["time_sec_ds"]       = data["time_sec_ds"].astype(object) 
       
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
     
        data.at[i, "zscore_ds"]     = tempdata["zscore_ds"]
        data.at[i, "time_sec_ds"]   = tempdata["time_sec_ds_zs"]
        
        # find DOB from the metadata file     
        for j in range(len(df_DOB)):
            if df_DOB.iloc[j]["mouseID"] == mouseID:
                DOB_temp  = str(df_DOB.iloc[j]["DOB"])
                break   
        
        try: 
            DOB_temp
        except:
            print ("**************%s has no DOB" %expID)
            continue
        
        expdate_temp = tempdata["date"]
        
        # find mouse age based on DOB and experiment date        
        DOB                   = date(int(DOB_temp[0:4]), int(DOB_temp[4:6]), int(DOB_temp[6:8])) 
        expdate               = date(int(expdate_temp[0:4]), int(expdate_temp[4:6]), int(expdate_temp[6:8])) 
        age_days              = (expdate - DOB).days  
    
        data.at[i, "age_wk"]      = math.floor(age_days/7)
        data.at[i, "age_mo"]      = math.floor(age_days*0.032876676)
        # print ("%s has age_wk %d and age_mo %d" %(expID, data.at[i, "age_wk"], data.at[i, "age_mo"]))
     
        #%##################################################################################
        #   trim z score data by pre-min and post-min 
        ######################################################################################  
        data.at[i, "zscore_ds"]      = data.at[i, "zscore_ds"][int((10-pre_min)*60*sampling_rate/ds_factor):int((10+post_min)*60*sampling_rate/ds_factor)] 
        data.at[i, "time_sec_ds"]    = data.at[i, "time_sec_ds"][int((10-pre_min)*60*sampling_rate/ds_factor):int((10+post_min)*60*sampling_rate/ds_factor)] 

        data.at[i, "zscore_ds_AUC"]     = simps(data.at[i, "zscore_ds"][frame_start_for_AUC:frame_end_for_AUC].clip(min=0), dx=1)
        if stim == "AngII":
            frame_start_for_peak        = int(pre_min*60*sampling_rate/ds_factor) # If stimulus is "AngII", count the very first peak
        else:
            frame_start_for_peak        = frame_start_for_AUC
        data.at[i, "zscore_ds_peak"]    = np.max(data.at[i, "zscore_ds"][frame_start_for_peak:frame_end_for_AUC])
    #     
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
    
            data.at[i, "zscore_ds_chow_g"]   = data.at[i, "zscore_ds"]/data.at[i, "chow_g"]
            # data.at[i, "zscore_ds_chow_BW"]  = data.at[i, "zscore_ds_chow_g"]/data.at[i, "BW"]    
    
            data.at[i, "zscore_ds_AUC_chow_g"]  = simps(data.at[i, "zscore_ds_chow_g"][frame_start_for_AUC:frame_end_for_AUC].clip(min=0), dx=1)
            # data.at[i, "zscore_ds_AUC_chow_BW"] = simps(data.at[i, "zscore_ds_chow_BW"][frame_start_for_AUC:frame_end_for_AUC].clip(min=0), dx=1)

         #%%##################################################################################
    #   #   Find time to peak value (rise), half time to baseline (drop)
    #   #####################################################################################    
        zscore_sm                         = smooth(data.at[i, "zscore_ds"], window_len = sm_factor)
        peak_zs                           = np.max(zscore_sm)         
        if ("PEG" in expID):
            peak_zs                       = np.max(zscore_sm[int((pre_min*1.5)*60*sampling_rate/ds_factor):])
            diff                          = np.absolute(zscore_sm[int((pre_min*1.5)*60*sampling_rate/ds_factor):int(pre_min*5*60*sampling_rate/ds_factor)]- peak_zs/2)            
            ind_T_halfrise_sm             = diff.argmin() +int((pre_min*1.5)*60*sampling_rate/ds_factor)            
    
        elif ("2MNaCl" in expID):
            diff                          = np.absolute(zscore_sm[int((pre_min*1.1)*60*sampling_rate/ds_factor):int(pre_min*1.35*60*sampling_rate/ds_factor)]- peak_zs/2)            
            ind_T_halfrise_sm             = diff.argmin() +int((pre_min*1.1)*60*sampling_rate/ds_factor)            
        elif ("AVP" in expID):
            diff                          = np.absolute(zscore_sm[int((pre_min*1.5)*60*sampling_rate/ds_factor):int(pre_min*3*60*sampling_rate/ds_factor)]- peak_zs/2)            
            ind_T_halfrise_sm             = diff.argmin() +int((pre_min*1.5)*60*sampling_rate/ds_factor)   
        else:          
            diff                          = np.absolute(zscore_sm[int((pre_min+1.1)*60*sampling_rate/ds_factor):int(pre_min*1.35*60*sampling_rate/ds_factor)]- peak_zs/2)            
            ind_T_halfrise_sm             = diff.argmin() +int((pre_min+1.1)*60*sampling_rate/ds_factor)
        data.at[i, "T_halfrise_sm"]       = data.at[i, "time_sec_ds"][ind_T_halfrise_sm]
 
        #%###################################################################################################
        #       FIGURE 4.  individual plot: z score, smoothened z score, time to peak
        ####################################################################################################
        
        fig, ax = plt.subplots(1, 1, figsize = (10, 6))

        ax0   = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)   
        ax0.axvline (x = 0, linewidth = 2, linestyle = "--", color = "gray")
        ax0.set_ylabel("z score", fontsize = 20)
        ax0.set_xlabel("time (min)")
        # ax0.set_ylim (min_zs_, max_zs_)

        ax1   = ax0.twinx()    
        ax1.axvline (x = data.at[i, "T_halfrise_sm"]/60 , linewidth = 2, linestyle = "--", color = "red")       
        # ax1.set_ylim (min_zs_, max_zs_)
        ax1.axes.get_yaxis().set_visible(False)

        suptitle = data.at[i,"mouseID"] + "_" + data.at[i,"expdate"] + "_" + str(int(data.at[i,"age_mo"])) + "mo"
        ax0.set_title(suptitle, fontsize = 20)
            
        suptitle = data.at[i,"mouseID"] + "_" + data.at[i,"expdate"] + "_" + str(int(data.at[i,"age_mo"])) + "mo" + "_zs_sm_rise"        
        plt.suptitle("%s" %suptitle, fontsize = 18) 
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])           
        save_file_path = analysis_path + "/" + suptitle + ".png"
        plt.savefig(save_file_path)

        # plt.show()
        plt.close()

    print ("files analyzed are: %s" %file_list)        
    
    
    
    
    #%%#############################################################################################
    #   PROCESS DATA   
    #   1) "Groupby" by "mouseID", save the mean zscore, mean AUC and age into the new dataframe
    #   2) Then "Groupby" by "age break points", save the mean and std into new variables 
    ################################################################################################
    
    #   1) "Groupby" by "mouseID", save the mean zscore, mean AUC and age into the new dataframe
    
    
    if chow == 1:  
        meandata   = pd.DataFrame({"zscore_ds_mean": data.groupby("mouseID")["zscore_ds"].apply(np.mean), 
                               # "zscore_ds_scaled_mean": data.groupby("mouseID")["zscore_ds_scaled"].apply(np.mean),
                               "zscore_ds_mean_chow_g": data.groupby("mouseID")["zscore_ds_chow_g"].apply(np.mean),
                               # "zscore_ds_mean_chow_BW": data.groupby("mouseID")["zscore_ds_chow_BW"].apply(np.mean),
                               "zscore_ds_AUC": data.groupby("mouseID")["zscore_ds_AUC"].apply(np.mean),
                               "zscore_ds_AUC_chow_g": data.groupby("mouseID")["zscore_ds_AUC_chow_g"].apply(np.mean),
                               # "zscore_ds_AUC_chow_BW": data.groupby("mouseID")["zscore_ds_AUC_chow_BW"].apply(np.mean),
                               "zscore_ds_peak": data.groupby("mouseID")["zscore_ds_peak"].apply(np.mean),
                               "chow_g": data.groupby("mouseID")["chow_g"].apply(np.mean),
                               "age_mo_min": data.groupby("mouseID")["age_mo"].apply(np.min),
                               "age_mo_max": data.groupby("mouseID")["age_mo"].apply(np.max),
                               "num_trials": data.groupby("mouseID").size(),
                               "T_halfrise_sm": data.groupby("mouseID")["T_halfrise_sm"].apply(np.mean)
                               }).reset_index()
        meanzscore_chow       = {}
        # meanzscore_chow_BW    = {}
    else:
        meandata   = pd.DataFrame({"zscore_ds_mean": data.groupby("mouseID")["zscore_ds"].apply(np.mean), 
                               # "zscore_ds_scaled_mean": data.groupby("mouseID")["zscore_ds_scaled"].apply(np.mean),
                               "zscore_ds_AUC": data.groupby("mouseID")["zscore_ds_AUC"].apply(np.mean),
                               "zscore_ds_peak": data.groupby("mouseID")["zscore_ds_peak"].apply(np.mean),
                               "age_mo_min": data.groupby("mouseID")["age_mo"].apply(np.min),
                               "age_mo_max": data.groupby("mouseID")["age_mo"].apply(np.max),
                               "num_trials": data.groupby("mouseID").size(),
                               "T_halfrise_sm": data.groupby("mouseID")["T_halfrise_sm"].apply(np.mean)
                               }).reset_index()   
        
    time_sec     = data.at[0, "time_sec_ds"]
    
    #   2) Then "Groupby" by "age break points", save the mean and std into new variables 
    # time_sec_raw
    ## Group the dataframe by age groups
    for i in range(len(age_breakpoints)):
        grouped[i]        = meandata.groupby(pd.cut(meandata["age_mo_min"], np.array(age_breakpoints[i])))
        group_key_list[i] = list(grouped[i].groups.keys())
    #
        meanzscore[i]        = grouped[i]["zscore_ds_mean"].apply(np.mean)
        if chow == 1:
            meanzscore_chow[i]     = grouped[i]["zscore_ds_mean_chow_g"].apply(np.mean)
            # meanzscore_chow_BW[i]  = grouped[i]["zscore_ds_mean_chow_BW"].apply(np.mean)
    

    ################################################################################
    #   Set max and mix for y-axis
    ################################################################################
    maxzscore_list  = [np.max(meandata["zscore_ds_mean"][x]) for x in range(len(meandata))]    
    minzscore_list  = [np.min(meandata["zscore_ds_mean"][x]) for x in range(len(meandata))]
    # maxAUC_list  = [np.max(meandata["zscore_ds_AUC"][x]) for x in range(len(meandata))]
    # maxpeak_list = [np.max(meandata["zscore_ds_peak"][x]) for x in range(len(meandata))]
    
    maxzscore       = np.max(maxzscore_list)
    minzscore       = np.min(minzscore_list)
    # maxAUC       = np.max(maxAUC_list)
    # maxpeak      = np.max(maxpeak_list)
    
    if fixed_ylim_ind == 1:
        minzscore = fixed_ylim_ind_range[0]
        maxzscore = fixed_ylim_ind_range[1]
    

    
    
    #%%###############################################################################
    #   FIGURE 1. zscore plot with mean + stderr by 1) three age groups 2) two age groups
    ################################################################################
    
    ################################################################################
    #   Set max and mix for y-axis
    ################################################################################
    maxzscore_list  = []
    
    for l in range(len(age_breakpoints)):
        ## age interval [0], [1], [2] needs to be called as [-3], [-2], [-1]
        ind_list   = [-i-1 for i in range(len(grouped[l]))]
        ind_list.sort(reverse = False)
        
        for i, k in enumerate(ind_list):   
            size        = grouped[l].size()[k]        
            zscore_values  = grouped[l].get_group(group_key_list[l][i])["zscore_ds_mean"].values    
            
            if size > 1:
                ste         = zscore_values.std()/math.sqrt(size-1)
            else:
                ste         = [np.nan]*len(zscore_values.mean())
    
            maxzscore_list.append(np.max(zscore_values.mean()) + np.max(ste))
        
    maxzscore       = np.nanmax(maxzscore_list)
    
    if fixed_ylim_mean == 1:
        minzscore = fixed_ylim_mean_range[0]
        maxzscore = fixed_ylim_mean_range[1]
    
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
               
            zscore_values  = grouped[l].get_group(group_key_list[l][i])["zscore_ds_mean"].values    
    
            if size > 1:
                ste     = zscore_values.std()/math.sqrt(size-1)
            else:
                ste     = [np.nan]*len(zscore_values.mean())
    
            color   = age_palette[i] 
            legend      = stim + " age: %d-%d mo (n=%d)" %(minage, maxage, size) 
        #    markers, caps, bars = plt.errorbar(time_sec/60, meanzscore[k], yerr = ste, fmt='-', linewidth = 4, ecolor = "dimgray", elinewidth = 1.5)       
        #    markers.set_color(color)
        #    [bar.set_alpha(0.5) for bar in bars]
            plt.plot(time_sec/60, zscore_values.mean(), linewidth = 2, color = color)
    
            plt.fill_between(time_sec/60, zscore_values.mean() - ste, zscore_values.mean() + ste, color = age_palette_ste[i], alpha = 0.5)
            plt.axvline(x = 0, linewidth = 2, linestyle = "--", color = "gray")
            plt.title(legend, fontsize = 20)
            plt.tick_params(axis = "both", which = "major", labelsize = 18)
            plt.xlabel("time (min)", fontsize = 18)
            plt.ylabel("zscore", fontsize = 18)
            plt.ylim(minzscore*1.05, maxzscore*1.05)
            plt.xlim(-pre_min, post_min)
            # plt.ylim(-15, 75)
            # plt.ylim(-3, 8)
            # plt.xlim(-pre_min, post_min)
            # plt.xlim(-5,30)
    # 
    suptitle = stim + "_" + add_to_csv + "_zscore" 
    plt.suptitle("%s" %suptitle, fontsize = 20)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])     
    if fixed_ylim_mean == 1:
        plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim.png")
        plt.savefig(analysis_path + "/" + suptitle + "_mean_ylim.pdf", dpi = 1200, format = 'pdf', bbox_inches = 'tight')   
    
    else:
        plt.savefig(analysis_path + "/" + suptitle + "_mean.png")
        plt.savefig(analysis_path + "/" + suptitle + "_mean.pdf", dpi = 1200, format = 'pdf', bbox_inches = 'tight')
       
    plt.show()    
    

    #%%##################################################################################
    #  HEATMAP   LICK AND BOUT Aligned to the water access. Mark first lick with a dotted line
    ####################################################################################

    ## Re-arrange the dataframe according to the peak value    
    pre_min_heatmap = pre_min
    post_min_heatmap = post_min
    heatmap_min = minzscore
    heatmap_max = maxzscore


    data                = data.sort_values("zscore_ds_peak", ascending = False).reset_index(drop = True)

    y_labels            = []
    df_heatmap          = pd.DataFrame(data.at[0, "zscore_ds"], columns = [data.at[0, "expID"]])
    for i in range (1, len(data)):
        df_heatmap[data.at[i, "expID"]] = pd.Series(data.at[i, "zscore_ds"])
        y_labels.append(data.at[i, "expID"].split("_")[0])
        
    plt.figure(figsize=(20, 10))
    ax        = plt.axes()
    sns.heatmap(df_heatmap.T, cmap = "jet", vmin = heatmap_min, vmax = heatmap_max, ax = ax) # xticklabels = True, yticklabels = True, ax = ax)
    
 
    
    x_labels  = np.arange(-pre_min_heatmap, post_min_heatmap+5, 5)
    x_ticks   = np.linspace(0, len(df_heatmap), len(x_labels))
    
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels)
    # ax.set_yticklabels(y_labels)
    ax.set_xlabel ("Time (min)")
    
  
    
    suptitle = add_to_csv + "_heatmap"              
    plt.suptitle(". first lick | first bout", size = 18)
    
    save_file_path = analysis_path + "/" + suptitle  + ".png"
    figure = ax.get_figure() 
    figure.savefig(save_file_path, dpi=600)
    save_file_path = analysis_path + "/" + suptitle  + ".pdf"
    figure = ax.get_figure() 
    figure.savefig(save_file_path, dpi=600)
    plt.show()

    
    #####################################################
    #   Save data and meandata as .csv files
    #####################################################
    if chow == 1:
        data_to_save     = data[["expID", "zscore_ds_AUC", "zscore_ds_AUC_chow_g", #"zscore_ds_AUC_chow_BW", 
                                 "zscore_ds_peak", "chow_g", "age_mo", "T_halfrise_sm"]]
    else:
        data_to_save     = data[["expID", "zscore_ds_AUC", "zscore_ds_peak", "age_mo", "T_halfrise_sm"]]    
    csv_path         = analysis_path + "/data_" + stim + ".csv"
    data_to_save.to_csv(csv_path, index = True)
    
    if chow == 1:
        meandata_to_save = meandata[["mouseID", "zscore_ds_AUC", "zscore_ds_AUC_chow_g", 
                                     #"zscore_ds_AUC_chow_BW", 
                                     "zscore_ds_peak", "chow_g", "age_mo_min", "age_mo_max", "num_trials", "T_halfrise_sm"]]
        
    else:
        meandata_to_save = meandata[["mouseID", "zscore_ds_AUC", "zscore_ds_peak", "age_mo_min", "age_mo_max", "num_trials", "T_halfrise_sm"]]    
           
    csv_path         = analysis_path + "/meandata_" + stim + "_" + add_to_csv + ".csv"
    meandata_to_save.to_csv(csv_path, index = True)


    
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
