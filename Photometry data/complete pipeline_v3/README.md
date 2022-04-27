Pipeline V3 for analyzing photometry data
Zack, 01/24/22

This is a more detailed version of text file “photom_pipeline_v2_explanation” that describes using the MATLAB pipeline written by James Grove.  

Overview
The pipeline contains scripts for converting Tucker Davis raw photometry data into MATLAB files, and then analyzing that data to generate various plots. The pipeline is designed to generate plots that show the mean response of multiple animals to the same stimulus (e.g. how a cohort of AgRP photometry mice respond to licking ensure). We will upload additional code that is designed to more efficiently generate plots for individual experiments in the near future. You should always scrutinize individual experiments before performing any grouped analyses.      

There are two main outputs of this pipeline:

1. licks_analysis_V3 aligns photometry data to TTL signals from a lickometer, finds lick bouts based on user defined parameters,and then generates an array of plots that show average responses across a cohort of mice. All of the analyses are plotted twice, once as dF/F and once as z-score. The main ouputs are PSTHs aligned to:  (1) first lick of the first bout,  (2) the last lick of the first bout, (3) the first lick of all bouts, and (4) the last lick of all bouts. In addition, the script generates plots showing a trace of the average response across time; the average lick rate across time; and the average cumulative number of licks across time. 

2.  IP_analysis_V3 is used for experiments that do not contain a lickometer (IP injections, solid food presentation, etc…). It generates mean traces for an experiment aligned to the manipulation time, expressed as both dF/F and z-score. It also generates a bar plot of the mean values during a time window that can be specified in the code. 

Detailed Protocol

A. Getting started. 

1. Install MATLAB on your computer. This is available for free via a UCSF site license. 

https://www.mathworks.com/academia/tah-portal/university-of-california-san-francisco-40716179.html

You will be prompted while running the pipeline to install additional toolboxes. These can be accessed by clicking on Apps, Get More Apps within the MATLAB interface.

2. Download the Photometry Pipeline folder from the lab Dropbox to your computer and unzip. 

3. Add the Photometry Pipeline folder and all its subfolders to the MATLAB search path. You can do this by clicking on Home, Set Path.

4. Transfer the raw photometry data from the lab computer to your computer and place in an appropriately labelled Experiment folder. The simplest way to organize the data is to create one folder per experiment (e.g.  AgRP neurons – Food Deprived – Ensure) and then place within that folder the various folders (“Blocks”) containing the individual photometry experiments (e.g. TL060_TL070-211210-102119).

	Note: By default, the TDT Synapse photometry software saves every experiment performed on a given day within a single folder that is named with the date. This folder is called “a Tank.” Within each Tank are subfolders that contain the data files from an individual experiment. Each subfolder within a Tank is called a “Block”. These Blocks are automically saved when the experiment ends with a name determined by the mouse ID and date (e.g. TL060_TL070-211210-102119). 

	Tank		photometryTTL-211230	
	Block		TL184_TL185-211230-114608
	Data		photometryTTL-211230_TL184_TL185-211230-114608.tsq

5. Navigate to the Experiment folder generated above so that it is the current folder in Matlab. 

B. Generating MATLAB files

1. Open writing_mat_files_fast_v3 in MATLAB.

2. Copy and paste the path for the Experiment folder where indicated in Line 6. 

3. Confirm again that the Experiment folder is the current folder in MATLAB.

4. Click Run in the Editor menu. You should see that one MATLAB file is generated for each Block and saved to the current folder.

C. Generate a metadata spreadsheet

The metadata spreadsheet is an Excel file that contains information about each experiment (described below). 

1. Open make_metadata_file_v2 in MATLAB. 

2. Copy and paste the path for the Experiment folder where indicated in Line 7, and confirm that the Experiment folder is the current folder in MATLAB.

3.  Click Run in the Editor menu. You should see that an Excel file called Metadata is generated and stored in the Experiment folder.

4. Open the Metadata file in Excel, enter all the relevant parameters for your experiment, and then save.

"exp" means the experiment name. This can be whatever you want except:

It must contain the word "licks" if a lickometer was used – otherwise compile_data will fail to load the TTL data and lick_analysis will fail. 

It cannot contain the word “intralipid” or any word with the letters “ip” – otherwise compile_data will classify this as an “ip” experiment and licks_analysis will fail.
 
"Z" means the id of the mouse in the top chamber. e.g. TL070

"K" means the id of the mouse in the bottom chamber

"Notes" is for anything you might want to add. This is usually blank.

"baseline_length" is the duration of time (in seconds) used for the baseline calculation and is assumed to end precisely at “manip_time.” 
Note that this is NOT the total amount of time the mouse was in chamber before the manipulation (which is often ~ 20 min and should always be longer than baseline_length due to the need for habituation). We typically use a baseline_length of 600 (i.e. 10 minutes).

"manip_time" is the timepoint (in seconds) when a manipulation (eg. Injection, lickometer access) was made. This is typically ~1200 (600 s habituation + 600 s baseline) and should be recorded precisely when you perform the experiment. 

"after_manip" is the duration of time (in seconds) the experiment lasts after manipulation was made. This is assumed to begin precisely at “manip_time” and is often 1800 (for a 30 minute experiment).

"downsample_to" is the rate (in Hz) you want to downsample to. This is usually 4 Hz.

Optional parameters (can be left blank; used for lickometer analysis):
"access_length" is the duration of time (in seconds) that access was given for
"access_time" is the timepoint (in seconds) when access was given

These two parameters can be useful when you want to analyze licking behavior for only a subset of the total time that the lickometer was available. For example, you may want to omit the first two minutes of licking data in order to eliminate any confounding effects of the AgRP neuron sensory drop. The file “Metadata example – ZK” in the Pipeline folder illustrates this. 

5. Additional things to note.
	If any of the times differ from mouse in Z vs. K, make separate lines for each

The following parameters must be the same across all the mice in a given experiment:
baseline_length, after_manip, downsample_to

5. The script make_metadata_file_v3 is identical to v2, except that it pre-populates the spreadsheet with some of the most common values for experiments in the lab. You must confirm that these are correct for your experiment and edit as necessary.  

D. Compile the data into a single MATLAB file.

1. Open compile_data_v3 in MATLAB.

2. Copy and paste the path for the Experiment folder where indicated in Line 6, and confirm that the Experiment folder is the current folder in MATLAB.

3. Enter the correct values for the TTL channels in lines 64 and 118 of the code. 
Explanation:  The TTL data is saved under specific channels (“Epocs”) that have different names depending on which photometry system you use. For the first floor systems, these channels are called “Ep2_” and “Ep4_” for systems Z and K, respectively. For the second floor systems, these channels are called “Epo1” and “Epo2” for Z and K.

4. Click Run in the Editor menu. You should see that a MATLAB file called all_photom_data is generated and stored in the Experiment folder.

Troubleshooting:

1. You may receive an error that “the polynomial is poorly conditioned…”. You can ignore this if the individual plots look reasonable.

2. If you get other errors, first check to make sure the metadata spreadsheet is correct. Usually the error is caused by a typo in the metadata spreadsheet, because it is saved in the wrong folder, etc… If necessary go back and regenerate the metadata spreadsheet from scratch.

3. The file all_photom_data contains a structure array named “d”. This array contains within it multiple fields that contain the data and parameters from the experiment. You can inspect this to investigate what was or was not compiled correctly. For example, double click on “d” in the workspace to view the individual fields in d. Double click on “mice” to see the names of the mice that were included. “TTL” contains a list of the TTL events for each mouse (the time in seconds when each lick was registered), etc… 

E. Run licks_analysis to analyze photometry experiments that use a lickometer. 

1. Open licks_analysis_v3 in MATLAB.

2. Confirm that the Experiment folder is the current folder in MATLAB.

3. Set the values in the code where indicated. 

Line 8 – access_time_by_lick 
	0 means that mouse received access for a set amount of time (e.g. 30 minutes)
	1 means that access was a set amount of time after the first lick
	This is usually 0 for most experiments in the lab

Line 12 – rebaseline_to_lick_beginning
	0 means that all bout PSTHs are baselined to the beginning of the entire experiment
	1 means that all bout PSTHs are re-baselined to the beginning of each lick bout
	This is usually 1 for most experiments in the lab.

Line 16 – length of bout
	This parameter (along with ili) determines what constitutes a bout. If this value is too large, there will be too few bouts and the PSTHs will be meaningless. In general 4 sec is a good place start. For other analyses you may want to consider only longer bouts (e.g. 6-10 sec or longer). You should test multiple values and record how many bouts are registered for each mouse to ensure the plots are meaningful.  To record how many bouts are detected, check the field “total_bouts” in “d” after running licks_analysis_v3.  Alternatively, use Licks_analysis_V4, which is the same as V3, except that it also outputs an Excel spreadsheet called “Statistics” that records some statistics for each mouse (including the number of licks and bouts).

4. Click Run in the Editor menu. You should see that a folder is generated called “lick graphs” that has an array of plots as .eps, .png., .html, and .mat files. You can open the .mat files in MATLAB to further edit the figures. 

5. It is a good idea to re-run licks_analysis using different bout sizes to test the robustness of the results to changes in parameters.  If you do this, first change the name of the “lick graphs” folder (e.g. to “lick graphs – 4 sec bouts”). Otherwise the plots will be overwritten with the new analysis.


F. Run IP_analysis to analyze photometry experiments that use a stimulus other than a lickometer (e.g. an IP injection or presentation of an object or food).  

1. Open IP_analysis_v3 in MATLAB.

2. Confirm that the Experiment folder is the current folder in MATLAB.

3. If you want to analyze the mean value during a specific window, set the values where indicated in lines 10-11. For example, if you want to know the mean photometry response from 1 to 10 minutes after the manipulation, enter 60 and 600.  Note that window_end cannot be longer than the total manipulation_length in the metadata spreadsheet. 

4. Click Run in the Editor menu.  A folder called graphs will be generated containing mean traces for the entire experiment (baseline to end of manipulation_length, as specified in the metadata spreadsheet) in 
\both F/F and z-score. It will also generate a bar plot of the mean values during the window indicated. 
