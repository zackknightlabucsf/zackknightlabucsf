%% Davis Analysis: with TTL adapter
%Don't run data through mdata_readmlist because will do z-score based on
%baseline before each stim. Will do that here
%Then run them through this script which will downsample the TTL file, find
%a single point for when the mouse licked (take the first one). 
%Next, cut the all_mrerank data around this TTL signal and re-baseline with
%z-score formula ((x-mean_b)/std_b).
%Concatenate data into a final storage cell called psth, where each cell is
%the neural activity for each presentation (all mice here). 
%Then can separate based on act/inhib/none and plot

%% Input stuffies
close all; clear all; 
txtfiles = { %txt files as you would with readmlist
    'TA120_fastedSucroseDavis2_1213.txt'
    'TA121_fastedSucroseDavis2_1213.txt'
    'TA124_fastedSucroseDavis2_1213.txt'
    'TA130_fastedSucroseDavis2_1213.txt'
    
    'TA120_fastedSucraloseDavis2_1215.txt'
    'TA121_fastedSucraloseDavis2_1215.txt'

};

infile = { %csv files
    '2021-12-13-11-30-04_gpio.csv'
    '2021-12-13-12-58-06_gpio.csv'
    '2021-12-13-14-17-10_gpio.csv'
    '2021-12-13-15-29-05_gpio.csv'
    
    '2021-12-15-11-12-10_gpio.csv'
    '2021-12-15-12-40-28_gpio.csv'

};

Davis_file = { %output file from Davis converted to csv
    'TA120_fastedSucrose2_1213.csv'
    'TA121_fastedSucrose2_1213.csv'
    'TA124_fastedSucrose2_1213.csv'
    'TA130_fastedSucrose2_1213.csv'
    
    'TA120_fastedSucralose2_1215.csv'
    'TA121_fastedSucralose2_1215.csv'
};

OpenTrials = { %which trials required a retry?
    [12]
    [30, 34]
    [34]
    [17]
    
    [8 24 27]
    [22 27 30 31]
    };
baseline = 30; %how many seconds pre and post tastant do you want to look at?

psth_all = []; psth_all_act=[]; 
for i = 1:length(infile)
    %% Read in TTLs
    Times = readtable(infile{i}); %GPIO starts at row 35
    Times = table2array(Times(35:end,1:2:3)); %just take the first column (timestamps)
    thr = max(Times(:,2))/2; %some files have small values, filter out noise
    cTimes = Times(find(Times(:,2)>thr),:); %cTimes = clean Times
    
    %now look for unique TTLs (the first in the series); should be 36 total
    %unless the mouse didn't lick one of them
    Lick_start = [cTimes(1)]; %storage vector starting with first TTL
    for j = 2:size(cTimes, 1)
        if cTimes(j) > cTimes(j-1)+30 %should be more than 30 seconds between tastants
            Lick_start = [Lick_start cTimes(j)];
        end
    end
    %% Downsample/filter raw data
    %use code from readmlist but don't z-score since we'll baseline to
    %pre-tastant for each of them
    zscoreyn = 0; rawyn = 1; txtfile = txtfiles{i};
    [ cC, mpsthnom,c_psthnom,all_mrerank,ind_mrerank,crerank,ind_crerank,c_corrm,para,prepostrange,cC_rerank] = cnmfe_yc_mdatalowpassfilter_readmtxt_decimateto1(txtfile,zscoreyn,rawyn);
    
    %% Take the PSTH and normalize the pre-lick baseline
    psth = []; %store each tastant as a separate cell
    tastant = []; %temp storage vector while cycle throught the neurons
    for k = 1:length(OpenTrials{i}) %go through each tastant
        for j = 1:size(all_mrerank,1) %go through every neuron
            try
                here = Lick_start(OpenTrials{i}(k))-60; %go before that lick
                a = all_mrerank(j, here-baseline:here+baseline); %cut out the times for that neuron
            catch
                here = Lick_start(OpenTrials{i}(k))-60;
                a = all_mrerank(j, here-baseline:end); %if not enough at end
            end
            za = (a-mean(a(1:baseline)))./std(a(1:baseline));
            tastant(j,:) = za; a = []; za = []; %store and reset
        end
        psth =[psth; tastant];
        tastant = []; %reset
    end
    psth_all = [psth_all; psth];
    %% Sort neurons based on act/inhib/none
    [psth_act]=Heatmap_quant_TA(psth_all, 'Title', baseline, baseline*2);

    %% Put into a storage vector that makes it easier to put into prism (put all mice together)
    psth_all_act = [psth_all_act; psth_act];
    
    clearvars -except infile txtfiles baseline psth_all_act i  ...
        psth_all Davis_file OpenTrials Davis_file
end
