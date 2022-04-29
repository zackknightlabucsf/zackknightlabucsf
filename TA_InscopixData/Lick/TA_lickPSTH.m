%% TA Inscopix Lick PSTH
% Some code taken from ttl2psth_for_inscopix_analysis_JG and TA_lickTTL
% Calculates the PSTH for each licking bout, separated by neurons that are
% Act/Inh/No change. Note that the error bars are not correct here!
close all; clear all; 
hw = 5; %bout threshold in sec 
baseline = 10; %how many seconds pre/post bout do you want to look at?
prism = 1; % Put 1 if need the mean and sem output for psth. 
exclude_first = 1; %put 1 if exclude first 3 minutes of session (based on water trial)
txtfiles = { %as you would in readmlist
    'TA68_fastedWater_011921.txt'
    'TA80_fastedWater_011921.txt'
    'TA81_fastedWater_011921.txt'
    'TA82_fastedWater_011921.txt'
    'TA83_fastedWater_011921.txt'
    
    %'TA81_fastedSucrose2_0211.txt'
    %'TA82_fastedSucrose2_0211.txt'
    %'TA83_fastedSucrose2_0211.txt'
    
    %'TA81_fastedSucraloseLow_0216.txt'
    %'TA82_fastedSucraloseLow_0216.txt'
    %'TA83_fastedSucraloseLow_0216.txt'
    
    %'TA59_fastedEnsure_0211.txt'
    %'TA81_fastedEnsure_1026.txt'
    %'TA82_fastedEnsure_1026.txt'
    %'TA83_fastedEnsure_1026.txt'
    
    %'TA68_Sucralose.txt'
    %'TA82_fastedSucraloseLick3_1201.txt'
    %'vLepR83_fastedSucraloseLick2.txt'
};

infiles = { %csv files
    '2021-01-19-16-40-26_video_gpio.csv'
    '2021-01-19-15-23-29_video_gpio.csv'
    '2021-01-19-14-00-39_video_gpio.csv'
    '2021-01-19-14-06-11_video_gpio.csv'
    '2021-01-19-15-22-32_video_gpio.csv'
    
    %'2021-02-11-14-40-14_video-gpio.csv'
    %'2021-02-11-14-46-13_video-gpio.csv'
    %'2021-02-11-15-44-18_video-gpio.csv'
    
    %'2021-02-16-13-38-18_video-gpio.csv'
    %'2021-02-16-13-41-20_video-gpio.csv'
    %'2021-02-16-14-40-04_video-gpio.csv'
    
    %'2020-02-11-10-50-03_licks.csv'
    %'TA81_EnsureLicks_1026.csv'
    %'TA82_EnsureLicks_1026.csv'
    %'TA83_fastedEnsureLick_csv.csv'
    
    %'TA68_Sucralose_0804.csv'
    %'TA82_SucraloseLicks3_1201.csv'
    %'TA83_fastedSucraloseLick2_1128_gpio.csv'
};
for i = 1:length(infiles)
%% Get TTLs from csv file
[ttl] = csv_reshape_for_inscopix_analysis_JG(infiles{i},'GPIO-1', 4, 9600);
temp=ttl;
%check for NaNs
positionNaN = isnan(temp);
temp(positionNaN) = 0;
    ttl_new= decimate(temp,4);
    ttl_logicNew = (ttl_new>2000);
    ttl_logic = (ttl>2000);
    ttl_logicNew = ttl_logicNew*20;
    for j = 1:length(ttl_logicNew)
        if ttl_logicNew(j) == 0
            ttl_logicNew(j) = -5;
        end
    end
%% Finding licking bouts

ttl_thr = (max(ttl_logicNew)+min(ttl_logicNew))/2;
ttl_time=find(ttl_logicNew>ttl_thr);
bouts={ttl_time(1)}; %first lick is first bout

%Now sorting individual TTLs into bouts
%half window is half of the graphing window size
for j=2:length(ttl_time)
    if ttl_time(j)-ttl_time(j-1)>hw*2 %if ttl on is further than hw*2 from last ttl
        bouts{length(bouts)+1}=ttl_time(j); %add this ttl
        %         keyboard();
    else
        if length(bouts)>=1
            bouts{length(bouts)}=[bouts{length(bouts)} ttl_time(j)]; 
        end
        %this ends with cell of frames for each bout
    end
end

%delete bouts under threshold
% blength=cellfun(@(x) x(end)-x(1),bouts); %find how long bouts last
% eind=find(blength<hw);
% bouts(eind)=[]; %remove those under threshold

%% Sort bouts based on time in session

if exclude_first == 1
    bouts_old = bouts; %store the old one because size will change
    bouts ={};
    for k = 1:size(bouts_old, 2)
        if bouts_old{k}(1) < 780 %remove bouts in first 3 minutes if choose to
            bouts{k} = [];
        else
            bouts{k} = bouts_old{k}; %if delete as I go, size will change
        end
    end
    bouts = bouts(~cellfun('isempty', bouts));
end
bout_period = []; %create a vector that stores which time period each bout happens in.
bout_length = []; %create a storage vector for the number of licks of each bout
for j = 1:length(bouts)
    if bouts{j}(1) < 1200
        bout_period(end+1) = 1; %1 for first ten minutes
        bout_length(end+1) = length(bouts{j});
    elseif bouts{j}(1) >=1200 && bouts{j}(1) < 1800
        bout_period(end+1) = 2; %second ten minutes
        bout_length(end+1) = length(bouts{j});
    else
        bout_period(end+1)=3; %last/third ten minutes
        bout_length(end+1) = length(bouts{j});
    end
end

%% Downsample and filter raw data
    %use code from readmlist but don't z-score since we'll baseline to
    %pre-bout for each of them
    zscoreyn = 0; rawyn = 1; txtfile = txtfiles{i};
    [ cC, mpsthnom,c_psthnom,all_mrerank,ind_mrerank,crerank,ind_crerank,c_corrm,para,prepostrange,cC_rerank] = cnmfe_yc_mdatalowpassfilter_readmtxt_decimateto1(txtfiles{i},zscoreyn,rawyn);

%% Cut neural activity around bout start and z-score to pre-bout
psth_all = {}; %store each bout as a separate cell
    for k = 1:length(bouts) %go through each tastant
        for j = 1:size(all_mrerank,1) %go through every neuron
            try
                a = all_mrerank(j, bouts{k}-baseline:bouts{k}+baseline); %cut out the times for that neuron
                za = (a-mean(a(1:baseline)))./std(a(1:baseline));
                psth_all{k}(j,:) = za; a = []; za = []; %store and reset
            catch
                a = all_mrerank(j, bouts{k}-baseline:end); %if not enough at end
                za = (a-mean(a(1:baseline)))./std(a(1:baseline));
                psth_all{k}(j,:) = za; a = []; za = []; %store and reset
            end
        end
    end

%% Sort two ways: act/inhib and based on bout time in session
% store psths in psth.first, psth.second, psth.third based on the bout_period indices
% determined earlier. Separately sort into psth.act psth.none psth.inhib.
% Compile all mice into here. 

if i == 1
    psth.first = [];psth.second=[];psth.third=[]; %if this is the first mouse, setup variables
    psth.act=[];psth.inh=[];psth.none=[];
end
%first sort based on bout period
for k = 1:length(bout_period)
   if bout_period(k) == 1
       psth.first = [psth.first; psth_all{k}];
   elseif bout_period(k) == 2
       psth.second = [psth.second; psth_all{k}];
   else
       try
        psth.third = [psth.third; psth_all{k}];
       catch %if dimensions don't agree (i.e. bout at very end)
        while size(psth_all{k},2)<baseline*2+1
            psth_all{k}(:,end+1) = NaN; %pad the end
        end
           psth.third = [psth.third; psth_all{k}];
       end
   end
end

%now sort based on activity 
act1 = []; none1=[]; inhib1=[];act=[];none=[];inhib=[];
for k = 1:size(psth_all,2)
        [act1, none1, inhib1]=Heatmap_quant_TAPSTH(psth_all{k}, 'Title', baseline, baseline*2);
        act = [act; act1]; none = [none; none1]; inhib = [inhib; inhib1];
end
psth.act = [psth.act; act]; psth.inh = [psth.inh; inhib]; psth.none = [psth.none; none];

%% Care most about activated neurons throughout session, so sort based on both
psth_act_sorted.first = []; psth_act_sorted.second = []; psth_act_sorted.third = [];
auc_act_sorted.first = []; auc_act_sorted.second = []; auc_act_sorted.third = [];
for k = 1:size(psth_all,2)
        [act2]=Heatmap_quant_TAPSTH(psth_all{k}, 'Title', baseline, baseline*2);
        if bout_period(k)==1
            psth_act_sorted.first = [psth_act_sorted.first; act2];
            temp_auc = trapz(act2(:,baseline:baseline*2),2)/bout_length(k); %normalize to number of licks in bout
            auc_act_sorted.first = [auc_act_sorted.first; temp_auc];
        elseif bout_period(k)==2
            psth_act_sorted.second = [psth_act_sorted.second; act2];
            temp_auc = trapz(act2(:,baseline:baseline*2),2)/bout_length(k);
            auc_act_sorted.second = [auc_act_sorted.second; temp_auc];
        else
            psth_act_sorted.third = [psth_act_sorted.third; act2];
            temp_auc = trapz(act2(:,baseline:baseline*2),2)/bout_length(k);
            auc_act_sorted.third = [auc_act_sorted.third; temp_auc];
        end
end

%% Reset for next mouse

clearvars -except exclude_first prism infiles txtfiles baseline hw psth auc_act_sorted psth_act_sorted
end

%% Calculate AUC while we're at it

psth.first_auc = [];psth.second_auc=[];psth.third_auc=[]; 
psth.act_auc=[];psth.inh_auc=[];psth.none_auc=[];

if psth.first ~=0
    psth.first_auc = trapz(psth.first(:, baseline:baseline*2),2);
end
if psth.second ~=0
    psth.second_auc = trapz(psth.second(:, baseline:baseline*2),2);
end
if psth.third ~=0
    psth.third_auc = trapz(psth.third(:, baseline:baseline*2),2);
end
if psth.act ~=0
    psth.act_auc = trapz(psth.act(:, baseline:baseline*2),2);
end
if psth.none ~=0
    psth.none_auc = trapz(psth.none(:, baseline:baseline*2),2);
end
if psth.inh ~=0
    psth.inh_auc = trapz(psth.inh(:, baseline:baseline*2),2);
end

%% Run this code for prism. PSTH is too long/large so calculate mean and SEM here
if prism ==1
    prism_SEM.first = [];prism_SEM.second=[];prism_SEM.third=[];prism_SEM.act=[];...
        prism_SEM.none=[];prism_SEM.inhib=[]; prism_SEM.firstAct = [];...
        prism_SEM.secondAct=[]; prism_SEM.thirdAct=[];%store all the SEMs
    
    prism_mean.first= [];prism_mean.second=[];prism_mean.third=[];prism_mean.act=[];...
        prism_mean.none=[];prism_mean.inhib=[];prism_mean.firstAct=[];...
        prism_mean.secondAct=[];prism_mean.thirdAct=[];%store all the means

    prism_SEM.firstAct = std(psth_act_sorted.first,0,1)./sqrt(size(psth_act_sorted.first,1));
    prism_mean.firstAct = mean(psth_act_sorted.first,1);
    
    prism_SEM.secondAct = std(psth_act_sorted.second,0,1)./sqrt(size(psth_act_sorted.second,1));
    prism_mean.secondAct = mean(psth_act_sorted.second,1);
    
    prism_SEM.thirdAct = std(psth_act_sorted.third,0,1)./sqrt(size(psth_act_sorted.third,1));
    prism_mean.thirdAct = mean(psth_act_sorted.third,1);
    
    prism_SEM.first = std(psth.first,0,1)./sqrt(size(psth.first,1));
    prism_mean.first = mean(psth.first,1);
    
    prism_SEM.second = std(psth.second,0,1)./sqrt(size(psth.second,1));
    prism_mean.second = mean(psth.second,1);
    
    prism_SEM.third = std(psth.third,0,1)./sqrt(size(psth.third,1));
    prism_mean.third = mean(psth.third,1);
    
    prism_SEM.act = std(psth.act,0,1)./sqrt(size(psth.act,1));
    prism_mean.act = mean(psth.act,1);
    
    prism_SEM.none = std(psth.none,0,1)./sqrt(size(psth.none,1));
    prism_mean.none = mean(psth.none,1);
    
    prism_SEM.inhib = std(psth.inh,0,1)./sqrt(size(psth.inh,1));
    prism_mean.inhib = mean(psth.inh,1);
end