close all
clear all
fclose all


% TA Lick Analysis no Bout Figures. 
% First uses code from James to find the lick TTLs, then I wrote code after
% it to quantify the area under the curve (AUC) during licking. Sorts the
% values based on if the neurons were activated/inhibited/nochange

%this is a list of csv files
list={
    %'licks_fastedEnsure_TA58.csv' %only change this
    %'licks_depWater_TA58'
    %'licks_depWater_TA59'
    %'licks_fedWater_TA58.csv'
    %'licks_fedWater_TA59'
    %'licks_fastedSucralose3_TA82.csv'
    %'licks_fastedChowWater_TA59.csv'
    %'licks_fedWater2_TA59.csv'
    %'licks_fastedLipid_TA81.csv'
    %'licks_fastedSiliconeOil_TA83.csv'
    %'licks_fastedWater_TA82.csv'
    %'licks_fastedSucrose2_TA83.csv'
    'licks_fastedEnsureSpaced_TA121.csv'
};

%if the following don't matter to you, ignore (the function won't call it)
zscoreyn=1;% zscoreyn: 0 - don't zscore; 1 - zscore to baselin; 2 - zscore to whole trial
rawyn=1; %normally set to 1
hw_sec=1; %bout threshold in secs
fr=4; %framerate
Title = 'TA121 EnsureSpaced Licks';
multi = 0; %1 = 2 lickometers, 0 = 1 lickometer
prestim = 600;
poststim = 1800;
n_sort = 0; %sort based on activated or inhibited? 1 = yes, 0 = no
boutthreshold_sec=hw_sec; %bout time in sec (nb. before was frame count)
firstbout=1; %0 - don't include; 1 - include; 2 - firstbout only
crossreffile=0; %0 - not a crossref file; 1 - crossref file
subtract_baseline=0; %0 - subtract baseline; 1 - subtract activity before first lick; 2-subtract activity before first and last lick separately

%% code from JG lick_analysis
crind={};
[cmap]=cbrewer('div', 'RdBu', 256); %need to get cbrewer function
cmap=flip(cmap);

if(zscoreyn==1)
   zscoreyn_words=char('zscored to prestim');
elseif(zscoreyn==2)
    zscoreyn_words=char('zscored to whole trial');
elseif(zscoreyn==0)
    zscoreyn_words=char('not zscored');
end

if(firstbout==1)
   firstbout_words=char('all bouts included');
elseif(firstbout==2)
    firstbout_words=char('first bout only');
elseif(firstbout==0)
    firstbout_words=char('all bouts but first');
end


for i=1:length(list)
    
    maxcv=20;
    mincv=-20;
    if zscoreyn
        maxcv=maxcv/5;
        mincv=mincv/5;
    end
    
%     %make sure the file ends in csv
%     if(endsWith(list{i},'.csv'))
%         csvfile_name=list{i}
%     else
%         csvfile_name=strcat(list{i}, '.csv');
%     end
    
    csvfile_name = list{i};
    %check if file exists
    fid=fopen(csvfile_name);
    if(fid==-1)
       fprintf(strcat(list{i}, ' could not be found. Double check spelling.\n'));
       return
    end
    clear fid
    
    [all_psth_start, all_psth_end,ttl_clean,cpsthrange, pretime, posttime, timestamps]=read_for_inscopix_analysis_JG( csvfile_name, rawyn, hw_sec, zscoreyn, boutthreshold_sec, firstbout, subtract_baseline);
    
    tot_start=[];
    tot_dff=[];
    tot_end=[];
    for(r=1:size(all_psth_start,1))
        tot_start=[tot_start; all_psth_start(r).avgpsth];
        tot_dff=[tot_dff; all_psth_start(r).dff];
        tot_end=[tot_end; all_psth_end(r).avgpsth];
    end
    
    [~, nind]=sort(tot_dff);
    tot_start=tot_start(nind,:);
    tot_end=tot_end(nind,:);
end

%% Analysis of AUC of licks (code by TA)

%for loop to step through in case there's multiple channels

% First, find where the licks are
Licks = find(ttl_clean);

%Find where the licks are consecutive
lick_bouts = {}; %make a cell array to store the bouts
temp = [Licks(1)]; %temporary storage variable to keep track of the bout
for i = 2:length(Licks)
    if Licks(i) == Licks(i-1)+1
        temp = [temp Licks(i)];
    elseif Licks(i) ~= Licks(i-1)+1
        lick_bouts{end+1} = temp;
        temp = [Licks(i)];
    end
    if i == length(Licks) % if at the last bout, include that too
        lick_bouts{end+1} = temp;
    end
end

%Make a cell array to store the calcium values at the lick times, then
%calculate the AUC. It's the same procedure for either sorting based on
%response type or not.

if n_sort == 1 %If you want to sort based on response type
    %First sort the neurons
    all_ca = [];
    [Neurons_act, Neurons_none, Neurons_inhib]=Heatmap_quant_TA(all_psth_start.zC, Title, prestim*fr, poststim*fr);
    %Store the neural responses during the lick times separately
    for i = 1:length(lick_bouts)
        all_ca.act{i} = Neurons_act(:, lick_bouts{i});
        all_ca.none{i} = Neurons_none(:, lick_bouts{i});
        all_ca.inhib{i} = Neurons_inhib(:, lick_bouts{i});
    end
    
    %Calculate the AUC
    AUC = []; %storage vectors for the values where each row is a different cell
    norm_AUC = []; %AUC per lick, so accounting for licks
    for i = 1:length(lick_bouts)
        if size(Neurons_act,1) ~=0 %check that it's not empty
            for j = 1:size(all_ca.act{i},1)
                if length(lick_bouts{i})==1 %if the bout is only 1 lick, just put the y value (area =1*y)
                    AUC.act(j,i) = all_ca.act{i}(j,:);
                    norm_AUC.act(j,i) = all_ca.act{i}(j,:); %since only 1 lick, no need to norm more
                else
                    AUC.act(j,i) = trapz(lick_bouts{i}, all_ca.act{i}(j,:));
                    norm_AUC.act(j,i) = trapz(lick_bouts{i}, all_ca.act{i}(j,:))/length(lick_bouts{i});
                end
            end
        end
        if size(Neurons_none,1) ~= 0
            for j = 1:size(all_ca.none{i},1)
                if length(lick_bouts{i})==1
                    AUC.none(j,i) = all_ca.none{i}(j,:);
                    norm_AUC.none(j,i) = all_ca.none{i}(j,:);
                else
                    AUC.none(j,i) = trapz(lick_bouts{i}, all_ca.none{i}(j,:));
                    norm_AUC.none(j,i) = trapz(lick_bouts{i}, all_ca.none{i}(j,:))/length(lick_bouts{i});
                end
            end
        end
        if size(Neurons_inhib,1) ~=0
            for j = 1:size(all_ca.inhib{i},1)
                if length(lick_bouts{i})==1
                    AUC.inhib(j,i) = all_ca.inhib{i}(j,:);
                    norm_AUC.inhib(j,i) = all_ca.inhib{i}(j,:);
                else
                    AUC.inhib(j,i) = trapz(lick_bouts{i}, all_ca.inhib{i}(j,:));
                    norm_AUC.inhib(j,i) = trapz(lick_bouts{i}, all_ca.inhib{i}(j,:))/length(lick_bouts{i});
                end
            end
        end
    end
        %Sort based on the mean time this occured.
    act_AUC.first=[]; act_AUC.second = []; act_AUC.third = [];
    none_AUC.first = []; none_AUC.second = []; none_AUC.third = [];
    inhib_AUC.first = []; inhib_AUC.second = []; inhib_AUC.third = [];
    
    nact_AUC.first=[]; nact_AUC.second = []; nact_AUC.third = [];
    nnone_AUC.first = []; nnone_AUC.second = []; nnone_AUC.third = [];
    ninhib_AUC.first = []; ninhib_AUC.second = []; ninhib_AUC.third = [];
    for i = 1:length(lick_bouts)
        %divide by 240 (60s*fr = 4)to convert to minutes
        value = mean(lick_bouts{i})/(60*fr);
        if value>=10 && value <20
            act_AUC.first = [act_AUC.first AUC.act(:,i)];
            none_AUC.first = [none_AUC.first AUC.none(:,i)];
            inhib_AUC.first = [inhib_AUC.first AUC.inhib(:,i)];
            
            nact_AUC.first = [nact_AUC.first norm_AUC.act(:,i)];
            nnone_AUC.first = [nnone_AUC.first norm_AUC.none(:,i)];
            ninhib_AUC.first = [ninhib_AUC.first norm_AUC.inhib(:,i)];
        elseif value>=20 && value <30
            act_AUC.second = [act_AUC.second AUC.act(:,i)];
            none_AUC.second = [none_AUC.second AUC.none(:,i)];
            inhib_AUC.second = [inhib_AUC.second AUC.inhib(:,i)];
            
            nact_AUC.second = [nact_AUC.second norm_AUC.act(:,i)];
            nnone_AUC.second = [nnone_AUC.second norm_AUC.none(:,i)];
            ninhib_AUC.second = [ninhib_AUC.second norm_AUC.inhib(:,i)];
        elseif value>=30 && value <40
            act_AUC.third = [act_AUC.third AUC.act(:,i)];
            none_AUC.third = [none_AUC.third AUC.none(:,i)];
            inhib_AUC.third = [inhib_AUC.third AUC.inhib(:,i)];
            
            nact_AUC.third = [nact_AUC.third norm_AUC.act(:,i)];
            nnone_AUC.third = [nnone_AUC.third norm_AUC.none(:,i)];
            ninhib_AUC.third = [ninhib_AUC.third norm_AUC.inhib(:,i)];
        end
    end
    %Take the mean across each cell within each time point. Make each
    %column a different time point (1st column is 1st 10 min, 2nd is 2nd,
    %etc.)
    mean_act = [mean(act_AUC.first,2) mean(act_AUC.second,2) mean(act_AUC.third,2)];
    mean_none = [mean(none_AUC.first,2) mean(none_AUC.second,2) mean(none_AUC.third,2)];
    mean_inhib = [mean(inhib_AUC.first,2) mean(inhib_AUC.second,2) mean(inhib_AUC.third,2)];
    
    %Take the means for when you look per lick
    norm_mean_act = [mean(nact_AUC.first,2) mean(nact_AUC.second,2) mean(nact_AUC.third,2)];
    norm_mean_none = [mean(nnone_AUC.first,2) mean(nnone_AUC.second,2) mean(nnone_AUC.third,2)];
    norm_mean_inhib = [mean(ninhib_AUC.first,2) mean(ninhib_AUC.second,2) mean(ninhib_AUC.third,2)];
    %Input into prism to plot
        
else %If you want to see all neurons together
    ca = {}; 
    for i = 1:length(lick_bouts)
        ca{i} = all_psth_start.zC(:,lick_bouts{i});
    end
    
    AUC = []; %storage vectors for the values where each row is a different cell
    for i = 1:length(lick_bouts)
        for j = 1:size(all_psth_start.zC,1)
            AUC(j,i) = trapz(lick_bouts{i}, ca{i}(j,:));
        end
    end
    %To help plot the results, sort the data in ten minute bins and average
    %over time
    
    %Sort based on the mean time this occured.
    act_AUC=[];
    for i = 1:length(lick_bouts)
        %divide by 240 (60s*fr = 4)to convert to minutes
        value = mean(lick_bouts{i})/(60*fr);
        if value>=10 && value <20
            act_AUC.first = [act_AUC.first AUC(:,i)];
        elseif value>=20 && value <30
            act_AUC.second = [act_AUC.second AUC(:,i)];
        elseif value>=30 && value <40
            act_AUC.third = [act_AUC.third AUC(:,i)];
        end
    end
    %Take the means across each cell to plot
    mean_first = mean(act_AUC.first,2);
    mean_second = mean(act_AUC.second,2);
    mean_third = mean(act_AUC.third,2);
end


%% Plot! Plotting overall, to plot based on response type, use Prism
if n_sort == 0
    %concatenate into one variable to make it easier
    mean_ca = [mean_first mean_second mean_third];
    X = ones(52,3); X(:,2) = 2; X(:,3) = 3; 

    figure; 
    scatter(X(:,1), mean_ca(:,1)); hold on; 
    scatter(X(:,2), mean_ca(:,2)); hold on; 
    scatter(X(:,3), mean_ca(:,3));
    set(gca, 'XTick', [0 1 2 3]); 
    xlim([0.5 3.5]); 
    hold on; 
    plot([0.5 3.5], [0 0], 'k--');
end


