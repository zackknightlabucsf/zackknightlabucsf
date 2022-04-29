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

%1 = Ensure
%2 = Polycal
%3 = Intralipid
%4 = Sucralose
%5 = Salt
%6 = Citric Acid
%7 = Silicone Oil
%8 = Sucrose
%9 = Quinine
%10 = Water

%% Input stuffies
close all; clear all; 
txtfiles = { %txt files as you would with readmlist
    %'TA120_fastedSucroseDavis2_1213.txt'
    %'TA121_fastedSucroseDavis2_1213.txt'
    %'TA124_fastedSucroseDavis2_1213.txt'
    %'TA130_fastedSucroseDavis2_1213.txt'
    
    %'TA120_fastedSucraloseDavis2_1215.txt'
    %'TA121_fastedSucraloseDavis2_1215.txt'
    %'TA124_fastedSucraloseDavis2_1215.txt'
    %'TA130_fastedSucraloseDavis2_1215.txt'
    
    %'TA121_fastedLipidCurve_0331.txt'
    %'TA124_fastedLipidCurve_0331.txt'
    %'TA175_fastedLipidCurve_0401.txt'
    %'TA176_fastedLipidCurve_0401.txt'
    %'TA186_fastedLipidCurve_0401.txt'
    
    'TA121_fastedTastePanel_0405.txt'
    'TA124_fastedTastePanel_0405.txt'
    'TA175_fastedTastePanel_0406.txt'
    'TA176_fastedTastePanel_0406.txt'
};

infile = { %csv files
    %'2021-12-13-11-30-04_gpio.csv'
    %'2021-12-13-12-58-06_gpio.csv'
    %'2021-12-13-14-17-10_gpio.csv'
    %'2021-12-13-15-29-05_gpio.csv'
    
    %'2021-12-15-11-12-10_gpio.csv'
    %'2021-12-15-12-40-28_gpio.csv'
    %'2021-12-15-14-07-20_gpio.csv'
    %'2021-12-15-15-34-28_gpio.csv'
    
    %'2022-03-31-14-05-48_gpio.csv'
    %'2022-03-31-15-25-49_gpio.csv'
    %'2022-04-01-12-00-11_gpio.csv'
    %'2022-04-01-13-17-14_gpio.csv'
    %'2022-04-01-14-30-05_gpio.csv'
    
    '2022-04-05-12-02-07_gpio.csv'
    '2022-04-05-13-20-18_video_green_gpio.csv'
    '2022-04-06-11-34-07_gpio.csv'
    '2022-04-06-12-52-27_gpio.csv'
};

Davis_file = { %output file from Davis converted to csv
    %'TA120_fastedSucrose2_1213.csv'
    %'TA121_fastedSucrose2_1213.csv'
    %'TA124_fastedSucrose2_1213.csv'
    %'TA130_fastedSucrose2_1213.csv'
    
    %'TA120_fastedSucralose2_1215.csv'
    %'TA121_fastedSucralose2_1215.csv'
    %'TA124_fastedSucralose2_1215.csv'
    %'TA130_fastedSucralose2_1215.csv'
};
baseline = 30; %how many seconds pre and post tastant do you want to look at?

perc_act = {}; psth_all = {}; psth_all_act={}; psth_all_none = {}; psth_all_inhib={};
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
    psth = {}; %store each tastant as a separate cell
    tastant = []; %temp storage vector while cycle throught the neurons
    for k = 1:size(Lick_start,2) %go through each tastant
        for j = 1:size(all_mrerank,1) %go through every neuron
            try
                a = all_mrerank(j, Lick_start(k)-baseline:Lick_start(k)+baseline); %cut out the times for that neuron
            catch
                a = all_mrerank(j, Lick_start(k)-baseline:end); %if not enough at end
            end
            za = (a-mean(a(1:baseline)))./std(a(1:baseline));
            tastant(j,:) = za; a = []; za = []; %store and reset
        end
        psth{k}= tastant;
        tastant = []; %reset
    end
    psth_all{i} = psth;
    %% Sort neurons based on act/inhib/none
    psth_act = {}; psth_none = {}; psth_inhib = {};
    num_act = []; total = size(all_mrerank,1);
    for k = 1:size(psth,2)
        [psth_act{k}, psth_none{k}, psth_inhib{k}]=Heatmap_quant_TA(psth{k}, 'Title', baseline, baseline);
        num_act(k) = size(psth_act{k},1);
    end
   
    %% Read in licks per presentation
%     davis = csvread(Davis_file{i}, 11, 8, [11 8 46 8]);
%     time = csvread(Davis_file{i},48,1);
%     I = find(time==0); time(I) =NaN;
    %% Calculate AUC for each neuron during each tastant
    %only look at 5s of licking (at least to start). 5s is too soon/short
    auc_a = {}; auc_i = {}; auc_n = {};
    for j = 1:size(psth_act,2) %just put the values in a new row below the psths
        if psth_act{j} ~=0 
            %             if isnan(time(j,1)) ~= 1 % ignore trials with only 1 lick
            %             auc_a{j} = trapz(psth_act{1,j}(:, baseline:baseline*2),2)./davis(j);
            %             end
            auc_a{j} = trapz(psth_act{1,j}(:, baseline:baseline*2),2);
        end
    end
    for j = 1:size(psth_none,2) %just put the values in a new row below the psths
        if psth_none{j} ~=0 
            %if isnan(time(j,1)) ~= 1
            %auc_n{j} = trapz(psth_none{1,j}(:, baseline:baseline*2),2)./davis(j);
            auc_n{j} = trapz(psth_none{1,j}(:, baseline:baseline*2),2);
            %end
        end
    end
    for j = 1:size(psth_inhib,2)  %just put the values in a new row below the psths
        if psth_inhib{j} ~=0 
            %if isnan(time(j,1)) ~= 1
            %auc_i{j} = trapz(psth_inhib{1,j}(:, baseline:baseline*2),2)./davis(j);
            %end
            auc_i{j} = trapz(psth_inhib{1,j}(:, baseline:baseline*2),2);
        end
    end
    %% Put into a storage vector that makes it easier to put into prism (put all mice together)
    psth_all_act{i} = psth_act;
    psth_all_none{i} = psth_none;
    psth_all_inhib{i} = psth_inhib;
    
    auc_all_a{i} = auc_a;
    auc_all_n{i} = auc_n;
    auc_all_i{i} = auc_i;
    
    perc_act{i} = [num_act'; total];
 
    clearvars -except infile txtfiles baseline psth_all_act psth_all_none...
        psth_all_inhib i auc_all_a auc_all_n auc_all_i psth_all Davis_file perc_act
end

%% Make it easier for putting into prism: put all the mice together
psth_all_prism = {}; psth_act_prism = {}; psth_inhib_prism = {}; psth_none_prism = {};
for j = 1: size(psth_all_act{1},2)
    tempa = []; tempi = []; tempn = []; tempall = [];
    for i = 1:length(infile)
     %tempa = [tempa; cell2mat(psth_all_act{i}(:,j))];    
     %tempi = [tempi; psth_all_inhib{i}(:,j)];
     %tempn = [tempn; psth_all_none{i}(:,j)];
     tempall = [tempall; cell2mat(psth_all{i}(:,j))];
    end
    %psth_act_prism{j} = tempa;
    %psth_inhib_prism{j} = tempi;
    %psth_none_prism{j} = tempn;
    psth_all_prism{j} = tempall; 
end   

%% Now take the average psth for each cell across solution
%Can only use this for situations where you have the same number of cells
%per tastant
psth_all_prism_1 = []; psth_all_prism_2 = []; psth_all_prism_3 = []; psth_all_prism_4 = [];
psth_all_prism_5 = []; psth_all_prism_6 = []; psth_all_prism_7 = []; psth_all_prism_8 = [];
psth_all_prism_9 = []; psth_all_prism_10 = [];
for i = 1:size(psth_all_prism{1},1)
    temp = [psth_all_prism{2}(i,:); psth_all_prism{14}(i,:); psth_all_prism{26}(i,:)];
    psth_all_prism_1 = [psth_all_prism_1; mean(temp,1)];
    
    temp = [psth_all_prism{3}(i,:); psth_all_prism{15}(i,:); psth_all_prism{27}(i,:)];
    psth_all_prism_2 = [psth_all_prism_2; mean(temp,1)];
    
    temp = [psth_all_prism{4}(i,:); psth_all_prism{16}(i,:); psth_all_prism{28}(i,:)];
    psth_all_prism_3 = [psth_all_prism_3; mean(temp,1)];
    
    temp = [psth_all_prism{5}(i,:); psth_all_prism{17}(i,:); psth_all_prism{29}(i,:)];
    psth_all_prism_4 = [psth_all_prism_4; mean(temp,1)];
    
    temp = [psth_all_prism{6}(i,:); psth_all_prism{18}(i,:); psth_all_prism{30}(i,:)];
    psth_all_prism_5 = [psth_all_prism_5; mean(temp,1)];
    
    temp = [psth_all_prism{7}(i,:); psth_all_prism{19}(i,:); psth_all_prism{31}(i,:)];
    psth_all_prism_6 = [psth_all_prism_6; mean(temp,1)];
    
    temp = [psth_all_prism{8}(i,:); psth_all_prism{20}(i,:); psth_all_prism{32}(i,:)];
    psth_all_prism_7 = [psth_all_prism_7; mean(temp,1)];
    
    temp = [psth_all_prism{10}(i,:); psth_all_prism{22}(i,:); psth_all_prism{34}(i,:)];
    psth_all_prism_8 = [psth_all_prism_8; mean(temp,1)];
    
    temp = [psth_all_prism{11}(i,:); psth_all_prism{23}(i,:); psth_all_prism{35}(i,:)];
    psth_all_prism_9 = [psth_all_prism_9; mean(temp,1)];
    
    temp = [psth_all_prism{1}(i,:); psth_all_prism{9}(i,:); psth_all_prism{12}(i,:); ...
        psth_all_prism{13}(i,:); psth_all_prism{21}(i,:); psth_all_prism{24}(i,:); ...
        psth_all_prism{25}(i,:); psth_all_prism{33}(i,:)]; %psth_all_prism{36}(i,:);]; exclude last presentation because session cut off
    psth_all_prism_10 = [psth_all_prism_10; mean(temp,1)];
    
end
    
psth_all_means_prism = {psth_all_prism_1 psth_all_prism_2 psth_all_prism_3 psth_all_prism_4 ...
    psth_all_prism_5 psth_all_prism_6 psth_all_prism_7 psth_all_prism_8 psth_all_prism_9 ...
    psth_all_prism_10};
clear psth_act_prism_1 psth_act_prism_2 psth_act_prism_3 psth_act_prism_4 ...
    psth_act_prism_5 psth_act_prism_6 psth_act_prism_7 psth_act_prism_8 psth_act_prism_9 ...
    psth_act_prism_10

%% Look at specific neurons across tastants

%Quantify are the neurons that respond to sucrose (#8) responding in the same
%direction during sucralose (#4)? 

%quantify direction of each neuron and store in a two column vector
Suc_signs = []; Suc_means = [];
num_same = 0; %counter for how many neurons respond same
for i = 1:size(psth_all_prism_4,1)
    Suclose = mean(psth_all_prism_8(i,31:end));
    Sucose = mean(psth_all_prism_4(i,31:end));
    Suc_means(i,1) = Suclose; Suc_means(i,2) = Sucose;
    if Suclose >1
        Suc_signs(i,1) = 1;
    elseif Suclose <-1
        Suc_signs(i,1) = -1; 
    else
        Suc_signs(i,1) = 0;
    end
    if Sucose >1
        Suc_signs(i,2) = 1;
    elseif Sucose<-1
        Suc_signs(i,2) = -1;
    else
        Suc_signs(i,2) = 0;
    end
    
    if Suc_signs(i,1) == 1 && Suc_signs(i,2)==1
        num_same = num_same+1;
    end
end

%% Now look at how neurons that respond to one tastant respond to other tastants
%Ensure is 2 14 26
%Quinine is 11 23 35
%Look at which neurons are activated during each of the three trials
[psth_a{1}, psth_none{1}, psth_inhib{1}, act{1}]=Heatmap_quant_TA(psth_all_prism{2}, 'Title', baseline, baseline);
[psth_a{2}, psth_none{2}, psth_inhib{2}, act{2}]=Heatmap_quant_TA(psth_all_prism{14}, 'Title', baseline, baseline);
[psth_a{3}, psth_none{3}, psth_inhib{3}, act{3}]=Heatmap_quant_TA(psth_all_prism{26}, 'Title', baseline, baseline);

%Pull out the neurons that are activated during all 3
these = intersect(intersect(act{1}, act{2}),act{3});

%Pull out those neurons from all the presentations
psth_EnsureAct = [];
psth_EnsureAct.Ensure = psth_all_prism_1(these,:);
psth_EnsureAct.Polycal = psth_all_prism_2(these,:);
psth_EnsureAct.Intralipid = psth_all_prism_3(these,:);
psth_EnsureAct.Sucralose = psth_all_prism_4(these,:);
psth_EnsureAct.Salt = psth_all_prism_5(these,:);
psth_EnsureAct.CitricAcid = psth_all_prism_6(these,:);
psth_EnsureAct.SiliconeOil = psth_all_prism_7(these,:);
psth_EnsureAct.Sucrose = psth_all_prism_8(these,:);
psth_EnsureAct.Quinine = psth_all_prism_9(these,:);
psth_EnsureAct.Water = psth_all_prism_10(these,:);

%Format the data into one matrix for kmeans clustering. Make each row the
%same neuron and append each tastant to make one long trace
for i = 1:size(psth_EnsureAct.Ensure,1)
    psth_EnsureAct.AllTastes(i,:) = [psth_EnsureAct.Ensure(i,:) psth_EnsureAct.Polycal(i,:)...
        psth_EnsureAct.Intralipid(i,:) psth_EnsureAct.Sucralose(i,:) psth_EnsureAct.Salt(i,:) ...
        psth_EnsureAct.CitricAcid(i,:) psth_EnsureAct.SiliconeOil(i,:) psth_EnsureAct.Sucrose(i,:) ...
        psth_EnsureAct.Quinine(i,:) psth_EnsureAct.Water(i,:)];
end

% normalize data from -1 to 1
norm_EnsureAct_all = [];
for i = 1:size(psth_EnsureAct.AllTastes,1)
    range = max(psth_EnsureAct.AllTastes(i,:))-min(psth_EnsureAct.AllTastes(i,:));
    temp = (psth_EnsureAct.AllTastes(i,:)-min(psth_EnsureAct.AllTastes(i,:)))/range;
    norm_EnsureAct_all(i,:) = (temp*2)-1;
end

%Do kmeans, unsupervised
klist=2:10;%the number of clusters you want to try
eva = evalclusters(norm_EnsureAct_all,'kmeans','Silhouette','klist',klist);

IDX=kmeans(norm_EnsureAct_all,eva.OptimalK, 'Replicates', 100);

%Separate Clusters and plot heatmaps
[colormap]=cbrewer('div', 'RdBu', 256);
colormap=flip(colormap);
figure;  
a = [30 91 152 213 274 335 396 457 518 579];
for i = 1:eva.OptimalK
   ToPlot = norm_EnsureAct_all(find(IDX==i),:);
   subplot(eva.OptimalK,1,i);
   h=heatmap_d(ToPlot,[],[],[],'MaxColorValue',1,'MinColorValue',-1,'Colormap',colormap,'Colorbar',1);
   hold on; 
   for j = 1:length(a) %lines to mark where taste was
    plot(zeros(1,50)+a(j), 1:50, 'k', 'LineWidth', 2);
   end
end

%% What are the neurons' responses to all solutions? Find neurons that respond to at least one

%Make a matrix where each row is a neuron and each column is the averaged
%response to a solution
activityAvg_allSolutions = [];
b = [0 1 2 3 4 5 6 7 8 9];
for i = 1:size(psth_allN_allT,1)
    for j = 1:length(b)
        stahrt = 31+b(j)*31; fin = stahrt+30;
        activityAvg_allSolutions(i,j) = mean(psth_allN_allT(i,stahrt:fin));
    end
end

IDX_responders = [];
for i = 1:size(activityAvg_allSolutions,1)
   if find(activityAvg_allSolutions(i,:)>1)
       IDX_responders = [IDX_responders; i];
   elseif find(activityAvg_allSolutions(i,:)<-1)
       IDX_responders = [IDX_responders; i];
   end
end

Responders_allT = psth_allN_allT(IDX_responders, :);
%% K Means for all the neurons

% normalize data from -1 to 1
norm_Resp_all = [];
for i = 1:size(Responders_allT,1)
    range = max(Responders_allT(i,:))-min(Responders_allT(i,:));
    temp = (Responders_allT(i,:)-min(Responders_allT(i,:)))/range;
    norm_Resp_all(i,:) = (temp*2)-1;
end

%Do kmeans, unsupervised
klist=2:10;%the number of clusters you want to try
eva = evalclusters(norm_Resp_all,'kmeans','Silhouette','klist',klist);

IDX=kmeans(norm_Resp_all,eva.OptimalK, 'Replicates', 100);

%Separate Clusters and plot heatmaps
[colormap]=cbrewer('div', 'RdBu', 256);
colormap=flip(colormap);
figure;  
a = [30 91 152 213 274 335 396 457 518 579];
for i = 1:eva.OptimalK
   ToPlot = norm_Resp_all(find(IDX==i),:);
   subplot(eva.OptimalK,1,i);
   h=heatmap_d(ToPlot,[],[],[],'MaxColorValue',1,'MinColorValue',-1,'Colormap',colormap,'Colorbar',1);
   hold on; 
   for j = 1:length(a) %lines to mark where taste was
    plot(zeros(1,50)+a(j), 1:50, 'k', 'LineWidth', 2);
   end
end

%% Take tallies of how many neurons respond to and across solutions
cell_totals = zeros(size(psth_all_prism_1,1), 10); %rows = neurons and columns = solution
cell_means = zeros(size(psth_all_prism_1,1), 10);
for i = 1:size(psth_all_prism_1,1)
    a = mean(psth_all_prism_1(i,31:61));
    if a>1 | a<-1
        %totals(1) = totals(1) + 1;
        cell_totals(i,1) = cell_totals(i)+1;
        cell_means(i,1) = a;
    end
    a = mean(psth_all_prism_2(i,31:61));
    if a >1 | a<-1
        %totals(2) = totals(2) + 1;
        cell_totals(i,2) = cell_totals(i)+1;
        cell_means(i,2) = a;
    end
    a = mean(psth_all_prism_3(i,31:61));
    if a >1 | a<-1
        %totals(3) = totals(3) + 1;
        cell_totals(i,3) = cell_totals(i)+1;
        cell_means(i,3) = a;
    end
    a = mean(psth_all_prism_4(i,31:61));
    if a >1 | a<-1
        %totals(4) = totals(4) + 1;
        cell_totals(i,4) = cell_totals(i)+1;
        cell_means(i,4) = a;
    end
    a = mean(psth_all_prism_5(i,31:61));
    if a >1 | a<-1
        %totals(5) = totals(5) + 1;
        cell_totals(i,5) = cell_totals(i)+1;
        cell_means(i,5) = a;
    end
    a = mean(psth_all_prism_6(i,31:61));
    if a >1 | a<-1
        %totals(6) = totals(6) + 1;
        cell_totals(i,6) = cell_totals(i)+1;
        cell_means(i,6) = a;
    end
    a = mean(psth_all_prism_7(i,31:61));
    if a >1 | a<-1
        %totals(7) = totals(7) + 1;
        cell_totals(i,7) = cell_totals(i)+1;
        cell_means(i,7) = a;
    end
    a = mean(psth_all_prism_8(i,31:61));
    if a >1 | a<-1
        %totals(8) = totals(8) + 1;
        cell_totals(i,8) = cell_totals(i)+1;
        cell_means(i,8) = a;
    end
    a = mean(psth_all_prism_9(i,31:61));
    if a >1 | a<-1
        %totals(9) = totals(9) + 1;
        cell_totals(i,9) = cell_totals(i)+1;
        cell_means(i,9) = a;
    end
    a = mean(psth_all_prism_10(i,31:61));
    if a >1 | a<-1
        %totals(10) = totals(10) + 1;
        cell_totals(i,10) = cell_totals(i)+1;
        cell_means(i,10) = a;
    end
        
end
    
% store the sign/direction of each cell for each tastant
cell_direcxns = zeros(size(psth_all_prism_1,1), 10);
for i = 1:length(cell_means)
    for j = 1:size(cell_means,2)
        if cell_means(i,j) >1
            cell_direcxns(i,j) = 1;
        elseif cell_means(i,j) <-1
            cell_direcxns(i,j) = -1;
        else
            cell_direcxns(i,j) = 0; %make it 2 so there's a marker if one cell wasn't tagged
        end
    end
end

% Now do some comparisons

%% plot heatmaps as one way to visualize all the neurons
% Titles = {'Ensure'; 'Polycal'; 'Intralipid'; 'Sucralose'; 'Salt'; 'Citric Acid'; 'Silicone Oil';...
%     'Sucrose'; 'Quinine'; 'Water'};
% 
% [colormap]=cbrewer('div', 'RdBu', 256);
% colormap=flip(colormap);
% 
% for i = 1:length(Titles)
%     figure;
%     h=heatmap_d(psth_all_means_prism{i},[],[],[],'MaxColorValue',4,'MinColorValue',-4,'Colormap',colormap,'Colorbar',1);
%     hold on;
%     Title = Titles{i};
%     title(Title);
% end
% title(Title);

%% Entropy Calculation

% normalize data from 0 to 1
% norm1_EnsureAct_all = [];
% for i = 1:size(psth_EnsureAct.AllTastes,1)
%     range = max(psth_EnsureAct.AllTastes(i,:))-min(psth_EnsureAct.AllTastes(i,:));
%     norm1_EnsureAct_all(i,:) = (psth_EnsureAct.AllTastes(i,:)-min(psth_EnsureAct.AllTastes(i,:)))/range;
% end

norm1_EnsureAct_tastes = [];
for i = 1:size(psth_EnsureAct_tastes,1)
    range = max(psth_EnsureAct_tastes(i,:))-min(psth_EnsureAct_tastes(i,:));
    norm1_EnsureAct_tastes(i,:) = (psth_EnsureAct_tastes(i,:)-min(psth_EnsureAct_tastes(i,:)))/range;
end

H = []; %vector to store entropy values for all neurons
% K = 1 for 10 stimuli because H = 1 for pi = 1/n (n=10)
% calculate entropy for each neuron using a sum of the activity 30s after
% first lick
a = [0 1 2 3 4 5]; c = a+1;
for i = 1:size(norm1_EnsureAct_tastes,1)
   %first do the summation of activity after first lick
   temp = [];
   for j = 1:length(a)
    stahrt = 31+a(j)*61; fin = stahrt+30;
    temp(j) = sum(norm1_EnsureAct_tastes(i,stahrt:fin));
   end
   b = sum(temp); temp2 = 0;
   for j = 1:length(c)
       p_i = temp(c(j))/b;
       temp2 = temp2+p_i*log10(p_i);
   end
   H(i) = -temp2;
end

%Sort entropy based on preferred taste
H_Pref.Sweet = H(find(TasteNum==1));
H_Pref.Fat = H(find(TasteNum==2));
H_Pref.Sour = H(find(TasteNum==3));
H_Pref.Salt = H(find(TasteNum==4));
H_Pref.Bitter = H(find(TasteNum==5));
H_Pref.Water = H(find(TasteNum==6));
%% Now calculate noise to signal ratio
%Take the highest response and divide by the second highest response.

%Make a matrix where each row is a neuron and each column is the averaged
%response to a solution
activityAvg = [];
a = [0 1 2 3 4 5];
for i = 1:size(psth_tastes,1)
    for j = 1:length(a) 
        stahrt = 31+a(j)*61; fin = stahrt+30;
        activityAvg(i,j) = mean(norm1_AllN_tastes(i,stahrt:fin));
    end
end

%Now go through each neuron and do the N2signal calculation
N2S = [];
for i = 1:size(psth_tastes,1)
    ascension = sort(activityAvg(i,:));
    first = ascension(end); second = ascension(end-1);
    N2S(i)= second/first;
end

% Mark which taste each neuron responded to the most
TasteNum = [];
for i = 1:size(psth_tastes,1)
    TasteNum(i) = find(activityAvg(i,:) == max(activityAvg(i,:)));
end

%Sort the N2S according to which taste each neuron preferred
N2S_Pref.Sweet = N2S(find(TasteNum==1));
N2S_Pref.Fat = N2S(find(TasteNum==2));
N2S_Pref.Sour = N2S(find(TasteNum==3));
N2S_Pref.Salt = N2S(find(TasteNum==4));
N2S_Pref.Bitter = N2S(find(TasteNum==5));
N2S_Pref.Water = N2S(find(TasteNum==6));
%% Entropy test only for tastes
%Use Sucralose (4), Silicone Oil (7), 
%Citric Acid(6), Salt(5), Quinine(9),
%Water(10)
%K = 1.2851 for 6 solutions (my math matches Brooke's paper)

%This one is for all neurons
psth_tastes = [psth_all_prism_4 psth_all_prism_7 psth_all_prism_6...
    psth_all_prism_5 psth_all_prism_9 psth_all_prism_10];

%This one is for neurons activated by Ensure
psth_EnsureAct_tastes = [psth_EnsureAct.Sucralose psth_EnsureAct.SiliconeOil ...
    psth_EnsureAct.CitricAcid psth_EnsureAct.Salt psth_EnsureAct.Quinine...
    psth_EnsureAct.Water];
norm1_AllN_tastes = [];
for i = 1:size(psth_tastes,1)
    range = max(psth_tastes(i,:))-min(psth_tastes(i,:));
    norm1_AllN_tastes(i,:) = (psth_tastes(i,:)-min(psth_tastes(i,:)))/range;
end

H = []; %vector to store entropy values for all neurons
a = [0 1 2 3 4 5]; c = a+1;
for i = 1:size(norm1_AllN_tastes,1)
   %first do the summation of activity after first lick
   temp = [];
   for j = 1:length(a)
    stahrt = 31+a(j)*61; fin = stahrt+30;
    temp(j) = sum(norm1_AllN_tastes(i,stahrt:fin));
   end
   b = sum(temp); temp2 = 0;
   for j = 1:length(c)
       p_i = temp(c(j))/b; %do the sum of the logs
       temp2 = temp2+p_i*log10(p_i);
   end
   H(i) = -1.2851*temp2;
end
%% Extra/old code
% auc_act_prism = {}; auc_inhib_prism = {}; auc_none_prism = {};
% for j = 1: size(auc_all_a,2)
%     tempa = []; tempi = []; tempn = [];
%     for i = 1:length(infile)
%      tempa = [tempa; cell2mat(auc_all_a{i}{j})];    
%      %tempi = [tempi; psth_all_inhib{i}(:,j)];
%      %tempn = [tempn; psth_all_none{i}(:,j)];
%     end
%     auc_act_prism{j} = tempa;
%     %auc_inhib_prism{j} = tempi;
%     %auc_none_prism{j} = tempn;
% end   

% temp_a = []; auc_act_prism = {};
% for i = 1:length(infile)
%     for j = 1:36
%         temp_a = [temp_a; cell2mat(mouse{j})];
%         auc_act_prism{j}=temp_a;
%     end
%     temp_a = [];
% end

% Suc32 = []; Suc16 = []; Suc8 = []; Suc4 = []; Suc2=[]; Suc1=[]; Suc0_5=[];
% Suc0_25=[]; Suc0_125=[]; Water=[];
% Suc32_perc = []; Suc16_perc = []; Suc8_perc = []; Suc4_perc = []; Suc2_perc=[]; Suc1_perc=[]; Suc0_5_perc=[];
% Suc0_25_perc=[]; Suc0_125_perc=[]; Water_perc=[];
% for j=1:36
%     for i = 1:length(infile)
%         mouse = auc_all_a{i};
%         now = perc_act{i};
%         %if size(mouse,2)>=j %stop the code in case there's nothing at the end
%             if j==9 || j==21 || j==33
%                 %Suc32 = [Suc32; mouse{j}];
%                 Suc32_perc = [Suc32_perc; now(j)/now(37)];
%             elseif j==1 || j==13 || j==25
%                 %Suc16 = [Suc16; mouse{j}];
%                 Suc16_perc = [Suc16_perc; now(j)/now(37)];
%             elseif j==2 || j==14 || j==26
%                 %Suc8 = [Suc8; mouse{j}];
%                 Suc8_perc = [Suc8_perc; now(j)/now(37)];
%             elseif j==4 || j==16 ||j==28
%                 %Suc4 = [Suc4; mouse{j}];
%                 Suc4_perc = [Suc4_perc; now(j)/now(37)];
%             elseif j==11 || j==23 || j==35
%                 %Suc2 = [Suc2; mouse{j}];
%                 Suc2_perc = [Suc2_perc; now(j)/now(37)];
%             elseif j==5 || j==17 || j==29
%                 %Suc1 = [Suc1; mouse{j}];
%                 Suc1_perc = [Suc1_perc; now(j)/now(37)];
%             elseif j==12 ||j==24 ||j==36
%                 %Suc0_5 = [Suc0_5; mouse{j}];
%                 Suc0_5_perc = [Suc0_5_perc; now(j)/now(37)];
%             elseif j==8 || j==20 || j==32
%                 %Suc0_25 = [Suc0_25; mouse{j}];
%                 Suc0_25_perc = [Suc0_25_perc; now(j)/now(37)];
%             elseif j==6||j==18||j==30
%                 %Suc0_125 = [Suc0_125; mouse{j}];
%                 Suc0_125_perc = [Suc0_125_perc; now(j)/now(37)];
%             else
%                 %Water = [Water; mouse{j}];
%                 Water_perc = [Water_perc; now(j)/now(37)];
%             end
%         end
%     end
% %end
% prism_perc = [Suc32_perc Suc16_perc Suc8_perc Suc4_perc Suc2_perc Suc1_perc Suc0_5_perc Suc0_25_perc Suc0_125_perc];
% auc.Suc32 = Suc32; auc.Suc16 = Suc16; auc.Suc8 = Suc8; auc.Suc4 = Suc4;
% auc.Suc2 = Suc2; auc.Suc1 = Suc1; auc.Suc0_5=Suc0_5; auc.Suc0_25 =Suc0_25; 
% auc.Suc0_125 =Suc0_125; auc.Water = Water; 
% perc_data = [Suc32_perc Suc16_perc Suc8_perc Suc4_perc Suc2_perc Suc1_perc Suc0_5_perc Suc0_25_perc Suc0_125_perc];
% perc_data = [perc_data [mean(Water_perc(1:3)); mean(Water_perc(4:6)); mean(Water_perc(7:9))]];
% 
% clear Suc0_125 Suc0_5 Suc0_25 Suc1 Suc2 Suc4 Suc8 Suc16 Suc32 Water;
