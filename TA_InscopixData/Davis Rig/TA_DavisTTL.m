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
    %'TA120_fastedSucroseDavis2_1213.txt'
    %'TA121_fastedSucroseDavis2_1213.txt'
    %'TA124_fastedSucroseDavis2_1213.txt'
    %'TA130_fastedSucroseDavis2_1213.txt'
    
    %'TA120_fastedSucraloseDavis2_1215.txt'
    %'TA121_fastedSucraloseDavis2_1215.txt'
    %'TA124_fastedSucraloseDavis2_1215.txt'
    %'TA130_fastedSucraloseDavis2_1215.txt'
    
    'TA121_fastedLipidCurve_0331.txt'
    'TA124_fastedLipidCurve_0331.txt'
    'TA175_fastedLipidCurve_0401.txt'
    'TA176_fastedLipidCurve_0401.txt'
    'TA186_fastedLipidCurve_0401.txt'
    
    %'TA121_fastedTastePanel_0405.txt'
    %'TA124_fastedTastePanel_0405.txt'
    %'TA175_fastedTastePanel_0406.txt'
    %'TA176_fastedTastePanel_0406.txt'
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
    
    '2022-03-31-14-05-48_gpio.csv'
    '2022-03-31-15-25-49_gpio.csv'
    '2022-04-01-12-00-11_gpio.csv'
    '2022-04-01-13-17-14_gpio.csv'
    '2022-04-01-14-30-05_gpio.csv'
    
    %'2022-04-05-12-02-07_gpio.csv'
    %'2022-04-05-13-20-18_video_green_gpio.csv'
    %'2022-04-06-11-34-07_gpio.csv'
    %'2022-04-06-12-52-27_gpio.csv'
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
     tempa = [tempa; cell2mat(psth_all_act{i}(:,j))];    
     %tempi = [tempi; psth_all_inhib{i}(:,j)];
     %tempn = [tempn; psth_all_none{i}(:,j)];
     tempall = [tempall; cell2mat(psth_all{i}(:,j))];
    end
    psth_act_prism{j} = tempa;
    %psth_inhib_prism{j} = tempi;
    %psth_none_prism{j} = tempn;
    psth_all_prism{j} = tempall; 
end   

% make a variable that has the means across different tastants


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
