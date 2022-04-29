%% Input stuffies
% This script takes the average activity after two different stimuli and
% outputs a heatmap showing the aligned traces, and the averages
% (mean_output, where each row is a neuron and column 1 is stim1, column2
% is for stim2). This is for sessions where you do both stimuli in the same
% session. This also re-baselines the traces before the second stimuli. You
% do this before hand by running "mdata_readmlist" but with the stim time
% at the second stimulus; This is what the "infile_suc" is. This script was
% originally written for a trial where a hormone was given before sucrose
% access, then modified for a trial where nVoke stimulation was given
% before Ensure access (I just want to explain the naming here). 
%
% At the end, the script quantifies how many were act/none/inhib by
% both/either/none of the stimuli
close all; clear all; 
infile_hormone = 'BNST_FastedSucrose_dff_wspace_noZ'; %output from readmlist centered around hormone
infile_suc = 'BNSTtoDMH Fasted Sucrose (ensure)_dff_wspace'; %output from readmlist centered around Ensure
Title = 'BNST2DMH_Sucralose';
Hormone = 'Sucrose'; %label for the plot on the left. I would put stim 
prestim = 600; %baseline time in seconds
Int = 300; %how many seconds after baseline do you want to average?
timeSuc = 1800; %what time was Ensure (in seconds)?

%% Setup stuff
load(infile_hormone);
Title = {Title}; Hormone = {Hormone};
ind_hormone = ind_mrerank; %store the neural number for each row
hormone = 1;
%% If doing hormone analysis, also quantify which neurons are modulated by licking

if hormone == 1
    mean_lick = []; mean_hormone = [];
    Int_end = timeSuc+600;
    [norm_rerank] = normab(cC, prestim, 1800, timestamp);
    for i = 1:size(all_mrerank,1) %step through each neuron
        %a = mean(all_mrerank(i, prestim:prestim+600));
        a = mean(all_mrerank(i, 600:900)); %was 1200
        mean_hormone = [mean_hormone; a];
    end
    
    all_mrerank_hormone = all_mrerank; 
    clear all_mrerank; clear ind_mrerank; 
    load(infile_suc);
    ind_suc = ind_mrerank;
    
    for i = 1:size(all_mrerank,1) %step through each neuron
        %a = mean(all_mrerank(i, prestim:prestim+600));
        a = mean(all_mrerank(i, 600:900));
        mean_lick = [mean_lick; a];
    end
    
else
    mean_hormone = 0;
    mean_lick = 0;
end

%% Match up the responses from 2 separate matrices into one

%lick first then hormone
mean_output = []; %the variable to put into prism, each row is a neuron
for i = 1:size(mean_lick,1)
    mean_output(i,1) = mean_lick(find(ind_suc == i));
    mean_output(i,2) = mean_hormone(find(ind_hormone == i));

end

%% Store the neuron traces based on reaction to stim

%Neurons_act = all_mrerank(act, :);
%Neurons_none = all_mrerank(none, :);
%Neurons_inhib = all_mrerank(inhib,:); 

%save(strcat(Title,'_quant'))
%save(Title);
%end

%% Plotting the two heatmaps next to eachother

%Ranking the sucrose neurons according to their order after hormone
all_mrerank_sucrose = [];
for i = 1:size(all_mrerank, 1)
    all_mrerank_sucrose(i,:) = all_mrerank(find(ind_mrerank == ind_hormone(i)), :);
end

figure; 
subplot(1,2,1)
h=heatmap_d(all_mrerank_hormone,[],[],[],'MaxColorValue',4,'MinColorValue',-4,'Colormap',colormap,'Colorbar',1);
hold on; title(Hormone);
y = 1:400; z = zeros(1,400)+600; x = zeros(1,400)+1800;
plot(z,y,'k', 'LineWidth', 2); hold on; plot(x,y, 'k', 'LineWidth', 2);
hold on;
yticks([1]); yticklabels([size(all_mrerank_hormone,1)]); hold on; 

subplot(1,2,2)
h=heatmap_d(all_mrerank_sucrose,[],[],[],'MaxColorValue',4,'MinColorValue',-4,'Colormap',colormap,'Colorbar',1);
hold on; title('IG Sucrose');
y = 1:400; z = zeros(1,400)+300;
plot(z,y,'k', 'LineWidth', 2); hold on;
yticks([1]); yticklabels([size(all_mrerank_hormone,1)]); 
%% Sorting categories
%e.g. AA = activated by stimulus 1 and 2; N = no change; I = inhibited
AA = 0; AN=0;AI=0;NA=0;NN=0;NI=0;IA=0;IN=0;II=0;
for i = 1:size(mean_output,1)
    if mean_output(i,1) >=1 
        if mean_output(i,2) >=1
            AA = AA+1;
        elseif mean_output(i,2) <=-1
            AI = AI+1;
        else
            AN = AN+1;
        end
    elseif mean_output(i,1) <1 && mean_output(i,1) >-1
        if mean_output(i,2) >=1
            NA = NA+1;
        elseif mean_output(i,2) <=-1
            NI = NI+1;
        else
            NN = NN+1;
        end
    elseif mean_output(i,1) <=-1
        if mean_output(i,2) >=1
            IA = IA+1;
        elseif mean_output(i,2) <=-1
            II = II+1;
        else
            IN = IN+1;
        end
    end
end
sort.AA = AA; sort.AN = AN; sort.AI = AI; sort.NA = NA; sort.NN = NN; ...
    sort.NI = NI; sort.IA = IA; sort.IN = IN; sort.II = II;
clear AA AN AI NA NN NI IA IN II; 