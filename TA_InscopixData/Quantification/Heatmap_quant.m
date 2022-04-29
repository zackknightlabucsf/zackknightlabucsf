%% Heatmap_Quant
% Use to sort neurons based on activity for some time after a stimulus
%% Input stuff
infile = 'BNSTtoDMH_FastedSucralose_dff_wspace'; %file name? This is the wspace file outputted from the heatmap function 
Title = 'Filler'; %title for the save file?
prestim = 600; %when was the stimulus given in seconds? 
Int = 600; %how many seconds post stimulus do you want to average over (e.g. 600=10min post lick access)
piechart = 0; %Do you want a pie chart plotted? 0 = no, 1 = yes


%% Load in stuff

load(infile, '-mat', 'all_mrerank'); %let's just load in all_mrerank to make the workspace cleaner. 

%% Quantify percentage of neurons activated and inhibited by stim

%setup storage variables
act = [];none=[];inhib = [];


for i = 1:size(all_mrerank,1) %step through each neuron
    a = mean(all_mrerank(i, prestim:prestim+Int));
    if a<-1
        inhib = [inhib i]; %store the neuron number; these are the indexes that we'll pull later
    elseif a>1
        act = [act i]; %store the neuron number
    else
        none = [none i];
    end
end

%% make a pie chart for how many neurons are act, inhib, or none
if piechart == 1
    X = [length(act), length(none), length(inhib)];
    
    figure;
    labels = {'Activated', 'No Change', 'Inhibited'};
    clear colormap
    pie(X);
    colormap([1 0 0; 1 1 1; 0 0 1]);
    legend(labels, 'Location', 'southoutside', 'Orientation', 'horizontal');
    title([Title]);
end

%% Store the neuron traces based on reaction to stim

%Each row is a neuron, and each column is a timepoint. Number of
%rows=number of neurons per category
Neurons_act = all_mrerank(act, :); %use the indexes we found above and pull out only those neurons
Neurons_none = all_mrerank(none, :);
Neurons_inhib = all_mrerank(inhib,:); 


save(strcat(Title,'_quant')) %save the file

