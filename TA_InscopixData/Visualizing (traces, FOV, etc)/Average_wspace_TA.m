%% Setting your parameters
%This code takes the file generated from making a heatmap of inscopix
%neurons, and plots a single line as an average for all of the neurons in
%the heatmap (plotting mean with SEM error bars).

%close all; 
clear all; %closes and clears everything to give a clean slate in your workspace

%the wspace file from the heatmap you generated. 
%Usually written as the name of the .txt file with _wspace tacked on at the end
infile = 'TL137_DavisTest_1027_dff_wspace'; 

max = 60; %how long was the session in minutes?
stim = 10; %optional, what's your stim time in minutes? Can't be equal to 0.
name = 'TL137 Davis Test'; %What do you want the graph to be titled?
load(infile); %this loads in the file

%% Generate the data to plot
%all_mrerank is the variable in the file with all of the neurons, ranked by
%response relative to baseline. Rows are different cells and columns are time

%the mean function takes the mean of all_mrerank across this first dimension (so collapsing the rows) 
neuron_means = mean(all_mrerank,1); 

%taking the standard error for our error bars at each time point
neurons_sem = std(all_mrerank,1)/sqrt(size(all_mrerank,1)); 

%generate the x data points or time.
x = 1:size(all_mrerank,2); %take the length of the trace
x = x/60; %converting seconds to minutes

%% Plot it

figure; %use this to generate a new window for the plot

%plot the error bars, I found this method from google but it basically
%makes a box out of the error bars then fills it in. Plot this first
%otherwise it obscures the mean line
fill([x';flipud(x')],[neuron_means'-neurons_sem';flipud(neuron_means'+neurons_sem')],[.9 .9 .9],'linestyle','none');
hold on;

%Plot the line
plot(x,neuron_means, 'k'); %plots a black ('k') line, you can look up this function for more options
hold on; %have to add this otherwise when you plot something else, it will write over it and erase it

%let's also plot a line at zero, so that it's easier to see changes.
%Comment this out if you don't want it.
a = 1:max;
b = zeros(1,length(a));
plot(a,b,'k--'); %'--' is the code for a dashed line
hold on;

%Plot a line when the stim happens, it's optional so use an if statement
YL = ylim; %get the y limits of the graph
if exist('stim') && stim ~=0 %check that stim exists and is not equal to 0
    a = YL(1):YL(2); %plot from top to bottom of graph
    c = zeros(1,length(a))+stim; 
    plot(c,a,'b--'); %plot a blue dashed line
end

%Adding labels to the figure
title(name);
xlabel('Time (minutes)');
ylabel('Z-Scored Change in Calcium Activity');
xlim([0 max]); %Setting the x axis limit