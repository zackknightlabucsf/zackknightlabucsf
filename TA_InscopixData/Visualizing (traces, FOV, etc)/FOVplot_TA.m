%%FOV plot: plot the neuron ROIs and color them according to their average
%%response after a stimulus

%% input stuffies
close all; clear all; 
infile2 = 'TA68_fastedEnsureWater2_wspace'; %small file (need neuron.Coor from here)
infile1 = '2020-06-26-14-22-06_video-pp-BP-MC_wspace_small_post'; %quant file
prestim = 600; %what time is stim given?
Int = 600; %what length of time interval post stim do you want to average activity?

%% load in code
load(infile1);
load(infile2);

poststim = prestim+Int;
mean_trace = []; %set up matrix to hold average values to color ROIs with
for i = 1:size(all_mrerank,1)
    mean_trace = [mean_trace; mean(all_mrerank(i,prestim:poststim))];
end

%find range of these mean values to setup color range
ran = round(range(mean_trace))+1;
min_value = min(mean_trace);
rgb = vals2colormap(mean_trace, [], [-4 4]);

figure; 
for i = 1:size(neuron.Coor,1)
    %ind_mrerank = indices in mrerank to match up to neuron.Coor
    %color_value = mean_trace(find(ind_mrerank==i)); %get the calcium value for this neuron
    %color_value = mean_trace(i);
    %can't do this because the center of the calcium data isn't the center
    %of the range.....
    %color_value = round(color_prevalue+abs(min_value))+1;%scale the calcium value to the color bar
    cont = medfilt1(neuron.Coor{i}')';
    %fill(cont(1,:), cont(2,:), rgb(i,:));
    plot(cont(1,:), cont(2,:), 'k');
    hold on;
end
colormap RdBu
colorbar 'RdBu';