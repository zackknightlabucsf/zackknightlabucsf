function [Number_pel_eaten, Time_pel_taken, y] = FeedingAnalysis(filename, number, expt)


%NOTE that this is old code from 2017, but might still be useful


%Function accepts csv files from freefeeding protocol (thus far only used
%for the 30prestim 60 free feeding with and without the laser on). 
%
%filename = the data file (.csv) to be read. 
%number = the number of data sets you'll analyze in sequence. So if you
%are planning to analyze 4 sets, and you're on the first one, put 1; on the
%second one, put 2, etc. 
%expt = description of experiment to be included in plot
%
%NOTE: file must already be imported and designated as "filename" before
%using the function
%Data output is peri stim, so all the feeding events are time stamped
%relative to when the stim turned off. 
%



%% Importing and formatting the data
Data = csvread(filename, 1, 7); 
%isolating when event occurs (column 1), current state (column 2), and
%when pellet beam is broken (Off1A3) (column 3)
Data = [Data(:, 1) Data(:, 2), Data(:, 11)]; 

%% Analysis
%Create a matrix that stores the time stamps when the beam breaks (pellet
%taken), and also adds to a counter
%n = size(Data_peri_stim,1); 
n = size(Data,1); 
Time_pel_taken = [];
Number_pel_eaten = 0;

% Add second conditional that the current state of the feeder (8) is also
% present in order to count it. 
%NOTE can't use this number notation in all situations, including the lickometer
for i = 1 : n
    if Data(i, 3) == 1 && Data(i-1, 2) == 3 %%situational
        Time_pel_taken = [Time_pel_taken; Data(i,1)];
        Number_pel_eaten = Number_pel_eaten + 1;
    else 
        Time_pel_taken = Time_pel_taken;
        Number_pel_eaten = Number_pel_eaten;
    end
end

%% Prep for plot, if desired
% Have y just add one at each time point in Time_pel_taken
y = 1:1:size(Time_pel_taken,1);

%Plot it!
figure(1);
if number == 1
    plot(Time_pel_taken, y, 'b.-');
    title('Pellets eaten over time'); 
    xlabel ('time peri stim (s)');
    %xlim([0,3600]);
    ylabel ('number pellets taken');
    text1 = ([num2str(Number_pel_eaten) ' total pellets eaten in ' num2str(expt)]);
    text(2000, 1,text1, 'Color', [0 0 1]); hold on; 
    x = zeros(1,70)+1800; y = 1:70; plot(x,y,'k'); 
    hold off
elseif number == 2
    plot(Time_pel_taken, y, 'g.-');
    title('Pellets eaten over time'); 
    xlabel ('time peri stim (s)');
    xlim([0,3600]);
    ylabel ('number pellets taken');
    text1 = ([num2str(Number_pel_eaten) ' total pellets eaten in ' num2str(expt)]);
    text(2000, 5,text1, 'Color', [0 1 0]);
    hold off
elseif number == 3
    plot(Time_pel_taken, y, 'm.-');
    title('Pellets eaten over time'); 
    xlabel ('time peri stim (s)');
    xlim([0,3600]);
    ylabel ('number pellets taken');
    text1 = ([num2str(Number_pel_eaten) ' total pellets eaten in ' num2str(expt)]);
    text(2000, 10,text1, 'Color', [1 0 1]);
    hold off
elseif number == 4
    plot(Time_pel_taken, y, 'k.-');
    title('Pellets eaten over time'); 
    xlabel ('time peri stim (s)');
    xlim([0,3600]);
    ylabel ('number pellets taken');
    text1 = ([num2str(Number_pel_eaten) ' total pellets eaten in ' num2str(expt)]);
    text(2000, 15,text1, 'Color', [0 0 0]);
    hold off
else 
    plot(Time_pel_taken, y, 'y.-');
    title('Pellets eaten over time'); 
    xlabel ('time peri stim (s)');
    xlim([0,3600]);
    ylabel ('number pellets taken');
    text1 = ([num2str(Number_pel_eaten) ' total pellets eaten in ' num2str(expt)]);
    text(2000, 20,text1, 'Color', [1 1 0]);
    hold off
end
