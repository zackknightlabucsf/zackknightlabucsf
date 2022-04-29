%%DavisRigAnalysis_TA
%
%For reading the data from the davis rig file and syncing with inscopix
%data (for when the TTL output adapter isn't working).
%Since there's a full minute between presentations and mice are only
%allowed to lick for 5s, treat each presentation as it's own bout. 
%
%Find the first lick in each bout and cut neural activity 30s peri this lick to
%create a PSTH for each presentation
%
%Note that this script is specific to one protocol!!

%% Inputs
IPI = 60000; %inter-presentation interval, in ms
infile = '1027TL137_DavisRig_nVista.ms8.csv'; %csv file converted from .txt. to .csv
%wspace_file = ''; % neural data output from CNMFE_readmlist
start_time = 661; %time when door opens (in s)
x = -30:30;

%% Read in the data
latency_1 = csvread(infile, 11, 9, [11 9 46 9]); %latency to first lick
ILIs = csvread(infile, 48, 0); %inter-lick times, column 1 is presentation #

%Which tastants were which presentations
Water = [1 9 12 13 21 24 25 33 36]; 
Ensure = [2 14 26];
Sucrose_low = [3 15 27];
Intralipid = [4 16 28];
Sucralose_high = [5 17 19];
Salt = [6 18 30];
Sucralose_low = [7 19 31];
SiliconeOil = [8 20 32];
Sucrose_high = [10 22 34];
quinine = [11 23 35];
%% Find the first lick in each bout
bout_start_ms = [latency_1(1)]; %storage vector. First bout starts at first latency

for i = 2:size(latency_1, 1)
    %time start of next bout is time of first lick in last bout + 5s +3s door close+ IPI + latency of first lick in second bout
    next_lick = bout_start_ms(i-1)+ 5000 +IPI+latency_1(i)+3000; 
    bout_start_ms = [bout_start_ms; next_lick];
    next_lick = 0; %reset for next bout
end

bout_start = round(bout_start_ms/1000); %convert from ms to s for neural recording
bout_start = bout_start +start_time;
%% Cut neural activity around this lick
psth = {}; %each cell is the activity for that presentation (1-36)
psth_norm = {}; %rebaselined
temp = []; temp2 = []; %temporary place holder for each neurons' activity

for i = 1:length(bout_start)-1
    temp = all_mrerank(:, bout_start(i)-30:bout_start(i)+30);
    psth{end+1} = temp;
    for j = 1:size(all_mrerank,1)
        %now re-baseline to pre-bout
        temp2 = [temp2; (temp(j, :)/mean(temp(j,1:30)))];
    end
    psth_norm{end+1} = temp2;
    temp = []; temp2 = [];
end

%have a special one for the last bout in case not enough end time
if bout_start(end)+30> size(all_mrerank, 2)
    temp = all_mrerank(:, bout_start(36)-30:end);
    psth{end+1} = temp;
    for j = 1:size(all_mrerank,1)
        temp2 = [temp2; (temp(j, :)/mean(temp(j,1:30)))];
    end
    psth_norm{end+1} = temp2;
else
    temp = all_mrerank(:, bout_start(36)-30:bout_start(36)+30);
    psth{end+1} = temp;
    for j = 1:size(all_mrerank,1)
        %now re-baseline to pre-bout
        temp2 = [temp2; (temp(j, :)/mean(temp(j,1:30)))];
    end
    psth_norm{end+1} = temp2;
end

%% Average the water traces into 1
water_psth_b1 = [psth{1}; psth{9}; psth{12}]; 
water_psth_b2 = [psth{13}; psth{21}; psth{24}]; 
water_psth_b3 = [psth{25}; psth{33}; psth{36}];
int = size(all_mrerank,1);

for i = 1:int
    a = [water_psth_b1(i,:); water_psth_b1(i+int,:); water_psth_b1(i+2*int, :)];
    water_psth_1(i,:) = mean(a, 1); 

    b = [water_psth_b2(i,:); water_psth_b2(i+int,:); water_psth_b2(i+2*int, :)];
    water_psth_2(i,:) = mean(b, 1);
    
    c = [water_psth_b3(i,:); water_psth_b3(i+int,:); water_psth_b3(i+2*int, :)];
    water_psth_3(i,:) = mean(c, 1);
end

%normalize the traces to baseline

water_psth_1n = []; water_psth_2n = []; water_psth_3n = [];
for j = 1:size(water_psth_1,1)
    water_psth_1n = [water_psth_1n; (water_psth_1(j, :)/mean(water_psth_1(j,1:30)))];
    water_psth_2n = [water_psth_2n; (water_psth_2(j, :)/mean(water_psth_2(j,1:30)))];
    water_psth_3n = [water_psth_3n; (water_psth_3(j, :)/mean(water_psth_3(j,1:30)))];
end

%% Plot
figure; 
%Water
plot(x, mean(water_psth_1,1)); hold on; 
plot(x, mean(water_psth_2, 1)); hold on; 
plot(x, mean(water_psth_3, 1)); hold on; 
title('Water'); legend ('First', 'Second', 'Third');

%Ensure
figure; 
plot(x, mean(psth{2}, 1)); hold on;
plot(x, mean(psth{14}, 1)); hold on;
plot(x, mean(psth{26}, 1)); hold on;
title('Ensure'); legend ('First', 'Second', 'Third');

%Sucrose High
figure; 
plot(x, mean(psth{10}, 1)); hold on;
plot(x, mean(psth{22}, 1)); hold on;
plot(x, mean(psth{34}, 1)); hold on;
title('Sucrose High'); legend ('First', 'Second', 'Third');

%Sucrose Low
figure; 
plot(x, mean(psth{3}, 1)); hold on;
plot(x, mean(psth{15}, 1)); hold on;
plot(x, mean(psth{27}, 1)); hold on;
title('Sucrose Low'); legend ('First', 'Second', 'Third');

%Sucralose High
figure; 
plot(x, mean(psth{5}, 1)); hold on;
plot(x, mean(psth{17}, 1)); hold on;
plot(x, mean(psth{29}, 1)); hold on;
title('Sucralose High'); legend ('First', 'Second', 'Third');

%Sucralose Low
figure; 
plot(x, mean(psth{7}, 1)); hold on;
plot(x, mean(psth{19}, 1)); hold on;
plot(x, mean(psth{31}, 1)); hold on;
title('Sucralose Low'); legend ('First', 'Second', 'Third');

%Intralipid
figure; 
plot(x, mean(psth{4}, 1)); hold on;
plot(x, mean(psth{16}, 1)); hold on;
plot(x, mean(psth{28}, 1)); hold on;
title('Intralipid'); legend ('First', 'Second', 'Third');

%Silicone Oil
figure; 
plot(x, mean(psth{8}, 1)); hold on;
plot(x, mean(psth{20}, 1)); hold on;
plot(x, mean(psth{32}, 1)); hold on;
title('Silicone Oil'); legend ('First', 'Second', 'Third');

%Salt
figure; 
plot(x, mean(psth{6}, 1)); hold on;
plot(x, mean(psth{18}, 1)); hold on;
plot(x, mean(psth{30}, 1)); hold on;
title('Salt'); legend ('First', 'Second', 'Third');

%Quinine
figure; 
plot(x, mean(psth{11}, 1)); hold on;
plot(x, mean(psth{23}, 1)); hold on;
plot(x, mean(psth{35}, 1)); hold on;
title('Quinine'); legend ('First', 'Second', 'Third');
