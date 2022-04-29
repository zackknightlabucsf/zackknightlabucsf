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

%% Notes to Self
%DONE:TA note 12/3/21: add in something to consider retries (an if statement)
%% Inputs
IPI = 60000; %inter-presentation interval, in ms
infile = 'TA121_fastedSucrose1115.csv'; %csv file converted from .txt. to .csv
%wspace_file = ''; % neural data output from CNMFE_readmlist
start_time = 799; %time when door opens (in s)
x = -30:30;

%% Read in the data
latency_1 = csvread(infile, 11, 9, [11 9 46 9]); %latency to first lick
ILIs = csvread(infile, 48, 0); %inter-lick times, column 1 is presentation #
Retries = csvread(infile, 11, 10, [11 10 46 10]);

%Which tastants were which presentations
%1-9 are high-low, and 10-12 are water
Order = [2 3 11 4 6 9 10 8 1 12 5 7 2 3 11 4 6 9 10 8 1 12 5 7 2 3 11 4 6 9 10 8 1 12 5 7];
Water = find(Order>9); 
Suc = [9 1 2 4 11 5 12 8 6];
%% Find the first lick in each bout
if Retries(1) == 0
    bout_start_ms = [latency_1(1)]; %storage vector. First bout starts at first latency
else 
    bout_start_ms = [latency_1(1)+8000*Retries(1)+IPI*Retries(1)];
end

for i = 2:size(latency_1, 1)
    %time start of next bout is time of first lick in last bout + 5s +3s door close+ IPI + latency of first lick in second bout
   if  Retries(i) ~= 0
       next_lick = bout_start_ms(i-1)+8000*Retries(i)+IPI*Retries(i)+latency_1(i);
    else
        next_lick = bout_start_ms(i-1)+ 5000 +IPI+latency_1(i)+3000; 
    end
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
    temp = Neurons_act(:, bout_start(i)-30:bout_start(i)+30);
    psth{end+1} = temp;
    for j = 1:size(Neurons_act,1)
        %now re-baseline to pre-bout
        temp2 = [temp2; (temp(j, :)/mean(temp(j,1:30)))];
    end
    psth_norm{end+1} = temp2;
    temp = []; temp2 = [];
end

%have a special one for the last bout in case not enough end time
if bout_start(end)+30> size(Neurons_act, 2)
    temp = Neurons_act(:, bout_start(36)-30:end);
    psth{end+1} = temp;
    for j = 1:size(Neurons_act,1)
        temp2 = [temp2; (temp(j, :)/mean(temp(j,1:30)))];
    end
    psth_norm{end+1} = temp2;
else
    temp = Neurons_act(:, bout_start(36)-30:bout_start(36)+30);
    psth{end+1} = temp;
    for j = 1:size(Neurons_act,1)
        %now re-baseline to pre-bout
        temp2 = [temp2; (temp(j, :)/mean(temp(j,1:30)))];
    end
    psth_norm{end+1} = temp2;
end

%% Average the water traces into 1
water_psth_b1 = [psth{Water(1)}; psth{Water(2)}; psth{Water(3)}]; 
water_psth_b2 = [psth{Water(4)}; psth{Water(5)}; psth{Water(6)}]; 
water_psth_b3 = [psth{Water(7)}; psth{Water(8)}; psth{Water(9)}];
int = size(Neurons_act,1);

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
%lightRed = [255 153 153];
%darkBlue = [0 0 139];
%darkRed = [83 0 0];
Blue = [0 0 1];
Red = [1 0 0];
redGRADIENT = @(n,nn) interp1([1/nn 1], [Red; Blue], n/nn);
%test = redGRADIENT((1:length(Suc)), length(Suc));

%First trial
figure; 
for i = 1:length(Suc)
    plot(x, mean(psth{Suc(i)}, 1), 'Color', [redGRADIENT(i,length(Suc))]); hold on;
end
plot(x, mean(water_psth_1,1), 'k'); hold on; 
title('First trial');  

%Second trial
figure; 
for i = 1:length(Suc)
    plot(x, mean(psth{Suc(i)+12}, 1), 'Color', [redGRADIENT(i,length(Suc))]); hold on;
end
plot(x, mean(water_psth_2,1), 'k'); hold on; 
title('Second trial'); 

%Third trial
figure; 
for i = 1:length(Suc)
    plot(x, mean(psth{Suc(i)+24},1), 'Color', [redGRADIENT(i,length(Suc))]); hold on;
end
plot(x, mean(water_psth_3,1), 'k'); hold on; 
title('Third trial'); 

