%adapted from Tara's LickAnalysis script
%filename = the data file (.csv) to be read. 
%% Importing and formatting the data
filename = '2022_02_02__09_57_02_TA137_220202.csv'; 
save_dir = 'TA133_FedEnsure_Inhibition2_112221';
bout_seperation = 1000; %in milliseconds 
%if two licks are more than x seconds apart they are considered as seperate bouts
Data = csvread(filename, 1, 7); 
%isolating timestamp (column 1) and licks (Column 7)
Data = [Data(:, 1), Data(:, 5)];
%turn all the seconds into milliseconds
Data(:,1) = Data(:,1)*1000;
%% Adding a lick counter
Counter = [0];
for i = 2:size(Data,1)
    if Data(i,2) == 0
        Counter = [Counter Counter(end)];
    else
        Counter = [Counter Counter(end)+1];
    end
end

CumLicks = [Data(:,1), Counter(:)];
%% Sort out the times when the mouse licks
Licks = []; % prep an empty vector to store the lick times
for i = 1:size(Data,1)
    if Data(i, 2) == 1
        Licks = [Licks Data(i,1)];
    end
end
%% Adding timestamps for each millisecond
tvec = 0:1:1800000;
TimeWhenLick = ismember(tvec, Licks);
LickData = [tvec(:), TimeWhenLick(:)];
%% Find the bouts
Bout_start = [Licks(1)]; index_bouts = [1];
for i = 2:length(Licks)
    if Licks(i)>=Licks(i-1)+ bout_seperation 
        Bout_start = [Bout_start Licks(i)];
        index_bouts = [index_bouts i];
    end
end
%% Calculate the length of the bouts (in time and number of licks)
times = []; number = []; 
for i = 2:length(index_bouts)
    times_temp = Licks(index_bouts(i)-1)-Licks(index_bouts(i-1));
    times = [times times_temp]; 
    
    number_temp = index_bouts(i)-index_bouts(i-1);
    number = [number number_temp];
end
%add statement to erase bouts with only one lick
bout.boutlength_time = times; bout.numberLicks = number; bout.licks = Licks; 
bout.boutnumber = length(index_bouts);
%% Lick rate during each bout
rate = []; 
for i = 1:length(times)
    if times(i) ~= 0
        rate = [rate number(i)/times(i)];
    end
end
bout.rate = rate;
%% Save in a file
%save(save_dir);



