%Count how long a mouse got opto stim
Data = csvread('2022_02_02__10_52_02_TA143_220202.csv', 1,7);
Times = Data(:,1); 
States = [];
for i = 1:size(Data,1)
    if Data(i,2) == 1 && Data(i,3) == 2 %if laser turns on
        States = [States 1]; %1=on
    elseif Data(i,2) == 2 && Data(i,3) == 1 %if laser turns off
        States = [States 0]; %0 is off
    else
        States = [States 2]; %2 is nothing
    end
end
On_Times = find(States == 1); Off_Times = find(States == 0);

Total_Times = 0; %Counter for how much stim time mouse got
for i = 1:length(On_Times)
    Total_Times = Total_Times + (Times(Off_Times(i))-Times(On_Times(i)));
end