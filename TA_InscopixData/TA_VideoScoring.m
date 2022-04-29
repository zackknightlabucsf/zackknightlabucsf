%% VideoScoring 
% This script translates an excel csv file from hand-scoring a behavior video into
% "TTLs" that can be aligned with calcium traces. 
%
% close all; clear all; 
% This is run right after running the mdata_readmlist
% Because need act/inhib data for each individual mouse
csvfile = readtable('TA83_0524.csv'); 
door_open = 601; %At how many seconds was the door opened IN THE VIDEO

%% Loading stuff
Scoring.filename = table2cell(csvfile(1,1));
Scoring.VidLength = table2array(csvfile(1,5));
Scoring.BiteStart = table2array(csvfile(:, 12)); %start of a bite
Scoring.BiteEnd = table2array(csvfile(:,13)); %end of a bite

if Scoring.BiteStart == Scoring.BiteEnd %Check to see if there are bites that take longer
    fprintf('No long bites!'); %maybe add in a conditional variable here if needed
    Scoring.Bites = Scoring.BiteStart;
else
    fprintf('There are long bites!');
end

if Scoring.Bites(1) < 600 %if the video started late, correct for the deficit
    diff = 600-door_open; %there should be 600s before door
    Scoring.BitesCorrected = Scoring.Bites +diff;
else
    Scoring.BitesCorrected = Scoring.Bites;
end

%% Turn the time stamps into a binary vector

Scoring.BitesRound = round(Scoring.BitesCorrected); %round to the nearest second

%output vector
Bites = zeros(1, size(all_mrerank,2));
Bites(Scoring.BitesRound) = 1; %put a 1 where a bite happened

save(strcat(list{1},'_BiteScored.mat'));


