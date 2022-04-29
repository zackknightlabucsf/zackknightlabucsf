%computing a continuous lick rate with sliding window of length K

%Note, the code is a bit rough here, but the idea is to compute a cross
%correlation between the calcium trace and the "lick trace" 

clear all; close all; 
%infile1 = 'TA68_fastedPB_0523_BiteScored'; %for licking use AUC file; bites use BiteScored
%infile2 = 'TA82_fastedWater_AUC'; %compare versus water control or sucralose b/c water has little activation

Bites = 1; 
Licks = 0;
plot = 0;
win = 2;

list2 = {
  %'TA81_fastedEnsure_AUC'
  %'TA82_fastedEnsure_AUC'
  %'TA83_fastedEnsure_AUC'
  
  %'TA59_fastedPB2_0228_BiteScored'
  %'TA68_fastedPB_0523_BiteScored'
  %'TA79_fastedPB_0921_BiteScored'
  %'TA80_fastedPB_0921_BiteScored'
  
  'TA68_fastedObj_0526_BiteScored'
  'TA79_fastedLego_0304_BiteScored'
};

if Licks == 1
    fr = 4;
    win = win*fr; %window length
elseif Bites == 1
    win = win; %don't include fr if doing bites b/c vector is diff length
end

%% Computations
%Act_XCorr.c = []; Act_XCorr.lags = [];
%Inhib_XCorr.c = []; Inhib_XCorr.lags = [];
for p = 1:length(list2)
    
    load(list2{p});
    if exist('Bites')
        ttl_clean = Bites;
    end

    Y = movmean(ttl_clean,win); %where X is the vector of licks, and the second number is the window length(K)

    x = mean(Neurons_act,1);
    if length(x) ~= length(Y)
        Y = Y(1:length(x)); %need both vectors to be same length; only look at bites during session (exclude end)
    end
    
    [c, lags] = xcorr(x,Y, 'coeff'); %where x and y are Calcium traces (Neuron_act) and Y from above
    Act_Xcorr{p, 1} = c; Act_Xcorr{p,2} = lags; 
    %Act_XCorr.lags = [Act_XCorr.lags; lags];
    
    x2 = mean(Neurons_inhib,1);
    [c2, lags2] = xcorr(x2, Y, 'coeff');
    Inhib_Xcorr{p, 1} = c2; Inhib_Xcorr{p,2} = lags2;
    %Inhib_XCorr.lags = [Inhib_XCorr.lags; lags2];
    
    clearvars -except list2 Act_Xcorr Inhib_Xcorr p fr win plot
end


%% Then plot r
if plot == 1
    figure;
    subplot(3,1,1);
    plot(1:length(Neurons_act), x, 'r'); hold on;
    plot(1:length(Neurons_inhib), x2, 'b'); hold on; title('Ca signal');
    
    subplot(3,1,2);
    plot(1:length(Y), Y, 'r'); hold on; title('bites'); hold on;
    plot(1:length(ttl_clean), ttl_clean, 'k'); hold on;
    
    subplot(3,1,3);
    plot(lags, c, 'r'); hold on; %plot([mean(lags) mean(lags)], [min(c) max(c)], 'r--', 'LineWidth', 2); hold on;
    plot(lags2, c2, 'b'); hold on; %plot([mean(lags2) mean(lags2)], [min(c2) max(c2)], 'b--', 'LineWidth', 2); hold on;
    legend({'Activated', 'Inhibited'}, 'FontSize', 12);
    title('Cross Correlation'); xlabel('lag (s)'); ylabel('cross-corr (a.u.)');%ylim([0 max([c, c2])]);
    
end

