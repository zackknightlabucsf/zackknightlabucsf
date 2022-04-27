clear all;
close all;
%% Importing Synapse files into MATLAB
addpath(genpath('E:\'));  %Add TDT scripts to path
%cd 'E:\'
%COPY and PASTE FOLDER NAME OF RECORDING BELOW
raw_data=TDTbin2mat('TL264_TL266-220322-144253','TYPE',{'epocs','streams'});  %Converts Synapse data folder into .mat data structure

%% Assembling data objects
%CHANGE THESE PARAMETERS FOR EACH PHOTOMETRY RECORDING
ds_rate=4;  %Downsample rate in Hz, 4 is preferred DONT CHANGE
recording_time=30; %define length of recording time after stimulus (in minutes)  
stimulus_z=1224/60; %define stimulus start time for z chamber  (input second to get minutes) 
stimulus_k=1229/60;   %define stimulus start time for k chamber  (input seconds to get minutes)
animal_1='TL264-FD-LickGlucoselow-032222';   %animal ID for top Z chamber
animal_2='TL266-FD-LickGlucoselow-032222';    %animal ID for bottom K chamber

%DONT CHANGE PARAMETERS BELOW
baseline_start_z=stimulus_z-10;    %define beginning of 10 minute baseline period
baseline_start_k=stimulus_k-10;    %define beginning of 10 minute baseline period
start_z=stimulus_z-20; %define start of time window for top z chamber(in minutes)
start_k=stimulus_k-20;    %define start of time window for top z chamber(in minutes)
end_time_z=stimulus_z+recording_time; %define end of recording for top z chamber (in minutes)
end_time_k=stimulus_k+recording_time;  %define end of recording for top z chamber (in minutes

 G_Z=downsample(raw_data.streams.G__Z.data,round(raw_data.streams.G__Z.fs/ds_rate));
 G_Z(:,1)=[];
 U_Z=downsample(raw_data.streams.U__Z.data,round(raw_data.streams.U__Z.fs/ds_rate));
 U_Z(:,1)=[];
 G_Z(:,round(end_time_z*60*4)+2:size(G_Z,2))=[];
 U_Z(:,round(end_time_z*60*4)+2:size(U_Z,2))=[];
   
 index_norm_z=round(start_z*60*4+1):round(stimulus_z*60*4);
 para_fit_z=polyfit(U_Z(index_norm_z),G_Z(index_norm_z),1); 
 Gfit_Z = para_fit_z(1).*U_Z + para_fit_z(2);
   
 trace_z=G_Z./Gfit_Z;   %traces for Z chamber 
 
 G_K=downsample(raw_data.streams.G__K.data,round(raw_data.streams.G__K.fs/ds_rate));
  G_K(:,1)=[];
 U_K=downsample(raw_data.streams.U__K.data,round(raw_data.streams.U__K.fs/ds_rate));
  U_K(:,1)=[];
 G_K(:,round(end_time_k*60*4)+2:size(G_K,2))=[];
 U_K(:,round(end_time_k*60*4)+2:size(U_K,2))=[];
   
 index_norm_k=round(start_k*60*4+1):round(stimulus_k*60*4);
 para_fit_k=polyfit(U_K(index_norm_k),G_K(index_norm_k),1); 
 Gfit_K = para_fit_k(1).*U_K + para_fit_k(2);
   
trace_k=G_K./Gfit_K;  %traces for K chamber

%Lick information
if isfield(raw_data.epocs,'Ep2_') == 1
 lick_data_1=raw_data.epocs.Ep2_.onset;     % Z lickometer
 lick_data_1_normalized=lick_data_1*ds_rate;    %normalized to 4Hz recording frequency of fluorescence trace 
end
if isfield(raw_data.epocs,'Ep4_') == 1
 lick_data_2=raw_data.epocs.Ep4_.onset;     % K  Lickometer
 lick_data_2_normalized=lick_data_2*ds_rate;    %normalized to 4Hz recording frequency of fluorescence trace
end
%% Normalizing Z and K traces for Delta F/F0
mean_baseline_trace_z=mean(trace_z(round(baseline_start_z*60*4):round(stimulus_z*60*4))); %mean baseline F0 between time -10 and 0 min
mean_baseline_trace_k=mean(trace_k(round(baseline_start_k*60*4):round(stimulus_k*60*4)));   %mean baseline F0 between time -10 and 0 min
deltaF_trace_z=(trace_z-mean_baseline_trace_z)/mean_baseline_trace_z*100;
deltaF_trace_k=(trace_k-mean_baseline_trace_k)/mean_baseline_trace_k*100;

%% Normalizing Z and K traces for Z-score
sd_baseline_trace_z=std(trace_z(round(baseline_start_z*60*4):round(stimulus_z*60*4)));
sd_baseline_trace_k=std(trace_k(round(baseline_start_k*60*4):round(stimulus_k*60*4))); 
z_score_trace_z=(trace_z-mean_baseline_trace_z)/sd_baseline_trace_z;
z_score_trace_k=(trace_k-mean_baseline_trace_k)/sd_baseline_trace_k;

%% Plotting figures
% Plot individual animals and save PDFs
time=[-10+1/(60*ds_rate):1/(60*ds_rate):recording_time];   %time between -10 minutes and 60 minutes 
y=deltaF_trace_z((round(baseline_start_z*60*4))+1:round((stimulus_z+recording_time)*60*4));      %Delta F/F for first animal from folder
y_zscore=z_score_trace_z((round(baseline_start_z*60*4))+1:round((stimulus_z+recording_time)*60*4)); 
y2=deltaF_trace_k((round(baseline_start_k*60*4))+1:round((stimulus_k+recording_time)*60*4));       %Delta F/F for second animal from folder
y2_zscore=z_score_trace_k((round(baseline_start_k*60*4))+1:round((stimulus_k+recording_time)*60*4));
time2=[0+1/(60*ds_rate):1/(60*ds_rate):recording_time]; %lick start time


figure(1)
%plot Gcamp Trace of first animal
subplot(2,1,1);
 if isfield(raw_data.epocs,'Ep2_') == 1
       yyaxis left
 end
plot(time,smooth(y,20),'b','LineWidth',1.5)
%patch([0 10 10 0],[-10 -10 25 25],'k','EdgeColor','none','FaceAlpha',[0.05]); 
set(gcf,'position',[100,100,700,600])
set(gca,'LineWidth',2,'FontSize',12,'TickDir','out');
box off
xlabel('Time (min)','FontWeight','bold');
ylabel('\DeltaF/F (%)','FontWeight','bold');
ylim([min(y)+(min(y)/4) max(y)+(max(y)/4)]);
%ylim([-30 200]);
xlim([-10 recording_time]);
hold on

%plot licks
 if isfield(raw_data.epocs,'Ep2_') == 1
    edges =[round(stimulus_z*60*4):1:round(end_time_z*60*4)];
    N = histcounts(lick_data_1_normalized,edges);
    smooth_licks=smooth(N);
    yyaxis right
    ylabel('Lick Rate (Hz)','FontWeight','bold'); 
    plot(time2,smooth_licks);
 end
title(animal_1);
%Plot z-score trace for first animal
subplot(2,1,2);
 if isfield(raw_data.epocs,'Ep2_') == 1
       yyaxis left
 end
plot(time,smooth(y_zscore,20),'b','LineWidth',1.5)
%patch([0 10 10 0],[-10 -10 25 25],'k','EdgeColor','none','FaceAlpha',[0.05]); 
set(gcf,'position',[100,100,700,600])
set(gca,'LineWidth',1.5,'FontSize',12,'TickDir','out');
box off
xlabel('Time (min)','FontWeight','bold');
ylabel('Z-score','FontWeight','bold');
ylim([min(y_zscore)+(min(y_zscore)/4) max(y_zscore)+(max(y_zscore)/4)]);
%ylim([-30 200]);
xlim([-10 recording_time]);
hold on

%plot licks
 if isfield(raw_data.epocs,'Ep2_') == 1
    edges =[round(stimulus_z*60*4):1:round(end_time_z*60*4)];
    N = histcounts(lick_data_1_normalized,edges);
    smooth_licks=smooth(N);
    yyaxis right
    ylabel('Lick Rate (Hz)','FontWeight','bold'); 
    plot(time2,smooth_licks);
 end
title(animal_1);
%annotation('textbox', [0.175 0.825 .1 .1],'String','Intragastric Infusion','FontWeight','bold','LineStyle','none');
hold off 
saveas(gcf,animal_1,'pdf');

figure(2)
%plot Gcamp Trace of second animal
subplot(2,1,1);
if isfield(raw_data.epocs,'Ep4_') == 1
    yyaxis left
end
plot(time,smooth(y2,20),'b','LineWidth',1.5)
%patch([0 10 10 0],[-10 -10 25 25],'k','EdgeColor','none','FaceAlpha',[0.05]); 
set(gcf,'position',[100,100,700,600])
set(gca,'LineWidth',1,'FontSize',12,'TickDir','out');
box off
xlabel('Time (min)','FontWeight','bold');
ylabel('\DeltaF/F (%)','FontWeight','bold');
xlim([-10 recording_time]);
%ylim([-30 200]);
ylim([min(y2)+(min(y2)/4) max(y2)+(max(y2)/4)]);
hold on

%plot licks
if isfield(raw_data.epocs,'Ep4_') == 1
    edges =[round(stimulus_k*60*4):1:round(end_time_k*60*4)];
    N = histcounts(lick_data_2_normalized,edges);
    smooth_licks=smooth(N);
    yyaxis right
    ylabel('Lick Rate (Hz)','FontWeight','bold'); 
    plot(time2,smooth_licks);
end
title(animal_2);

%Plot z-score trace of second animal
subplot(2,1,2);
if isfield(raw_data.epocs,'Ep4_') == 1
    yyaxis left
end
plot(time,smooth(y2_zscore,20),'b','LineWidth',1.5)
%patch([0 10 10 0],[-10 -10 25 25],'k','EdgeColor','none','FaceAlpha',[0.05]); 
set(gcf,'position',[100,100,700,600])
set(gca,'LineWidth',1,'FontSize',12,'TickDir','out');
box off
xlabel('Time (min)','FontWeight','bold');
ylabel('Z-score','FontWeight','bold');
xlim([-10 recording_time]);
%ylim([-10 20]);
ylim([min(y2_zscore)+(min(y2_zscore)/4) max(y2_zscore)+(max(y2_zscore)/4)]);
hold on

%plot licks
if isfield(raw_data.epocs,'Ep4_') == 1
    edges =[round(stimulus_k*60*4):1:round(end_time_k*60*4)];
    N = histcounts(lick_data_2_normalized,edges);
    smooth_licks=smooth(N);
    yyaxis right
    ylabel('Lick Rate (Hz)','FontWeight','bold'); 
    plot(time2,smooth_licks);
end
title(animal_2);
%annotation('textbox', [0.175 0.825 .1 .1],'String','Intragastric Infusion','FontWeight','bold','LineStyle','none');
hold off
saveas(gcf,animal_2,'pdf');


% % Plot average trace from multiple animals
% all_fluorescence_traces=cat(1,y,y2);  %matrix with fluorescence trace of all animals
% average_fluorescence_trace=mean(all_fluorescence_traces);  %mean fluorescence trace of all animals
% sem_average_fluorescence=std(all_fluorescence_traces)/size(all_fluorescence_traces,1)^0.5;
% x1=x;
% x2=fliplr(x);
% Y1=average_fluorescence_trace-sem_average_fluorescence;
% Y2=average_fluorescence_trace+sem_average_fluorescence;
% Y3=fliplr(Y2);
% 
% figure(3)
% patch([x1 x2],[Y1 Y3],'b','EdgeColor','none','FaceAlpha',[0.05]);
% hold on 
% plot(x,average_fluorescence_trace,'b','LineWidth',1);
% hold on
% patch([0 10 10 0],[-10 -10 25 25],'k','EdgeColor','none','FaceAlpha',[0.05]); 
% set(gca,'LineWidth',1,'FontSize',12,'TickDir','out');
% box off
% xlabel('Time (min)','FontWeight','bold');
% ylabel('\DeltaF/F (%)','FontWeight','bold');
% xlim([-10 60]);
% ylim([-10 25]);
% annotation('textbox', [0.175 0.825 .1 .1],'String','Intragastric Infusion','FontWeight','bold','LineStyle','none');
% hold off
