%making figures for licks files
%for these traces, repeats are not combined

clear all; close all
load('all_photom_data.mat');
ind = find(contains(d.cat,'licks'));

access_time_by_lick=0;
%0-access is given a set amount of time
%1-access is given a set amount of time after first lick

rebaseline_to_lick_beginning=1;
%0-data is df/f and z-scored to baseline
%1-data is df/f and z-scored to the beginning of the lick bout

bt=4; %length of a bout (seconds)
hw=15; %half-window for graphing (seconds) and also baseline for rebaselining to beginning of lick bout
ili=2; %licks within ili are in a bout (seconds)

ds_psth=4; %amount to further downsample by for graphing the psth
ds_whole_trace=4; %amount to further downsample by for graphing the whole trace


for(i =ind)
    bouts={};
    
    access_start=d.access_time{i};
    ttl_time=d.ttl{i}(d.ttl{i}>access_start);
    if(access_time_by_lick)
        access_end=ttl_time(1)+d.access_length{i};
    else
        access_end=d.access_time{i}+d.access_length{i};
    end
    ttl_time(ttl_time>access_end)=[];
    d.ttl{i}=ttl_time;
    
    if(size(ttl_time,1)>0)
        bouts={ttl_time(1)};
        for i2=2:length(ttl_time);
            if ttl_time(i2)-ttl_time(i2-1)>ili; %if ttl on is further than ili from last ttl
                bouts{length(bouts)+1}=[ttl_time(i2)]; %add this ttl to a new bout
                
            else %if ttl on is closer than ili from last ttl
                bouts{length(bouts)}=[bouts{length(bouts)} ttl_time(i2)];%add this ttl to the old bout
            end
        end
        %all bouts are now at least one ili apart
        
        bouts2={};
        for i2=1:size(bouts,2)
            if((bouts{i2}(length(bouts{i2}))-bouts{i2}(1))>bt)
                bouts2{size(bouts2,2)+1}=bouts{i2};
            end
        end
        %all bouts are now at least one ili apart and at least bt long
        
        bouts3={};
        for i2=1:size(bouts2,2)
            if(bouts2{i2}(length(bouts2{i2}))<(access_end-hw))
                bouts3{size(bouts3,2)+1}=bouts2{i2};
            end
        end
        %all bouts end before inaccess-hw
        
        d.total_licks{i}=size(ttl_time,1);
        
        d.bouts{i}=bouts3;
        d.total_bouts{i}=size(bouts3,2);
    else
        d.bouts{i}=[];
    end
end

%now have bouts

licks=struct;

%must use string()
licks.mouse=[];
licks.experiment=[];
licks.full_trace=[];
licks.first_lick_start=[];
licks.first_lick_end=[];
licks.all_licks_start=[];
licks.all_licks_end=[];
licks.ttl={};
licks.firstlick=[];
licks.first_lick_start_z=[];
licks.first_lick_end_z=[];
licks.all_licks_start_z=[];
licks.all_licks_end_z=[];
licks.rate=[];
licks.access=[];
licks.repeat=[];

ind = find(contains(d.cat,'licks')& cellfun('size',d.bouts,2)>0);

for(i =ind)
    licks.mouse=[licks.mouse,string(d.mice{i})];
    licks.experiment=[licks.experiment,string(d.experiments{i})];
    licks.rate=[licks.rate,d.rate{i}];
    first_lick=d.bouts{1,i}{1}(1);
    licks.firstlick=[licks.firstlick;first_lick];
    licks.ttl{length(licks.ttl)+1}=d.ttl{i};
    
    licks.access=[licks.access,d.access_time{i}];
    
    rate=d.rate{i};
    licks.full_trace=[licks.full_trace,d.traces{i,1}];
    
    starts=[];
    ends=[];
    zstarts=[];
    zends=[];
    for(i2=1:size(d.bouts{1,i},2))
        
        
        mean_start=d.traces{i}(round((d.bouts{1,i}{i2}(1)-hw)*rate): round((d.bouts{1,i}{i2}(1)+hw)*rate));
        lick_beginning=mean(d.traces{i}(round((d.bouts{1,i}{i2}(1)-hw)*rate): round((d.bouts{1,i}{i2}(1))*rate)));
        lick_beginning_std=std(d.traces{i}(round((d.bouts{1,i}{i2}(1)-hw)*rate): round((d.bouts{1,i}{i2}(1))*rate)));
        mean_end=d.traces{i}(round((d.bouts{1,i}{i2}(size(d.bouts{1,i}{i2},2))-hw)*rate): round((d.bouts{1,i}{i2}(size(d.bouts{1,i}{i2},2))+hw)*rate));
        
        baseline_std=std(d.traces{i}(1:d.manipulation{i}));
        
        if(rebaseline_to_lick_beginning)
            starts=[starts,mean_start-lick_beginning];
            zstarts=[zstarts,(mean_start-lick_beginning)/lick_beginning_std];
            ends=[ends,mean_end-lick_beginning];
            zends=[zends,(mean_end-lick_beginning)/lick_beginning_std];
        else
            starts=[starts,mean_start];
            zstarts=[zstarts,mean_start/baseline_std];
            ends=[ends,mean_end];
            zends=[zends,mean_end/baseline_std];
        end
    end
    licks.first_lick_start=[licks.first_lick_start;starts(:,1)'];
    licks.first_lick_end=[licks.first_lick_end;ends(:,1)'];
    licks.all_licks_start=[licks.all_licks_start;mean(starts,2)'];
    licks.all_licks_end=[licks.all_licks_end;mean(ends,2)'];
    
    licks.first_lick_start_z=[licks.first_lick_start_z;zstarts(:,1)'];
    licks.first_lick_end_z=[licks.first_lick_end_z;zends(:,1)'];
    licks.all_licks_start_z=[licks.all_licks_start_z;mean(zstarts,2)'];
    licks.all_licks_end_z=[licks.all_licks_end_z;mean(zends,2)'];
    
    if(sum(contains(licks.mouse,string(d.mice{i}))&contains(licks.experiment,string(d.experiments{i})))>1)
        licks.repeat=[licks.repeat,1];
    else
        licks.repeat=[licks.repeat,0];
    end
end

exps=unique(licks.experiment);
close all
for(e=1:size(exps,2))
    ind=find(contains(licks.experiment, exps(e)) & licks.repeat==0);
    if(length(ind)>0)
        [x,ny1,ndy1]=JG_downsample_v3(licks.first_lick_start(ind,:),ds_psth,ds_psth);
        x=(x-1)*ds_psth/licks.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([mean(x) mean(x)], ylim)
        title(strcat(exps(e), " ","first lick bout start: dff trace"))
        xlabel('sec')
        ylabel('dff')
        
        [x,ny1,ndy1]=JG_downsample_v3(licks.first_lick_start_z(ind,:),ds_psth,ds_psth);
        x=(x-1)*ds_psth/licks.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([mean(x) mean(x)], ylim)
        title(strcat(exps(e), " ","first lick bout start: z trace"))
        xlabel('sec')
        ylabel('z')
        
        [x,ny1,ndy1]=JG_downsample_v3(licks.first_lick_end(ind,:),ds_psth,ds_psth);
        x=(x-1)*ds_psth/licks.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([mean(x) mean(x)], ylim)
        title(strcat(exps(e), " ","first lick bout end: dff trace"))
        xlabel('sec')
        ylabel('dff')
        
        [x,ny1,ndy1]=JG_downsample_v3(licks.first_lick_end_z(ind,:),ds_psth,ds_psth);
        x=(x-1)*ds_psth/licks.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([mean(x) mean(x)], ylim)
        title(strcat(exps(e), " ","first lick bout end: z trace"))
        xlabel('sec')
        ylabel('z')
        
        [x,ny1,ndy1]=JG_downsample_v3(licks.all_licks_start(ind,:),ds_psth,ds_psth);
        x=(x-1)*ds_psth/licks.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([mean(x) mean(x)], ylim)
        title(strcat(exps(e), " ","all lick bouts start: dff trace"))
        xlabel('sec')
        ylabel('dff')
        
        [x,ny1,ndy1]=JG_downsample_v3(licks.all_licks_start_z(ind,:),ds_psth,ds_psth);
        x=(x-1)*ds_psth/licks.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([mean(x) mean(x)], ylim)
        title(strcat(exps(e), " ","all lick bouts start: z trace"))
        xlabel('sec')
        ylabel('z')
        
        [x,ny1,ndy1]=JG_downsample_v3(licks.all_licks_end(ind,:),ds_psth,ds_psth);
        x=(x-1)*ds_psth/licks.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([mean(x) mean(x)], ylim)
        title(strcat(exps(e), " ","all lick bouts end: dff trace"))
        xlabel('sec')
        ylabel('dff')
        
        [x,ny1,ndy1]=JG_downsample_v3(licks.all_licks_end_z(ind,:),ds_psth,ds_psth);
        x=(x-1)*ds_psth/licks.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([mean(x) mean(x)], ylim)
        title(strcat(exps(e), " ","all lick bouts end: z trace"))
        xlabel('sec')
        ylabel('z')
        
        [x,ny1,ndy1]=JG_downsample_v3(licks.all_licks_start(ind,:),ds_psth,ds_psth);
        x=(x-1)*ds_psth/licks.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([mean(x) mean(x)], ylim)
        title(strcat(exps(e), " ","all lick bouts start: dff trace"))
        xlabel('sec')
        ylabel('dff')
        
        [x,ny1,ndy1]=JG_downsample_v3(licks.all_licks_start_z(ind,:),ds_psth,ds_psth);
        x=(x-1)*ds_psth/licks.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([mean(x) mean(x)], ylim)
        title(strcat(exps(e), " ","all lick bouts start: z trace"))
        xlabel('sec')
        ylabel('z')
        
        [x,ny1,ndy1]=JG_downsample_v3(licks.all_licks_end(ind,:),ds_psth,ds_psth);
        x=(x-1)*ds_psth/licks.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([mean(x) mean(x)], ylim)
        title(strcat(exps(e), " ","all lick bouts end: dff trace"))
        xlabel('sec')
        ylabel('dff')
        
        [x,ny1,ndy1]=JG_downsample_v3(licks.full_trace(:,ind)',ds_whole_trace,ds_whole_trace);
        x=x*ds_whole_trace/licks.rate(ind(1))/60;
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([licks.access(ind(1))/60 licks.access(ind(1))/60], ylim)
        title(strcat(exps(e), " ","whole trace: dff trace to baseline: access shown"))
        xlabel('min')
        ylabel('dff')
        
        y=licks.full_trace(:,ind)';
        y=y./std(y(:,1:(licks.access(ind(1))*licks.rate(ind(1))))')';
        [x,ny1,ndy1]=JG_downsample_v3(y,ds_whole_trace,ds_whole_trace);
        x=x*ds_whole_trace/licks.rate(ind(1))/60;
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([licks.access(ind(1))/60 licks.access(ind(1))/60], ylim)
        title(strcat(exps(e), " ","whole trace: z-scored to baseline: access shown"))
        xlabel('min')
        ylabel('dff')
        
        
        ys=[];
        x=[];
        for(i=ind)
            [x,y]=times_to_csum(licks.ttl{i},size(licks.full_trace(:,ind)',2)/licks.rate(i),1,0);
            ys=[ys;y];
        end
        
        [x,ny1,ndy1]=JG_downsample_v3(ys,1,1);
        x=x/60;
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([licks.access(ind(1))/60 licks.access(ind(1))/60], ylim)
        title(strcat(exps(e), " ","cumulative licks: access shown"))
        xlabel('min')
        ylabel('licks')
        
        
        num_licks2=[];
        for(i=ind)
            ttls=licks.ttl{i};
            num_licks=[];
            for(i2=(1:(size(licks.full_trace(:,ind)',2)/licks.rate(i))))
                num_licks=[num_licks,length(find(ttls<i2&ttls>i2-1))];
            end
            num_licks2=[num_licks2;num_licks];
        end
        
        [x,ny1,ndy1]=JG_downsample_v3(num_licks2,ds_whole_trace,ds_whole_trace);
        x=x*ds_whole_trace/60;
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        xlim([min(x) max(x)])
        hold on
        plot([licks.access(ind(1))/60 licks.access(ind(1))/60], ylim)
        title(strcat(exps(e), " ","binned licks: access shown"))
        xlabel('min')
        ylabel('hz')
        
    end
end

try
    mkdir('licks graphs')
end
home=cd;
cd('licks graphs');

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  set(0, 'CurrentFigure', FigHandle)
%   set(gcf,'renderer','Painters')
  FigName=get(get(gca, 'Title'),'string');
  FigName=regexprep(FigName, ' +', '_');
  FigName=regexprep(FigName, ',+', '');
  FigName=regexprep(FigName, ':+', '');
  savefig(strcat(FigName, '.fig'));
  saveas(gcf, FigName, 'svg');
  saveas(gcf, FigName, 'png');
  saveas(gcf, FigName, 'eps');
end

%Create a table with lick statistics for each mouse
Row_labels={'Mouse';'Experiment';'Access_Length';'Total_Licks';'Total_Bouts'};
Stats_alone=[d.mice;d.experiments;d.access_length;d.total_licks; d.total_bouts];
Stats_with_labels=[Row_labels Stats_alone];
writetable(cell2table(Stats_with_labels), 'Statistics.xlsx');

cd(home)