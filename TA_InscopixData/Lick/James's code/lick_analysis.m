close all
clear all
fclose all

%nb. this script requires suplabel.m

%addpath('/Users/jamesgrove/Documents/MATLAB/cbrewer') %comment out if using windows

%this is a list of csv files
list={
    %'licks_fastedEnsure_TA58.csv' %only change this
    %'licks_depWater_TA58'
    %'licks_depWater_TA59'
    %'licks_fedWater_TA58.csv'
    %'licks_fedWater_TA59'
    %'licks_fastedEnsure_TA59.csv'
    %'licks_fedWater2_TA58.csv'
    %'licks_fastedDry_TA68.csv'
    %'licks_fastedSucralose_TA68.csv'
    %'licks_fastedEnsureSpaced2_TA124.csv'
    'licks_fastedSucraloseSpaced_TA120.csv'
};

%if the following don't matter to you, ignore (the function won't call it)
zscoreyn=1;% zscoreyn: 0 - don't zscore; 1 - zscore to baselin; 2 - zscore to whole trial
rawyn=1; %normally set to 1
hw_sec=10; %bout threshold in secs
fr=4; %framerate
boutthreshold_sec=hw_sec; %bout time in sec (nb. before was frame count)
firstbout=1; %0 - don't include; 1 - include; 2 - firstbout only
crossreffile=0; %0 - not a crossref file; 1 - crossref file
subtract_baseline=1; %0 - subtract baseline; 1 - subtract activity before first lick; 2-subtract activity before first and last lick separately

%%
crind={};
[cmap]=cbrewer('div', 'RdBu', 256); %need to get cbrewer function
cmap=flip(cmap);

if(zscoreyn==1)
   zscoreyn_words=char('zscored to prestim');
elseif(zscoreyn==2)
    zscoreyn_words=char('zscored to whole trial');
elseif(zscoreyn==0)
    zscoreyn_words=char('not zscored');
end

if(firstbout==1)
   firstbout_words=char('all bouts included')
elseif(firstbout==2)
    firstbout_words=char('first bout only')
elseif(firstbout==0)
    firstbout_words=char('all bouts but first')
end


for i=1:length(list)
    
    maxcv=20;
    mincv=-20;
    if zscoreyn
        maxcv=maxcv/5;
        mincv=mincv/5;
    end
    
    %make sure the file ends in csv
    if(endsWith(list{i},'.csv'))
        csvfile_name=list{i}
    else
        csvfile_name=strcat(list{i}, '.csv');
    end
    
    csvfile_name = list{i};
    %check if file exists
    fid=fopen(csvfile_name);
    if(fid==-1)
       fprintf(strcat(list{i}, ' could not be found. Double check spelling.\n'));
       return
    end
    clear fid
    
    [all_psth_start, all_psth_end,ttl_clean,cpsthrange, pretime, posttime, timestamps]=read_for_inscopix_analysis_JG( csvfile_name, rawyn, hw_sec, zscoreyn, boutthreshold_sec, firstbout, subtract_baseline);
    
    tot_start=[];
    tot_dff=[];
    tot_end=[];
    for(r=1:size(all_psth_start,1))
        tot_start=[tot_start; all_psth_start(r).avgpsth];
        tot_dff=[tot_dff; all_psth_start(r).dff];
        tot_end=[tot_end; all_psth_end(r).avgpsth];
    end
    
    [~, nind]=sort(tot_dff);
    tot_start=tot_start(nind,:);
    tot_end=tot_end(nind,:);
    
    figure
    h=heatmap_d(tot_start,[],[],[],'MaxColorValue',maxcv/2,'MinColorValue',mincv/2,'Colormap',cmap,'Colorbar',1);
    xticks([0, 5,10,15,20,25,30,35,40]*4+1)
    xticklabels([0, 5,10,15,20,25,30,35,40])
    yticks([1])
    yticklabels([size(tot_end,2)])
    hold on
    plot([size(tot_end,2)/2-1 size(tot_end,2)/2-1], [0 size(tot_end,1)], 'color', 'black')
    hold off
    title1=('Beginning of lick bout');
    if(firstbout==1)
        title1={title1; 'first bout included'};
    elseif(firstbout==2)
        title1={title1; 'first bout only'};
    else
        title1={title1; 'first bout not included'};
    end
    
    if(subtract_baseline==0)
        title1={title1; 'baseline activity subtracted'};
    else
        title1={title1; 'pre-bout activity subtracted'};
    end
    title('Beginning of lick bout');

    figure
    h=heatmap_d(tot_end,[],[],[],'MaxColorValue',maxcv/2,'MinColorValue',mincv/2,'Colormap',cmap,'Colorbar',1);
    xticks([0, 5,10,15,20,25,30,35,40]*4+1)
    xticklabels([0, 5,10,15,20,25,30,35,40])
    yticks([1])
    yticklabels([size(tot_end,2)])
    hold on
    plot([size(tot_end,2)/2-1 size(tot_end,2)/2-1], [0 size(tot_end,1)], 'color', 'black')
    hold off
    title1=('End of lick bout');
    if(firstbout==1)
        title1={title1; 'first bout included'};
    elseif(firstbout==2)
        title1={title1; 'first bout only'};
    else
        title1={title1; 'first bout not included'};
    end
    
    if(subtract_baseline==0)
        title1={title1; 'baseline activity subtracted'};
    elseif(subtract_baseline==1)
        title1={title1; 'pre-bout activity subtracted'};
    else
        title1=[title1; 'bout activity subtracted'];
    end
    title('End of lick bout');
    
    tot_zC=[];
    zC_dff=[];
    for(r=1:size(all_psth_start,1))
        ts=timestamps(r);
        pt2=posttime(r)+ts;
        for(n=1:size(all_psth_start(r).zC,1))
            zC_dff=[zC_dff;mean(all_psth_start(r).zC(n,ts:pt2))];
        end
        delta=ts-min(timestamps)+1;
        pt2=ts+min(posttime);
        tot_zC=[tot_zC; all_psth_start(r).zC(:,delta:pt2)];
    end
    
    [~, nind2]=sort(zC_dff);
    tot_zC2=tot_zC(nind2,:);
    figure
    h=heatmap_d(tot_zC2,[],[],[],'MaxColorValue',maxcv,'MinColorValue',mincv,'Colormap',cmap,'Colorbar',1);
    xticks([0, 10, 20, 30 ,40, 50, 60]*4*60)
    xticklabels([0, 10, 20, 30 ,40, 50, 60])
    yticks([1])
    yticklabels([size(tot_end,2)])
    hold on
    plot([min(timestamps) min(timestamps)], [0 size(tot_end,1)], 'color', 'black')
    hold off
    title('Activity following access, whole trial, sorted by activity after access')
    
    [~, nind]=sort(tot_dff);
    tot_zC2=tot_zC(nind,:);
    figure
    h=heatmap_d(tot_zC2,[],[],[],'MaxColorValue',maxcv,'MinColorValue',mincv,'Colormap',cmap,'Colorbar',1);
    xticks([0, 10, 20, 30 ,40, 50, 60]*4*60)
    xticklabels([0, 10, 20, 30 ,40, 50, 60])
    yticks([1])
    yticklabels([size(tot_end,2)])
    hold on
    plot([min(timestamps) min(timestamps)], [0 size(tot_end,1)], 'color', 'black')
    hold off
    title1=('Activity following access, whole trial,');
    if(firstbout==1)
        title1={title1; 'sorted by activity during all lick bouts'};
    elseif(firstbout==2)
        title1={title1; 'sorted by activity during first lick bout'};
    else
        title1={title1; 'sorted by activity during all lick bouts except the first'};
    end
    title('Activity following access, whole trial,')
    
    

end
