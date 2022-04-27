%Lickometer analysis
%By: James Grove
%Adapted from scripts by: Yiming Chen

close all;clear all;fclose all

%nb. this script requires suplabel.m, cbrewer.m, JG_cmap_v1, suplabel, and
%CNMFe master folder

%this script will look at 2 timescales: activity changes during a lick bout
%and activity changes after access is removed


%this is a list of csv files containing metadata
list={
    'metadata.csv'
    };


zscoreyn=1;%0 - don't zscore; 1 - zscore to pre-stim; 2 - zscore to whole trial
baseline=0 %0-baseline licks to prestim; 1-baseline licks to pre-first lick
%what it is baselined to, is what dff is calculated using
hw_sec=60; %psth window size in seconds
fr=4; %framerate in Hz
crossreffile=0 %1 - contains crossregistered neurons; 0 - doesn't
firstbout=2; %0 - don't include; 1 - include; 2 - firstbout only

%a bout is defined as a series of licks >boutthreshold_sec with pauses <ili_sec
boutthreshold_sec=10; %bout time threshold (in sec)
ili_sec=2 %maximum time between licks in a bout (in sec)

%don't change the script after this line

%%
crind={};
[cmap]=cbrewer('div', 'RdBu', 256);
cmap=flip(cmap);
rawyn=1

if(zscoreyn==1)
    zscoreyn_words=char("zscored to prestim")
elseif(zscoreyn==2)
    zscoreyn_words=char("zscored to whole trial")
elseif(zscoreyn==0)
    zscoreyn_words=char("not zscored")
end

if(firstbout==1)
    firstbout_words=char("all bouts included")
elseif(firstbout==2)
    firstbout_words=char("first bout only")
elseif(firstbout==0)
    firstbout_words=char("all bouts but first")
end

if(baseline==0)
    baseline_words=char("baselined to prestim")
elseif(baseline==1)
    baseline_words=char("baselined to pre-first-lick")
end

for i=1:length(list)
    
    maxcv=20;
    mincv=-20;
    if zscoreyn
        maxcv=4;
        mincv=-4;
    end
    
    %make sure the file ends in csv
    if(endsWith(list{i},'.csv'))
        csvfile_name=list{i}
    else
        csvfile_name=strcat(list{i}, '.csv');
    end
    
    %check if file exists
    fid=fopen(csvfile_name);
    if(fid==-1)
        fprintf(strcat(list{i}, ' could not be found. Double check spelling.\n'));
        return
    end
    clear fid
    
    [all_psth_start, all_psth_end,ttl_clean,cpsthrange,all_ttl]=read_for_inscopix_analysis_JG( csvfile_name, rawyn, hw_sec, zscoreyn, boutthreshold_sec, ili_sec, firstbout, baseline);
    
    num_licks2=[];
    for(i=1:size(all_ttl,1))
        this=[];
        for t=0:1800
            this=[this,nansum(all_ttl(i,1:(4*t+1)))];
        end
        num_licks2=[num_licks2;this];
    end
    y=mean(num_licks2,1);
    dy=std(num_licks2,1)/sqrt(size(num_licks2,1));
    x=1:size(y,2);
    x=x/60;
    figure
    y1=y+dy;
    y2=y-dy;
    X=[x,fliplr(x)];              %#create continuous x value array for plotting
    Y=[y1,fliplr(y2)];              %#create y values for out and then back
    fill(X,Y,[.5 .5 .5],'linestyle','none');
    hold on
    plot(x,y, 'Color', [0 0 0])
    title("cumulative licks, water")
    xlim([0 x(size(y,2))])
    ylim([0 400])
    
    
    
    if(crossreffile==0)
        avgpsth_start=[];
        avgpsth_end=[];
        dff=[];
        dff_total=[];
        zC_timed_to_first_lick=[];
        zS_timed_to_first_lick=[];
        C_timed_to_first_lick=[];
        S_timed_to_first_lick=[];
        S_whole_trial=[];
        C_whole_trial=[];
        coor=vertcat(all_psth_start.coor);
        zC=[];
        for(r=1:size(all_psth_start,1))
            avgpsth_start=[avgpsth_start;all_psth_start(r).avgpsth];
            avgpsth_end=[avgpsth_end;all_psth_end(r).avgpsth];
            dff=[dff;all_psth_start(r).dff];
            dff_total=[dff_total;all_psth_start(r).dff_total];
            zC=[zC;all_psth_start(r).zC];
            zC_timed_to_first_lick=[zC_timed_to_first_lick;all_psth_start(r).zC_timed_to_first_lick];
            C_timed_to_first_lick=[C_timed_to_first_lick;all_psth_start(r).C_timed_to_first_lick];
            
            S_timed_to_first_lick=[S_timed_to_first_lick;all_psth_start(r).S_timed_to_first_lick];
            S_whole_trial=[S_whole_trial;all_psth_start(r).S_whole_trial(:,1:10*4*60)];
                        C_whole_trial=[C_whole_trial;all_psth_start(r).C_whole_trial(:,1:10*4*60)];

            zS_timed_to_first_lick=[zS_timed_to_first_lick;all_psth_start(r).zS_timed_to_first_lick];
            
        end
        
        [~,ind]=unique(avgpsth_start,'rows');
        avgpsth_start=avgpsth_start(ind,:);
        avgpsth_end=avgpsth_end(ind,:);
        dff=dff(ind,:);
        
     
        
        
        [~,ind]=unique(zC,'rows');
        zC=zC(ind,:);
        zC_timed_to_first_lick=zC_timed_to_first_lick(ind,:);
        dff_total=dff_total(ind,:);
        
        [~,rerank1]=sort(dff);
        [~,rerank2]=sort(dff_total);
        
        avgpsth_start=avgpsth_start(rerank1,:);
        avgpsth_end=avgpsth_end(rerank1,:);
        
        zC=zC(rerank2,:);
        zC_timed_to_first_lick_old=zC_timed_to_first_lick;
        zC_timed_to_first_lick=zC_timed_to_first_lick(rerank2,:);
        %% 
        %Create pie charts of activity changes during lick bout and after
        %access is removed
        
        figure
        pie([length(find(dff>0.5)) length(find(dff<-0.5)) length(find(dff<0.5 & dff>-0.5))])
        colormap([1 0 0;0 0 1; 0.5 0.5 0.5])
        title(strcat(all_psth_start(1).experiment, ", licks effect: ", zscoreyn_words, ', ', firstbout_words, ', ', baseline_words))

        figure
        pie([length(find(dff_total>0.5)) length(find(dff_total<-0.5)) length(find(dff_total<0.5 & dff_total>-0.5))])
        colormap([1 0 0;0 0 1; 0.5 0.5 0.5])
        title(strcat(all_psth_start(1).experiment, ", long-term effect: ", zscoreyn_words, ', ', firstbout_words, ', ', baseline_words))

        %%
        %Making heatmaps of lick bout activity
        
        zmax=maxcv;
        coloron=zmax;
        coloroff=1;
        [cmap]=JG_cmap_v1(zmax,coloron,coloroff);

        figure()
        subplot(121)
        h=heatmap_d(avgpsth_start,[],[],[],'MaxColorValue',maxcv,'MinColorValue',mincv,'Colormap',cmap,'Colorbar',1);
        title('start')
        yticks([1 size(avgpsth_start,1)])
        yticklabels([1 size(avgpsth_start,1)])
        subplot(122)
        h=heatmap_d(avgpsth_end,[],[],[],'MaxColorValue',maxcv,'MinColorValue',mincv,'Colormap',cmap,'Colorbar',1);
        title('end')
        [ax,h3]=suplabel([all_psth_start(1).experiment, ', ranked by start, ' zscoreyn_words, ', ', firstbout_words, ', ', baseline_words]  ,'t');
        yticks([1 size(avgpsth_end,1)])
        yticklabels([1 size(avgpsth_end,1)])
        xticks([1 ((size(avgpsth_start,2)-1)/2)+1 size(avgpsth_start,2)])
        xticklabels([-((size(avgpsth_start,2)-1)/2/4) 0 ((size(avgpsth_start,2)-1)/2/4)])
        subplot(122)
        yticks([1 size(avgpsth_end,1)])
        yticklabels([1 size(avgpsth_end,1)])
        xticks([1 ((size(avgpsth_start,2)-1)/2)+1 size(avgpsth_start,2)])
        xticklabels([-((size(avgpsth_start,2)-1)/2/4) 0 ((size(avgpsth_start,2)-1)/2/4)])
        
        %%
        %Making average time series of lick bout activity
        
        yrange=[-1 3];
        yts=(-1:2);
        y=avgpsth_start;
        filterby=4;
        ds=4; %false
        [x,ny1,ndy1]=JG_downsample_v3(y,filterby,ds);
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.2 .9 .9],'linestyle','none');
        hold on
        plot(ny1);set(gca, 'YLim', yrange)
        title(strcat(all_psth_start(1).experiment, ", licks start, filtered, 1s, ", zscoreyn_words, ', ', firstbout_words, ', ', baseline_words))
        xlim([1 size(x,1)])
        xticks([ds ((size(avgpsth_start,2)-1)/2)+1 size(avgpsth_start,2)]/ds)
        xticklabels([-((size(avgpsth_start,2)-1)/2/4) 0 ((size(avgpsth_start,2)-1)/2/4)])
        yticks(yts)
        
        y=avgpsth_end;
        filterby=4;
        ds=4; %false
        [x,ny1,ndy1]=JG_downsample_v3(y,filterby,ds);
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.2 .9 .9],'linestyle','none');
        hold on
        plot(ny1);set(gca, 'YLim', yrange)
        title(strcat(all_psth_start(1).experiment, ", licks end, filtered, 1s, ", zscoreyn_words, ', ', firstbout_words, ', ', baseline_words))
        xlim([1 size(x,1)])
        xticks([ds ((size(avgpsth_start,2)-1)/2)+1 size(avgpsth_start,2)]/ds)
        xticklabels([-((size(avgpsth_start,2)-1)/2/4) 0 ((size(avgpsth_start,2)-1)/2/4)])
        yticks(yts)
        
        %%
        %Making heatmaps of whole trial activity
        
        figure()
        h=heatmap_d(zC,[],[],[],'MaxColorValue',maxcv,'MinColorValue',mincv,'Colormap',cmap,'Colorbar',1);
        title(strcat(all_psth_start(1).experiment, ', whole trial, timed to access'))
        yticks([1 size(zC,1)])
        yticklabels([1 size(zC,1)])
        
        figure()
        h=heatmap_d(zC_timed_to_first_lick,[],[],[],'MaxColorValue',maxcv,'MinColorValue',mincv,'Colormap',cmap,'Colorbar',1);
        title(strcat(all_psth_start(1).experiment, ', whole trial, timed to first lick'))
        yticks([1 size(zC_timed_to_first_lick,1)])
        yticklabels([1 size(zC_timed_to_first_lick,1)])
        
        %%
        %Making average time series of whole trial activity

        y=zC;
        filterby=20;
        ds=20; %true
        [x,ny1,ndy1]=JG_downsample_v3(y,filterby,ds);
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.2 .9 .9],'linestyle','none');
        hold on
        plot(ny1);set(gca, 'YLim', yrange)
        title(strcat(all_psth_start(1).experiment, ", whole trial, timed to access, downsampled & filtered, 5s"))
        xlim([1 size(x,1)])
        xticks([5,10,15,20,25,30,35,40]*60*4/ds)
        xticklabels([5,10,15,20,25,30,35,40])
        yticks(yts)
        ylim([floor(min(ny1+ndy1)) ceil(max(ny1+ndy1))])
        
        y=zC_timed_to_first_lick;
        filterby=20;
        ds=20; %true
        [x,ny1,ndy1]=JG_downsample_v3(y,filterby,ds);
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.2 .9 .9],'linestyle','none');
        hold on
        plot(ny1);set(gca, 'YLim', yrange)
        title(strcat(all_psth_start(1).experiment, ", whole trial, timed to first lick bout, downsampled & filtered, 5s"))
        xlim([1 size(x,1)])
        xticks([5,10,15,20,25,30,35,40]*60*4/ds)
        xticklabels([5,10,15,20,25,30,35,40])
        yticks(yts)
        ylim([floor(min(ny1+ndy1)) ceil(max(ny1+ndy1))])
        
        %%
        %Making average time series of neurons inhibited and excited over
        %long-term
        
        figure
        y=zC_timed_to_first_lick_old(find(dff_total>0.5),:);
        filterby=20;
        ds=4
        [nx,ny,ndy]=JG_downsample_v3(y,filterby,ds);
        fill([nx;flipud(nx)],[ny-ndy;flipud(ny+ndy)],[.9 .2 .2],'linestyle','none');
        hold on
        plot(ny)
        xticks([5,10,15,20,25,30,35,40,45]*60*4/ds)
        xticklabels([5,10,15,20,25,30,35,40,45])
        act=size(y,1);
        hold on
        y=zC_timed_to_first_lick_old(find(dff_total<-0.5),:);
        filterby=20;
        ds=4
        [nx,ny,ndy]=JG_downsample_v3(y,filterby,ds);
        fill([nx;flipud(nx)],[ny-ndy;flipud(ny+ndy)],[.2 .2 .9],'linestyle','none');
        hold on
        plot(ny)
        xticks([5,10,15,20,25,30,35,40,45]*60*4/ds)
        xticklabels([5,10,15,20,25,30,35,40,45])
        inhib=size(y,1);
        hold on
        y=zC_timed_to_first_lick_old(find(dff_total>-0.5&dff_total<0.5),:);
        filterby=20;
        ds=4
        [nx,ny,ndy]=JG_downsample_v3(y,filterby,ds);
        fill([nx;flipud(nx)],[ny-ndy;flipud(ny+ndy)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(ny)
        xticks([5,10,15,20,25,30,35,40,45]*60*4/ds)
        xticklabels([5,10,15,20,25,30,35,40,45])
        non=size(y,1);;
        title([strcat(all_psth_start(1).experiment,", mean activated cells: ",string(round(act/length(dff_total)*100)), " percent")
            strcat(all_psth_start(1).experiment,", mean inhibited cells: ",string(round(inhib/length(dff_total)*100)), " percent"); 
            strcat(all_psth_start(1).experiment,", mean nonresponsive cells: ",string(round(non/length(dff_total)*100)), " percent")])
        xlim([1 max(nx)])
        
      
        %% 
        %Compare activity during lick bout and after lick bout
       
        yrange=[-10 12];
        xrange=[-4 10];
        fl=fitlm(dff,dff_total) %p=0.262; r2=0.00743
        figure
        h=plot(fl);
        hx=get(h,'Xdata');
        hy=get(h,'Ydata');
        x=hx{3}';
        xs=hx{1}';
        ys=hy{1}';
        y=hy{2}';
        upper=hy{4}';
        lower=hy{3}';
        fill([x; flipud(x)], [lower; flipud(upper)], [.2 .9 .9],'linestyle','none')
        hold on
        plot(x,y)
        hold on
        scatter(xs,ys, 'k')
        ylim(yrange)
        xlim(xrange)
        [~, p]=ttest(dff,dff_total) %p=0.1140, paired t-test, welch's
        signrank(dff,dff_total) %p=0.2735, wilcoxon signed rank
        title(['activated during lick vs after access, p=', char(string([anova(fl).pValue(1)])), ', r2=', char(string(fl.Rsquared.Ordinary))])
        xlabel([all_psth_start(1).experiment,' activity during lick '])
        ylabel([all_psth_start(1).experiment, ' activity long-term'])
        
        A=[length(find(dff>.5)) length(find(dff_total>.5)) size(dff,1)];
        I=[length(find(dff>.5&dff_total>.5)) length(find(dff>.5)) length(find(dff_total>.5)) size(dff,1)];
        figure
        I(4)=I(4)*2 %stupid work around
        A(3)=A(3)*2 %stupid work around
        venn(A,I)
        title(['activated during lick: ', char(string(A(1))), '; vs after access: ', char(string(A(2))), '; both=', char(string(I(1))), '; all=', char(string(size(dff,1)))])

        
        %% making contour maps
        % making contour map
        
        for (mouse=1:size(all_psth_start,1))
            for dfs=1:2
                figure
                zmax=2;
                coloron=zmax;
                coloroff=.5;
                [cmap]=JG_cmap_v1_redder(zmax,coloron,coloroff);
                zmap=-zmax:(zmax*2)/255:zmax;
                coor=all_psth_start(mouse).coor;
                dff1=all_psth_start(mouse).dff;
                dff2=all_psth_start(mouse).dff_total;
                
                for c=1:size(coor)
                    if(dfs==1)
                        this_z=dff(c);
                    elseif(dfs==2)
                        this_z=dff_total(c);
                    end
                    if(this_z>zmax)
                        this_z=zmax;
                    elseif(this_z<-zmax)
                        this_z=-zmax;
                    end
                    p =coor{c}';
                    p=unique(p,'rows','stable');
                    p=p*0.8125; %adjustment to microns
                    
                    %removing the smallest area around large point deviations
                    p=[p; p(1,:)];
                    p_change=[sqrt(diff((p(:,1))).^2 + diff((p(:,2))).^2)];
                    three_std=std(diff(p_change))*3+mean(diff(p_change));
                    
                    if(size(find(abs(diff(p_change))>three_std),1)>0)
                        these=find(abs(diff(p_change))>three_std);
                        area1=polyarea(p(these(1):these(length(these)),1),p(these(1):these(length(these)),2));
                        p2=p;
                        p2(these(1):these(length(these)),:)=[];
                        area2=polyarea(p2(:,1),p2(:,2));
                        if(area2>area1)
                            p(these(1):(these(length(these))+1),:)=[];
                        else
                            p=p(these(1):these(length(these)),:);
                        end
                        p_change=[sqrt(diff((p(:,1))).^2 + diff((p(:,2))).^2)];
                    end
                    
                    p=flip(p);
                    %p(1,:)=[];
                    p=[p; p(1,:)];
                    
                    p_change=[sqrt(diff((p(:,1))).^2 + diff((p(:,2))).^2)];
                    three_std=std(diff(p_change))*3+mean(diff(p_change));
                    
                    if(size(find(abs(diff(p_change))>three_std),1)>0)
                        these=find(abs(diff(p_change))>three_std);
                        area1=polyarea(p(these(1):these(length(these)),1),p(these(1):these(length(these)),2));
                        p2=p;
                        p2(these(1):these(length(these)),:)=[];
                        area2=polyarea(p2(:,1),p2(:,2));
                        if(area2>area1)
                            p(these(1):(these(length(these))+1),:)=[];
                        else
                            p=p(these(1):these(length(these)),:);
                        end
                        p_change=[sqrt(diff((p(:,1))).^2 + diff((p(:,2))).^2)];
                    end
                    %
                    
                    n = size(p,1); n1 = n-1;
                    for    i=0:1:n1
                        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values
                    end
                    l=[];
                    UB=[];
                    for u=0:0.002:1
                        for d=1:n
                            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
                        end
                        l=cat(1,l,UB);                                      %catenation
                    end
                    P=l*p;
                    fill(P(:,1),P(:,2), [interp1(zmap,cmap(:,1),this_z) interp1(zmap,cmap(:,2),this_z) interp1(zmap,cmap(:,3),this_z)])
                    hold on
                end
                if(dfs==2)
                    dfstr='long-term effect';
                else
                    dfstr='lick effect';
                end
                title(['mouse ' num2str(mouse) ', ' dfstr, ', axes in microns'])
                
                maxs=max(cat(2,coor{:}),[],2)*0.8125;
                mins=min(cat(2,coor{:}),[],2)*0.8125;
                xlim([mins(1) maxs(1)])
                ylim([mins(2) maxs(2)])
                
                xlim([mean(xlim)-150 mean(xlim)+150]);ylim([mean(ylim)-150 mean(ylim)+150])
                
                circ=circle(mean(xlim), mean(ylim),270/2);
                %plot(circ(:,1),circ(:,2), 'color', [0 0 0])
                hold off
                
                colormap(cmap)
                c=colorbar('Ticks',[0 1],'TickLabels',[-zmax,zmax]);
                c.Label.String = 'z-score';
            end
        end
        
    
    
    if(crossreffile==1)
        %make graphs sorted by A start
        sorted_by= "A start"
        [A_combined_psth_start, A_combined_psth_end, B_combined_psth_start, B_combined_psth_end] = rerank_for_inscopix_analysis_JG(all_psth_start, all_psth_end, sorted_by);
        
        fl=fitlm(A_combined_psth_start.dff,B_combined_psth_start.dff) %p=0.0759; r2=0.259
        figure
        plot(fl)
        ylim([-4 4]) %change this
        [~, p]=ttest(A_combined_psth_start.dff,B_combined_psth_start.dff) %p=0.7324, paired t-test, welch's
        signrank(A_combined_psth_start.dff,B_combined_psth_start.dff) %p=0.4548, wilcoxon signed rank
        title(['lick effect for ', A_combined_psth_start.experiment(1,:), ' versus ', B_combined_psth_start.experiment(1,:)])
        xlabel(['lick effect for ', A_combined_psth_start.experiment(1,:)])
        ylabel(['lick effect for ', B_combined_psth_start.experiment(1,:)])
        
        
        fl=fitlm(A_combined_psth_start.dff_total,B_combined_psth_start.dff_total) %p=0.0222; r2=0.391
        figure
        plot(fl)
        ylim([-4 4]) %change this
        [~, p]=ttest(A_combined_psth_start.dff_total,B_combined_psth_start.dff_total) %p=0.4195, paired t-test, welch's
        signrank(A_combined_psth_start.dff_total,B_combined_psth_start.dff_total) %p=0.2439, wilcoxon signed rank
        title(['long-term effect for ', A_combined_psth_start.experiment(1,:), ' versus ', B_combined_psth_start.experiment(1,:)])
        xlabel(['long-term effect for ', A_combined_psth_start.experiment(1,:)])
        ylabel(['long-term effect for ', B_combined_psth_start.experiment(1,:)])
        
        figure()
        subplot(121)
        h=heatmap_d(A_combined_psth_start.avgpsth,[],[],[],'MaxColorValue',maxcv/4,'MinColorValue',mincv/4,'Colormap',cmap,'Colorbar',1);
        title(A_combined_psth_start.experiment(1,:))
        yticks([1 size(A_combined_psth_start.avgpsth,1)])
       yticklabels([1 size(A_combined_psth_start.avgpsth,1)])
        subplot(122)
        h=heatmap_d(B_combined_psth_start.avgpsth,[],[],[],'MaxColorValue',maxcv/4,'MinColorValue',mincv/4,'Colormap',cmap,'Colorbar',1);
        title(B_combined_psth_start.experiment(1,:))
       yticks([1 size(B_combined_psth_start.avgpsth,1)])
       yticklabels([1 size(B_combined_psth_start.avgpsth,1)])
        [ax,h3]=suplabel(['Bout start: sorted by bout start in ', A_combined_psth_start.experiment(1,:), ', ' zscoreyn_words, ', ', firstbout_words, ', ', baseline_words]  ,'t');
        
        figure()
        subplot(121)
        h=heatmap_d(A_combined_psth_end.avgpsth,[],[],[],'MaxColorValue',maxcv/4,'MinColorValue',mincv/4,'Colormap',cmap,'Colorbar',1);
        title(A_combined_psth_start.experiment(1,:))
        yticks([1 size(A_combined_psth_end.avgpsth,1)])
       yticklabels([1 size(A_combined_psth_end.avgpsth,1)])
        subplot(122)
        h=heatmap_d(B_combined_psth_end.avgpsth,[],[],[],'MaxColorValue',maxcv/4,'MinColorValue',mincv/4,'Colormap',cmap,'Colorbar',1);
        title(B_combined_psth_start.experiment(1,:))
       yticks([1 size(B_combined_psth_end.avgpsth,1)])
       yticklabels([1 size(B_combined_psth_end.avgpsth,1)])
        [ax,h3]=suplabel(['Bout end: sorted by bout start in ', A_combined_psth_start.experiment(1,:), ', ' zscoreyn_words, ', ', firstbout_words, ', ', baseline_words]  ,'t');
        
        sorted_by= "A total"
        [A_combined_psth_start, A_combined_psth_end, B_combined_psth_start, B_combined_psth_end] = rerank_for_inscopix_analysis_JG(all_psth_start, all_psth_end, sorted_by);
        
        figure()
        subplot(121)
        h=heatmap_d(A_combined_psth_end.zC,[],[],[],'MaxColorValue',maxcv/4,'MinColorValue',mincv/4,'Colormap',cmap);
        title(A_combined_psth_start.experiment(1,:))
        yticks([1 size(A_combined_psth_end.zC,1)])
       yticklabels([1 size(A_combined_psth_end.zC,1)])
        subplot(122)
        h=heatmap_d(B_combined_psth_end.zC,[],[],[],'MaxColorValue',maxcv/4,'MinColorValue',mincv/4,'Colormap',cmap);
        title(B_combined_psth_start.experiment(1,:))
        yticks([1 size(B_combined_psth_end.zC,1)])
       yticklabels([1 size(B_combined_psth_end.zC,1)])
        [ax,h3]=suplabel(['Whole trial: sorted by long term effect in ', A_combined_psth_start.experiment(1,:), ', ' zscoreyn_words, ', ', firstbout_words]  ,'t');
        
        
    end
end
end
clear all

try
    mkdir('graphs');
end


old_dir=cd;
FolderName = strcat(cd,'/graphs');   % Your destination folder
cd(FolderName)
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  set(0, 'CurrentFigure', FigHandle)
  FigName=get(get(gca, 'Title'),'string');
  if(iscell(FigName))
  FigName=FigName{1};
  end
  FigName=regexprep(FigName, ' +', '_');
  FigName=regexprep(FigName, ',+', '');
  savefig(strcat(FigName, '.fig'));
  saveas(gcf, FigName, 'epsc');
end
cd(old_dir)
