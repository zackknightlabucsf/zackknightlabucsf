function [psth_start,psth_end, ttl_clean,bouts] = ttl2psth_for_inscopix_analysis_JG(zC, ttl,hw, boutthreshold, firstbout,crossreg, mouse_id, experiment, subtract_baseline)

%UNTITLED4 C:NxT matrix of activity trace; ttl: 1xT array of TTL; hw: (unit: array postiion) size of halfwindow width defining bouts
%   C: psth data of whole trial; ttl: timestamp of events; hw: half window
% firstbout: 0 - don't include; 1 - include; 2 - firstbout only
%   of psth; boutthreshold: number of events needed to be considered as a bout-
%   extra ones will be excluded from the ttl;
psth_start=struct('psth',{[]},'avgpsth',[0],'ttl_bout',ttl*0,'time_ttl_bout',[0],'dff',[],'boutsize',[], 'crossreg',[], 'mouseid',[], 'zC', [], 'experiment', []);
psth_end=struct('psth',{[]},'avgpsth',[0],'ttl_bout',ttl*0,'time_ttl_bout',[0],'dff',[],'boutsize',[], 'crossreg', [], 'mouseid',[], 'zC', [], 'experiment', []);

ttlthr=(max(ttl)+min(ttl))/2;
ttl_time=find(ttl>ttlthr);
bouts={ttl_time(1)};

%% attribute TTL into individual bouts based on half window
%half window is half of the graphing window size
for i=2:length(ttl_time)
    if ttl_time(i)-ttl_time(i-1)>hw*2 %if ttl on is further than hw*2 from last ttl
        bouts{length(bouts)+1}=[ttl_time(i)]; %add this ttl
        %         keyboard();
    else
        if length(bouts)>=1
            bouts{length(bouts)}=[bouts{length(bouts)} ttl_time(i)]; 
        end
        %this ends with cell of frames for each bout
    end
end
%delete bouts under threshold

blength=cellfun(@(x) x(end)-x(1),bouts); %find how long bouts last
eind=find(blength<boutthreshold);
bouts(eind)=[]; %remove those under threshold

if firstbout==1 %include all
else if firstbout==2 %first bout only
       bouts=bouts(1); 
    else %remove first bout
        bouts(1)=[];
    end
end


%% find bout start/end
ttl_time=[]; %all counted bouts
n=1; %all removed bouts
eind=[];
for i=1:length(bouts);
    if bouts{i}(end)<size(zC,2)-hw %if bout ends before end of zC
        psth_start.time_ttl_bout(n)=bouts{i}(1);
        psth_end.time_ttl_bout(n)=bouts{i}(end);
        %         keyboard();
        %
        ttl_time=[ttl_time bouts{i}];
    else
        eind=[eind i];
    end
    n=n+1;
    %
end
bouts(eind)=[];

blength=cellfun(@(x) x(end)-x(1),bouts); %find how long bouts last

psth_start.boutsize=blength;
psth_end.boutsize=blength;

%% get ttl trace that exclude the licks in bouts below threshold
ttl_temp=ttl*0;
% keyboard();
ttl_temp(psth_start.time_ttl_bout)=1;
psth_start.ttl_bout=ttl_temp;ttl_temp=ttl*0;

ttl_temp=ttl*0;
ttl_temp(psth_end.time_ttl_bout)=1;
psth_end.ttl_bout=ttl_temp;

ttl_temp=ttl*0;
ttl_temp(ttl_time)=1;
ttl_clean=ttl_temp;

%% analyze psth
for i=1:size(zC,1)
    %% psth at the start of bout
    mpsth=[];
    tempC=zC(i,:);
    mean_start=[];
    for i2=1:length(psth_start.time_ttl_bout)
        %         try
        psth=tempC(psth_start.time_ttl_bout(i2)-hw:psth_start.time_ttl_bout(i2)+hw);
        %         catch
        %             keyboard();
        %         end
        mean_start=[mean_start, mean(psth(1:hw))];
        if(subtract_baseline~=0)
            psth=psth-mean(psth(1:hw)); %NOT re-zscored, but average of before is subtracted
        end
        mpsth=[mpsth;psth];
    end
    psth_start.psth{i}=mpsth; %every psth
    psth_start.dff=[psth_start.dff; mean(mean(mpsth(:,1:hw)')-mean(mpsth(:,hw+1:2*hw)'))]; %diff between before and after
    %     keyboard();
    if isempty(mpsth) %if no ttl signal
        psth_start.avgpsth(i,1:hw*2+1)=zeros(1,hw*2+1);
    else
        psth_start.avgpsth(i,1:hw*2+1)=mean(mpsth,1); %averaged psth
    end
    %% psth end
    mpsth=[];
    tempC=zC(i,:);
    for i2=1:length(psth_end.time_ttl_bout)
%         try
            psth=tempC(psth_end.time_ttl_bout(i2)-hw:psth_end.time_ttl_bout(i2)+hw);
            psths=tempC(psth_start.time_ttl_bout(i2)-hw:psth_start.time_ttl_bout(i2)+hw);
%         catch
%             keyboard();
%         end
        if(subtract_baseline==2)
            psth=psth-mean(psth(1:hw)); %NOT re-zscored, but average of before is subtracted
        elseif(subtract_baseline==1)
            psth=psth-mean_start(i2);
        end
        mpsth=[mpsth;psth];
    end
    psth_end.psth{i}=mpsth;%every psth
    psth_end.dff=[psth_end.dff; mean(mean(mpsth(:,1:hw)')-mean(mpsth(:,hw+1:2*hw)'))]; %diff between before and after
    
    if isempty(mpsth)%if no ttl signal
        psth_end.avgpsth(i,1:hw*2+1)=zeros(1,hw*2+1);
    else
        psth_end.avgpsth(i,1:hw*2+1)=mean(mpsth,1); %averaged psth
    end
end
psth_start.crossreg=[psth_start.crossreg; crossreg];
psth_end.crossreg=[psth_end.crossreg; crossreg];
psth_start.mouseid=[psth_start.mouseid; mouse_id];
psth_end.mouseid=[psth_end.mouseid; mouse_id];

psth_start.zC=[psth_start.zC; zC];
psth_end.zC=[psth_end.zC; zC];

psth_start.experiment=[psth_start.experiment; experiment];
psth_end.experiment=[psth_end.experiment; experiment];
end