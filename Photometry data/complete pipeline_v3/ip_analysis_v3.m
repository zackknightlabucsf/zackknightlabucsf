%making figures for non-licking photometry data

clear all; close all
load('all_photom_data.mat');

ds=4; %amount to further downsample by for graphing purposes

%for analysis, we look at effects over certain time windows after the start of a manipulation
%here, we look at an effect occuring between "window_start" to "window_end"
window_start=300; %beginning of time window analyzed (seconds after manipulation)
window_end=2100; %end of time window analyzed (seconds after manipulation)
%p-values come from two-sided permutation tests

ind=1:size(d.cat,2);
%ind = find(~contains(d.cat,'licks')); %removing licking photometry data

exps=struct;%contains repeats
exps.mouse=[];
exps.experiment=[];
exps.full_trace=[];
exps.z_trace=[];
exps.repeat=[];
exps.rate=[];
exps.manip=[];

for(i =ind)
    exps.mouse=[exps.mouse,string(d.mice{i})];
    exps.experiment=[exps.experiment,string(d.experiments{i})];
    exps.rate=[exps.rate,d.rate{i}];
    t=d.traces{i};
    m=d.manipulation{i};
    z=(t-mean(t(1:m)))/std(t(1:m));
    exps.manip=[exps.manip,m];
    exps.full_trace=[exps.full_trace;reshape(t,1,length(t))];
    exps.z_trace=[exps.z_trace;reshape(z,1,length(z))];
    
    if(sum(contains(exps.mouse,string(d.mice{i}))&contains(exps.experiment,string(d.experiments{i})))>1)
        exps.repeat=[exps.repeat,1];
    else
        exps.repeat=[exps.repeat,0];
    end
end

ip=struct; %repeats are averaged
ip.mouse=[];
ip.experiment=[];
ip.full_trace=[];
ip.z_trace=[];
ip.rate=[];
ip.manip=[];
for(m=unique(exps.mouse))
    for(e=unique(exps.experiment))
        ind = find(contains(exps.experiment,e) & contains(exps.mouse,m));
        if(length(ind)>0)
        ind2=ind(1);
        m=exps.manip(ind2);
        ip.rate=[ip.rate,exps.rate(ind2)];
        ip.mouse=[ip.mouse, m];
        ip.experiment=[ip.experiment,e];
        ip.full_trace=[ip.full_trace; mean(exps.full_trace(ind,:),1)];
        ip.z_trace=[ip.z_trace; mean(exps.z_trace(ind,:),1)];
        ip.manip=[ip.manip,m];
        end
    end
end


exps=unique(ip.experiment);
dffs={}; %cell array of effects (in df/f)
zs={}; %cell array of effects (in z-score)
close all
for(e=1:size(exps,2))
    ind=find(contains(ip.experiment, exps(e)));
    if(length(ind)>0)
        y=ip.full_trace(ind,:);
        [x,ny1,ndy1]=JG_downsample_v3(y,ds,ds);
        x=x/60*ds/ip.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        title(strcat(exps(e),": dff trace"))
        ylabel('df/f')
        xlabel('min')
        hold on
        plot([(ip.manip(ind(1))-1)/ip.rate(ind(1))/60 (ip.manip(ind(1))-1)/ip.rate(ind(1))/60], ylim)
        
        figure
        ws=window_start*ip.rate(ind(1))+ip.manip(ind(1));
        we=window_end*ip.rate(ind(1))+ip.manip(ind(1));
        dff=mean(y(:,ws:we),2);
        bar(mean(dff))
        hold on
        er=errorbar(1,mean(dff),std(dff)/sqrt(length(dff)));er.Color = [0 0 0];er.LineStyle = 'none';er.LineWidth=1;er.CapSize=50
        plot(dff*0+1,dff,'o','MarkerFaceColor','black', 'Color', 'black')
        ylabel('dff')
        title([strcat(exps(e),": dff barplot"); strcat("p=", p_to_star(permutationTest(dff,dff*0,10000,'exact',1)),"; n=", string(length(dff)))])
        dffs{size(dffs,2)+1}=dff;
        
        y=ip.z_trace(ind,:);
        [x,ny1,ndy1]=JG_downsample_v3(y,ds,ds);
        x=x/60*ds/ip.rate(ind(1));
        figure
        fill([x;flipud(x)],[ny1-ndy1;flipud(ny1+ndy1)],[.4 .4 .4],'linestyle','none');
        hold on
        plot(x,ny1,'k')
        title(strcat(exps(e),": zscored trace"))
        ylabel('z-score')
        xlabel('min')
        hold on
        plot([(ip.manip(ind(1))-1)/ip.rate(ind(1))/60 (ip.manip(ind(1))-1)/ip.rate(ind(1))/60], ylim)
        
        figure
        ws=window_start*ip.rate(ind(1))+ip.manip(ind(1));
        we=window_end*ip.rate(ind(1))+ip.manip(ind(1));
        z=mean(y(:,ws:we),2);
        bar(mean(z))
        hold on
        er=errorbar(1,mean(z),std(z)/sqrt(length(z)));er.Color = [0 0 0];er.LineStyle = 'none';er.LineWidth=1;er.CapSize=50
        plot(z*0+1,z,'o','MarkerFaceColor','black', 'Color', 'black')
        ylabel('z-score')
        title([strcat(exps(e),": zscored barplot"); strcat("p=", p_to_star(permutationTest(z,z*0,10000,'exact',1)),"; n=", string(length(z)))])
        zs{size(zs,2)+1}=z;
    end
end


try
    mkdir('graphs')
end
home=cd;
cd('graphs');

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  set(0, 'CurrentFigure', FigHandle)
%   set(gcf,'renderer','Painters')
  FigName=get(get(gca, 'Title'),'string');
  if(iscell(FigName))FigName=FigName{1};end
  FigName=regexprep(FigName, ' +', '_');
  FigName=regexprep(FigName, ',+', '');
  FigName=regexprep(FigName, ':+', '');
  savefig(strcat(FigName, '.fig'));
  saveas(gcf, FigName, 'png');
  saveas(gcf, FigName, 'eps');
end
cd(home)