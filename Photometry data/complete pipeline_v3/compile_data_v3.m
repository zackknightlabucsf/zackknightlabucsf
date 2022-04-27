%this script is used to organize all photometry data into one file

clear all;close all;

f=readtable('metadata.xlsx'); %f=files
data_dir='C:\Users\ZK\Desktop\PRLH Analysis\FD Empty Bottle - First 2 min omitted';

%Note: this script assumes that Epo1+Epo2 are  the TTL channels for Z+K
%If they are not, change the code on the noted lines

%%%%%%%%%
home=cd;
cd(data_dir)
traces={};
rawtraces={};
mice={};
experiments={};
cat={};
notes={};
problems={};
ttl={};
access_length={};
access_time={};
rate={};
manipulation={};

for nf=1:size(f,1)
    try load(string(f.files{nf}));
        %
        manip_time=f.manip_time(nf);if(iscell(manip_time));manip_time=str2double(manip_time);end
        after_manip=f.after_manip(nf);if(iscell(after_manip));after_manip=str2double(after_manip);end
        baseline_length=f.baseline_length(nf);if(iscell(baseline_length));baseline_length=str2double(baseline_length);end
        downsample_to=f.downsample_to(nf);if(iscell(downsample_to));downsample_to=str2double(downsample_to);end
        %
        if(length(f.Z{nf})>1)
            [~,G_Z]=JG_downsample_v3(d.streams.G__Z.data,downsample_to,round(d.streams.G__Z.fs/downsample_to));
            [~,U_Z]=JG_downsample_v3(d.streams.U__Z.data,downsample_to,round(d.streams.U__Z.fs/downsample_to));
            index_norm=((manip_time-baseline_length)*downsample_to+1):(manip_time*downsample_to);
            para_fit=polyfit(U_Z(index_norm),G_Z(index_norm),1);
            Gfit_Z = para_fit(1).*U_Z + para_fit(2);
            y=G_Z./Gfit_Z;
            y(1)=[];
            y=(y-1)*100;
            
            new_beginning=(manip_time-baseline_length)*downsample_to+1;
            new_end=(manip_time+after_manip)*downsample_to;
            y=y(new_beginning:new_end); %removing stuff before baseline and after after_manip
            G_Z=G_Z(new_beginning:new_end); %removing stuff before baseline and after after_manip
            traces{size(traces,1)+1,1}=y;
            rawtraces{size(rawtraces,1)+1,1}=G_Z;
            rate{size(rate,2)+1}=downsample_to;
            manipulation{size(manipulation,2)+1}=baseline_length*downsample_to+1;
            
            mice{size(mice,2)+1}=f.Z{nf};
            %notes{size(notes,2)+1}=f.Notes{nf};
            notes{size(notes,2)+1}=NaN;
            experiments{size(experiments,2)+1}=lower(f.exp{nf});
            if(contains(lower(f.exp{nf}), 'ip'))
                cat{size(cat,2)+1}='ip';
                ttl{size(ttl,1)+1,1}=[];
            elseif(contains(lower(f.exp{nf}), 'licks'))
                cat{size(cat,2)+1}='licks';
                try
                    ttl{size(ttl,1)+1,1}=d.epocs.Ep2_.onset-(new_beginning-1)/downsample_to; %change to ttl channel for Z
                catch
                    ttl{size(ttl,1)+1,1}=[];
                end
            elseif(contains(lower(f.exp{nf}), 'ig'))
                cat{size(cat,2)+1}='ig';
                ttl{size(ttl,1)+1,1}=[];
            else
                cat{size(cat,2)+1}='misc';
                ttl{size(ttl,1)+1,1}=[];
            end
            
        if(iscell(f.access_length(nf)))
            access_length{size(access_length,2)+1}=0;
        else
            access_length{size(access_length,2)+1}=f.access_length(nf);
        end
        
        if(iscell(f.access_time(nf)))
            access_time{size(access_time,2)+1}=0;
        else
            access_time{size(access_time,2)+1}=f.access_time(nf)-((new_beginning-1)/downsample_to);
        end
        end
        
        if(length(f.K{nf})>1)
            [~,G_K]=JG_downsample_v3(d.streams.G__K.data,downsample_to,round(d.streams.G__K.fs/downsample_to));
            [~,U_K]=JG_downsample_v3(d.streams.U__K.data,downsample_to,round(d.streams.U__K.fs/downsample_to));
            index_norm=((manip_time-baseline_length)*downsample_to+1):(manip_time*downsample_to);
            para_fit=polyfit(U_K(index_norm),G_K(index_norm),1);
            Gfit_K = para_fit(1).*U_K + para_fit(2);
            y=G_K./Gfit_K;
            y(1)=[];
            y=(y-1)*100;
            new_beginning=(manip_time-baseline_length)*downsample_to+1;
            new_end=(manip_time+after_manip)*downsample_to;
            y=y(new_beginning:new_end); %removing stuff before baseline and after after_manip
            G_K=G_K(new_beginning:new_end); %removing stuff before baseline and after after_manip
            traces{size(traces,1)+1,1}=y;
            rawtraces{size(rawtraces,1)+1,1}=G_K;
            rate{size(rate,2)+1}=downsample_to;
            manipulation{size(manipulation,2)+1}=baseline_length*downsample_to+1;
            
            
            mice{size(mice,2)+1}=f.K{nf};
            %notes{size(notes,2)+1}=f.Notes{nf};
            notes{size(notes,2)+1}=NaN;
            experiments{size(experiments,2)+1}=lower(f.exp{nf});
            if(contains(lower(f.exp{nf}), 'ip'))
                cat{size(cat,2)+1}='ip';
                ttl{size(ttl,1)+1,1}=[];
            elseif(contains(lower(f.exp{nf}), 'licks'))
                cat{size(cat,2)+1}='licks';
                try
                    ttl{size(ttl,1)+1,1}=d.epocs.Ep4_.onset-(new_beginning-1)/downsample_to; %change to ttl channel for K
                catch
                    ttl{size(ttl,1)+1,1}=[];
                end
            elseif(contains(lower(f.exp{nf}), 'ig'))
                cat{size(cat,2)+1}='ig';
                ttl{size(ttl,1)+1,1}=[];
            else
                cat{size(cat,2)+1}='misc';
                ttl{size(ttl,1)+1,1}=[];
            end
            
        if(iscell(f.access_length(nf)))
            access_length{size(access_length,2)+1}=0;
        else
            access_length{size(access_length,2)+1}=f.access_length(nf);
        end
        
        if(iscell(f.access_time(nf)))
            access_time{size(access_time,2)+1}=0;
        else
            access_time{size(access_time,2)+1}=f.access_time(nf)-((new_beginning-1)/downsample_to);
        end
        end
    catch
        fprintf('\n')
        fprintf('problem with')
        problems{size(problems,1)+1,1}=string(f.files(nf));
    end
    fprintf('\n')
    fprintf(string(nf))
end
cd(home)

d=struct;
d.traces=traces;
d.rawtraces=rawtraces;
d.mice=mice;
d.experiments=experiments;
d.cat=cat;
d.notes=notes;
d.ttl=ttl;
d.access_length=access_length;
d.access_time=access_time;
d.rate=rate;
d.manipulation=manipulation;


save('all_photom_data.mat','d')