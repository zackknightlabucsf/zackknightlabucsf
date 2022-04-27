function [ttl] = csv_reshape_for_inscopix_analysis_JG( lmfile,lmfile_channel, matfile_fr, lastframe)

gpiofile=readtable([lmfile]);

these=ismember(gpiofile.ChannelName,lmfile_channel);
gpiofile=gpiofile(these,:);

signal=gpiofile.Value;
ttl=[];

for i=1:lastframe;
    first_time=i-1;
    last_time=i+1;
    these=gpiofile.Value(gpiofile.Time_s_<last_time/matfile_fr & gpiofile.Time_s_>first_time/matfile_fr);
    
    if(length(these)==0)
        ttl(i)=min(gpiofile.Value);
    else
        ttl(i)=mean(these);
    end
end

clear these
clear first_time
clear last_time
clear signal
clear gpiofile
clear signal

end
