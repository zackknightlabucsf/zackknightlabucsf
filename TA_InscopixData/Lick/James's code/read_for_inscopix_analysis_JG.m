function [all_psth_start, all_psth_end,ttl_clean,cpsthrange, pretimes, posttimes, timestamps]=read_for_inscopix_analysis_JG( csvfile_name, rawyn, hw_sec, zscoreyn, boutthreshold_sec, firstbout, subtract_baseline)

csvfile=readtable(csvfile_name);

%make sure all these variables are actually used
all_psth_start=[];
all_psth_end=[];
avg_psth_start=[];
avg_psth_end=[];
cpsthrange=[];
ttl_clean=[];
pretimes=[];
posttimes=[];
timestamps=[];


for csv_row = 1:size(csvfile,1)
    
    %% read csvfile
    mouse_id=table2cell(csvfile(csv_row,1));mouse_id=mouse_id{1};
    matfile=table2cell(csvfile(csv_row,2));matfile=matfile{1};
    crossreg=table2cell(csvfile(csv_row,3));crossreg=crossreg{1};
    experiment=table2cell(csvfile(csv_row,4));experiment=experiment{1};
    matfile_fr=table2cell(csvfile(csv_row,5));matfile_fr=matfile_fr{1};
    vidfile=table2cell(csvfile(csv_row,6));vidfile=vidfile{1};
    vidfile_fr=table2cell(csvfile(csv_row,7));vidfile_fr=vidfile_fr{1};
    lmfile=table2cell(csvfile(csv_row,8));lmfile=lmfile{1};
    lmfile_channel=table2cell(csvfile(csv_row,9));lmfile_channel=lmfile_channel{1};
    timestamp=table2cell(csvfile(csv_row,10));timestamp=timestamp{1};
    pretime=table2cell(csvfile(csv_row,11));pretime=pretime{1};
    posttime=table2cell(csvfile(csv_row,12));posttime=posttime{1};
    
    sind=strfind(timestamp,':');
    minu=str2double(timestamp(1:sind-1));
    seco=str2double(timestamp(sind+1:end));
    timestamp=(minu*60+seco)*matfile_fr;
    timestamps=[timestamps, timestamp];
    
    sind=strfind(pretime,':');
    minu=str2double(pretime(1:sind-1));
    seco=str2double(pretime(sind+1:end));
    pretime=(minu*60+seco)*matfile_fr;
    pretimes=[pretimes, pretime];
    
    sind=strfind(posttime,':');
    minu=str2double(posttime(1:sind-1));
    seco=str2double(posttime(sind+1:end));
    posttime=(minu*60+seco)*matfile_fr;
    posttimes=[posttimes, posttime];
    
    hw=hw_sec*matfile_fr; %now hw is in frames, not secs
    boutthreshold=boutthreshold_sec*matfile_fr; %now boutthreshold is in frames, not secs
    
    %% load right C file (depends on crossreg)
    
    if(crossreg=='NA')
        load(matfile,'neuron')
        if rawyn
            C=neuron.C_raw;
        else
            C=neuron.C;
        end
    elseif(crossreg=='A')
        load(matfile,'mt_neurons')
        if rawyn
            C=mt_neurons.C_raw_A;
        else
            C=mt_neurons.C_A;
        end
    elseif(crossreg=='B')
        load(matfile,'mt_neurons')
        if rawyn
            C=mt_neurons.C_raw_B;
        else
            C=mt_neurons.C_B;
        end
    else
        ['Error: check row ', csv_row, ' in ', csvfile_name]
        return
    end
    
    %% create z-scored C, or zC
    zC=mat2cell(C,ones(size(C,1),1),size(C,2));
    clear C
    
    if zscoreyn==1 %zscore to prestim
        zC=cellfun(@(x) x/std(x(timestamp-pretime:timestamp)), zC, 'UniformOutput', false);
    else if zscoreyn==2 %zscore to whole trial
            zC=cellfun(@(x) x/std(x), zC, 'UniformOutput', false);
        else if zscoreyn==0 %don't zscore
%                 zC=cellfun(@(x) (x+100), zC, 'UniformOutput', false);
%                 this above line existed in the original YC file but I
%                 don't understand why
                zC=cellfun(@(x) 100*x/mean(x(timestamp-pretime:timestamp)), zC, 'UniformOutput', false);
                %returns dF/F in percent (i.e. 90m= 90%)
            end
            
        end
    end
    
    zC=cellfun(@(x) x-mean(x(timestamp-pretime:timestamp)), zC, 'UniformOutput', false); %subtract mean of pretime
    zC=cell2mat(zC);
    
    lastframe=size(zC,2);
    ttl=[];
    [ttl] = csv_reshape_for_inscopix_analysis_JG( lmfile,lmfile_channel, matfile_fr, lastframe);
    
    ttl(1:timestamp)=min(ttl); %ignore ttl signal from before timestamp
    
    [psth_start,psth_end, ttl_clean,bouts] = ttl2psth_for_inscopix_analysis_JG(zC, ttl,hw, boutthreshold, firstbout, crossreg, mouse_id, experiment, subtract_baseline);
    
    %put everything into a large file
    all_psth_start=[all_psth_start;psth_start];
    all_psth_end=[all_psth_end;psth_end];
    ttl_clean=ttl_clean;
    cpsthrange=[timestamp-pretime:timestamp+posttime];
end
end