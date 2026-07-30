close all
clear all
fclose all
filename = 'FILENAME'
list={
    % 'CSV' %Lists individual recordings to be combined (same stimulus)
    };

        % Initialize output matrices.
neurons_raw = [];       % Baseline-subtracted only (raw)
neurons_zscored = [];   % Baseline-subtracted then z-scored

for i=1:length(list)
    txtfile=strcat(list{i},'.txt'); %combines all the neurons/columns from each mouse in the file listed in the text document
    
    %% initiate txtfile
% txtfile='.txt';
fid=fopen(txtfile);
t=fscanf(fid,'%s',[1]);
expname=t;
t=fscanf(fid,'%s',[1]);

% Read first token from text file
while t
    %% read info
    inmatfile=t
    fprintf(inmatfile)
    t=fscanf(fid,'%s',[1]);
    fr=str2double(t);
    
    t=fscanf(fid,'%s',[1]);
    timestamp=t;
    sind=strfind(timestamp,':');
    minu=str2double(timestamp(1:sind-1));
    seco=str2double(timestamp(sind+1:end));
    timestamp=(minu*60+seco)*fr;
    
    t=fscanf(fid,'%s',[1]);
    pretime=str2double(t)*fr;
    
    t=fscanf(fid,'%s',[1]);
    posttime=str2double(t)*fr;
    
    %% process data
    traces = csvread(inmatfile,2,1);
    C = traces';
    
    % Fixed window lengths (in samples)
    pre_samp  = 600 * fr;   % always 2400 samples of baseline
    post_samp = 1800 * fr;  % always 7200 samples post

    % Extract fixed window: timestamp is the infusion onset frame
    seg_start = timestamp - pre_samp + 1;
    seg_end   = timestamp + post_samp;

    if seg_end <= size(C, 2)
        psth = C(:, seg_start : seg_end);
    else
        psth = C(:, seg_start : end);
        padSize = seg_end - size(C, 2);
        psth(:, end+1 : end+padSize) = NaN;
    end

    % Z-score using the fixed 600s baseline (first pre_samp frames of psth)
    raw_segment    = zeros(size(psth));
    zscore_segment = zeros(size(psth));

    for j = 1:size(psth, 1)
        baseline   = mean(psth(j, 1:pre_samp));
        std_val    = std(psth(j,  1:pre_samp));
        raw_segment(j, :)    = psth(j, :) - baseline;
        zscore_segment(j, :) = raw_segment(j, :) / std_val;
    end

    % Append the processed data for this recording to the overall arrays.
    neurons_raw = [neurons_raw; raw_segment];
    neurons_zscored = [neurons_zscored; zscore_segment];
    
    % Read the next file token.
    t = fscanf(fid, '%s', 1);
end

fclose(fid);
end

% Save both output arrays to a MAT-file for future use.
save(filename, 'neurons_raw', 'neurons_zscored');
