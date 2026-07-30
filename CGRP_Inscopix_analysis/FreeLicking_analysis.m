%% BCJ Inscopix Analysis to free licking for solutions
% Adapted from TA code Licks_PSTH_cnmfecells_activatedcellsonly
%01/31/25

clear; close all;

%% Input Parameters (Explicitly in SECONDS)
neural_sampling_rate = 4; % Neural data at 4 Hz
min_ili = 0.1; % Minimum inter-lick interval (50 ms, defines individual licks)
inter_bout_interval = 10; % **5 second IBI to define separate bouts**
min_licks_per_bout = 3; % Minimum number of licks per bout
min_bout_duration = 10; %need this at least at 1 sec otherwise loses resolution to match to a 4hz recording

hw=10; % baseline for rebaselining to beginning of lick bout (second)
hw2 = 20; %Time after first lick bout to be used to determine if the neuron is lick-activated (half-window for graphing PSTH around licking, in seconds)

prestim = 600;
poststim = 1800;
access_end = prestim + poststim; %at what sec does access end? This is to make sure that the last bout is a hw before it

zs_threshold = 3;
bout_threshold = 0.6; % 60% threshold

filename = 'Ensure'
graph_data = 0; %put 1 if want to graph things
exclude_first = 0; %put 1 if exclude exclude_start seconds of session
exclude_start = prestim + 300;

txtfiles = { 
     'TXT_FILE.txt' %See "Sample_TXT_File"

};

infiles = { %csv files
     'CSV_of_licking.csv' %Output GPIO from inscopix IDPS software


 };

%%
num_files = length(txtfiles);
if length(infiles) ~= num_files
   error('txtfiles and infiles must have the same length!');
end

%Initialize the "results" array of structions ad global aggregator
% We'll store the “per-file” outputs in 'results' (array of structs).
results = struct(); %will hold per-file details
allNeurons = struct('act', [], 'inhib', [], 'none', [], 'act_neg', [], 'all', []);
categories = {'act','inhib','none'};


% We also want a structure that collects "raw_mean", "licknorm", etc.
allCatRes = struct();
results(num_files,1) = struct();
allFilesPSTH       = [];   % will be [total_bouts × window_length]
allFilesPSTH_end   = [];
allLicks_Act    = [];
allZ_Act        = [];
allLicks_Inhib  = [];
allZ_Inhib      = [];

for c = 1:length(categories)
    cat_name = categories{c};
    % For each category, store arrays that are [#neurons x #windows].
    allCatRes.(cat_name).raw_mean  = [];   % will vstack neurons from each file
    allCatRes.(cat_name).licknorm  = [];   % likewise
    allCatRes.(cat_name).R_raw  = [];   % likewise
    allCatRes.(cat_name).P_raw  = [];   % likewise
   
%     allCatRes.(cat_name).bout_avgZS   = [];  % (#neurons_cat x 1), e.g. "mean_bout_average_zs"
%     allCatRes.(cat_name).zscore_lickM = [];  % (#neurons_cat x 1), e.g. "zscore_per_lick_mean"
    allCatRes.(cat_name).per_neuron_boutmean_z    = [];
    allCatRes.(cat_name).per_neuron_boutmean_auc  = [];
    allCatRes.(cat_name).per_neuron_boutnorm_z    = [];
    allCatRes.(cat_name).mean_bout_subwin     = [];
    allCatRes.(cat_name).per_bout_z           = [];
    allCatRes.(cat_name).per_bout_auc         = [];
    allCatRes.(cat_name).per_bout_norm_z      = [];
    
    allCatRes.(cat_name).avg_nonBout_all = [];
    allCatRes.(cat_name).avg_nonBout_windows = [];
    %For PSTH metric (from all bouts)
    allCatRes.(cat_name).psth_perBout     = [];   % will be [total_bouts_cat × window]
    allCatRes.(cat_name).psth_end_perBout = [];
    allCatRes.(cat_name).PSTH_act = []; %each neuron's responses to all bouts will be averaged
    allCatRes.(cat_name).PSTH_means_smoothed = [];
    allCatRes.(cat_name).PSTH_end_neurons = [];
    allCatRes.(cat_name).PSTH_mean_cat = [];
    allCatRes.(cat_name).PSTH_peak = [];
     allCatRes.(cat_name).PSTH_peak_smoothed = [];
    allCatRes.(cat_name).PSTH_auc = [];
    allCatRes.(cat_name).PSTH_mean_peak = [];
    allCatRes.(cat_name).PSTH_mean_auc = [];
    allCatRes.(cat_name).PSTH_time_peak = [];
    allCatRes.(cat_name).PSTH_time_peak_smoothed = [];
    allCatRes.(cat_name).PSTH_end_time_peak_smoothed = [];
    
    %For first bout (separate)
    allCatRes.(cat_name).PSTH_first       = [];
    allCatRes.(cat_name).first_peak       = [];
    allCatRes.(cat_name).first_auc        = [];
    allCatRes.(cat_name).first_time_peak  = [];
    
    % For subsequent bouts (all except first)
    allCatRes.(cat_name).PSTH_subsequent_avg = [];
    allCatRes.(cat_name).subsequent_peak     = [];
    allCatRes.(cat_name).subsequent_auc      = [];
    allCatRes.(cat_name).subsequent_time_peak= [];%
    
    %for calculating bouts vs. licks
    allCatRes.(cat_name).collapsed_bout_means = [];
    allCatRes.(cat_name).collapsed_bout_licks = [];
    allCatRes.(cat_name).collapsed_bout_auc   = [];
end

% Initialize "all" category to store metrics for all neurons (not subdivided)
allCatRes.all.raw_mean  = [];
allCatRes.all.licknorm  = [];
allCatRes.all.R_raw     = [];
allCatRes.all.P_raw     = [];
% allCatRes.all.bout_avgZS = [];
% allCatRes.all.zscore_lickM = [];

allCatRes.all.per_neuron_boutmean_z    = [];
allCatRes.all.per_neuron_boutmean_auc  = [];
allCatRes.all.per_neuron_boutnorm_z    = [];
allCatRes.all.mean_bout_subwin     = [];
allCatRes.all.per_bout_z           = [];
allCatRes.all.per_bout_auc         = [];
allCatRes.all.per_bout_norm_z      = [];
allCatRes.all.avg_nonBout_all = [];
allCatRes.all.avg_nonBout_windows = [];
allCatRes.all.PSTH_act = [];
allCatRes.all.PSTH_means_smoothed = [];
allCatRes.all.PSTH_end_neurons = [];
allCatRes.all.PSTH_mean_cat = [];
allCatRes.all.PSTH_peak = [];
allCatRes.all.PSTH_auc = [];
allCatRes.all.PSTH_mean_peak = [];
allCatRes.all.PSTH_mean_auc = [];
allCatRes.all.PSTH_time_peak = [];
allCatRes.all.PSTH_first = [];
allCatRes.all.first_peak = [];
allCatRes.all.first_auc = [];
allCatRes.all.first_time_peak = [];
allCatRes.all.PSTH_subsequent_avg = [];
allCatRes.all.subsequent_peak = [];
allCatRes.all.subsequent_auc = [];
allCatRes.all.subsequent_time_peak = [];
allCatRes.all.collapsed_bout_means = [];
allCatRes.all.collapsed_bout_licks = [];
allCatRes.all.collapsed_bout_auc   = [];

%% %%%%%%%%%% FILE LOOP %%%%%%%%%%%%%%
%%Pull variables and neuron traces from txt files

for idx = 1:num_files
    thisTxt = txtfiles{idx};
    thisCSV = infiles{idx};
    fprintf('Processing File %d: %s, %s\n', idx, thisTxt, thisCSV);
    
    %% Read metadata from txt file
    % 1) Read the text file
    fileID = fopen(thisTxt, 'r');
    data   = textscan(fileID, '%s %d %s %d %d', 'Delimiter', '\t');
    fclose(fileID);

    % Extract stim_time from the txt file
    timeStr = data{3}{2};  % Extract the time string from the 3rd column, 2nd line, which contains time like '12:30'
    timeParts = textscan(timeStr, '%d:%d'); % Convert time string to minutes and seconds
    minutes = timeParts{1};
    seconds = timeParts{2};
    stim_time = minutes * 60 + seconds; % Calculate total seconds
    matfile_fr = double(data{2}(2)); % Extract imaging frequency

    %% Extract neural traces from the csv file
    raw_traces = [];
    for row_idx = 2:length(data{1})
        %Extract csv filename from the first column
        csvFileName = strtrim(data{1}{row_idx}); % Remove any whitespace
        if exist(csvFileName, 'file') ~= 2
            error('The CSV file %s does not exist.', csvFileName);
        end
    
        opts = detectImportOptions(csvFileName);
        opts.DataLines = [3, Inf]; % Start from the third row to the end
        variableNames = opts.VariableNames;
        opts.SelectedVariableNames = variableNames(2:end); % Select all columns except the first

        % Read the data into a table
        tracesTable = readtable(csvFileName, opts);

        % Convert the table to a numeric matrix
        C = table2array(tracesTable);
        C = C';
        raw_traces = [raw_traces; C];
    end

    %% Adjust the length of traces based on stim_time and access_end and create a z-scored version
    % Preallocate a matrix to store segments (assuming all traces are the same length)
    num_traces = size(raw_traces,1);
    start_idx = (stim_time-prestim) *matfile_fr;
    end_idx = (stim_time + poststim) *matfile_fr-1;

    % % Extract segments
    traces = raw_traces(:, start_idx:end_idx);
    % 
    % %Z-score traces based on prestim period
    prestim_length = prestim *matfile_fr;

    prestim_mean    = mean(traces(:, 1:prestim_length), 2);
    prestim_std     = std(traces(:,  1:prestim_length), 0, 2);

    % Z-score the traces
    trace_zs = (traces - prestim_mean) ./ prestim_std;
    [num_neurons, num_frames] = size(trace_zs);

    % The offset in absolute session time for frame #1 in 'trace_zs':
    time_offset = stim_time - prestim;

    %% 2) ============Read the lick CSV => "filtered_licks", "bouts" ================
    lickData = readtable(thisCSV);
    gpio1data = lickData(strcmp(lickData.ChannelName, 'GPIO-1'), :); %extracts just the gpio-1 channel into a table format, includes timestamp and ttl value

    % Extract Lick Events Using Rising Edges
    ttl_values = gpio1data.Value;
    time_values = gpio1data.Time_s_; % Ensure time is in seconds

    % Compute dynamic TTL threshold
    ttl_thr = (max(ttl_values) + min(ttl_values)) / 2;

    % Detect rising edges (first transition from low to high)
    rising_edges = find(diff(ttl_values > ttl_thr) == 1) + 1;
    lick_event_times = time_values(rising_edges);

    % Exclude any lick time > (stim_time + poststim)
    end_time = stim_time + poststim;
    lick_event_times = lick_event_times(lick_event_times <= end_time);
    lick_event_times = lick_event_times(lick_event_times >= stim_time);

    % Filter Out Licks with a Physiologically Impossible ILI (in SECONDS)
    filtered_licks = lick_event_times(1); % Initialize with first lick
    for i = 2:length(lick_event_times)
        if (lick_event_times(i) - filtered_licks(end)) > min_ili
            filtered_licks = [filtered_licks; lick_event_times(i)]; % Add valid lick
        end
    end

    %% =============Detect Lick Bouts Using Inter-lick interval criteria========================
    bouts = {}; % Initialize bout storage
    current_bout = filtered_licks(1); % Start first bout

    for i = 2:length(filtered_licks)
        if (filtered_licks(i) - filtered_licks(i-1)) > inter_bout_interval
            % We have a gap => close out the current_bout if it meets criteria
            if length(current_bout) >= min_licks_per_bout
                bout_duration = current_bout(end) - current_bout(1);
                if bout_duration >= min_bout_duration
                    bouts{end+1} = current_bout;  % Save completed bout
                end
            end
            % Start a new bout
            current_bout = filtered_licks(i);
        else
            % Still part of the same bout
            current_bout = [current_bout; filtered_licks(i)];
        end
    end

    % Add the last bout if it meets the lick count requirement
    % Check the very last bout
    if length(current_bout) >= min_licks_per_bout
        bout_duration = current_bout(end) - current_bout(1);
        if bout_duration >= min_bout_duration
            bouts{end+1} = current_bout;
        end
    end

    %% Compute Licks Per Bout
    licks_per_bout = cellfun(@length, bouts); % Count licks per bout
    bout_durations = cellfun(@(x) x(end) - x(1), bouts); % Compute bout durations in seconds
    avg_bout_duration = mean(bout_durations);
    total_bouts = length(bouts);
    total_filtered_licks = length(filtered_licks);
    
    %% Output Results
    disp(['Total Licks Detected: ', num2str(total_filtered_licks)]);
    disp(['Total Lick Bouts: ', num2str(total_bouts)]);
    disp(['Average Bout Duration (s): ', num2str(avg_bout_duration)]);

    %% Plot Lick Events and Bouts
%     figure;
%     histogram(filtered_licks, 'BinWidth', 1/neural_sampling_rate, 'EdgeColor', 'k', 'FaceAlpha', 0.7);
%     xlabel('Time (s)');
%     ylabel('Lick Count');
%     title('Binned Lick Events (Aligned to 4 Hz Neural Data)');
% 
%     % Plot number of licks per bout
%     figure;
%     bar(1:total_bouts, licks_per_bout);
%     xlabel('Bout Number');
%     ylabel('Licks per Bout');
%     title('Licks Per Bout');
% 
%     % Plot bout durations
%     figure;
%     bar(1:total_bouts, bout_durations);
%     xlabel('Bout Number');
%     ylabel('Bout Duration (s)');
%     title('Bout Durations');

%% Use Exact Bout Start Times for PSTH Calculations
% Extract bout start times in seconds (preserving full resolution)
bout_start_times = cellfun(@(x) x(1), bouts);

if exclude_first == 1 %keep only bouts AFTER "exclude_start"
    bout_indices = (bout_start_times > exclude_start);
else
    %keep all bouts that occur AFTER prestim
    bout_indices = (bout_start_times > prestim);
end

% Subset the bouts to keep for PSTH calculation
bout_delayed = bouts(bout_indices); %these get PSTH

%% Cut neural activity around bout, PSTH around start/end
threshold_std = 1e-3; % Threshold for standard deviation

psth_all      = cell(1, num_traces);
psth_all_end  = cell(1, num_traces);

% Initialize each cell as empty matrices
for j = 1:num_traces
    psth_all{j} = [];
    psth_all_end{j} = [];
end

log_skipped_bouts = struct('neuron', {}, 'bout', {}, 'reason', {});

for k = 1:length(bout_delayed)  % bout_delayed are the final chosen bouts
    % The absolute time (in seconds) of this bout's start
    bout_start_sec = bout_delayed{k}(1);
    bout_end_sec   = bout_delayed{k}(end);

    % Convert to frame indices in your truncated data
    %  (1) subtract time_offset, (2) multiply by fr, (3) round
    lick_idx_start = round( (bout_start_sec - time_offset) * matfile_fr );
    lick_idx_end   = round( (bout_end_sec   - time_offset) * matfile_fr );


    %% ----------------- PSTH around bout start ------------------ 
  for j = 1:num_traces
      start_idx = lick_idx_start - hw * matfile_fr;
      end_idx   = lick_idx_start + hw2 * matfile_fr;
        
        % Check array bounds
        if start_idx < 1 || end_idx > size(traces, 2)
            continue;
        end
% Extract neural segment
        seg_start = traces(j, start_idx:end_idx);
        pre_bout  = seg_start(1 : hw * matfile_fr);
        mean_pre = mean(pre_bout, 'omitnan');
        std_pre  = std(pre_bout, 0, 'omitnan');
        za = (seg_start - mean_pre) ./ std_pre;
        mean_post = mean(za(hw2*matfile_fr+1 : end), 'omitnan');
        if mean_post > 60
            log_skipped_bouts(end+1) = struct('neuron', j, 'bout', k, ...
               'reason','Abnormal z-score');
            continue;
        end
        
        % Store
        psth_all{j}(end+1, :) = za;
        
 %% ----------------- PSTH around bout end -------------------- %%
        start_idx_end = lick_idx_end - hw * matfile_fr;
        end_idx_end   = lick_idx_end + hw2 * matfile_fr;
        
        if start_idx_end < 1 || end_idx_end > size(traces, 2)
            continue;
        end
        
        seg_end    = traces(j, start_idx_end:end_idx_end);
        pre_bout_e = seg_end(1 : hw * matfile_fr);
        mean_pre_e = mean(pre_bout_e, 'omitnan');
        std_pre_e  = std(pre_bout_e, 0, 'omitnan');
        za_end = (seg_end - mean_pre_e) ./ std_pre_e;
        mean_post_end = mean(za_end(hw2*matfile_fr+1 : end), 'omitnan');
        if mean_post_end > 60
            log_skipped_bouts(end+1) = struct('neuron', j, 'bout', k, ...
               'reason','Abnormal z-score (end)');
            continue;
        end
        
        psth_all_end{j}(end+1, :) = za_end;
        
       
    end
end

%% flatten per-neuron cells into a big matrix for this file, so can
        % save responses per bout across all files, rather than just per neuron
        %— 1) flatten & save “all neurons” per‐bout PSTH, once per file ——
    thisPSTH      = vertcat( psth_all{:} );        % [sum_j bouts_j × window]
    thisPSTH_end  = vertcat( psth_all_end{:} );

    results(idx).psth_perBout     = thisPSTH;
    results(idx).psth_end_perBout = thisPSTH_end;

    allFilesPSTH     = [allFilesPSTH;     thisPSTH];
    allFilesPSTH_end = [allFilesPSTH_end; thisPSTH_end];


%% --- PSTH Metrics for All Bouts Combined ---
    % Average PSTH across all delayed bouts (from psth_all)
    PSTH_means = zeros(num_traces, hw2*matfile_fr + hw*matfile_fr + 1);
    PSTH_end_means = zeros(num_traces, hw2*matfile_fr + hw*matfile_fr + 1);
    for j = 1:num_traces
        if ~isempty(psth_all{j})
            PSTH_means(j,:) = mean(psth_all{j}, 1);
        end
    end
    for j = 1:num_traces
        if ~isempty(psth_all_end{j})
            PSTH_end_means(j,:) = mean(psth_all_end{j}, 1);
        end
    end
    
    nBins        = size(PSTH_means,2);
    baselineF    = hw  * matfile_fr;       % #frames before bout
    postF        = hw2 * matfile_fr;       % #frames after bout
    frame_range  = -baselineF : postF;      % e.g. -60:60 for hw=15s at 4Hz
    bout_t       = frame_range / matfile_fr;  % from -hw … +hw2 in seconds
    
    evoked_idx        = bout_t >= 0;       % logical mask, length = nBins
    evoked_time_vector = bout_t(evoked_idx);% times from 0 … +hw2
    
     % Compute the area under the curve (AUC) for each neuron using the trapezoidal rule:
    auc_response = trapz(evoked_time_vector, PSTH_means(:, evoked_idx), 2);

    % Compute the peak response for each neuron:
    peak_response = max(PSTH_means(:, evoked_idx), [], 2);
    
    % pick a smoothing window of ~33% of the evoked period
    spanFrames = round(0.33 * (postF+1));
    % movmean smooths along columns for each row (neuron)
    smooth_PSTH_means = movmean(PSTH_means, spanFrames, 2);
    peak_response_smoothed = max(smooth_PSTH_means(:, evoked_idx), [], 2);

    % Calculate time to peak response within the evoked window
    time_to_peak    = nan(size(PSTH_means,1),1);
    time_to_peak_sm = nan(size(smooth_PSTH_means,1),1);
    
    for n = 1:size(PSTH_means,1)
         evoked_data = PSTH_means(n, evoked_idx);
        if isempty(evoked_data) || all(isnan(evoked_data))
            % Instead of skipping, assign a default value (e.g., 0 or the overall mean)
            time_to_peak(n) = 0;  
        else
            [~, idxPeak] = max(evoked_data);
            time_to_peak(n) = evoked_time_vector(idxPeak);
        end
    end
    for n = 1:size(smooth_PSTH_means,1)
         evoked_data = smooth_PSTH_means(n, evoked_idx);
        if isempty(evoked_data) || all(isnan(evoked_data))
            % Instead of skipping, assign a default value (e.g., 0 or the overall mean)
            time_to_peak_sm(n) = 0;  
        else
            [~, idxPeak] = max(evoked_data);
            time_to_peak_sm(n) = evoked_time_vector(idxPeak);
        end
    end
    
 % --- compute for PSTH_end_means (trough instead of peak) ---

% 1) raw trough response
min_response_end             = min( PSTH_end_means(:, evoked_idx), [], 2 );

% 2) smoothed trough response
smooth_PSTH_end_means        = movmean(PSTH_end_means, spanFrames, 2);
min_response_smoothed_end    = min( smooth_PSTH_end_means(:, evoked_idx), [], 2 );

% 3) time to raw trough
time_to_min_end      = nan( size(PSTH_end_means,1), 1 );
for n = 1:size(PSTH_end_means,1)
    ed = PSTH_end_means(n, evoked_idx);
    if all(isnan(ed))
        time_to_min_end(n) = 0;
    else
        [~, idxMin] = min(ed);
        time_to_min_end(n) = evoked_time_vector(idxMin);
    end
end

% 4) time to smoothed trough
time_to_min_sm_end   = nan( size(smooth_PSTH_end_means,1), 1 );
for n = 1:size(smooth_PSTH_end_means,1)
    ed = smooth_PSTH_end_means(n, evoked_idx);
    if all(isnan(ed))
        time_to_min_sm_end(n) = 0;
    else
        [~, idxMin] = min(ed);
        time_to_min_sm_end(n) = evoked_time_vector(idxMin);
    end
end

%% Extract and analyze just the the first bout

for i = 1:num_neurons;
    psth_firstBout(i,:) = psth_all{i}(1,:);
end

for n = 1:size(psth_firstBout,1)
         evoked_data = psth_firstBout(n, evoked_idx);
        if isempty(evoked_data) || all(isnan(evoked_data))
            % Instead of skipping, assign a default value (e.g., 0 or the overall mean)
            time_to_peak(n) = 0;  
        else
            [~, idxPeak] = max(evoked_data);
            first_time_peak_all(n) = evoked_time_vector(idxPeak);
        end
end


first_peak_all = max(psth_firstBout(:, evoked_idx), [], 2);
% first_auc_all  = trapz(time_vector_firstBout(evoked_idx_first), psth_firstBout(:, evoked_idx_first), 2);
first_time_peak_all = NaN(num_neurons, 1);
for n = 1:num_neurons
    [~, pIdx] = max( psth_firstBout(n, evoked_idx), [], 2 );
    first_time_peak_all(n) = evoked_time_vector(pIdx);
end


% Optionally, store the result in your results structure:
    results(idx).psth_firstBout = psth_firstBout;
    results(idx).first_peak_all = first_peak_all;
%     results(idx).first_auc_all  = first_auc_all;
    results(idx).first_time_peak_all = first_time_peak_all;

    allCatRes.all.PSTH_first      = [allCatRes.all.PSTH_first;      psth_firstBout];
    allCatRes.all.first_peak      = [allCatRes.all.first_peak;      first_peak_all];
    allCatRes.all.first_time_peak = [allCatRes.all.first_time_peak; first_time_peak_all];
  %% --- PSTH Metrics for Subsequent Bouts (All But First) ---
    % For subsequent bouts, average the PSTH segments from bouts 2:end.
    % 1) collect each neuron’s mean‐over‐subsequent‐bouts into a cell
    temp = cellfun(@(x) mean(x(2:end,:),1,'omitnan'), psth_all,'UniformOutput', false);
    psth_subsequent_avg = vertcat(temp{:});   % gives [num_neurons × nBins]

%     % Compute metrics for subsequent bouts (using psth_subsequent_avg)
    auc_subsequent = trapz(evoked_time_vector, psth_subsequent_avg(:, evoked_idx), 2);
    peak_subsequent = max(psth_subsequent_avg(:, evoked_idx), [], 2);
    time_to_peak_subsequent = zeros(num_traces, 1);
    for n = 1:size(psth_subsequent_avg,1)
        data_sub = psth_subsequent_avg(n, evoked_idx);
        if isempty(data_sub) || all(isnan(data_sub))
            time_to_peak_subsequent(n) = 0;
        else
            [~, idxPeak] = max(data_sub);
            time_to_peak_subsequent(n) = evoked_time_vector(idxPeak);
        end
    end

    %% Store Combined PSTH Metrics for All Bouts (Combined & Subsequent)
    results(idx).PSTH_means = PSTH_means;
    results(idx).auc_response = auc_response;
    results(idx).peak_response = peak_response;
    results(idx).time_to_peak = time_to_peak;
    results(idx).psth_subsequent_avg = psth_subsequent_avg;
    results(idx).auc_subsequent = auc_subsequent;
    results(idx).peak_subsequent = peak_subsequent;
    results(idx).time_to_peak_subsequent = time_to_peak_subsequent;
    
%% ========Sort activated vs. non-activated cells (activated, inhibited, and none)================
    act = []; 
    inhib = [];
    none = [];

    prestim_idx = prestim * matfile_fr;
    idx_1200 = 1200 * matfile_fr;
    idx_900 = 900*matfile_fr;
    idx_1800 = 1800 * matfile_fr;
    idx_2100 = 2100 * matfile_fr;
    idx_2400 = 2397 * matfile_fr;

    for v = 1:size(trace_zs,1) %go through every neuron and calculate mean responses in specified windows
            a = mean(trace_zs(v,prestim_idx:idx_1200));
            b = mean(trace_zs(v,prestim_idx:idx_900));
            c = mean(trace_zs(v,idx_1200:idx_1800));
            d = mean(trace_zs(v,idx_1800:idx_2400));

            if a > 1 || b > 1
                act = [act v];

            elseif a< -1 %|| c < -1 || d < -1
                inhib = [inhib v];
    %             inhib_indices = [inhib_indices, remain_neuron_indices(v)];
            else
                none = [none v];
    %             none_indices = [none_indices, remain_neuron_indices(v)];
            end
    end

    Neurons_act = trace_zs(act, :);
    Neurons_inhib = trace_zs(inhib, :);
    Neurons_none = trace_zs(none, :);


    disp(['Number activated neurons:', num2str(length(act))]);
    disp(['Number inhibited neurons:', num2str(length(inhib))]);
    disp(['Number none neurons:', num2str(length(none))]);
%% ==============Sort for lick-activated, secondary classification================

    act_lick = [];
    act_non_lick = [];

    for j = 1:length(act)
        neuron_idx = act(j);
        if ~isempty(psth_all{neuron_idx}) && size(psth_all{neuron_idx}, 1) > 0

            bout_responses = psth_all{neuron_idx}(:, hw*matfile_fr+1:end); % Size: [num_bouts x hw]

            % Find the maximum z-score in the post-bout period for each bout
            max_responses = max(bout_responses, [], 2); % Max along the time axis for each bout

            % Count the number of bouts where the neuron responds (z-score > threshold)
            num_bouts = size(bout_responses, 1);
            num_responsive_bouts = sum(max_responses > zs_threshold);

            % Calculate the percentage of responsive bouts
            responsive_percentage = num_responsive_bouts / num_bouts;

            % Classify neuron based on the responsive percentage
            if responsive_percentage > bout_threshold
                act_lick = [act_lick, neuron_idx];
            else
                act_non_lick = [act_non_lick, neuron_idx];
            end
        else
            % If no data is available, classify as 'remain' or handle accordingly
            act_non_lick = [act_non_lick, neuron_idx];
        end
    end

    %% Store z-scored traces for cells identified as lick_activated and other
    Neurons_lick_act = trace_zs(act_lick, :);
    Neurons_act_non_lick = trace_zs(act_non_lick, :);

%% Split out the psth per but to categorize

% now loop over your categories using the same act/inhib/none variables:
for c = 1:numel(categories)
  cat     = categories{c};       % 'act','inhib','none'
  neurons = eval(cat);           % now eval('act') → your act vector

  catPSTH     = vertcat( psth_all{neurons} );
  catPSTH_end = vertcat( psth_all_end{neurons} );

  results(idx).(cat).psth_perBout     = catPSTH;
  results(idx).(cat).psth_end_perBout = catPSTH_end;

  allCatRes.(cat).psth_perBout     = [allCatRes.(cat).psth_perBout;     catPSTH];
  allCatRes.(cat).psth_end_perBout = [allCatRes.(cat).psth_end_perBout; catPSTH_end];
  allCatRes.(cat).PSTH_first      = [allCatRes.(cat).PSTH_first;      psth_firstBout(neurons,:)];
  allCatRes.(cat).first_peak      = [allCatRes.(cat).first_peak;      first_peak_all(neurons)];
  allCatRes.(cat).first_time_peak = [allCatRes.(cat).first_time_peak; first_time_peak_all(neurons)];
end


    
    %% ===========Caclulate averaged z-score responses over different time periods for each neuron category, also normalized to licks=================

% 1) Time‑window metrics (10–20, 20–30, 10–30 min)
time_windows = [ ...
    600,  1200;   % A: 10–20 min
    1201, 1800;   % B: 20–30 min
    1801, 600+poststim % C: 30–40 min
    600, 600+poststim];  % D: 10–40 min

time_window_indices = round(time_windows * matfile_fr);
n_win = size(time_window_indices,1);

lick_counts_window = zeros(n_win,1);
for w=1:n_win
  fs = time_window_indices(w,1);
  fe = time_window_indices(w,2);
  t0 = (fs-1)/matfile_fr + time_offset;
  t1 = (fe-1)/matfile_fr + time_offset;
  lick_counts_window(w) = sum(filtered_licks>=t0 & filtered_licks<=t1);
end

% --- 2) Compute z-score over the different time windows, and also normalized by licks ---
% Pre-allocate for all neurons
raw_mean_all   = NaN(num_neurons,n_win);
licknorm_all   = NaN(num_neurons,n_win);

for w=1:n_win
  fs = time_window_indices(w,1);
  fe = time_window_indices(w,2);
  seg = trace_zs(:,fs:fe);
  raw_mean_all(:,w) = mean(seg,2,'omitnan');
  if lick_counts_window(w)>0
    licknorm_all(:,w) = raw_mean_all(:,w)/lick_counts_window(w);
  end
end

results(idx).raw_mean   = raw_mean_all;
results(idx).licknorm   = licknorm_all;

% 2a) Accumulate into allCatRes
allCatRes.all.raw_mean = [allCatRes.all.raw_mean; raw_mean_all];
allCatRes.all.licknorm = [allCatRes.all.licknorm; licknorm_all];

% 3) Slice out each category for raw_mean and licknorm
    for c = 1:numel(categories)
        cat = categories{c};
        idx = eval(cat);  % e.g. act, none, or inhib
        
        % append them into your global structure
        allCatRes.(cat).raw_mean   = [allCatRes.(cat).raw_mean;   raw_mean_all(idx,:)];
        allCatRes.(cat).licknorm   = [allCatRes.(cat).licknorm;  licknorm_all(idx,:)];
    end


%% ========Calculate Average Z-Score During Each BOUT, also normalized to licks================
% 4) Bout Metrics
% Convert sub-windows to frame indices
sub_windows = [600,1200; 1201,1800; 600,1800];
sub_windows = double(sub_windows);        % make sure it’s doubles
time_offset = double(time_offset);
sub_windows_idx     = round((sub_windows - time_offset) * matfile_fr);
sub_windows_idx(:,1)= max(1, sub_windows_idx(:,1));
sub_windows_idx(:,2)= min(num_frames, sub_windows_idx(:,2));
n_sub       = size(sub_windows_idx,1);

% Pre‑allocate matrices for ALL neurons × bouts
num_bouts = numel(bouts);
mean_z_per_bout_all    = NaN(num_neurons,num_bouts);
auc_per_bout_all       = NaN(num_neurons,num_bouts);
mean_z_per_bout_subwin = NaN(num_neurons,num_bouts,n_sub);

% Loop over bouts once, compute for every neuron
for b=1:num_bouts
  bs = round((bouts{b}(1)   - time_offset)*matfile_fr);
  be = round((bouts{b}(end) - time_offset)*matfile_fr);
  bs = max(1,bs);  be = min(num_frames,be);
  if bs<be
    seg = trace_zs(:,bs:be);
    mean_z_per_bout_all(:,b) = mean(seg,2,'omitnan');
    auc_per_bout_all(:,b)    = trapz(seg,2) * (1/matfile_fr);
    for sw=1:n_sub
      os = max(bs,sub_windows_idx(sw,1));
      oe = min(be,sub_windows_idx(sw,2));
      if os<oe
        subseg = trace_zs(:,os:oe);
        mean_z_per_bout_subwin(:,b,sw) = mean(subseg,2,'omitnan');
      end
    end
  end
end

% 4) Compute per‑neuron averages across bouts
mean_bout_average_zs = mean(mean_z_per_bout_all,2,'omitnan');
mean_bout_AUC_zs     = mean(auc_per_bout_all,2,'omitnan');
mean_bout_subwin     = squeeze(mean(mean_z_per_bout_subwin,2,'omitnan'));

% --- Normalize per bout by lick count ---
zscore_per_lick_all = NaN(num_neurons,num_bouts);
for b=1:num_bouts
  nlick=licks_per_bout(b);
  if nlick>0
    zscore_per_lick_all(:,b) = mean_z_per_bout_all(:,b)/nlick;
  end
end
per_neuron_norm_z = mean(zscore_per_lick_all,2,'omitnan');

% ---- ALL neurons ----
allCatRes.all.per_neuron_boutmean_z    = [allCatRes.all.per_neuron_boutmean_z;   mean_bout_average_zs];
allCatRes.all.per_neuron_boutmean_auc  = [allCatRes.all.per_neuron_boutmean_auc; mean_bout_AUC_zs];
allCatRes.all.per_neuron_boutnorm_z    = [allCatRes.all.per_neuron_boutnorm_z;   per_neuron_norm_z];
allCatRes.all.mean_bout_subwin     = [allCatRes.all.mean_bout_subwin;   mean_bout_subwin];

allCatRes.all.per_bout_z = [allCatRes.all.per_bout_z; mean_z_per_bout_all(:)];
allCatRes.all.per_bout_auc       = [allCatRes.all.per_bout_auc;   auc_per_bout_all(:)];
allCatRes.all.per_bout_norm_z    = [allCatRes.all.per_bout_norm_z; zscore_per_lick_all(:)];

% 5) ---- Slice into act/inhib/none ----
for iCat = 1:numel(categories)
    cat     = categories{iCat};
    idx_neur = eval(cat);

    % per‑neuron
    allCatRes.(cat).per_neuron_boutmean_z   = [allCatRes.(cat).per_neuron_boutmean_z;   mean_bout_average_zs(idx_neur)];
    allCatRes.(cat).per_neuron_boutmean_auc = [allCatRes.(cat).per_neuron_boutmean_auc; mean_bout_AUC_zs(idx_neur)];
    allCatRes.(cat).per_neuron_boutnorm_z  = [allCatRes.(cat).per_neuron_boutnorm_z;   per_neuron_norm_z(idx_neur)];
    allCatRes.(cat).mean_bout_subwin  = [allCatRes.(cat).mean_bout_subwin;  mean_bout_subwin(idx,:)];

    % per‑bout
    temp = mean_z_per_bout_all(idx_neur,:);    % N_cat × B bouts
    allCatRes.(cat).per_bout_z = [allCatRes.(cat).per_bout_z;temp(:)];   % flatten to (N_cat*B)×1                          
    tempA = auc_per_bout_all(idx_neur,:);
    allCatRes.(cat).per_bout_auc        = [allCatRes.(cat).per_bout_auc; tempA(:)];
    tempN = zscore_per_lick_all(idx_neur,:);
    allCatRes.(cat).per_bout_norm_z     = [allCatRes.(cat).per_bout_norm_z; tempN(:)];
%     allCatRes.(cat).per_bout_z          = [allCatRes.(cat).per_bout_z;     mean_z_per_bout_all(idx_neur,:)];
%     allCatRes.(cat).per_bout_auc        = [allCatRes.(cat).per_bout_auc;   auc_per_bout_all(idx_neur,:)];
%     allCatRes.(cat).per_bout_norm_z     = [allCatRes.(cat).per_bout_norm_z; zscore_per_lick_all(idx_neur,:)];
end

% === Store in results if desired ===
% results(idx).mean_bout_average_zs = mean_bout_average_zs;
% results(idx).mean_bout_AUC_zs     = mean_bout_AUC_zs;
% results(idx).mean_bout_subwin     = mean_bout_subwin;
% results(idx).zscore_per_lick_all  = zscore_per_lick_all;


%% 4) ================== Correlation with # of Licks =================

% 3) Correlation with # of licks (per‑bout raw mean)
        R_raw = NaN(num_neurons,1);
        P_raw = NaN(num_neurons,1);
for n=1:num_neurons
  y = mean_z_per_bout_all(n,:);
  x = licks_per_bout;
  v = ~isnan(y)&~isnan(x);
  if sum(v)>1
    [R_raw(n),P_raw(n)] = corr(y(v)', x(v)', 'Type','Pearson');
  end
end
% results(idx).R_raw = R_raw;
% results(idx).P_raw = P_raw;

    for n = 1:num_neurons
        data_y   = mean_z_per_bout_all(n,:);
        data_x   = licks_per_bout;  % # of licks in each bout
        valid    = ~isnan(data_y) & ~isnan(data_x);   % for raw

        if sum(valid) > 1
            [R_raw(n), P_raw(n)] = corr(data_y(valid)', data_x(valid)', 'Type','Pearson');
        end
    end
    
    % 1) Identify the subset of "act" neurons that have R_raw < 0
    act_neg_idx = act(R_raw(act) < 0); 
    act_corr = act(R_raw(act) <.05);
    % 2) Extract their traces
    Neurons_act_neg = trace_zs(act_neg_idx, :);
    Neurons_correlated = trace_zs(act_corr, :);
    
    % 3) Store them in 'results' so you can quickly retrieve
%     results(idx).Neurons_act_neg = Neurons_act_neg;
%     results(idx).act_neg_idx     = act_neg_idx;

% R_raw(n), P_raw(n) => correlation & p-value for raw average z-score for
% each neuron across bouts

% Create a structure to hold the aggregated data for each category
for icat = 1:length(categories)
    cat = categories{icat};  % e.g., 'act', 'inhib', or 'none'
    idx=eval(cat);
          z    = mean_z_per_bout_all(idx,:);
          a    = auc_per_bout_all(idx,:);
          
      ltk  = repmat(licks_per_bout,numel(idx),1);
      % Flatten the matrix into one vector (each element represents a bout from one neuron)
      allCatRes.(cat).collapsed_bout_means = [allCatRes.(cat).collapsed_bout_means; z(:)];
      allCatRes.(cat).collapsed_bout_auc   = [allCatRes.(cat).collapsed_bout_auc;   a(:)];
      allCatRes.(cat).collapsed_bout_licks = [allCatRes.(cat).collapsed_bout_licks; ltk(:)];
  
end

z_all = mean_z_per_bout_all(:);
a_all = auc_per_bout_all(:);
l_all = repmat(licks_per_bout,num_neurons,1);
allCatRes.all.collapsed_bout_means = [allCatRes.all.collapsed_bout_means; z_all];
allCatRes.all.collapsed_bout_auc   = [allCatRes.all.collapsed_bout_auc;   a_all];
allCatRes.all.collapsed_bout_licks = [allCatRes.all.collapsed_bout_licks; l_all(:)];

%% ==========Determine z-score when NOT licking/between bouts=======================
inBoutMask = false(1, num_frames);
for b = 1:length(bouts)
    bout_start_sec = bouts{b}(1);
    bout_end_sec   = bouts{b}(end);

    sdx = round((bout_start_sec - time_offset)*matfile_fr);
    edx = round((bout_end_sec   - time_offset)*matfile_fr);
    sdx = max(sdx,1); 
    edx = min(edx,num_frames);
    if sdx < edx
        inBoutMask(sdx:edx) = true;
    end
end
nonBoutMask = ~inBoutMask;

% 2) Entire session average "non-bout" z-score
avg_nonBout_all = NaN(num_neurons,1);
frames_nonbout = find(nonBoutMask);
for n = 1:num_neurons
    avg_nonBout_all(n) = mean(trace_zs(n, frames_nonbout), 'omitnan');
end
% results(idx).avg_zscore_nonbouts_all = avg_nonBout_all;

% 3) Time windows for the 0–10, 20–30 min
time_windows_nb = [ ...
    prestim,   1200;  % 0–10 min
    1201,1800;  % 10–20
    1801,prestim+poststim; % 20–30
];
num_wins = size(time_windows_nb,1);
avg_nonBout_windows = NaN(num_neurons,num_wins);

for w = 1:num_wins
    start_sec = time_windows_nb(w,1);
    end_sec   = time_windows_nb(w,2);

    start_idx_win = round(start_sec * matfile_fr +1);
    end_idx_win   = round(end_sec   * matfile_fr);
    start_idx_win = max(start_idx_win,1);
    end_idx_win   = min(end_idx_win,num_frames);

    if start_idx_win < end_idx_win
        validFrames = start_idx_win:end_idx_win;
        validFrames = validFrames(nonBoutMask(validFrames));
        for n = 1:num_neurons
            if ~isempty(validFrames)
                avg_nonBout_windows(n,w) = mean(trace_zs(n,validFrames), 'omitnan');
            end
        end
    end
end

%————— compute session‐wide lick count —————
sessLicks = total_filtered_licks;

% %————— compute avg non-bout z-score over last 10 min, this is so can run licking_non_bout_bin after —————
nFrames     = size(trace_zs,2);
last10Win   = max(1, nFrames - 600*matfile_fr + 1) : nFrames;
validFrames = last10Win(nonBoutMask(last10Win));          % only non-bout
% per‐neuron
avgPerNeuron = mean(trace_zs(:, validFrames), 2, 'omitnan' );
% activated neurons
allLicks_Act   = [allLicks_Act;   repmat(sessLicks,numel(act),1)];
allZ_Act       = [allZ_Act;       avgPerNeuron(act)];

% inhibited neurons
allLicks_Inhib = [allLicks_Inhib; repmat(sessLicks,numel(inhib),1)];
allZ_Inhib     = [allZ_Inhib;     avgPerNeuron(inhib)];

% results(idx).avg_zscore_nonbouts_windows = avg_nonBout_windows;
% columns => time windows, rows => neurons

%% ========Graph the neurons and make sure the lick bouts are properly aligned=============
 % Time vector for the neural data (0 = first frame of extracted segment)
    time_vector = (0 : size(trace_zs,2)-1) / matfile_fr; 
if graph_data == 1
    
    % 11) Plot Example: All Cell Traces with Lick Events Overlaid
    figure; hold on;
    offsetIncrement = 20;  % vertical offset for each neuron’s trace
    num_neurons = size(trace_zs,1);
    yMin = 0;
    yMax = num_neurons * offsetIncrement;

    % Time vector for the neural data (0 = first frame of extracted segment)
    time_vector = (0 : size(trace_zs,2)-1) / matfile_fr; 

    % Plot vertical lines for each lick time 
    % (Shift them if you want the PSTH window to be zeroed at “stim_time - prestim”)
    for i = 1:length(filtered_licks)
        xLine = filtered_licks(i) - time_offset; 
        if xLine >= 0 && xLine <= max(time_vector)
            line([xLine xLine], [yMin yMax], 'Color',[0.8 0.8 0.8], 'LineWidth',1.2);
        end
    end

    % Plot each neuron’s z-scored trace, offset by i*offsetIncrement
    for i2 = 1:num_neurons
        plot(time_vector, trace_zs(i2,:) + (i2-1)*offsetIncrement);
    end
    xlabel('Time (s)');
    ylabel('Cell Activity (Z-score + offset)');
    title('All Cell Traces with Lick Events');
    hold off;

    %Plot the lick-activated neurons
    figure; hold on;
    yMax =  size(Neurons_lick_act,1) * offsetIncrement;
    % Plot vertical lines for each lick time 
    % (Shift them if you want the PSTH window to be zeroed at “stim_time - prestim”)
    for i = 1:length(filtered_licks)
        xLine = filtered_licks(i) - time_offset; 
        if xLine >= 0 && xLine <= max(time_vector)
            line([xLine xLine], [yMin yMax], 'Color',[0.8 0.8 0.8], 'LineWidth',1.2);
        end
    end

    % Plot each neuron’s z-scored trace, offset by i*offsetIncrement
    for i2 = 1:size(Neurons_lick_act,1)
        plot(time_vector, Neurons_lick_act(i2,:) + (i2-1)*offsetIncrement);
    end
    xlabel('Time (s)');
    ylabel('Cell Activity (Z-score + offset)');
    title('Lick-Activated Neurons');
    hold off;

    %Plot the lick-activated neurons
    figure; hold on;
    yMax =  size(Neurons_act_non_lick,1) * offsetIncrement;
    % Plot vertical lines for each lick time 
    % (Shift them if you want the PSTH window to be zeroed at “stim_time - prestim”)
    for i = 1:length(filtered_licks)
        xLine = filtered_licks(i) - time_offset; 
        if xLine >= 0 && xLine <= max(time_vector)
            line([xLine xLine], [yMin yMax], 'Color',[0.8 0.8 0.8], 'LineWidth',1.2);
        end
    end

    % Plot each neuron’s z-scored trace, offset by i*offsetIncrement
    for i2 = 1:size(Neurons_act_non_lick,1)
        plot(time_vector, Neurons_act_non_lick(i2,:) + (i2-1)*offsetIncrement);
    end
    xlabel('Time (s)');
    ylabel('Cell Activity (Z-score + offset)');
    title('Non-Lick-Activated Neurons');
    hold off;
end

   %% % Compute category-specific metrics by indexing the overall metrics
% results(idx).first_peak_act = first_peak_all(act);
% % results(idx).first_auc_act  = first_auc_all(act);
% results(idx).first_time_peak_act = first_time_peak_all(act);
% 
% results(idx).first_peak_inhib = first_peak_all(inhib);
% % results(idx).first_auc_inhib  = first_auc_all(inhib);
% results(idx).first_time_peak_inhib = first_time_peak_all(inhib);
% 
% results(idx).first_peak_none = first_peak_all(none);
% % results(idx).first_auc_none  = first_auc_all(none);
% results(idx).first_time_peak_none = first_time_peak_all(none);
%     

%% ========== Saving results ==============
results_by_category = struct();  % will store all category-based outputs
cat_indices.act   = act;   
cat_indices.inhib = inhib; 
cat_indices.none  = none; 

for c = 1:length(categories)
    cat_name  = categories{c};        % e.g. 'act'
    idx_cat   = cat_indices.(cat_name);   % e.g. act's indices
    
     % 1) Per-bout raw average z‑score for these neurons
%     mean_z_cat = mean_bout_average_zs(idx_cat, :);   % [#neuron_in_cat x #bouts]
%     % 2) Per-bout for normalized z‑score (per lick) for these neurons
%     z_lick_cat = zscore_per_lick_bout_mean(idx_cat, :);   % same size
    % 3) Correlations
%     r_raw_cat  = R_raw(idx_cat);
%     p_raw_cat  = P_raw(idx_cat);
%     % Store them in a struct for each category
%     % Store them in a struct under 'results_by_category.(cat_name)'
%     % Or a single average across all bouts & neurons:
%     results_by_category.(cat_name).mean_z_cat  = mean_z_cat;
%     results_by_category.(cat_name).z_lick_cat  = z_lick_cat;
%     results_by_category.(cat_name).R_raw   = r_raw_cat;
%     results_by_category.(cat_name).P_raw   = p_raw_cat;
end
% results(idx).by_category = results_by_category;
%% 12) Package Results for This Session
%%Store Results for the Current File Pair
% results(idx).filename = thisTxt;
% results(idx).Neurons_inhib = Neurons_inhib;
% results(idx).Neurons_none = Neurons_none;
% results(idx).Neurons_lick_act = Neurons_lick_act;
% results(idx).Neurons_act_non_lick = Neurons_act_non_lick;
% results(idx).Neurons_act = Neurons_act;
% results(idx).Neurons_act_neg = Neurons_act_neg;
% 
% results(idx).PSTH_inhib = PSTH_means(inhib, :);
% results(idx).PSTH_none = PSTH_means(none, :);
% results(idx).PSTH_lick_act = PSTH_means(act_lick,:);
% results(idx).PSTH_act_non_lick = PSTH_means(act_non_lick, :);
% 
% results(idx).PSTHend_inhib = PSTH_end_means(inhib, :);
% results(idx).PSTHend_none = PSTH_end_means(none, :);
% results(idx).PSTHend_lick_act = PSTH_end_means(act_lick,:);
% results(idx).PSTHend_act_non_lick = PSTH_end_means(act_non_lick, :);

% Per-bout data for all neurons
% results(idx).mean_z_per_bout     = mean_z_per_bout;  % [num_neurons x num_bouts]
% results(idx).zscore_per_lick     = zscore_per_lick;  % [num_neurons x num_bouts]
% results(idx).licks_per_bout      = licks_per_bout;   % [1 x num_bouts]
% results(idx).bout_durations        = bout_durations;       % e.g. cellfun(@(x) x(end)-x(1), bouts)
% results(idx).avg_bout_duration     = avg_bout_duration;     % e.g. mean(bout_durations)
% results(idx).total_bouts           = total_bouts;           % e.g. length(bouts)
% results(idx).total_filtered_licks  = total_filtered_licks;  % e.g. length(filtered_licks)
% results(idx).R_raw               = R_raw;
% results(idx).P_raw               = P_raw;
% results(idx).meanR_activated = meanRact;
% results(idx).time_to_peak = time_to_peak; % for all bouts combined
  
% Now include the category-based breakdown:
% results(idx).by_category         = results_by_category;

%% Save the file
save([thisTxt '_results.mat'],'results','-v7.3'); 

%% Append the category-based arrays to "allCatRes"
for c = 1:length(categories)
        cat_name = categories{c};
        cat_idx  = cat_indices.(cat_name); % e.g. [2,5,7,...]
        
        % 1) "raw_mean" from mean_responses.(cat_name).raw_mean => [#neurons_cat x #windows]
%         raw_mean_cat  = mean_responses.(cat_name).raw_mean;
%         licknorm_cat  = mean_responses.(cat_name).licknorm;

        % 2) "R_raw" etc. We subset them by cat_idx
        R_cat = R_raw(cat_idx);
        P_cat = P_raw(cat_idx);

        % 3) "mean_bout_average_zs" (size = [num_neurons x 1]), also subset
%         bout_avgZS_cat = mean_bout_average_zs(cat_idx);
        
        % 4) "zscore_per_lick_mean" likewise
%         zscore_lickM_cat = zscore_per_lick_bout_mean(cat_idx);
        
        %5) 
        zscore_nonbouts_all_cat = avg_nonBout_all(cat_idx);
        zscore_nonbouts_windows_cat = avg_nonBout_windows(cat_idx,:);
        
        %6)
        PSTH_act = PSTH_means(cat_idx,:);
        PSTH_means_smoothed = smooth_PSTH_means(cat_idx,:);
        PSTH_end_neurons = PSTH_end_means(cat_idx,:);
        PSTH_mean_cat = mean(PSTH_means(cat_idx,:), 1, 'omitnan');
        % For the activated group:
        PSTH_peak = peak_response(cat_idx);
        PSTH_peak_smoothed = peak_response_smoothed(cat_idx);
        PSTH_auc = auc_response(cat_idx,:);
        PSTH_mean_peak = mean(PSTH_peak, 1, 'omitnan');
        PSTH_mean_auc = mean(PSTH_auc, 1, 'omitnan');
        PSTH_time_peak = time_to_peak(cat_idx);
        PSTH_time_peak_smoothed = time_to_peak_sm(cat_idx);
        PSTH_end_time_peak_smoothed = time_to_min_sm_end(cat_idx)
        
        % For subsequent bouts:
        % (Assuming psth_subsequent_avg is [num_neurons x #bins])
        PSTH_sub_cat = psth_subsequent_avg(cat_idx, :);
        subsequent_peak_cat = max(PSTH_sub_cat(:, evoked_idx), [], 2);
        subsequent_auc_cat = trapz(time_vector(evoked_idx), PSTH_sub_cat(:, evoked_idx), 2);
        subsequent_time_peak = NaN(size(PSTH_sub_cat,1),1);
        for n = 1:size(PSTH_sub_cat,1)
            [~, pIdx] = max(PSTH_sub_cat(n, evoked_idx));
            subsequent_time_peak(n) = evoked_time_vector(pIdx);
        end
        % Now we vertically concatenate:

%         allCatRes.(cat_name).raw_mean = [allCatRes.(cat_name).raw_mean ; raw_mean_cat];
%         allCatRes.(cat_name).licknorm = [allCatRes.(cat_name).licknorm; licknorm_cat];
        allCatRes.(cat_name).R_raw = [allCatRes.(cat_name).R_raw; R_cat];
        allCatRes.(cat_name).P_raw = [allCatRes.(cat_name).P_raw; P_cat];
%         allCatRes.(cat_name).bout_avgZS = [allCatRes.(cat_name).bout_avgZS; bout_avgZS_cat];
%         allCatRes.(cat_name).zscore_lickM = [allCatRes.(cat_name).zscore_lickM;zscore_lickM_cat];
        allCatRes.(cat_name).avg_nonBout_all = [allCatRes.(cat_name).avg_nonBout_all;zscore_nonbouts_all_cat];
        allCatRes.(cat_name).avg_nonBout_windows = [allCatRes.(cat_name).avg_nonBout_windows;zscore_nonbouts_windows_cat];
        allCatRes.(cat_name).PSTH_act = [allCatRes.(cat_name).PSTH_act;PSTH_act];
        allCatRes.(cat_name).PSTH_means_smoothed = [allCatRes.(cat_name).PSTH_means_smoothed;PSTH_means_smoothed];
        allCatRes.(cat_name).PSTH_end_neurons = [allCatRes.(cat_name).PSTH_end_neurons;PSTH_end_neurons];
        allCatRes.(cat_name).PSTH_mean_cat = [allCatRes.(cat_name).PSTH_mean_cat;PSTH_mean_cat];
        allCatRes.(cat_name).PSTH_peak = [allCatRes.(cat_name).PSTH_peak;PSTH_peak];
        allCatRes.(cat_name).PSTH_peak_smoothed = [allCatRes.(cat_name).PSTH_peak_smoothed;PSTH_peak_smoothed];
        allCatRes.(cat_name).PSTH_auc = [allCatRes.(cat_name).PSTH_auc;PSTH_auc];
        allCatRes.(cat_name).PSTH_mean_peak = [allCatRes.(cat_name).PSTH_mean_peak;PSTH_mean_peak];
        allCatRes.(cat_name).PSTH_mean_auc = [allCatRes.(cat_name).PSTH_mean_auc;PSTH_mean_auc];
        allCatRes.(cat_name).PSTH_time_peak = [allCatRes.(cat_name).PSTH_time_peak;PSTH_time_peak];
        allCatRes.(cat_name).PSTH_time_peak_smoothed = [allCatRes.(cat_name).PSTH_time_peak_smoothed;PSTH_time_peak_smoothed];
        allCatRes.(cat_name).PSTH_end_time_peak_smoothed = [allCatRes.(cat_name).PSTH_end_time_peak_smoothed;PSTH_end_time_peak_smoothed];
        allCatRes.(cat_name).PSTH_subsequent_avg = [allCatRes.(cat_name).PSTH_subsequent_avg; PSTH_sub_cat];
        allCatRes.(cat_name).subsequent_peak = [allCatRes.(cat_name).subsequent_peak; subsequent_peak_cat];
        allCatRes.(cat_name).subsequent_auc = [allCatRes.(cat_name).subsequent_auc; subsequent_auc_cat];
        allCatRes.(cat_name).subsequent_time_peak = [allCatRes.(cat_name).subsequent_time_peak; subsequent_time_peak];
  
end
%% ====== Update "allNeurons" to Concatenate Data ======
    allNeurons.act       = [allNeurons.act;       trace_zs(act, :)];
    allNeurons.inhib     = [allNeurons.inhib;     trace_zs(inhib, :)];
    allNeurons.none      = [allNeurons.none;      trace_zs(none, :)];
    allNeurons.act_neg = [allNeurons.act_neg;       trace_zs(act_neg_idx, :)];
    allNeurons.all = [allNeurons.all; trace_zs];

end

save('allFiles_concatenated.mat','allNeurons','allCatRes','results','allFilesPSTH','allFilesPSTH_end','-v7.3');
disp('All files processed. Combined data saved in allFiles_concatenated.mat');