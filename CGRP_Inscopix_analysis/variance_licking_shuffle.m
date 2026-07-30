%% Chow Feeding: Lick-Based Neural Variance Partitioning
%  Asks how much variance in single-cell z-scored activity is explained by:
%    1. bout_active      — binary lick-bout state (1 during bout, 0 during IBI)
%    2. time_from_access — linear time elapsed since food access (s)
%                          captures slow neural drift, arousal decay, and
%                          satiety-related ramping independent of lick rate
%
%  Neural traces are lightly Gaussian-smoothed before regression to reduce
%  high-frequency noise that inflates residuals and deflates R².
%  Smoothing sigma is set in USER PARAMETERS (set to 0 to disable).
%
%  Bout definition: no minimum bout duration — only min lick count applies,
%  so short but valid licking episodes are included in the bout-active vector.
%
%  Cumulative lick curves are computed and plotted as a reference but are
%  NOT used as a predictor in this version.
%
%  Inputs:
%    - txtfiles{}  : one .txt metadata file per animal
%                    [csvFile  fr  stim_mm:ss  prestim_s  end_s  videoOffset_mm:ss]
%                    (same format as FreeLicking_analysis_v3 / Chow.txt)
%    - infiles{}   : one lick GPIO CSV per animal (Inscopix GPIO export)
%                    Must contain 'ChannelName' column with 'GPIO-1' rows
%                    and columns 'Value' + 'Time_s_'
%
%  BCJ 2026

clear; close all;

%% ========== USER PARAMETERS ==========
% --- Lick filtering ---
min_ili             = 0.1;      % Min inter-lick interval (s) — remove double-detections

% --- Bout parameters ---
inter_bout_interval = 10;       % IBI gap (s) that splits licks into bouts
min_licks_per_bout  = 3;        % Min licks for a valid bout
%  No minimum bout duration — short bouts are included

% --- Smoothing ---
smooth_sigma_s      = 2;        % Gaussian kernel SD in seconds (2s = 8 frames at 4 Hz)
%  Set to 0 to disable smoothing. Applied after z-scoring, per neuron.

% --- Recording structure ---
prestim             = 600;      % Pre-access baseline (s) — used for z-scoring

% --- Cumulative lick plot (reference only, not a predictor) ---
cumul_plot_dur_s    = 3600;     % Duration to plot (s from food access)
cumul_bin_s         = 60;       % Bin size for group cumulative plot (s)

graph_data          = 1;        % 1 = produce figures

% --- Shuffle control ---
n_shuffles          = 500;      % circular shifts per neuron
% Circular shift preserves the temporal autocorrelation of the bout
% regressor while breaking its alignment with neural activity.

%% ========== FILE LISTS ==========
% One txt file per animal paired 1-to-1 with one lick GPIO CSV.
% txtfiles use Chow.txt-style format (see Chow.txt in project).
% infiles are Inscopix GPIO exports (see FreeLicking_analysis_v3 infiles{}).

txtfiles = {
    % ---- Chow cohort (uncomment/edit as needed) ----
     'FastedEnsure_BCJ1.txt'
    'FastedEnsure_AR3.txt'
    'FastedEnsure_AR6.txt'
    'FastedEnsure_CGRP1-6.txt'
    'FastedEnsure2_CGRP1-7.txt'
    'FastedEnsure_CGRP1-8.txt'
      'FastedEnsure_BCJ7.txt'
     'FastedEnsure_CGRP5.txt'
    'FastedEnsure_CGRP1-17.txt'
    'FastedEnsure_CGRP1-18.txt'
};

infiles = {
    % ---- Lick GPIO CSVs (must match txtfiles row-for-row) ----
        '2022-03-30-13-40-13_BCJ1_FastedEnsure.csv'
    '2021-12-16-16-16-24_AR3Ensure_licks.csv'
    '2021-12-16-13-30-08_AR6Ensure-licks.csv'
    'BCJ_CGRP1-6_20231222_lickingEnsure-gpio.csv'
    'BCJ_CGRP1-7_20231227_lickingEnsure2-gpio.csv'
    'BCJ_CGRP1-8_20231222_lickingEnsure-gpio.csv'
      '20220817_BCJ7_Ensure_licks.csv'
    'CGRP5_ensure_licks.csv'
    'CGRP1-17_ensure_licks.csv'
    'CGRP1-18_ensure_licks.csv'
};

num_files = length(txtfiles);
if length(infiles) ~= num_files
    error('txtfiles and infiles must have the same number of entries!');
end

%% ========== OUTPUT STRUCTURES ==========
categories = {'act','inhib','none'};
cat_colors = [0.85 0.33 0.10;   % act   = orange-red
              0.00 0.45 0.74;   % inhib = blue
              0.50 0.50 0.50];  % none  = grey

empty_cat = struct('act',[],'inhib',[],'none',[]);

% Individual R2 (each predictor alone, OLS)
allR2_bout = empty_cat;   % Lick-bout binary (1 = licking, 0 = IBI)
allR2_time = empty_cat;   % Time-from-access (normalised linear ramp)

% Partial R2: unique variance explained above the other predictor
allPartR2_bout = empty_cat;
allPartR2_time = empty_cat;

% Full model R2 (bout_active + time_from_access together)
allR2_full = empty_cat;

% Shuffle null distributions — circular shift, mean R2 per neuron
allR2_bout_shuf = empty_cat;   % mean shuffle R2 for bout predictor
allR2_time_shuf = empty_cat;   % mean shuffle R2 for time predictor

% Cumulative lick curves (reference / plotting only — not a predictor)
cumul_time_axis = (0 : cumul_bin_s : cumul_plot_dur_s)';
allCumulCurves  = [];   % [nAnimals x nTimeBins]

allCounts = struct('act',0,'inhib',0,'none',0);

%% ========== FILE LOOP ==========
for idx = 1:num_files

    thisTxt = txtfiles{idx};
    thisCSV = infiles{idx};
    fprintf('\n[%d/%d] %s  |  %s\n', idx, num_files, thisTxt, thisCSV);

    %% ----------------------------------------------------------------
    %%  STEP A — READ METADATA FROM TXT FILE
    %% ----------------------------------------------------------------
    % Supports two txt formats:
    %   (a) Chow.txt-style: one header row, data rows below (tab-delimited)
    %       [csvFile  fr  stim_mm:ss  prestim_s  end_s  videoOffset_mm:ss]
    %   (b) FreeLicking-style: no explicit header, data starts at row 2
    % Tries (a) first; falls back to (b) if empty.

    fileID = fopen(thisTxt,'r');
    data   = textscan(fileID, '%s %d %s %d %d %s', 'Delimiter','\t', 'HeaderLines',1);
    fclose(fileID);

    if isempty(data{1})
        fileID = fopen(thisTxt,'r');
        data   = textscan(fileID, '%s %d %s %d %d %s', 'Delimiter','\t', 'HeaderLines',0);
        fclose(fileID);
    end

    row = 1;
    while row <= length(data{1}) && isempty(strtrim(data{1}{row}))
        row = row + 1;
    end
    if row > length(data{1})
        warning('No valid data row in %s — skipping.', thisTxt);
        continue
    end

    csvFileName = strtrim(data{1}{row});
    matfile_fr  = double(data{2}(row));

    tparts    = strsplit(strtrim(data{3}{row}), ':');
    stim_time = str2double(tparts{1})*60 + str2double(tparts{2});

    prestim_s = double(data{4}(row));
    end_s     = double(data{5}(row));
    poststim  = end_s - stim_time;

    if ~isempty(data{6}) && row <= length(data{6}) && ~isempty(strtrim(data{6}{row}))
        vparts = strsplit(strtrim(data{6}{row}), ':');
        video_neural_offset = str2double(vparts{1})*60 + str2double(vparts{2});
    else
        video_neural_offset = 0;
    end

    fprintf('  stim=%.0fs  prestim=%.0fs  poststim=%.0fs  offset=%.0fs\n', ...
        stim_time, prestim_s, poststim, video_neural_offset);

    %% ----------------------------------------------------------------
    %%  STEP B — LOAD NEURAL TRACES, Z-SCORE, AND SMOOTH
    %% ----------------------------------------------------------------
    if ~exist(csvFileName,'file')
        warning('Neural CSV not found: %s — skipping.', csvFileName);
        continue
    end

    opts = detectImportOptions(csvFileName);
    opts.DataLines             = [3, Inf];                   % skip cell-ID row (1) and status row (2)
    opts.SelectedVariableNames = opts.VariableNames(2:end);  % drop time column
    raw_traces = table2array(readtable(csvFileName, opts))'; % [neurons x frames_full]

    start_idx_rec = round((stim_time - prestim_s) * matfile_fr) + 1;
    end_idx_rec   = min(round((stim_time + poststim) * matfile_fr), size(raw_traces,2));
    traces        = raw_traces(:, start_idx_rec:end_idx_rec);

    % Z-score using pre-stimulus baseline
    prestim_length = prestim_s * matfile_fr;
    prestim_mean   = mean(traces(:, 1:prestim_length), 2, 'omitnan');
    prestim_std    = std( traces(:, 1:prestim_length), 0, 2, 'omitnan');
    prestim_std(prestim_std < 1e-6) = 1e-6;
    trace_zs = (traces - prestim_mean) ./ prestim_std;  % [neurons x frames]

    % Gaussian smoothing — applied after z-scoring so baseline stats are unaffected.
    % A 2 s sigma (8 frames at 4 Hz) preserves bout-scale dynamics while
    % attenuating frame-to-frame noise that would otherwise inflate residuals.
    if smooth_sigma_s > 0
        sigma_frames = smooth_sigma_s * matfile_fr;
        trace_zs     = gausssmooth_rows(trace_zs, sigma_frames);
        fprintf('  Smoothed: sigma=%.1fs (%.0f frames)\n', smooth_sigma_s, sigma_frames);
    end

    [num_neurons, num_frames] = size(trace_zs);
    time_offset = stim_time - prestim_s;
    time_vec    = time_offset + (0:num_frames-1)' / matfile_fr;  % [frames x 1]

    fprintf('  Neurons: %d   Frames: %d\n', num_neurons, num_frames);

    %% ----------------------------------------------------------------
    %%  STEP C — DETECT LICKS AND LICK BOUTS FROM GPIO CSV
    %%  Rising-edge logic identical to FreeLicking_analysis_v3.m
    %% ----------------------------------------------------------------
    if ~exist(thisCSV,'file')
        warning('Lick CSV not found: %s — skipping.', thisCSV);
        continue
    end

    lickData  = readtable(thisCSV);
    gpio1data = lickData(strcmp(lickData.ChannelName, 'GPIO-1'), :);

    ttl_values  = gpio1data.Value;
    time_values = gpio1data.Time_s_;

    ttl_thr      = (max(ttl_values) + min(ttl_values)) / 2;
    rising_edges = find(diff(ttl_values > ttl_thr) == 1) + 1;
    lick_event_times = time_values(rising_edges);

    access_end = stim_time + poststim;
    lick_event_times = lick_event_times( ...
        lick_event_times >= stim_time & lick_event_times <= access_end);

    if isempty(lick_event_times)
        warning('No lick events found in %s — skipping.', thisCSV);
        continue
    end

    % Remove impossibly short inter-lick intervals
    filtered_licks = lick_event_times(1);
    for i = 2:length(lick_event_times)
        if (lick_event_times(i) - filtered_licks(end)) > min_ili
            filtered_licks = [filtered_licks; lick_event_times(i)]; %#ok<AGROW>
        end
    end

    % Group licks into bouts — NO minimum bout duration
    bouts        = {};
    current_bout = filtered_licks(1);
    for i = 2:length(filtered_licks)
        if (filtered_licks(i) - filtered_licks(i-1)) > inter_bout_interval
            if length(current_bout) >= min_licks_per_bout
                bouts{end+1} = current_bout; %#ok<AGROW>
            end
            current_bout = filtered_licks(i);
        else
            current_bout = [current_bout; filtered_licks(i)];
        end
    end
    if length(current_bout) >= min_licks_per_bout
        bouts{end+1} = current_bout;
    end

    if isempty(bouts)
        warning('No valid lick bouts in %s — skipping.', thisTxt);
        continue
    end

    total_bouts      = length(bouts);
    licks_per_bout   = cellfun(@length, bouts);
    bout_start_times = cellfun(@(x) x(1),   bouts);
    bout_end_times   = cellfun(@(x) x(end), bouts);

    fprintf('  Licks: %d   Valid bouts: %d   Mean licks/bout: %.1f\n', ...
        length(filtered_licks), total_bouts, mean(licks_per_bout));

    %% ----------------------------------------------------------------
    %%  STEP D — BUILD FRAME-ALIGNED PREDICTOR VECTORS
    %%
    %%  bout_active      : 1 during a valid lick bout, 0 during IBI
    %%  time_from_access : seconds elapsed since food access onset
    %%                     (0 during prestim, ramps linearly post-access)
    %%
    %%  time_from_access is normalised to [0,1] before regression so that
    %%  regression coefficients remain numerically comparable across sessions.
    %% ----------------------------------------------------------------
    bout_active      = false(num_frames, 1);
    time_from_access = zeros(num_frames, 1);

    for f = 1:num_frames
        t = time_vec(f);
        if any(t >= bout_start_times & t <= bout_end_times)
            bout_active(f) = true;
        end
        time_from_access(f) = max(0, t - stim_time);
    end

    % Restrict regression to post-stimulus frames only
    post_idx = time_vec >= stim_time;
    n_post   = sum(post_idx);

    bout_v      = double(bout_active(post_idx));
    time_v      = time_from_access(post_idx);
    time_v_norm = time_v / max(time_v);   % normalise to [0,1]

    fprintf('  Bout-active frames: %d / %d (%.1f%%)\n', ...
        sum(bout_v), n_post, 100*mean(bout_v));

    %% ----------------------------------------------------------------
    %%  STEP E — CUMULATIVE LICK CURVE (reference / plotting only)
    %% ----------------------------------------------------------------
    cumul_interp = zeros(1, length(cumul_time_axis));
    for b = 1:length(cumul_time_axis)
        cumul_interp(b) = sum(filtered_licks <= stim_time + cumul_time_axis(b));
    end
    allCumulCurves = [allCumulCurves; cumul_interp]; %#ok<AGROW>

    %% ----------------------------------------------------------------
    %%  STEP F — CLASSIFY NEURONS (act / inhib / none)
    %%  Mean z-score over first 600 s post-access:
    %%    > 1  → activated
    %%    < -1 → inhibited
    %%  Classification is done on smoothed traces (same signal as regression).
    %% ----------------------------------------------------------------
    pidx      = round(prestim_s * matfile_fr);
    idx_10min = min(round((prestim_s + 600) * matfile_fr), num_frames);
    idx_5min = min(round((prestim_s + 300) * matfile_fr), num_frames);

    act = []; inhib = []; none_idx = [];
    for v = 1:num_neurons
        a = mean(trace_zs(v, pidx:idx_10min), 'omitnan');
        b = mean(trace_zs(v, pidx:idx_5min), 'omitnan');

        if     a >  1 || b > 1; act      = [act      v]; %#ok<AGROW>
        elseif a < -1; inhib    = [inhib    v]; %#ok<AGROW>
        else;          none_idx = [none_idx v]; %#ok<AGROW>
        end
    end

    fprintf('  Act: %d  |  Inhib: %d  |  None: %d\n', ...
        numel(act), numel(inhib), numel(none_idx));

    allCounts.act   = allCounts.act   + numel(act);
    allCounts.inhib = allCounts.inhib + numel(inhib);
    allCounts.none  = allCounts.none  + numel(none_idx);

    %% ----------------------------------------------------------------
    %%  STEP G — PER-NEURON VARIANCE PARTITIONING + SHUFFLE CONTROL
    %%
    %%  Two predictors:
    %%    bout_v      — binary lick-bout ON/OFF
    %%                  answers: how much variance tracks consummatory behavior?
    %%    time_v_norm — normalised time-from-access [0,1]
    %%                  answers: how much variance tracks slow session-level change?
    %%
    %%  Design matrices:
    %%    X_full       = [1  bout_v  time_v_norm]
    %%    X_bout_only  = [1  bout_v]
    %%    X_time_only  = [1  time_v_norm]
    %%
    %%  Computed per neuron:
    %%    R2_bout     = individual OLS R2 for bout alone
    %%    R2_time     = individual OLS R2 for time alone
    %%    R2_full     = OLS R2 for both together
    %%    partR2_bout = R2_full - R2_time   (unique variance from bout)
    %%    partR2_time = R2_full - R2_bout   (unique variance from time)
    %%
    %%  Shuffle: bout_v and time_v_norm are each circularly shifted by a
    %%  random offset (min 30 s) to generate a null distribution.
    %%  Mean shuffle R2 per neuron is stored for paired comparison.
    %% ----------------------------------------------------------------
    X_full      = [ones(n_post,1), bout_v,  time_v_norm];
    X_bout_only = [ones(n_post,1), bout_v];
    X_time_only = [ones(n_post,1), time_v_norm];

    cat_indices = {act, inhib, none_idx};

    % Pre-generate shift offsets (shared across categories for this session)
    min_shift = round(30 * matfile_fr);
    max_shift = n_post - min_shift;
    if max_shift <= min_shift
        shift_offsets = randi(n_post, n_shuffles, 1);
    else
        shift_offsets = randi([min_shift, max_shift], n_shuffles, 1);
    end

    for ci = 1:3
        cat   = categories{ci};
        nrns  = cat_indices{ci};
        n_cat = numel(nrns);
        if n_cat == 0, continue; end

        r2b_c       = NaN(n_cat,1);
        r2t_c       = NaN(n_cat,1);
        r2f_c       = NaN(n_cat,1);
        pr2b_c      = NaN(n_cat,1);
        pr2t_c      = NaN(n_cat,1);
        r2b_shuf_c  = NaN(n_cat,1);
        r2t_shuf_c  = NaN(n_cat,1);

        for ni = 1:n_cat
            n = nrns(ni);
            y = trace_zs(n, post_idx)';
            if sum(isfinite(y)) < 20, continue; end

            r2b = compute_R2(X_bout_only, y);
            r2t = compute_R2(X_time_only, y);
            r2f = compute_R2(X_full,      y);

            r2b_c(ni)  = r2b;
            r2t_c(ni)  = r2t;
            r2f_c(ni)  = r2f;
            pr2b_c(ni) = max(0, r2f - r2t);
            pr2t_c(ni) = max(0, r2f - r2b);

            % Shuffle: circularly shift bout_v and time_v_norm independently
            shuf_r2b = NaN(n_shuffles,1);
            shuf_r2t = NaN(n_shuffles,1);
            for s = 1:n_shuffles
                bout_shuf = circshift(bout_v,      shift_offsets(s));
                time_shuf = circshift(time_v_norm, shift_offsets(s));
                shuf_r2b(s) = compute_R2([ones(n_post,1), bout_shuf], y);
                shuf_r2t(s) = compute_R2([ones(n_post,1), time_shuf], y);
            end
            r2b_shuf_c(ni) = mean(shuf_r2b, 'omitnan');
            r2t_shuf_c(ni) = mean(shuf_r2t, 'omitnan');
        end

        allR2_bout.(cat)      = [allR2_bout.(cat);      r2b_c];
        allR2_time.(cat)      = [allR2_time.(cat);      r2t_c];
        allR2_full.(cat)      = [allR2_full.(cat);      r2f_c];
        allPartR2_bout.(cat)  = [allPartR2_bout.(cat);  pr2b_c];
        allPartR2_time.(cat)  = [allPartR2_time.(cat);  pr2t_c];
        allR2_bout_shuf.(cat) = [allR2_bout_shuf.(cat); r2b_shuf_c];
        allR2_time_shuf.(cat) = [allR2_time_shuf.(cat); r2t_shuf_c];
    end

end  % file loop

%% ========== POPULATION-LEVEL STATISTICS ==========
anal_cats   = {'act','inhib'};
cat_labels  = {'ACTIVATED','INHIBITED'};
pred_labels = {'Bout (licking ON/OFF)', 'Time-from-access'};

fprintf('\n\n===== NEURON COUNTS =====\n');
fprintf('  Act: %d  |  Inhib: %d  |  None (excluded from stats): %d\n', ...
    allCounts.act, allCounts.inhib, allCounts.none);

r2_structs   = {allR2_bout,      allR2_time};
pr2_structs  = {allPartR2_bout,  allPartR2_time};
shuf_structs = {allR2_bout_shuf, allR2_time_shuf};

for ci = 1:2
    cat = anal_cats{ci};
    fprintf('\n----- %s (n=%d neurons) -----\n', cat_labels{ci}, numel(allR2_full.(cat)));

    fprintf('  Individual R2 and unique partial R2:\n');
    for k = 1:2
        vals  = r2_structs{k}.(cat);  vals  = vals(isfinite(vals));
        pvals = pr2_structs{k}.(cat); pvals = pvals(isfinite(pvals));
        fprintf('    %-25s  R2=%.4f (SD=%.4f)  Partial R2=%.4f (SD=%.4f)\n', ...
            pred_labels{k}, mean(vals), std(vals), mean(pvals), std(pvals));
    end

    rf = allR2_full.(cat); rf = rf(isfinite(rf));
    fprintf('    %-25s  R2=%.4f (SD=%.4f)\n', 'Full model (bout+time)', mean(rf), std(rf));

    % Paired t-tests: are bout and time predictors different from each other?
    fprintf('  Paired t-tests (bout vs time, per neuron):\n');
    va = allR2_bout.(cat);     vb = allR2_time.(cat);
    ok = isfinite(va) & isfinite(vb);
    if sum(ok) > 2
        [~,pv,~,st] = ttest(va(ok), vb(ok));
        fprintf('    Individual R2:   t(%.0f)=%.3f  p=%.4f\n', st.df, st.tstat, pv);
    end
    pa = allPartR2_bout.(cat); pb = allPartR2_time.(cat);
    ok2 = isfinite(pa) & isfinite(pb);
    if sum(ok2) > 2
        [~,pv2,~,st2] = ttest(pa(ok2), pb(ok2));
        fprintf('    Partial R2:      t(%.0f)=%.3f  p=%.4f\n', st2.df, st2.tstat, pv2);
    end

    % Shuffle comparison: observed R2 vs mean circular-shift shuffle R2
    fprintf('  Shuffle control (circular shift, n=%d):\n', n_shuffles);
    for k = 1:2
        obs  = r2_structs{k}.(cat);
        shuf = shuf_structs{k}.(cat);
        ok_s = isfinite(obs) & isfinite(shuf);
        if sum(ok_s) > 2
            [~,ps,~,sts] = ttest(obs(ok_s), shuf(ok_s), 'Tail','right');
            fprintf('    %-25s  Observed=%.4f  Shuffle=%.4f  t(%.0f)=%.3f  p=%.4f\n', ...
                pred_labels{k}, mean(obs(ok_s)), mean(shuf(ok_s)), sts.df, sts.tstat, ps);
        end
    end
end

fprintf('\n===== CUMULATIVE LICKS — reference (n=%d animals) =====\n', size(allCumulCurves,1));
fprintf('  Mean total licks at %.0f min: %.1f ± %.1f SEM\n', ...
    cumul_plot_dur_s/60, ...
    mean(allCumulCurves(:,end)), std(allCumulCurves(:,end))/sqrt(size(allCumulCurves,1)));

%% ========== FIGURES ==========
if graph_data

    %% Figure 1 — Individual R2: bout vs time-from-access, act vs inhib
    figure('Name','R2: Bout vs Time-from-Access','Position',[50 100 700 500]);
    hold on;
    bw      = 0.3;
    offsets = [-0.2, 0.2];
    leg_h   = gobjects(2,1);
    for ci = 1:2
        cat = anal_cats{ci};
        col = cat_colors(ci,:);
        r2s = {allR2_bout.(cat), allR2_time.(cat)};
        for k = 1:2
            vals = r2s{k}; vals = vals(isfinite(vals));
            if isempty(vals), continue; end
            xp = k + offsets(ci);
            b  = bar(xp, mean(vals), bw, 'FaceColor',col, 'EdgeColor','k','FaceAlpha',0.75);
            if k == 1, leg_h(ci) = b; end
            errorbar(xp, mean(vals), std(vals)/sqrt(numel(vals)), ...
                'k','LineWidth',1.5,'CapSize',6,'LineStyle','none');
            jit = (rand(size(vals))-0.5)*bw*0.7;
            scatter(xp+jit, vals, 14, col*0.6, 'filled','MarkerFaceAlpha',0.4);
        end
    end
    xlim([0.5 2.5]); xticks(1:2); xticklabels(pred_labels);
    ylabel('R^2');
    title('Variance Explained: Lick-Bout State vs Time-from-Access');
    legend(leg_h, {'Activated','Inhibited'}, 'Location','northeast');
    box on;

    %% Figure 2 — Partial R2: unique variance per predictor
    figure('Name','Partial R2: Unique Variance','Position',[200 100 700 500]);
    hold on;
    pr2s = {allPartR2_bout, allPartR2_time};
    for ci = 1:2
        cat = anal_cats{ci};
        col = cat_colors(ci,:);
        for k = 1:2
            vals = pr2s{k}.(cat); vals = vals(isfinite(vals));
            if isempty(vals), continue; end
            xp = k + offsets(ci);
            bar(xp, mean(vals), bw, 'FaceColor',col, 'EdgeColor','k','FaceAlpha',0.75);
            errorbar(xp, mean(vals), std(vals)/sqrt(numel(vals)), ...
                'k','LineWidth',1.5,'CapSize',6,'LineStyle','none');
            jit = (rand(size(vals))-0.5)*bw*0.7;
            scatter(xp+jit, vals, 14, col*0.6, 'filled','MarkerFaceAlpha',0.4);
        end
    end
    xlim([0.5 2.5]); xticks(1:2); xticklabels(pred_labels);
    ylabel('Partial R^2 (unique)');
    title('Unique Variance: Lick-Bout State vs Time-from-Access');
    legend(leg_h, {'Activated','Inhibited'}, 'Location','northeast');
    box on;

    %% Figure 3 — Per-neuron partial R2 scatter: bout vs time
    figure('Name','Per-neuron Partial R2 Scatter','Position',[350 100 520 500]);
    hold on;
    all_pts = [];
    for ci = 1:2
        cat = anal_cats{ci};
        col = cat_colors(ci,:);
        pb  = allPartR2_bout.(cat);
        pt  = allPartR2_time.(cat);
        ok  = isfinite(pb) & isfinite(pt);
        scatter(pb(ok), pt(ok), 25, col, 'filled','MarkerFaceAlpha',0.45, ...
            'DisplayName', sprintf('%s (n=%d)', cat_labels{ci}, sum(ok)));
        all_pts = [all_pts; pb(ok), pt(ok)]; %#ok<AGROW>
    end
    if ~isempty(all_pts)
        ax_max = max(all_pts(:)) * 1.1;
        plot([0 ax_max],[0 ax_max],'k--','LineWidth',1.2);
    end
    xlabel('Partial R^2 (Lick-Bout)');
    ylabel('Partial R^2 (Time-from-Access)');
    title('Unique Variance per Neuron: Bout vs Time');
    axis equal; box on; legend('Location','northwest');

    %% Figure 4 — Shuffle control: observed vs shuffle scatter per neuron
    shuf_structs_fig = {allR2_bout_shuf, allR2_time_shuf};
    figure('Name','Shuffle Control','Position',[500 100 950 450]);
    sp = 0;
    for ci = 1:2
        cat = anal_cats{ci};
        col = cat_colors(ci,:);
        for k = 1:2
            sp = sp + 1;
            subplot(2,2,sp); hold on;
            obs  = r2_structs{k}.(cat);
            shuf = shuf_structs_fig{k}.(cat);
            ok_s = isfinite(obs) & isfinite(shuf);
            scatter(shuf(ok_s), obs(ok_s), 18, col, 'filled','MarkerFaceAlpha',0.4);
            ax_max = max([shuf(ok_s); obs(ok_s); 0.01]) * 1.1;
            plot([0 ax_max],[0 ax_max],'k--','LineWidth',1.2);
            xlabel('Shuffle R^2 (mean)'); ylabel('Observed R^2');
            title(sprintf('%s — %s', cat_labels{ci}, pred_labels{k}));
            frac = mean(obs(ok_s) > shuf(ok_s));
            text(ax_max*0.05, ax_max*0.9, sprintf('%.0f%% above shuffle', frac*100), ...
                'FontSize',8,'Color',col*0.7);
            axis([0 ax_max 0 ax_max]); axis square; box on;
        end
    end

    %% Figure 5 — Observed vs shuffle R2: bar comparison (mean ± SEM)
    figure('Name','Observed vs Shuffle R2','Position',[550 150 800 480]);
    bw_s = 0.35;
    for ci = 1:2
        cat = anal_cats{ci};
        col = cat_colors(ci,:);
        subplot(1,2,ci); hold on;
        obs_means  = [mean(allR2_bout.(cat),'omitnan'),      mean(allR2_time.(cat),'omitnan')];
        obs_sems   = [std(allR2_bout.(cat),'omitnan')/sqrt(sum(isfinite(allR2_bout.(cat)))), ...
                      std(allR2_time.(cat),'omitnan')/sqrt(sum(isfinite(allR2_time.(cat))))];
        shuf_means = [mean(allR2_bout_shuf.(cat),'omitnan'), mean(allR2_time_shuf.(cat),'omitnan')];
        shuf_sems  = [std(allR2_bout_shuf.(cat),'omitnan')/sqrt(sum(isfinite(allR2_bout_shuf.(cat)))), ...
                      std(allR2_time_shuf.(cat),'omitnan')/sqrt(sum(isfinite(allR2_time_shuf.(cat))))];
        for k = 1:2
            bar(k-0.2, obs_means(k),  bw_s, 'FaceColor', col,         'EdgeColor','k','FaceAlpha',0.8);
            bar(k+0.2, shuf_means(k), bw_s, 'FaceColor', col*0.4+0.6, 'EdgeColor','k','FaceAlpha',0.8,'LineStyle','--');
            errorbar(k-0.2, obs_means(k),  obs_sems(k),  'k','LineWidth',1.5,'CapSize',6,'LineStyle','none');
            errorbar(k+0.2, shuf_means(k), shuf_sems(k), 'k','LineWidth',1.5,'CapSize',6,'LineStyle','none');
        end
        xlim([0.5 2.5]); xticks(1:2); xticklabels(pred_labels);
        ylabel('R^2'); title(cat_labels{ci});
        b1 = bar(nan, nan, 'FaceColor', col,         'EdgeColor','k','FaceAlpha',0.8);
        b2 = bar(nan, nan, 'FaceColor', col*0.4+0.6, 'EdgeColor','k','FaceAlpha',0.8,'LineStyle','--');
        legend([b1 b2], {'Observed','Shuffle'}, 'Location','northeast');
        yline(0,'k-'); box on;
    end

    %% Figure 6 — Cumulative lick curves (reference / not a predictor)
    figure('Name','Cumulative Licks — Group Reference','Position',[500 100 750 500]);
    hold on;
    n_anim   = size(allCumulCurves,1);
    grp_mean = mean(allCumulCurves, 1, 'omitnan');
    grp_sem  = std(allCumulCurves,  0, 1, 'omitnan') / sqrt(n_anim);
    t_min    = cumul_time_axis / 60;

    for a = 1:n_anim
        plot(t_min, allCumulCurves(a,:), '-', ...
            'Color',[0.4 0.4 0.4 0.35], 'LineWidth',0.8);
    end
    fill([t_min; flipud(t_min)], ...
         [grp_mean - grp_sem, fliplr(grp_mean + grp_sem)]', ...
         [0.2 0.5 0.8], 'FaceAlpha',0.25, 'EdgeColor','none');
    plot(t_min, grp_mean, '-', 'Color',[0.2 0.5 0.8], 'LineWidth',2.5);

    xlabel('Time from food access (min)');
    ylabel('Cumulative licks');
    title(sprintf('Cumulative Lick Intake — reference (n=%d animals)', n_anim));
    xlim([0 cumul_plot_dur_s/60]); box on;
    text(cumul_plot_dur_s/60*0.98, grp_mean(end), ...
        sprintf('  %.0f ± %.0f', grp_mean(end), grp_sem(end)), ...
        'HorizontalAlignment','left', 'FontSize',9, 'Color',[0.2 0.5 0.8]);

end

%% ========== SAVE RESULTS ==========
save('lick_variance_results.mat', ...
    'allR2_bout','allR2_time','allR2_full', ...
    'allPartR2_bout','allPartR2_time', ...
    'allR2_bout_shuf','allR2_time_shuf', ...
    'allCumulCurves','cumul_time_axis', ...
    'allCounts','categories','cat_colors', ...
    'smooth_sigma_s','inter_bout_interval','min_licks_per_bout','n_shuffles');
fprintf('\nResults saved to lick_variance_results.mat\n');

%% ========== LOCAL FUNCTIONS ==========

function out = gausssmooth_rows(M, sigma_frames)
% Gaussian-smooth each row of M along columns (time axis).
% sigma_frames : SD of Gaussian kernel in frames.
% Zero-padded by 3*sigma on each side to reduce edge effects.
    half_w = ceil(3 * sigma_frames);
    x      = -half_w : half_w;
    kernel = exp(-x.^2 / (2 * sigma_frames^2));
    kernel = kernel / sum(kernel);
    out    = zeros(size(M));
    pad    = zeros(size(M,1), half_w);
    M_pad  = [pad, M, pad];
    for n = 1:size(M,1)
        smoothed  = conv(M_pad(n,:), kernel, 'valid');
        out(n,:)  = smoothed(1:size(M,2));
    end
end

function r2 = compute_R2(X, y)
% OLS R2 = 1 - SSres/SStot
    valid = isfinite(y);
    X = X(valid,:); y = y(valid);
    if size(X,1) < size(X,2) + 2
        r2 = NaN; return
    end
    b      = X \ y;
    yhat   = X * b;
    SS_res = sum((y - yhat).^2);
    SS_tot = sum((y - mean(y)).^2);
    if SS_tot < 1e-10
        r2 = NaN; return
    end
    r2 = max(0, 1 - SS_res/SS_tot);
end