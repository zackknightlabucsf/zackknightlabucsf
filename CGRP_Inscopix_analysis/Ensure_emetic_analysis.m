%% Ensure + LiCl Analysis v3
% Classifies neurons as Ensure / LiCl / Both / None across all files in
% EnsureLiCl_combo.txt.
%
% Ensure (stim 1): first filter only – mean z-score > elevated_threshold
%                  over Int1 seconds post-stimulus.
%
% LiCl   (stim 2): two-filter approach (adapted from IG_sort_activationTimes):
%   Filter 1 – mean z-score > elevated_threshold over Int2 seconds post-stim.
%   Filter 2 – sliding-window check: within any window_size_licl-second
%              window after the first threshold crossing, the neuron must
%              spend >= min_duration_licl seconds above threshold.         
% Average z-scores are stored per neuron, per stimulus, per category for
% downstream group comparisons.

clear; close all;

%% ── Input Parameters ────────────────────────────────────────────────────
prestim1   = 600;   % Pre-stim baseline duration, Ensure  (s)
prestim2   = 480;   % Pre-stim baseline duration, LiCl/DON   (s)
poststim1  = 3000;  % Recording window after Ensure stim =3000, DON = 1800  (s)
poststim2  = 1750;  % Recording window after LiCl stim = 1750, DON = 300 (s)
Int1       = 600;   % Post-stim averaging window, Ensure, LiCl = 600, DON = 300 (s)
Int2       = 900;   % Post-stim averaging window, LiCl = 900, DON = 300   (s)
Int3       = 900;   % Post-stim averaging window, Ensure  (s)
piechart   = 1;     % 1 = generate pie charts

% Shared threshold
elevated_threshold = 1;  % z-score

% ── LiCl sliding-window parameters (from IG_sort_activationTimes) ────────
window_size_licl   = 180;   % Sliding window length (s)
min_duration_licl  = 108;   % Min seconds above threshold within that window (s)
                            % (42 s ≈ 70 % of a 60-s window at 4 Hz)

%% ── Load File List ───────────────────────────────────────────────────────
txt_file = 'IGEnsure_LiCl_combo.txt';
fileID   = fopen(txt_file, 'r');
data     = textscan(fileID, '%s %d %s %s %d', 'Delimiter', '\t');
fclose(fileID);

csv_files   = data{1};
frequencies = data{2};
stim_times1 = data{3};
stim_times2 = data{4};

%% ── Initialise Cumulative Storage ────────────────────────────────────────
% Traces aligned to Ensure stim, split by category
cumulative_traces_both   = [];
cumulative_traces_ensure = [];
cumulative_traces_licl   = [];
cumulative_traces_none   = [];

% Traces aligned to LiCl stim, split by category
cumulative_traces_both_licl   = [];
cumulative_traces_ensure_licl = [];
cumulative_traces_licl_licl   = [];
cumulative_traces_none_licl   = [];

cumulative_traces_all = [];

% Per-neuron average z-scores for each stimulus × category combination.
% Each row = one neuron; columns depend on the category vector populated
% inside the loop. Collected here as flat arrays for group-level stats.
cumAvgZ = struct( ...
    'ensure_stim', struct('both',[],'ensure',[],'licl',[],'none',[]), ...
    'licl_stim',   struct('both',[],'ensure',[],'licl',[],'none',[]));

% Cumulative means (legacy – kept for backwards compatibility)
cumulative_means = struct( ...
    'ensure', struct('both',[],'ensure',[],'licl',[],'none',[]), ...
    'licl',   struct('both',[],'ensure',[],'licl',[],'none',[]));

results         = struct();
results_indices = struct();

%% ── Per-File Loop ────────────────────────────────────────────────────────
for idx = 1:length(csv_files)

    csv_file    = csv_files{idx};
    matfile_fr  = double(frequencies(idx));   % Hz

    % Parse MM:SS stim times → seconds
    p1 = textscan(stim_times1{idx}, '%d:%d');
    stim_time1 = double(p1{1})*60 + double(p1{2});

    p2 = textscan(stim_times2{idx}, '%d:%d');
    stim_time2 = double(p2{1})*60 + double(p2{2});

    % ── Load CSV ─────────────────────────────────────────────────────────
    opts = detectImportOptions(csv_file);
    opts.DataLines = [3, Inf];
    opts.SelectedVariableNames = opts.VariableNames(2:end);
    raw_traces = table2array(readtable(csv_file, opts))';   % [neurons × samples]
    num_traces = size(raw_traces, 1);

    % ── Z-score relative to Ensure baseline ──────────────────────────────
    start_idx1     = (stim_time1 - prestim1) * matfile_fr + 1;
    end_idx1       = (stim_time1 + poststim1) * matfile_fr - 1;
    traces1        = raw_traces(:, start_idx1:end_idx1);
    pre1_samp      = prestim1 * matfile_fr;
    prestim_mean1  = mean(traces1(:, 1:pre1_samp), 2);
    prestim_std1   = std( traces1(:, 1:pre1_samp), 0, 2);
    traces_zs_Ensure = (traces1 - prestim_mean1) ./ prestim_std1;

    % ── Z-score relative to LiCl baseline ────────────────────────────────
    start_idx2    = (stim_time2 - prestim2) * matfile_fr;
    end_idx2      = (stim_time2 + poststim2) * matfile_fr - 1;
    traces2       = raw_traces(:, start_idx2:end_idx2);
    pre2_samp     = prestim2 * matfile_fr;
    prestim_mean2 = mean(traces2(:, 1:pre2_samp), 2);
    prestim_std2  = std( traces2(:, 1:pre2_samp), 0, 2);
    traces_zs_LiCl = (traces2 - prestim_mean2) ./ prestim_std2;

    time_step = 1 / matfile_fr;   % seconds per sample

    % ════════════════════════════════════════════════════════════════════
    %  FILTER 1 – mean z-score screen (applied to BOTH stimuli)
    % ════════════════════════════════════════════════════════════════════
    ensure_post_start = pre1_samp + 1;
    ensure_post_end1   = pre1_samp + Int1 * matfile_fr;
    ensure_post_end2   = pre1_samp + Int3 * matfile_fr;

    licl_post_start   = pre2_samp;
    licl_post_end     = pre2_samp + Int2 * matfile_fr;

    ensure_pass_f1 = false(num_traces, 1);
    licl_pass_f1   = false(num_traces, 1);

    for i = 1:num_traces
        a = mean(traces_zs_Ensure(i, ensure_post_start:ensure_post_end1));
        b = mean(traces_zs_LiCl(  i, licl_post_start :licl_post_end  ));
        c = mean(traces_zs_Ensure(i, ensure_post_start:ensure_post_end2));
        ensure_pass_f1(i) = (a > elevated_threshold || c > elevated_threshold);
        licl_pass_f1(i)   = (b > elevated_threshold);
    end

    % ════════════════════════════════════════════════════════════════════
    %  FILTER 2 – sliding window (LiCl only)
    %  Applied only to neurons that passed filter 1 for LiCl.
    %  Mirrors the logic in IG_sort_activationTimes_v6_test.m.
    % ════════════════════════════════════════════════════════════════════
    win_samp_licl = window_size_licl * matfile_fr;  % samples in window

    licl_pass_f2        = false(num_traces, 1);
    first_licl_act_time = zeros(num_traces, 1);   % seconds after LiCl stim onset

    for i = 1:num_traces
        if ~licl_pass_f1(i)
            continue   % only run filter 2 on filter-1 candidates
        end

        post_trace      = traces_zs_LiCl(i, licl_post_start:end);
        total_post      = length(post_trace);
        above           = post_trace > elevated_threshold;

        % Rising-edge crossings (below → above threshold)
        crosses = [above(1), diff(above) == 1];
        cross_idx = find(crosses);

        for c_idx = cross_idx
            win_end  = min(c_idx + win_samp_licl - 1, total_post);
            act_secs = sum(above(c_idx:win_end)) * time_step;

            if act_secs >= min_duration_licl
                licl_pass_f2(i)        = true;
                first_licl_act_time(i) = c_idx * time_step;   % s after LiCl onset
                break
            end
        end
    end

    % ── Final LiCl classification: must pass both filters ────────────────
    licl_responds = licl_pass_f1 & licl_pass_f2;

    % ════════════════════════════════════════════════════════════════════
    %  ASSIGN CATEGORIES
    % ════════════════════════════════════════════════════════════════════
    both   = find( ensure_pass_f1 &  licl_responds)';
    ensure = find( ensure_pass_f1 & ~licl_responds)';
    licl   = find(~ensure_pass_f1 &  licl_responds)';
    none   = find(~ensure_pass_f1 & ~licl_responds)';

    % ════════════════════════════════════════════════════════════════════
    %  PER-NEURON AVERAGE Z-SCORES
    %  For each neuron, average z-score is computed over Int1 (Ensure)
    %  and Int2 (LiCl) post-stim windows and stored by category.
    % ════════════════════════════════════════════════════════════════════

    % Helper: average z over a fixed post-stim window for an index set
    avgZ_ensure = @(idx_set) mean(traces_zs_Ensure(idx_set, ...
                      ensure_post_start:ensure_post_end1), 2);   % [n×1]
    avgZ_licl   = @(idx_set) mean(traces_zs_LiCl(  idx_set, ...
                      licl_post_start  :licl_post_end  ), 2);   % [n×1]

    % Per-neuron means, per stimulus, per category
    pnAvgZ_ensure_both   = avgZ_ensure(both);
    pnAvgZ_ensure_ensure = avgZ_ensure(ensure);
    pnAvgZ_ensure_licl   = avgZ_ensure(licl);
    pnAvgZ_ensure_none   = avgZ_ensure(none);

    pnAvgZ_licl_both   = avgZ_licl(both);
    pnAvgZ_licl_ensure = avgZ_licl(ensure);
    pnAvgZ_licl_licl   = avgZ_licl(licl);
    pnAvgZ_licl_none   = avgZ_licl(none);

    % Accumulate across files
    cumAvgZ.ensure_stim.both   = [cumAvgZ.ensure_stim.both;   pnAvgZ_ensure_both];
    cumAvgZ.ensure_stim.ensure = [cumAvgZ.ensure_stim.ensure; pnAvgZ_ensure_ensure];
    cumAvgZ.ensure_stim.licl   = [cumAvgZ.ensure_stim.licl;   pnAvgZ_ensure_licl];
    cumAvgZ.ensure_stim.none   = [cumAvgZ.ensure_stim.none;   pnAvgZ_ensure_none];

    cumAvgZ.licl_stim.both   = [cumAvgZ.licl_stim.both;   pnAvgZ_licl_both];
    cumAvgZ.licl_stim.ensure = [cumAvgZ.licl_stim.ensure; pnAvgZ_licl_ensure];
    cumAvgZ.licl_stim.licl   = [cumAvgZ.licl_stim.licl;   pnAvgZ_licl_licl];
    cumAvgZ.licl_stim.none   = [cumAvgZ.licl_stim.none;   pnAvgZ_licl_none];

    % Legacy cumulative_means (row vector per file, for backwards compat.)
    cumulative_means.ensure.both   = [cumulative_means.ensure.both;   pnAvgZ_ensure_both];
    cumulative_means.ensure.ensure = [cumulative_means.ensure.ensure; pnAvgZ_ensure_ensure];
    cumulative_means.ensure.licl   = [cumulative_means.ensure.licl;   pnAvgZ_ensure_licl];
    cumulative_means.ensure.none   = [cumulative_means.ensure.none;   pnAvgZ_ensure_none];

    cumulative_means.licl.both   = [cumulative_means.licl.both;   pnAvgZ_licl_both];
    cumulative_means.licl.ensure = [cumulative_means.licl.ensure; pnAvgZ_licl_ensure];
    cumulative_means.licl.licl   = [cumulative_means.licl.licl;   pnAvgZ_licl_licl];
    cumulative_means.licl.none   = [cumulative_means.licl.none;   pnAvgZ_licl_none];

    % ════════════════════════════════════════════════════════════════════
    %  STORE TRACES BY CATEGORY
    % ════════════════════════════════════════════════════════════════════
    Neurons_both   = traces_zs_Ensure(both,   :);
    Neurons_ensure = traces_zs_Ensure(ensure, :);
    Neurons_licl   = traces_zs_Ensure(licl,   :);
    Neurons_none   = traces_zs_Ensure(none,   :);

    Neuron_both_licl   = traces_zs_LiCl(both,   :);
    Neuron_ensure_licl = traces_zs_LiCl(ensure, :);
    Neuron_licl_licl   = traces_zs_LiCl(licl,   :);
    Neuron_none_licl   = traces_zs_LiCl(none,   :);

    cumulative_traces_both         = [cumulative_traces_both;         Neurons_both];
    cumulative_traces_ensure       = [cumulative_traces_ensure;       Neurons_ensure];
    cumulative_traces_licl         = [cumulative_traces_licl;         Neurons_licl];
    cumulative_traces_none         = [cumulative_traces_none;         Neurons_none];
    cumulative_traces_both_licl    = [cumulative_traces_both_licl;    Neuron_both_licl];
    cumulative_traces_ensure_licl  = [cumulative_traces_ensure_licl;  Neuron_ensure_licl];
    cumulative_traces_licl_licl    = [cumulative_traces_licl_licl;    Neuron_licl_licl];
    cumulative_traces_none_licl    = [cumulative_traces_none_licl;    Neuron_none_licl];
    cumulative_traces_all          = [cumulative_traces_all;          traces_zs_Ensure];

    % ── Per-file results struct ───────────────────────────────────────────
    results_indices(idx).csv_file    = csv_file;
    results_indices(idx).num_neurons = num_traces;
    results_indices(idx).both        = both;
    results_indices(idx).ensure      = ensure;
    results_indices(idx).licl        = licl;
    results_indices(idx).none        = none;

    results(idx).csv_file          = csv_file;
    results(idx).Neurons_both      = Neurons_both;
    results(idx).Neurons_ensure    = Neurons_ensure;
    results(idx).Neurons_licl      = Neurons_licl;
    results(idx).Neurons_none      = Neurons_none;
    results(idx).first_licl_act_time = first_licl_act_time; % s after LiCl stim; 0 = no sustained window

    % Per-neuron average z-scores stored in results for this file
    results(idx).avgZ_ensure_stim.both   = pnAvgZ_ensure_both;
    results(idx).avgZ_ensure_stim.ensure = pnAvgZ_ensure_ensure;
    results(idx).avgZ_ensure_stim.licl   = pnAvgZ_ensure_licl;
    results(idx).avgZ_ensure_stim.none   = pnAvgZ_ensure_none;
    results(idx).avgZ_licl_stim.both     = pnAvgZ_licl_both;
    results(idx).avgZ_licl_stim.ensure   = pnAvgZ_licl_ensure;
    results(idx).avgZ_licl_stim.licl     = pnAvgZ_licl_licl;
    results(idx).avgZ_licl_stim.none     = pnAvgZ_licl_none;

    % ── Per-animal printout ───────────────────────────────────────────────
    n_b = numel(both);   n_e = numel(ensure);
    n_l = numel(licl);   n_nn = numel(none);
    fprintf('\n────────────────────────────────────────────────\n');
    fprintf('File %d / %d : %s\n', idx, length(csv_files), csv_file);
    fprintf('────────────────────────────────────────────────\n');
    fprintf('  %-14s  %5s  %6s\n', 'Category', 'n', '%');
    fprintf('  %-14s  %5d  %5.1f%%\n', 'Both',        n_b,  100*n_b /num_traces);
    fprintf('  %-14s  %5d  %5.1f%%\n', 'Ensure only', n_e,  100*n_e /num_traces);
    fprintf('  %-14s  %5d  %5.1f%%\n', 'LiCl only',   n_l,  100*n_l /num_traces);
    fprintf('  %-14s  %5d  %5.1f%%\n', 'None',         n_nn, 100*n_nn/num_traces);
    fprintf('  %-14s  %5d\n',          'Total',         num_traces);
end

%% ── Consolidated Summary ─────────────────────────────────────────────────
n_both   = size(cumulative_traces_both,   1);
n_ensure = size(cumulative_traces_ensure, 1);
n_licl   = size(cumulative_traces_licl,   1);
n_none   = size(cumulative_traces_none,   1);
n_total  = n_both + n_ensure + n_licl + n_none;

fprintf('\n════════════════════════════════════════\n');
fprintf('CONSOLIDATED SUMMARY  (all files)\n');
fprintf('════════════════════════════════════════\n');
fprintf('  Both    : %d  (%.1f%%)\n', n_both,   100*n_both/n_total);
fprintf('  Ensure  : %d  (%.1f%%)\n', n_ensure, 100*n_ensure/n_total);
fprintf('  LiCl    : %d  (%.1f%%)\n', n_licl,   100*n_licl/n_total);
fprintf('  None    : %d  (%.1f%%)\n', n_none,   100*n_none/n_total);
fprintf('  Total   : %d\n', n_total);

%% ── Average Z-Score Summary (across all files, per category) ─────────────
fprintf('\n── Mean avg-z per category ────────────────────────────────────\n');
fprintf('                 Ensure stim         LiCl stim\n');
fprintf('  Both    : %8.3f ± %.3f    %8.3f ± %.3f\n', ...
    mean(cumAvgZ.ensure_stim.both),   std(cumAvgZ.ensure_stim.both),   ...
    mean(cumAvgZ.licl_stim.both),     std(cumAvgZ.licl_stim.both));
fprintf('  Ensure  : %8.3f ± %.3f    %8.3f ± %.3f\n', ...
    mean(cumAvgZ.ensure_stim.ensure), std(cumAvgZ.ensure_stim.ensure), ...
    mean(cumAvgZ.licl_stim.ensure),   std(cumAvgZ.licl_stim.ensure));
fprintf('  LiCl    : %8.3f ± %.3f    %8.3f ± %.3f\n', ...
    mean(cumAvgZ.ensure_stim.licl),   std(cumAvgZ.ensure_stim.licl),   ...
    mean(cumAvgZ.licl_stim.licl),     std(cumAvgZ.licl_stim.licl));
fprintf('  None    : %8.3f ± %.3f    %8.3f ± %.3f\n', ...
    mean(cumAvgZ.ensure_stim.none),   std(cumAvgZ.ensure_stim.none),   ...
    mean(cumAvgZ.licl_stim.none),     std(cumAvgZ.licl_stim.none));

%% ── Pie Chart ────────────────────────────────────────────────────────────
if piechart == 1
    X      = [n_both, n_ensure, n_licl, n_none];
    labels = {'Both', 'Ensure Only', 'LiCl Only', 'None'};
    figure;
    pie(X);
    colormap([1 0 1; 1 1 0; 0 1 0; 0.85 0.85 0.85]);
    legend(labels, 'Location', 'southoutside', 'Orientation', 'horizontal');
    title('Consolidated Classification Across All Files');
end

%% ── Save ─────────────────────────────────────────────────────────────────
save('Cumulative_Traces_v4.mat', ...
    'cumulative_traces_both',   'cumulative_traces_ensure', ...
    'cumulative_traces_licl',   'cumulative_traces_none', ...
    'cumulative_traces_both_licl',   'cumulative_traces_ensure_licl', ...
    'cumulative_traces_licl_licl',   'cumulative_traces_none_licl');

save('Consolidated_Results_v4.mat', ...
    'results', 'results_indices', 'cumulative_means', 'cumAvgZ', ...
    'cumulative_traces_both',   'cumulative_traces_ensure', ...
    'cumulative_traces_licl',   'cumulative_traces_none');


%% ── Outlier Diagnostic: find neurons with extreme LiCl z-scores ──────────
% Identifies which file and neuron row produced high z-score, and
% reports the likely cause (near-zero baseline std, artifact spike, etc.)
 
z_outlier_thresh = 30;   % flag any neuron whose peak |z| exceeds this
 
fprintf('\n════════════════════════════════════════════════════════════\n');
fprintf('OUTLIER DIAGNOSTIC  (|peak z| > %g in LiCl-aligned traces)\n', z_outlier_thresh);
fprintf('════════════════════════════════════════════════════════════\n');
 
any_found = false;
 
for idx = 1:length(csv_files)
    csv_file   = csv_files{idx};
    matfile_fr = double(frequencies(idx));
 
    % Re-derive LiCl baseline indices (same logic as main loop)
    p2 = textscan(stim_times2{idx}, '%d:%d');
    stim_time2 = double(p2{1})*60 + double(p2{2});
 
    start_idx2 = (stim_time2 - prestim2) * matfile_fr;
    end_idx2   = (stim_time2 + poststim2) * matfile_fr - 1;
 
    opts = detectImportOptions(csv_file);
    opts.DataLines = [3, Inf];
    opts.SelectedVariableNames = opts.VariableNames(2:end);
    raw_traces = table2array(readtable(csv_file, opts))';
 
    traces2       = raw_traces(:, start_idx2:end_idx2);
    pre2_samp     = prestim2 * matfile_fr;
    prestim_mean2 = mean(traces2(:, 1:pre2_samp), 2);
    prestim_std2  = std( traces2(:, 1:pre2_samp), 0, 2);
    traces_zs_LiCl = (traces2 - prestim_mean2) ./ prestim_std2;
 
    % Find neurons with extreme peak z-score
    peak_z = max(abs(traces_zs_LiCl), [], 2);
    bad_neurons = find(peak_z > z_outlier_thresh);
 
    for i = 1:numel(bad_neurons)
        n = bad_neurons(i);
        any_found = true;
 
        bl_mean = prestim_mean2(n);
        bl_std  = prestim_std2(n);
        pk_z    = peak_z(n);
        raw_bl  = traces2(n, 1:pre2_samp);
        raw_max = max(abs(traces2(n, :)));
 
        % Guess at root cause
        if bl_std < 0.01
            cause = 'NEAR-ZERO BASELINE STD (flat/dead baseline)';
        elseif bl_std < 0.1
            cause = 'Very low baseline std (nearly flat trace)';
        elseif max(abs(raw_bl)) > 10 * median(abs(raw_bl))
            cause = 'Artifact spike IN BASELINE inflating std asymmetrically';
        elseif raw_max > 5 * (bl_std * z_outlier_thresh + abs(bl_mean))
            cause = 'Massive post-stim spike in raw trace';
        else
            cause = 'Unknown — inspect raw trace manually';
        end
 
        fprintf('\nFile %d: %s\n', idx, csv_file);
        fprintf('  Neuron row      : %d\n',   n);
        fprintf('  Peak |z|        : %.1f\n', pk_z);
        fprintf('  Baseline mean   : %.4f\n', bl_mean);
        fprintf('  Baseline std    : %.6f\n', bl_std);
        fprintf('  Raw trace max   : %.4f\n', raw_max);
        fprintf('  Likely cause    : %s\n',   cause);
 
        % Plot raw baseline and full LiCl-window trace for visual inspection
        t_bl  = (1:pre2_samp) / matfile_fr;
        t_all = (1:size(traces2,2)) / matfile_fr;
 
        figure('Name', sprintf('Outlier: File %d Neuron %d', idx, n));
        subplot(3,1,1);
        plot(t_bl, raw_bl, 'k');
        xline(prestim2, '--r', 'Stim onset');
        xlabel('Time (s)'); ylabel('Raw F');
        title(sprintf('File %d | Neuron %d — Baseline window (%.4f ± %.6f)', ...
            idx, n, bl_mean, bl_std));
 
        subplot(3,1,2);
        plot(t_all, traces2(n,:), 'b');
        xline(prestim2, '--r', 'LiCl');
        xlabel('Time (s)'); ylabel('Raw F');
        title('Full LiCl-window raw trace');
 
        subplot(3,1,3);
        plot(t_all, traces_zs_LiCl(n,:), 'm');
        xline(prestim2, '--r', 'LiCl');
        yline(z_outlier_thresh,  '--k', sprintf('+%g threshold', z_outlier_thresh));
        yline(-z_outlier_thresh, '--k');
        xlabel('Time (s)'); ylabel('Z-score');
        title('Z-scored LiCl trace');
 
        sgtitle(strrep(sprintf('%s  —  %s', csv_file, cause), '_', '\_'), ...
            'Interpreter', 'tex', 'FontWeight', 'bold');
    end
end
 
if ~any_found
    fprintf('  No neurons exceeded |z| = %g. Try lowering z_outlier_thresh.\n', ...
        z_outlier_thresh);
end