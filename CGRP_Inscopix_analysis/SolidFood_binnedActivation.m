%% Chow Feeding: Binned Activation Analysis
%  Calculates the percentage of activated (z>1) and inhibited (z<-1)
%  neurons across food access in 5-minute bins, for:
%    1. Whole population
%    2. Rapidly activated neurons (z>1 mean in first 5 OR first 10 min)
%    3. Late-activated neurons (NOT rapid, but z>1 mean in last 10 min)
%
%  Inputs:
%    - Chow.txt   : [csvFile  fr  stim_mm:ss  prestim_s  end_s  videoOffset_mm:ss]
%    - *_chow.csv : calcium traces (rows=timepoints, cols=neurons;
%                   row 1 = cell IDs, row 2 = status, row 3+ = data)


clear; close all;

%% ========== USER PARAMETERS ==========
chow_txt        = 'Chow.txt';      % Master file list
prestim_s       = 600;             % Pre-stimulus baseline duration (s)
bin_size_s      = 300;             % Bin size (5 min = 300 s)
analysis_dur_s  = 3600;            % Post-access duration to analyse (s).
                                   % Set to the longest session you want to
                                   % include (e.g. 3600 = 60 min).
                                   % Animals whose recording ends earlier will
                                   % have NaN in the remaining bins.
act_threshold   = 1;               % Z-score threshold for activation
inhib_threshold = -1;              % Z-score threshold for inhibition

% Classification windows (seconds from food access onset)
rapid_win1_s    = [0, 300];        % First 5 min
rapid_win2_s    = [0, 600];        % First 10 min
late_win_s      = [-600, 0];       % Last 10 min (offset from session end)

%% ========== READ FILE LIST ==========
fileID = fopen(chow_txt, 'r');
header = fgetl(fileID);  % skip header line
filelines = {};
while ~feof(fileID)
    line = strtrim(fgetl(fileID));
    if ~isempty(line)
        filelines{end+1} = line;
    end
end
fclose(fileID);

num_files = numel(filelines);
fprintf('Found %d animal(s) in %s\n', num_files, chow_txt);

%% ========== OUTPUT ACCUMULATORS ==========
% We store per-bin data across all animals, then pool
% Each cell array entry: [nAnimals x nBins] of fraction activated/inhibited

all_bin_edges          = {};   % bin edges per animal (for alignment)
% Activated subgroups (classified by z > act_threshold in early/late windows)
all_frac_act_all       = {};   % % of all neurons activated  — whole population
all_frac_inhib_all     = {};   % % of all neurons inhibited  — whole population
all_frac_act_rap       = {};   % % of all neurons — rapidly activated subgroup
all_frac_inhib_rap     = {};   % % of all neurons — rapidly inhibited subgroup
all_frac_act_late      = {};   % % of all neurons — late activated subgroup
all_frac_inhib_late    = {};   % % of all neurons — late inhibited subgroup

% Pooled z-score matrices: all neurons concatenated across animals
% Each cell: [nNeurons_animal x nBins] mean-z per neuron per bin
pool_binZ_all       = {};   % whole population
pool_binZ_rap_act   = {};   % rapidly activated neurons
pool_binZ_late_act  = {};   % late activated neurons
pool_binZ_rap_inh   = {};   % rapidly inhibited neurons
pool_binZ_late_inh  = {};   % late inhibited neurons
pool_nbins          = {};

% Neuron counts per category
total_n_all       = 0;
total_n_rap_act   = 0;
total_n_late_act  = 0;
total_n_rap_inh   = 0;
total_n_late_inh  = 0;

%% ========== MAIN LOOP ==========
for idx = 1:num_files

    parts       = strsplit(strtrim(filelines{idx}), '\t');
    csvFileName = strtrim(parts{1});
    matfile_fr  = str2double(parts{2});

    tparts    = strsplit(strtrim(parts{3}), ':');
    stim_time = str2double(tparts{1})*60 + str2double(tparts{2});

    % end_s from Chow.txt is not used for windowing — we always take
    % exactly prestim_s before stim_time and analysis_dur_s after it.
    % This ensures all animals are aligned on the same grid regardless
    % of when in the session food access occurred.

    fprintf('\n[%d/%d] %s  stim=%.0fs  pre=%.0fs  post=%.0fs\n', ...
        idx, num_files, csvFileName, stim_time, prestim_s, analysis_dur_s);

    if ~exist(csvFileName, 'file')
        warning('File not found: %s — skipping.', csvFileName);
        continue
    end

    %% --- Load and z-score traces ---
    opts = detectImportOptions(csvFileName);
    opts.DataLines             = [3, Inf];
    opts.SelectedVariableNames = opts.VariableNames(2:end);  % drop time col
    raw_traces = table2array(readtable(csvFileName, opts))'; % [neurons x frames_full]

    % Window: [stim_time - prestim_s,  stim_time + analysis_dur_s]
    % Clamp end to actual file length — animals whose recording ends
    % before analysis_dur_s will have NaN-filled trailing bins.
    start_idx_rec = round((stim_time - prestim_s) * matfile_fr) + 1;
    end_idx_rec   = min(round((stim_time + analysis_dur_s) * matfile_fr), size(raw_traces, 2));
    traces        = raw_traces(:, start_idx_rec:end_idx_rec);

    prestim_len  = round(prestim_s * matfile_fr);
    prestim_mean = mean(traces(:, 1:prestim_len), 2, 'omitnan');
    prestim_std  = std( traces(:, 1:prestim_len), 0, 2, 'omitnan');
    prestim_std(prestim_std < 1e-6) = 1e-6;
    trace_zs     = (traces - prestim_mean) ./ prestim_std;  % [neurons x frames]

    [num_neurons, num_frames] = size(trace_zs);
    time_offset = stim_time - prestim_s;
    time_vec    = time_offset + (0:num_frames-1) / matfile_fr;  % absolute time

    % Post-access frames (t >= stim_time)
    post_mask = time_vec >= stim_time;
    post_idx  = find(post_mask);
    t_post    = time_vec(post_idx) - stim_time;   % seconds from food access

    fprintf('  Neurons: %d  |  Post-access frames: %d  (%.1f min)\n', ...
        num_neurons, numel(post_idx), max(t_post)/60);

    %% --- 5-min bin structure (fixed length from analysis_dur_s) ---
    % n_bins is the SAME for every animal; bins beyond the recording end
    % will be NaN for that animal (handled inside compute_fractions /
    % neuron_bin_means when no frames fall in that bin).
    n_bins      = floor(analysis_dur_s / bin_size_s);
    bin_edges   = (0:n_bins) * bin_size_s;         % edges in s from access
    bin_centers = bin_edges(1:end-1) + bin_size_s/2;

    fprintf('  Bins: %d x %.0f s  |  actual frames cover %.1f min\n', ...
        n_bins, bin_size_s, max(t_post)/60);

    %% --- Neuron classification ---
    % Helper: mean z for a neuron in a time window [t_start, t_end] (s from access)
    mean_z_window = @(nrn, t_start, t_end) ...
        mean(trace_zs(nrn, post_idx(t_post >= t_start & t_post < t_end)), 'omitnan');

    % --- Precompute window means for every neuron (vectorised) ---
    mz5_all  = arrayfun(@(n) mean_z_window(n, rapid_win1_s(1), rapid_win1_s(2)), 1:num_neurons)';
    mz10_all = arrayfun(@(n) mean_z_window(n, rapid_win2_s(1), rapid_win2_s(2)), 1:num_neurons)';

    t_session_end = max(t_post);
    late_start    = max(0, t_session_end + late_win_s(1));
    mz_late_all   = arrayfun(@(n) mean_z_window(n, late_start, t_session_end+1), 1:num_neurons)';

    % --- Five mutually exclusive groups, assigned in priority order ---
    % Each neuron lands in exactly one group; no neuron appears in two.
    %
    %   Priority 1: Rapid-activated  — mean z > act_threshold  in 0–5 OR 0–10 min
    %   Priority 2: Rapid-inhibited  — mean z < inhib_threshold in 0–5 OR 0–10 min
    %                                  (only if NOT rapid-activated)
    %   Priority 3: Late-activated   — mean z > act_threshold  in last 10 min
    %                                  (only if NOT rapid)
    %   Priority 4: Late-inhibited   — mean z < inhib_threshold in last 10 min
    %                                  (only if NOT rapid and NOT late-activated)
    %   Priority 5: Other            — everything else

    is_rapid_act = (mz5_all > act_threshold)   | (mz10_all > act_threshold);
    is_rapid_inh = ~is_rapid_act & ...
                   ((mz5_all < inhib_threshold) | (mz10_all < inhib_threshold));
    is_rapid     = is_rapid_act | is_rapid_inh;
    is_late_act  = ~is_rapid & (mz_late_all > act_threshold);
    is_late_inh  = ~is_rapid & ~is_late_act & (mz_late_all < inhib_threshold);

    rapid_act_idx = find(is_rapid_act);
    rapid_inh_idx = find(is_rapid_inh);
    late_act_idx  = find(is_late_act);
    late_inh_idx  = find(is_late_inh);
    all_idx       = 1:num_neurons;

    n_other = num_neurons - numel(rapid_act_idx) - numel(rapid_inh_idx) ...
                          - numel(late_act_idx)  - numel(late_inh_idx);
    fprintf('  Rapid-act: %d  Rapid-inh: %d  Late-act: %d  Late-inh: %d  Other: %d\n', ...
        numel(rapid_act_idx), numel(rapid_inh_idx), ...
        numel(late_act_idx),  numel(late_inh_idx), n_other);

    total_n_all      = total_n_all      + num_neurons;
    total_n_rap_act  = total_n_rap_act  + numel(rapid_act_idx);
    total_n_rap_inh  = total_n_rap_inh  + numel(rapid_inh_idx);
    total_n_late_act = total_n_late_act + numel(late_act_idx);
    total_n_late_inh = total_n_late_inh + numel(late_inh_idx);

    %% --- Compute per-bin fractions ---
    % Denominator is always num_neurons (whole animal) so all values are
    % "% of all neurons" and groups are directly comparable / additive.
    compute_bin_fracs = @(neuron_list) compute_fractions( ...
        trace_zs, neuron_list, post_idx, t_post, ...
        n_bins, bin_edges, act_threshold, inhib_threshold, num_neurons);

    [frac_act_all,  frac_inhib_all]  = compute_bin_fracs(all_idx);
    [frac_act_rap,  frac_inhib_rap]  = compute_bin_fracs(rapid_act_idx);
    [frac_act_late, frac_inhib_late] = compute_bin_fracs(late_act_idx);
    % For inhibited subgroups we only need the inhibited output
    [~, frac_inhib_rap_inh]  = compute_bin_fracs(rapid_inh_idx);
    [~, frac_inhib_late_inh] = compute_bin_fracs(late_inh_idx);

    %% --- Per-neuron mean-z per bin (for pooled calculation) ---
    nbm = @(idx) neuron_bin_means(trace_zs, idx, post_idx, t_post, n_bins, bin_edges);
    pool_binZ_all{end+1}      = nbm(all_idx);        
    pool_binZ_rap_act{end+1}  = nbm(rapid_act_idx);  
    pool_binZ_late_act{end+1} = nbm(late_act_idx);   
    pool_binZ_rap_inh{end+1}  = nbm(rapid_inh_idx); 
    pool_binZ_late_inh{end+1} = nbm(late_inh_idx);   
    pool_nbins{end+1}         = n_bins;               

    %% --- Store per-animal fractions ---
    all_bin_edges{end+1}        = bin_edges;           
    all_frac_act_all{end+1}     = frac_act_all;        
    all_frac_inhib_all{end+1}   = frac_inhib_all;      
    all_frac_act_rap{end+1}     = frac_act_rap;        
    all_frac_inhib_rap{end+1}   = frac_inhib_rap_inh;    % rapidly inhibited
    all_frac_act_late{end+1}    = frac_act_late;      
    all_frac_inhib_late{end+1}  = frac_inhib_late_inh;   % late inhibited

end  % file loop

%% ========== AGGREGATE ACROSS ANIMALS ==========
% All animals share the same n_bins (derived from analysis_dur_s).
% Animals with shorter recordings contribute NaN in late bins; mean/SEM
% are computed with 'omitnan' so those bins are simply averaged over
% however many animals actually had data there.
n_bins_fixed    = floor(analysis_dur_s / bin_size_s);
fprintf('\nFixed bin count: %d  (%.0f min total)\n', n_bins_fixed, analysis_dur_s/60);

% Per-animal matrices [n_animals x n_bins_fixed], all values % of total neurons
mat_act_all      = cell2mat(all_frac_act_all')    * 100;
mat_inhib_all    = cell2mat(all_frac_inhib_all')  * 100;
mat_act_rap      = cell2mat(all_frac_act_rap')    * 100;  % rapidly activated subgroup
mat_inhib_rap    = cell2mat(all_frac_inhib_rap')  * 100;  % rapidly inhibited subgroup
mat_act_late     = cell2mat(all_frac_act_late')   * 100;  % late activated subgroup
mat_inhib_late   = cell2mat(all_frac_inhib_late') * 100;  % late inhibited subgroup

bin_centers_min = ((1:n_bins_fixed) - 0.5) * (bin_size_s/60);  % bin centers in minutes

% Group mean ± SEM across animals
mean_sem = @(M) deal(mean(M,1,'omitnan'), std(M,0,1,'omitnan') ./ sqrt(sum(~isnan(M),1)));

[m_act_all,    s_act_all]    = mean_sem(mat_act_all);
[m_inhib_all,  s_inhib_all]  = mean_sem(mat_inhib_all);
[m_act_rap,    s_act_rap]    = mean_sem(mat_act_rap);
[m_inhib_rap,  s_inhib_rap]  = mean_sem(mat_inhib_rap);
[m_act_late,   s_act_late]   = mean_sem(mat_act_late);
[m_inhib_late, s_inhib_late] = mean_sem(mat_inhib_late);

%% ========== POOLED POPULATION (all neurons concatenated) ==========
poolZ_all      = cell2mat(pool_binZ_all');      % [total_n_all      x n_bins_fixed]
poolZ_rap_act  = cell2mat(pool_binZ_rap_act');  % [total_n_rap_act  x n_bins_fixed]
poolZ_late_act = cell2mat(pool_binZ_late_act'); % [total_n_late_act x n_bins_fixed]
poolZ_rap_inh  = cell2mat(pool_binZ_rap_inh');  % [total_n_rap_inh  x n_bins_fixed]
poolZ_late_inh = cell2mat(pool_binZ_late_inh'); % [total_n_late_inh x n_bins_fixed]

% All fractions expressed as % of total_n_all — directly comparable across panels
pool_frac_act_all    = sum(poolZ_all      > act_threshold,   1,'omitnan') / total_n_all * 100;
pool_frac_inhib_all  = sum(poolZ_all      < inhib_threshold, 1,'omitnan') / total_n_all * 100;
pool_frac_act_rap    = sum(poolZ_rap_act  > act_threshold,   1,'omitnan') / total_n_all * 100;
pool_frac_inhib_rap  = sum(poolZ_rap_inh  < inhib_threshold, 1,'omitnan') / total_n_all * 100;
pool_frac_act_late   = sum(poolZ_late_act > act_threshold,   1,'omitnan') / total_n_all * 100;
pool_frac_inhib_late = sum(poolZ_late_inh < inhib_threshold, 1,'omitnan') / total_n_all * 100;

%% ========== PRINT SUMMARY TABLE ==========
fprintf('\n\n====================================================\n');
fprintf('NEURON COUNTS  (pooled across %d animals)\n', num_files);
fprintf('  Total neurons    : %d\n', total_n_all);
fprintf('  Rapid activated  : %d (%.1f%%)\n', total_n_rap_act,  100*total_n_rap_act/total_n_all);
fprintf('  Late activated   : %d (%.1f%%)\n', total_n_late_act, 100*total_n_late_act/total_n_all);
fprintf('  Rapid inhibited  : %d (%.1f%%)\n', total_n_rap_inh,  100*total_n_rap_inh/total_n_all);
fprintf('  Late inhibited   : %d (%.1f%%)\n', total_n_late_inh, 100*total_n_late_inh/total_n_all);

fprintf('\n--- SANITY CHECK: rapid + late + other should sum to whole population ---\n');
fprintf('  (pooled, bin 1)  whole=%.1f%%  rapid=%.1f%%  late=%.1f%%  sum_sub=%.1f%%\n', ...
    pool_frac_act_all(1), pool_frac_act_rap(1), pool_frac_act_late(1), ...
    pool_frac_act_rap(1)+pool_frac_act_late(1));

fprintf('\n--- WHOLE POPULATION: %% Activated per 5-min bin ---\n');
fprintf('Bin (min)\t%%Act (mean±SEM)\t%%Inhib (mean±SEM)\n');
for b = 1:n_bins_fixed
    fprintf('  %4.0f–%4.0f\t%.1f ± %.1f\t\t%.1f ± %.1f\n', ...
        (b-1)*5, b*5, ...
        m_act_all(b), s_act_all(b), ...
        m_inhib_all(b), s_inhib_all(b));
end

fprintf('\n--- RAPID ACT NEURONS: %% per 5-min bin (of all neurons) ---\n');
fprintf('Bin (min)\t%%Act (mean±SEM)\n');
for b = 1:n_bins_fixed
    fprintf('  %4.0f–%4.0f\t%.1f ± %.1f\n', (b-1)*5, b*5, m_act_rap(b), s_act_rap(b));
end

fprintf('\n--- RAPID INHIB NEURONS: %% per 5-min bin (of all neurons) ---\n');
fprintf('Bin (min)\t%%Inhib (mean±SEM)\n');
for b = 1:n_bins_fixed
    fprintf('  %4.0f–%4.0f\t%.1f ± %.1f\n', (b-1)*5, b*5, m_inhib_rap(b), s_inhib_rap(b));
end

fprintf('\n--- LATE ACT NEURONS: %% per 5-min bin (of all neurons) ---\n');
fprintf('Bin (min)\t%%Act (mean±SEM)\n');
for b = 1:n_bins_fixed
    fprintf('  %4.0f–%4.0f\t%.1f ± %.1f\n', (b-1)*5, b*5, m_act_late(b), s_act_late(b));
end

fprintf('\n--- LATE INHIB NEURONS: %% per 5-min bin (of all neurons) ---\n');
fprintf('Bin (min)\t%%Inhib (mean±SEM)\n');
for b = 1:n_bins_fixed
    fprintf('  %4.0f–%4.0f\t%.1f ± %.1f\n', (b-1)*5, b*5, m_inhib_late(b), s_inhib_late(b));
end

fprintf('\n--- POOLED: Whole Population ---\n');
fprintf('Bin (min)\t%%Act\t%%Inhib\n');
for b = 1:n_bins_fixed
    fprintf('  %4.0f–%4.0f\t%.1f\t%.1f\n', (b-1)*5, b*5, pool_frac_act_all(b), pool_frac_inhib_all(b));
end

fprintf('\n--- POOLED: Rapid Activated / Rapid Inhibited ---\n');
fprintf('Bin (min)\t%%RapAct\t%%RapInhib\n');
for b = 1:n_bins_fixed
    fprintf('  %4.0f–%4.0f\t%.1f\t%.1f\n', (b-1)*5, b*5, pool_frac_act_rap(b), pool_frac_inhib_rap(b));
end

fprintf('\n--- POOLED: Late Activated / Late Inhibited ---\n');
fprintf('Bin (min)\t%%LateAct\t%%LateInhib\n');
for b = 1:n_bins_fixed
    fprintf('  %4.0f–%4.0f\t%.1f\t%.1f\n', (b-1)*5, b*5, pool_frac_act_late(b), pool_frac_inhib_late(b));
end

%% ========== FIGURE 1: Per-animal mean ± SEM ==========
% -----------------------------------------------------------------------
% Subplot 1 — Whole Population
%   Activated line : m_act_all   ± s_act_all    — ALL neurons, binwise z>1  (mat_act_all)
%   Inhibited line : m_inhib_all ± s_inhib_all  — ALL neurons, binwise z<-1 (mat_inhib_all)
%   No pre-selection; every neuron tested fresh each bin
%
% Subplot 2 — Rapid subgroups (classified first, take priority)
%   Activated line : m_act_rap   ± s_act_rap    — neurons with mean z>1 in 0–5 OR 0–10 min
%                    Source: mat_act_rap  [n_animals x n_bins_fixed]
%                    Classifier: is_rapid_act = (mz5>1)|(mz10>1)  → rapid_act_idx
%   Inhibited line : m_inhib_rap ± s_inhib_rap  — neurons with mean z<-1 in 0–5 OR 0–10 min
%                    Source: mat_inhib_rap [n_animals x n_bins_fixed]
%                    Classifier: is_rapid_inh = ~is_rapid_act & (mz5<-1|mz10<-1) → rapid_inh_idx
%
% Subplot 3 — Late subgroups (only neurons NOT classified as rapid)
%   Activated line : m_act_late   ± s_act_late   — NOT rapid, mean z>1 in last 10 min
%                    Source: mat_act_late  [n_animals x n_bins_fixed]
%                    Classifier: is_late_act = ~is_rapid & mz_late>1  → late_act_idx
%   Inhibited line : m_inhib_late ± s_inhib_late — NOT rapid AND NOT late-act, mean z<-1 in last 10 min
%                    Source: mat_inhib_late [n_animals x n_bins_fixed]
%                    Classifier: is_late_inh = ~is_rapid & ~is_late_act & mz_late<-1 → late_inh_idx
%
% X-axis: bin_centers_min [1 x n_bins_fixed] (min from food access)
% Y-axis: % of all neurons (denominator = num_neurons per animal)
% Shading: ±1 SEM across animals
% -----------------------------------------------------------------------
figure('Name','Binned Activation — Per-animal mean±SEM','Position',[50 50 1300 500]);
colors = struct( ...
    'act',   [0.85 0.33 0.10], ...   % orange-red
    'inhib', [0.00 0.45 0.74], ...   % blue
    'rap',   [0.93 0.69 0.13], ...   % gold
    'late',  [0.47 0.67 0.19]);      % green

subplot_titles = {'Whole Population', ...
                  'Rapid (act=gold, inhib=blue)', ...
                  'Late (act=green, inhib=blue)'};
mA  = {m_act_all,   m_act_rap,   m_act_late};
sA  = {s_act_all,   s_act_rap,   s_act_late};
mI  = {m_inhib_all, m_inhib_rap, m_inhib_late};
sI  = {s_inhib_all, s_inhib_rap, s_inhib_late};
act_col = {colors.act, colors.rap, colors.late};

for p = 1:3
    subplot(1,3,p); hold on;

    fill_x = [bin_centers_min, fliplr(bin_centers_min)];
    fill(fill_x, [mA{p}+sA{p}, fliplr(mA{p}-sA{p})], act_col{p},   'FaceAlpha',0.25, 'EdgeColor','none');
    fill(fill_x, [mI{p}+sI{p}, fliplr(mI{p}-sI{p})], colors.inhib, 'FaceAlpha',0.25, 'EdgeColor','none');

    plot(bin_centers_min, mA{p}, '-o', 'Color', act_col{p},   'LineWidth',2, 'MarkerSize',5);
    plot(bin_centers_min, mI{p}, '-s', 'Color', colors.inhib, 'LineWidth',2, 'MarkerSize',5);

    yline(0, 'k--', 'LineWidth',0.8);
    xlabel('Time from food access (min)');
    ylabel('% of all neurons');
    title(subplot_titles{p});
    legend({'Activated (z>1)','Inhibited (z<-1)'}, 'Location','best');
    xlim([0 max(bin_centers_min)+bin_size_s/60/2]);
    ylim([0 100]);
    set(gca,'TickDir','out','Box','off','FontSize',11);
end

sgtitle(sprintf('5-min Binned Activity — Mean \pm SEM across animals  (n=%d)', num_files), ...
    'FontSize',13,'FontWeight','bold');

%% ========== FIGURE 2: Pooled population (single line per group) ==========
% -----------------------------------------------------------------------
% Subplot 1 — Whole Population (pooled)
%   Activated line : pool_frac_act_all   — sum(poolZ_all>1) / total_n_all * 100
%   Inhibited line : pool_frac_inhib_all — sum(poolZ_all<-1) / total_n_all * 100
%   Source: poolZ_all [total_n_all x n_bins_fixed]
%
% Subplot 2 — Rapid subgroups (pooled, take classification priority)
%   Activated line : pool_frac_act_rap   — sum(poolZ_rap_act>1)  / total_n_all * 100
%   Inhibited line : pool_frac_inhib_rap — sum(poolZ_rap_inh<-1) / total_n_all * 100
%   Sources: poolZ_rap_act  [total_n_rap_act x n_bins_fixed]  (rapid_act_idx neurons)
%            poolZ_rap_inh  [total_n_rap_inh x n_bins_fixed]  (rapid_inh_idx: ~rapid_act & early z<-1)
%
% Subplot 3 — Late subgroups (pooled, only neurons NOT classified as rapid)
%   Activated line : pool_frac_act_late   — sum(poolZ_late_act>1)  / total_n_all * 100
%   Inhibited line : pool_frac_inhib_late — sum(poolZ_late_inh<-1) / total_n_all * 100
%   Sources: poolZ_late_act [total_n_late_act x n_bins_fixed]  (~rapid & late z>1)
%            poolZ_late_inh [total_n_late_inh x n_bins_fixed]  (~rapid & ~late_act & late z<-1)
%
% X-axis: bin_centers_min [1 x n_bins_fixed]
% Y-axis: % of total_n_all — all panels share the same denominator
% No SEM — single value per bin across entire pooled dataset
% -----------------------------------------------------------------------
figure('Name','Binned Activation — Pooled neurons','Position',[50 580 1300 420]);

pool_act_dat   = {pool_frac_act_all,  pool_frac_act_rap,  pool_frac_act_late};
pool_inhib_dat = {pool_frac_inhib_all,pool_frac_inhib_rap,pool_frac_inhib_late};
pool_titles    = { ...
    sprintf('Whole Population  (N=%d)', total_n_all), ...
    sprintf('Rapid  (act n=%d, inh n=%d)', total_n_rap_act, total_n_rap_inh), ...
    sprintf('Late  (act n=%d, inh n=%d)',  total_n_late_act, total_n_late_inh)};

for p = 1:3
    subplot(1,3,p); hold on;
    plot(bin_centers_min, pool_act_dat{p},   '-o', 'Color', act_col{p},   'LineWidth',2, 'MarkerSize',5);
    plot(bin_centers_min, pool_inhib_dat{p}, '-s', 'Color', colors.inhib, 'LineWidth',2, 'MarkerSize',5);
    yline(0,'k--','LineWidth',0.8);
    xlabel('Time from food access (min)');
    ylabel('% of all neurons');
    title(pool_titles{p});
    legend({'Activated (z>1)','Inhibited (z<-1)'},'Location','best');
    xlim([0 max(bin_centers_min)+bin_size_s/60/2]);
    ylim([0 100]);
    set(gca,'TickDir','out','Box','off','FontSize',11);
end

sgtitle(sprintf('5-min Binned Activity — Pooled neurons  (%d animals)', num_files), ...
    'FontSize',13,'FontWeight','bold');

%% ========== SAVE ==========
save('BinnedActivation_results.mat', ...
    'mat_act_all','mat_inhib_all', ...
    'mat_act_rap','mat_inhib_rap', ...
    'mat_act_late','mat_inhib_late', ...
    'pool_frac_act_all','pool_frac_inhib_all', ...
    'pool_frac_act_rap','pool_frac_inhib_rap', ...
    'pool_frac_act_late','pool_frac_inhib_late', ...
    'bin_centers_min','bin_size_s', ...
    'm_act_all','s_act_all','m_inhib_all','s_inhib_all', ...
    'm_act_rap','s_act_rap','m_inhib_rap','s_inhib_rap', ...
    'm_act_late','s_act_late','m_inhib_late','s_inhib_late', ...
    'total_n_all','total_n_rap_act','total_n_late_act', ...
    'total_n_rap_inh','total_n_late_inh', ...
    'pool_binZ_all','pool_binZ_rap_act','pool_binZ_rap_inh', ...
    'pool_binZ_late_act','pool_binZ_late_inh');

fprintf('\nResults saved to BinnedActivation_results.mat\n');
% =========================================================================
%% LOCAL FUNCTION 2
% =========================================================================
function binZ = neuron_bin_means( ...
    trace_zs, neuron_list, post_idx, t_post, n_bins, bin_edges)
% Returns mean z-score per neuron per bin — used for pooled population fractions.
%
% Output:
%   binZ  [n_neurons x n_bins]  mean z-score of each neuron in each bin

    n_neurons = numel(neuron_list);
    binZ = NaN(n_neurons, n_bins);

    if n_neurons == 0
        return
    end

    for b = 1:n_bins
        bin_mask   = t_post >= bin_edges(b) & t_post < bin_edges(b+1);
        bin_frames = post_idx(bin_mask);
        if isempty(bin_frames), continue; end
        binZ(:,b) = mean(trace_zs(neuron_list, bin_frames), 2, 'omitnan');
    end
end
function [frac_act, frac_inhib] = compute_fractions( ...
    trace_zs, neuron_list, post_idx, t_post, ...
    n_bins, bin_edges, act_thr, inhib_thr, n_total)
% Returns fraction of neurons in neuron_list that are activated / inhibited
% in each 5-min bin, expressed as a fraction of n_total (the whole-session
% neuron count). Pass n_total = numel(neuron_list) for within-group fractions,
% or n_total = num_neurons (all cells) to get subgroup as % of whole population.
%
% Inputs:
%   trace_zs    [nNeurons x nFrames]  full z-scored trace
%   neuron_list  vector of neuron indices to evaluate
%   post_idx     indices into trace_zs columns that are post-access
%   t_post       time in seconds from food access for each post_idx frame
%   n_bins       number of bins
%   bin_edges    [1 x n_bins+1] bin edges in seconds from access
%   act_thr      activation threshold
%   inhib_thr    inhibition threshold
%   n_total      denominator — use num_neurons (all cells) so subgroup
%                fractions read as "% of whole population"
%
% Outputs:
%   frac_act    [1 x n_bins] fraction (0-1) of n_total that are activated
%   frac_inhib  [1 x n_bins] fraction (0-1) of n_total that are inhibited

    n_sub  = numel(neuron_list);
    frac_act   = NaN(1, n_bins);
    frac_inhib = NaN(1, n_bins);

    if n_sub == 0 || n_total == 0
        return
    end

    for b = 1:n_bins
        t_lo = bin_edges(b);
        t_hi = bin_edges(b+1);
        bin_mask = t_post >= t_lo & t_post < t_hi;
        bin_frames = post_idx(bin_mask);

        if isempty(bin_frames)
            continue
        end

        % Mean z-score per neuron in this bin
        bin_z = mean(trace_zs(neuron_list, bin_frames), 2, 'omitnan');  % [n_sub x 1]

        % Divide by n_total so rapid/late values are "% of all neurons"
        frac_act(b)   = sum(bin_z >  act_thr)  / n_total;
        frac_inhib(b) = sum(bin_z <  inhib_thr) / n_total;
    end
end