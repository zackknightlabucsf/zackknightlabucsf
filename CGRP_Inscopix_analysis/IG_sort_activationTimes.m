%% Quant IG and IP responses (or compare two time points during the trial)
% Use to sort neurons based on activity for some time after a stimulus
%% Input parameters
clear;

infile = '1-5mLIGEnsureCGRP1-6_zs'; %file name? This is the wspace file outputted from the heatmap function 
Title = '1-5mLIGEnsureCGRP1-6_zs_0528266_v6'
sampling_rate = 4;  % Hz
prestim1 = 600 * sampling_rate;     %when was the first stimulus given in seconds? 
prestim2 = 1500 * sampling_rate;    %when was the second stimulus given in seconds? 800 for IP, 1200 for 1 mL IG, 1500 for 1.5 mL IG
prestim3 = 600 * sampling_rate;      %Time of stimulus in seconds (600 for IG infusions)
Int = 600 * sampling_rate;          %how many seconds post stimulus do you want to average over (e.g. 600=10min post lick access), 600 1 mL IG, 900 1.5 mL IG
Int2 = 900 * sampling_rate;          %how many seconds post stimulus do you want to average over for the second stimulus (e.g. 600=10min post lick access)600 IP, 600 IG
Int3 = 1799 * sampling_rate;          %100 for IG infusions 600 for IP

window_size = 60;               % Sliding window size (seconds) 60
min_duration = 42 ;  % Minimum duration in consecutive data points for the second activation screen (60 for 60 seconds) 42
average_duration = 600; % Duration for averaging after activation (seconds)
piechart = 1; %Do you want a pie chart plotted? 0 = no, 1 = yes

prestim_baseline_start = 500;          % sample index – well before stimulus
prestim_baseline_end   = prestim1;     % up to (but not including) stimulus onset
prestim_elevation_thresh = 1;          % z-score; same as activation threshold
%% Load Data
load(infile, '-mat', 'neurons_zscored'); %just load in neurons_zscored to make the workspace cleaner. 

% Average traces, smoothed for graphing
graph_all_neurons = [];
 for i = 1:size(neurons_zscored,1)
     graph = smooth(neurons_zscored(i,:),20);
     graph_all_neurons(i,:) = graph;
 end
 
%% First Filter: Classify Neurons based on mean activity
%Initialize storage variables
act = [];
none=[];
inhib = [];

for i = 1:size(neurons_zscored,1) %step through each neuron
    a = mean(neurons_zscored(i, prestim1:prestim1+Int));
    b = mean(neurons_zscored(i, prestim2:prestim2+Int2));
    c = mean(neurons_zscored(i, prestim1:prestim1+300*sampling_rate));
    
    if a> 1 || b> 1
        act = [act i]; %Neurons are activated at any timepoint
       
     elseif a<-1 || b<-1
        inhib = [inhib i]; %Store the neuron number for neurons that are classified as "inhibited"--index we'll pull later
       
    else
        none = [none i]; %Neurons show no significant change
    end
end

%% Store Neuron Traces Based on Initial Classification
%Each row is a neuron, and each column is a timepoint. Number of rows=number of neurons per category 
Neurons_screen = neurons_zscored(act, :); %use the indexes we found above and pull out only those neurons
Neurons_none = neurons_zscored(none, :);
Neurons_inhib = neurons_zscored(inhib,:); 

%% Second Filter: Identify Neurons with Sustained Activation (excludes large single events). Determine activation time when threshold after stimulus > 1 z-score
% Sampling rate and time step
time_step = 1 / sampling_rate; % 1 second
cells = Neurons_screen;
threshold = 1; % Threshold for activation (z-score)


% Convert window size to number of samples
window_size_samples = window_size * sampling_rate;

% Total number of samples in the trace
total_samples = size(cells, 2);

% Preallocate array for first activation times
first_act_time = zeros(size(cells, 1), 1); % 0 for neurons with no activation

% Time vector
time_vector = (1:total_samples) * time_step;

%% Loop over each neuron in Neurons_screen
% Rising-edge crossings (below → above) within the post-stim trace
   %% Loop over each neuron in Neurons_screen
for i = 1:size(cells, 1)
    full_trace = cells(i, :);


pre_stim_mean = mean(full_trace(prestim_baseline_start:prestim_baseline_end));
    if pre_stim_mean > prestim_elevation_thresh
        first_act_time(i) = 0;  % pre-existing elevation – exclude
        continue;
    end
 
    % Post-stimulus trace (starting at sample prestim1+1)
    post_stim_start_idx = prestim1 + 1;
    post_stim_trace     = full_trace(post_stim_start_idx:end);
    total_post_stim     = length(post_stim_trace);
    above_threshold = post_stim_trace > threshold;


    crosses_from_below = [above_threshold(1), diff(above_threshold) == 1];
    crossing_indices   = find(crosses_from_below);  % indices in post_stim_trace
 
    found = false;
 
    for idx = crossing_indices
        % Define window in post-stim coordinates
        window_end_idx  = min(idx + window_size_samples - 1, total_post_stim);
        window_indices  = idx:window_end_idx;
 
        % Total seconds above threshold within this window
        total_active_duration = sum(above_threshold(window_indices)) * time_step;
 
        if total_active_duration >= min_duration
            % Convert post-stim index to time relative to stimulus onset
            % idx is 1-based within post_stim_trace; sample (prestim1 + idx)
            % Time relative to stim = idx * time_step  (since post_stim_trace(1) = t=time_step after stim)
            first_act_time(i)     = idx * time_step;   % seconds AFTER prestim1
            found = true;
            break;
        end
    end
 
    if ~found
        first_act_time(i) = 0;
    end
end
    

%% Identify Neurons that pass the second filter
active_cells_idx = find(first_act_time > 0); % Indices in Neurons_screen
Neurons_act = cells(active_cells_idx, :); % Array of neurons that became active

% first_act_time_passed is already in seconds relative to stimulus onset
first_act_time_passed = first_act_time(active_cells_idx);   % seconds after prestim1

%% Update 'none' and 'act' Category with Neurons Failing Second Filter
% Map indices back to original neuron indices in neurons_zscored
act_neurons_pass_both_filters = act(active_cells_idx); % Neuron indices that pass both filters
act_neurons_fail_second_filter = setdiff(act, act_neurons_pass_both_filters); % Neurons that fail second filter

% Update 'none' category
none = unique([none, act_neurons_fail_second_filter]);  % Update 'none' category

% Update indices
Neurons_none = neurons_zscored(none, :);
act = act_neurons_pass_both_filters;

%% Response time for Inhibited Neurons
first_inhib_time = zeros(length(inhib), 1);
 
for i = 1:length(inhib)
    full_trace      = neurons_zscored(inhib(i), :);
    post_stim_trace = full_trace(prestim1+1:end);
    total_post_stim = length(post_stim_trace);
    below_threshold = post_stim_trace < -threshold;
 
    % Falling-edge crossings (include first sample if already below)
    crosses       = [below_threshold(1), diff(below_threshold) == 1];
    crossing_idxs = find(crosses);
 
    for idx = crossing_idxs
        win_end         = min(idx + window_size_samples - 1, total_post_stim);
        suppressed_secs = sum(below_threshold(idx:win_end)) * time_step;
        if suppressed_secs >= min_duration
            first_inhib_time(i) = idx * time_step;  % seconds after prestim1
            break;
        end
    end
    % if no window found, first_inhib_time(i) remains 0
end
 
first_inhib_time_passed = first_inhib_time(first_inhib_time > 0);

%% Verify Categories are Mutually Exclusive and Exhaustive
% Total neurons
total_cells = size(neurons_zscored, 1);

% Verify that the sum of neurons in each category equals the total number of neurons
assert(length(act_neurons_pass_both_filters) + length(none) + length(inhib) == total_cells, ...
    'The categories do not sum up to the total number of neurons.');

%% Calculate the average z-score response during the IG infusion or during the post-infusion

% Define the intervals
interval1_start = prestim1;   % Start of Interval 1
interval1_end = prestim1 + Int;    % End of Interval 1

interval2_start = prestim2;  % Start of Interval 2
interval2_end = prestim2 + Int2;    % End of Interval 2

interval3_start = prestim3;  % Start of Interval 3
interval3_end = prestim3 + Int3;    % End of Interval 3

% % Preallocate arrays to store mean activities
% mean_activity_IG = zeros(size(neurons_zscored, 1), 1);
% mean_activity_postIG_10 = zeros(size(neurons_zscored, 1), 1);

% Compute mean activities during each interval
mean_activity_IG = mean(neurons_zscored(:, interval1_start:interval1_end), 2);
mean_activity_postIG_10 = mean(neurons_zscored(:, interval2_start:interval2_end), 2);

%% Curve of when neurons are activated
% exclude the zeros for the "activated only" curve
rel_times = first_act_time_passed(first_act_time_passed>0);

% define a uniform time axis from 0 to the max activation time
t_end     = ceil(max(rel_times));       % whole seconds
t_vec     = 0:time_step:t_end;          % sec
t_vec_min = t_vec/60;

% preallocate
num_all      = size(neurons_zscored,1);
num_activated = numel(rel_times);

prop_all     = zeros(size(t_vec));
prop_active  = zeros(size(t_vec));

% fill in with vectorized cumsum via histcounts
edges      = [-time_step/2, t_vec + time_step/2];  
counts     = histcounts(rel_times, edges);        % how many neurons first activated in each bin
cum_counts = cumsum(counts);                     

% percent of the *entire* pop that's activated by each t
prop_all    = cum_counts / num_all;              

% percent of the *ever*–activated neurons by each t
prop_active = cum_counts / num_activated;        

%% now plot it
figure;
plot(t_vec_min, prop_all,   'LineWidth',2); hold on;
plot(t_vec_min, prop_active,'--','LineWidth',2);
xlabel('Time since infusion (s)');
ylabel('Proportion of neurons activated');
legend('of all neurons','of neurons that ever activate','Location','southeast');
grid on;

%% Calculate Average Z-Score relative to response time for activated and inhibited neurons

% Initialize vector to store average z-scores
average_z_scores = zeros(length(active_cells_idx), 1);
total_cells = size(neurons_zscored, 1);
total_samples_all   = size(neurons_zscored, 2);

% For each neuron, calculate the average z-score response
for idx = 1:length(active_cells_idx)
    neuron_row      = act(idx);                              % row in neurons_zscored
    rel_onset_sec   = first_act_time(active_cells_idx(idx)); % seconds after prestim1
 
    start_idx = prestim1 + round(rel_onset_sec * sampling_rate);
    end_idx   = start_idx + round(average_duration * sampling_rate);
    start_idx = max(start_idx, 1);
    end_idx   = min(end_idx, total_samples_all);
 
    average_z_scores_act(idx) = mean(neurons_zscored(neuron_row, start_idx:end_idx));
end

%% Average Z-Score Response – INHIBITED neurons
% Computed from suppression onset; for neurons with no detected onset (time=0),
% averages from prestim1 instead.
average_z_scores_inhib = zeros(length(inhib), 1);
 
for idx = 1:length(inhib)
    neuron_row    = inhib(idx);
    rel_onset_sec = first_inhib_time(idx);  % 0 if no sustained window found
 
    if rel_onset_sec > 0
        start_idx = prestim1 + round(rel_onset_sec * sampling_rate);
    else
        start_idx = prestim1;  % fall back to infusion onset
    end
    end_idx   = start_idx + round(average_duration * sampling_rate);
    start_idx = max(start_idx, 1);
    end_idx   = min(end_idx, total_samples_all);
 
    average_z_scores_inhib(idx) = mean(neurons_zscored(neuron_row, start_idx:end_idx));
end

%% Calculate means for groups
ig_mean_act = mean_activity_IG(act, :);
ig_mean_inhib = mean_activity_IG(inhib, :);
ig_mean_none = mean_activity_IG(none, :);

postIG10_mean_act = mean_activity_postIG_10(act, :);
postIG10_mean_inhib = mean_activity_postIG_10(inhib, :);
postIG10_mean_none = mean_activity_postIG_10(none, :);

Allmean_IG_act = mean(ig_mean_act);
Allmean_IG_inhib = mean(ig_mean_inhib);
Allmean_PostIG_act = mean(postIG10_mean_act);
Allmean_PostIG_inhib = mean(postIG10_mean_inhib);

% %Calculate populated weighted values
fraction_act = length(ig_mean_act)/length(mean_activity_IG);
fraction_inhib = length(ig_mean_inhib)/length(mean_activity_IG);
pop_weight_IG_act = ig_mean_act * fraction_act;
pop_weight_IG_inhib = ig_mean_inhib * fraction_inhib;
pop_weight_postIG_act = postIG10_mean_act * fraction_act;
pop_weight_postIG_inhib = postIG10_mean_inhib * fraction_inhib;

% Sort the data in ascending order for CDF graphs
sortedData_IG = sort(mean_activity_IG);
n = length(sortedData_IG);

sortedData_PostIG10 = sort(mean_activity_postIG_10);
m = length (sortedData_PostIG10);

% Calculate the cumulative probability for each data point
cdfValues_IG = (1:n) / n;
cdfValues_postIG10 = (1:m)/m;
% cdfValues_postIG20 = (1:o)/o;

% Plot the CDF
figure;
plot(sortedData_IG, cdfValues_IG, 'LineWidth', 2);
xlabel('Activity Values');
ylabel('Cumulative Probability');
title('Empirical CDF of Mean Activity IG');
xlim([-5, 20]);
grid on;
%% Plot Pie Chart
if piechart == 1
    % Prepare data for pie chart
    X = [length(act_neurons_pass_both_filters), length(none), length(inhib)];
    labels = {'Activated', 'No Change', 'Inhibited'};
    
    % Plot pie chart
    figure;
    pie(X);
    colormap([0 1 0; 0.8 0.8 0.8; 1 0 0]); % Green, Gray, Red
    legend(labels, 'Location', 'southoutside', 'Orientation', 'horizontal');
    title(Title);
end

%% Calculate AUC from prestim3 across Int3 for all neurons
% Window: prestim3 to prestim3+Int3 (samples), e.g. 600s to 2400s
auc_window = neurons_zscored(:, interval3_start:interval3_end);  % [neurons x samples]

% AUC via trapezoidal integration; time_step converts samples to seconds
auc_all = trapz(auc_window, 2) * time_step;  % [neurons x 1], units = z-score * seconds

% Split AUC by final category
auc_act   = auc_all(act,   :);   % activated neurons (pass both filters)
auc_inhib = auc_all(inhib, :);   % inhibited neurons
auc_none  = auc_all(none,  :);   % non-responsive neurons

% Summary statistics
auc_mean_act   = mean(auc_act);
auc_mean_inhib = mean(auc_inhib);
auc_mean_none  = mean(auc_none);

fprintf('\n--- AUC Summary (prestim3 to prestim3+Int3) ---\n');
fprintf('  Activated  (n=%d): mean AUC = %.2f z*s\n', length(act),   auc_mean_act);
fprintf('  Inhibited  (n=%d): mean AUC = %.2f z*s\n', length(inhib), auc_mean_inhib);
fprintf('  No Change  (n=%d): mean AUC = %.2f z*s\n', length(none),  auc_mean_none);

%% Plot AUC sorted by category (Activated | Inhibited)
% Sort within each group by descending AUC
[auc_act_sorted,   sort_idx_act]   = sort(auc_act,   'descend');
[auc_inhib_sorted, sort_idx_inhib] = sort(auc_inhib, 'descend');

% Concatenate: activated first, then inhibited (none excluded for clarity)
auc_combined  = [auc_act_sorted;   auc_inhib_sorted];
n_act_plot    = length(auc_act_sorted);
n_inhib_plot  = length(auc_inhib_sorted);
x_combined    = 1:length(auc_combined);

figure;
hold on;

% Bar plot, colored by group
bar(1:n_act_plot,                       auc_act_sorted,   'FaceColor', [0.20 0.63 0.17], 'EdgeColor', 'none');
bar((n_act_plot+1):(n_act_plot+n_inhib_plot), auc_inhib_sorted, 'FaceColor', [0.84 0.10 0.11], 'EdgeColor', 'none');

% Dividing line between groups
xline(n_act_plot + 0.5, '--k', 'LineWidth', 1.2);

% Zero reference
yline(0, '-k', 'LineWidth', 0.8);

% Labels
xlabel('Neuron (sorted within group)');
ylabel('AUC (z-score \times s)');
title({['AUC: ' strrep(Title,'_','\_')], ...
       sprintf('Window: prestim3 to +%ds  |  Act n=%d, Inhib n=%d', ...
               Int3/sampling_rate, n_act_plot, n_inhib_plot)});

% Bracket labels above plot
ax = gca;
yl = ylim;
text(n_act_plot/2,          yl(2)*0.97, 'Activated', 'HorizontalAlignment','center', ...
     'FontWeight','bold', 'Color',[0.20 0.63 0.17], 'FontSize', 10);
text(n_act_plot + n_inhib_plot/2, yl(2)*0.97, 'Inhibited', 'HorizontalAlignment','center', ...
     'FontWeight','bold', 'Color',[0.84 0.10 0.11], 'FontSize', 10);

legend({'Activated','Inhibited'}, 'Location','northeast');
grid on;
hold off;

%% Peak Z-Score and Time to Peak During IG Infusion Window
% Window: infusion onset (prestim1) to infusion end (prestim1 + Int)
% Computed for ALL neurons; then split by category.

infusion_frames = interval1_start : interval1_end;   % sample indices

% Peak z-score (max absolute positive deflection) per neuron
peak_z_all        = max(neurons_zscored(:, infusion_frames), [], 2);   % [nNeurons x 1]

% Time to peak (seconds relative to infusion onset)
[~, peak_frame_rel] = max(neurons_zscored(:, infusion_frames), [], 2); % frame index within window
time_to_peak_all    = (peak_frame_rel - 1) / sampling_rate;            % seconds after prestim1

% Split by category
peak_z_act   = peak_z_all(act,   :);
peak_z_inhib = peak_z_all(inhib, :);
peak_z_none  = peak_z_all(none,  :);

time_to_peak_act   = time_to_peak_all(act,   :);
time_to_peak_inhib = time_to_peak_all(inhib, :);
time_to_peak_none  = time_to_peak_all(none,  :);

fprintf('\n--- Peak Z-Score During IG Infusion ---\n');
fprintf('  Activated (n=%d): mean peak z = %.2f,  mean time-to-peak = %.1f s\n', ...
    length(act),   mean(peak_z_act),   mean(time_to_peak_act));
fprintf('  Inhibited (n=%d): mean peak z = %.2f,  mean time-to-peak = %.1f s\n', ...
    length(inhib), mean(peak_z_inhib), mean(time_to_peak_inhib));
fprintf('  No Change (n=%d): mean peak z = %.2f,  mean time-to-peak = %.1f s\n', ...
    length(none),  mean(peak_z_none),  mean(time_to_peak_none));

% --- Figure: Time-to-Peak histogram per category ---
figure;
hold on;
histogram(time_to_peak_act   / 60, 'BinWidth', 1, 'FaceColor', [0.20 0.63 0.17], ...
    'FaceAlpha', 0.6, 'DisplayName', 'Activated');
histogram(time_to_peak_inhib / 60, 'BinWidth', 1, 'FaceColor', [0.84 0.10 0.11], ...
    'FaceAlpha', 0.6, 'DisplayName', 'Inhibited');
xlabel('Time to Peak Z (min after infusion onset)');
ylabel('Number of Neurons');
title({['Time to Peak Z During IG Infusion: ' strrep(Title,'_','\_')], ...
       sprintf('Window: 0–%d min', Int/sampling_rate/60)});
legend('Location','northeast');
grid on; hold off;

%% 5-Minute Epoch Mean Z-Score (Across Full Recording)
% Divides the entire recording into non-overlapping 5-min bins.
% Each neuron gets a mean z-score for each epoch -> compare Ensure vs Saline
% at matched time points using the exported epoch_means matrix.

epoch_dur_s      = 300;                            % 5 min in seconds
epoch_dur_samp   = epoch_dur_s * sampling_rate;    % samples per epoch
total_samp       = size(neurons_zscored, 2);
n_epochs         = floor(total_samp / epoch_dur_samp);

% Epoch means: [neurons x epochs]
epoch_means = NaN(size(neurons_zscored, 1), n_epochs);
for ep = 1:n_epochs
    ep_start = (ep - 1) * epoch_dur_samp + 1;
    ep_end   =  ep      * epoch_dur_samp;
    epoch_means(:, ep) = mean(neurons_zscored(:, ep_start:ep_end), 2);
end

% Epoch time axis (center of each bin, in minutes from recording start)
epoch_centers_samp = ((1:n_epochs) - 0.5) * epoch_dur_samp;
epoch_time_min     = epoch_centers_samp / sampling_rate / 60;

% Mark which epochs overlap the IG infusion window
infusion_start_min = prestim1 / sampling_rate / 60;
infusion_end_min   = (prestim1 + Int) / sampling_rate / 60;

% Split epoch_means by category
epoch_means_act   = epoch_means(act,   :);   % [n_act   x n_epochs]
epoch_means_inhib = epoch_means(inhib, :);   % [n_inhib x n_epochs]
epoch_means_none  = epoch_means(none,  :);   % [n_none  x n_epochs]

% Population mean ± SEM per epoch
epoch_pop_mean = mean(epoch_means, 1, 'omitnan');
epoch_pop_sem  = std(epoch_means,  0, 1, 'omitnan') / sqrt(size(epoch_means, 1));

epoch_act_mean = mean(epoch_means_act,   1, 'omitnan');
epoch_act_sem  = std(epoch_means_act,    0, 1, 'omitnan') / sqrt(size(epoch_means_act, 1));

epoch_inhib_mean = mean(epoch_means_inhib, 1, 'omitnan');
epoch_inhib_sem  = std(epoch_means_inhib,  0, 1, 'omitnan') / sqrt(size(epoch_means_inhib, 1));

fprintf('\n--- 5-Min Epoch Means: %d epochs of %d s ---\n', n_epochs, epoch_dur_s);
fprintf('  Epoch time axis (min): '); fprintf('%.1f  ', epoch_time_min); fprintf('\n');

% --- Figure: Population mean z per epoch (all neurons) ---
figure;
hold on;
% Shade infusion window
patch([infusion_start_min infusion_end_min infusion_end_min infusion_start_min], ...
      [-10 -10 10 10], [0.8 0.9 1], 'EdgeColor','none', 'FaceAlpha',0.4, ...
      'DisplayName','IG infusion');

% Plot mean ± SEM for each category
errorbar(epoch_time_min, epoch_act_mean,   epoch_act_sem,   '-o', ...
    'Color',[0.20 0.63 0.17], 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','Activated');
errorbar(epoch_time_min, epoch_inhib_mean, epoch_inhib_sem, '-s', ...
    'Color',[0.84 0.10 0.11], 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','Inhibited');
errorbar(epoch_time_min, epoch_pop_mean,   epoch_pop_sem,   '-^', ...
    'Color',[0.3 0.3 0.3], 'LineWidth',1.2, 'MarkerSize',4, 'DisplayName','All neurons');

xline(infusion_start_min, '--b', 'Infusion ON',  'LabelVerticalAlignment','bottom');
xline(infusion_end_min,   '--r', 'Infusion OFF', 'LabelVerticalAlignment','bottom');

xlabel('Time (min from recording start)');
ylabel('Mean Z-Score (\pm SEM)');
title({['5-Min Epoch Mean Z: ' strrep(Title,'_','\_')]});
legend('Location','northwest');
ylim([-3 8]); grid on; hold off;

%% Store the neuron traces based on reaction to stim

save(strcat(Title, '_quant'))