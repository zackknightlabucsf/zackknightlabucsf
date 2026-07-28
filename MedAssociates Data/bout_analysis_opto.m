clear all; close all;

% this code sets time zero as the time of lickometer access (which it will ask you), if annimal licks
% beforehand it will only count after access (this is rare) but will flag
% the animals in the excel sheet

% === SETTINGS ===
time_windows = [0 1800; 0 3600; 0 7200]; % 0–30,0-60, 0–120, and all
win_labels = {'0to30','0to60','All'};  % must match rows in time_windows
assert(numel(win_labels) == size(time_windows,1), 'win_labels must match time_windows rows');

ili = 5; %If the time gap between two licks exceeds the ILI threshold , the code treats that as the end of a bout and starts a new one: 
min_bout_length = 3;
bin_size = 1; %for cumulative licks
time_bins = 0:bin_size:7200;
imi = 5*60; %inter-meal interval
min_meal_length = 5;

% === ASK FOR ACCESS TIME ===
access_min = inputdlg('When did the animal get access to the lickometer (min)?', ...
                      'Lickometer Access Time', [1 50], {'10'});
if isempty(access_min), error('❌ No access time provided.'); end
access_time = str2double(access_min{1}) * 60; 
disp(['✅ Using access time: ' num2str(access_time) ' seconds']);

% === SELECT FOLDER ===
folder = uigetdir();
fileList = dir(fullfile(folder, '*.txt'));

% === SORT FILES BY SUBJECT NUMBER ===
subjectNums = nan(length(fileList),1);
for i = 1:length(fileList)
    fname = fileList(i).name;
    token = regexp(fname, 'Subject\s*[A-Za-z]*?(\d+)', 'tokens');
    if ~isempty(token), subjectNums(i) = str2double(token{1}{1}); end
end
[~, sort_idx] = sort(subjectNums, 'ascend', 'MissingPlacement', 'last');
fileList = fileList(sort_idx);

% === CREATE RESULTS FOLDER ===
[~, folder_name] = fileparts(folder);
results_dir = fullfile(folder, ['results_' folder_name]);
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

% === COLUMN HEADERS ===
metric_names = { ...
    'Total_Licks', ...
    'Number_Bouts', ...
    'Mean_Bout_Size', ...
    'Mean_Bout_Dur_s', ...
    'Median_Bout_Dur_s', ...
    'Number_Meals', ...
    'Mean_Meal_Size', ...
    'Mean_Meal_Dur_s', ...
    'Median_Meal_Dur_s'};

headers = {'Subject'};
for w = 1:numel(win_labels)
    for m = 1:numel(metric_names)
        headers{end+1} = sprintf('%s_%s', metric_names{m}, win_labels{w});
    end
end

headers = [headers, {'First_Meal_Size', 'First_Meal_Dur_s'}];
headers = [headers, {'Latency_to_Lick','Flag_Early_Licks'}];

results = {};
all_cum_licks = {};
subject_labels = {};
animal_ids = {};
flag_list = {};
trial_timestamps = {};  
protected_zero_rows = [];

nCols = numel(headers);
% helper to create an empty row with default values
empty_row = repmat({NaN}, 1, nCols);
empty_row{1} = '';  % subject column default
nW = size(time_windows,1);
nM = numel(metric_names);   % 9 metrics per window

% === MAIN LOOP ===
for f = 1:length(fileList)
    filename = fileList(f).name;
    file_stamp = regexp(filename, '\d{4}-\d{2}-\d{2}_\d{2}h\d{2}m', 'match', 'once');
        if isempty(file_stamp)
            file_stamp = sprintf('file%03d', f); % fallback
        end
    fullpath = fullfile(folder, filename);
    filetext = fileread(fullpath);
    lines = splitlines(filetext);

    subj_name_match = regexp(filename, 'Subject\s*(\w+)', 'tokens');
    if isempty(subj_name_match), continue; end
    subj_name = upper(subj_name_match{1}{1});

    subject_lines = find(contains(lines, 'Subject'));
    for s = 1:length(subject_lines)
        start_idx = subject_lines(s);
        if s < length(subject_lines)
            end_idx = subject_lines(s+1) - 1;
        else
            end_idx = numel(lines);
        end
        block = lines(start_idx:end_idx);

        for lick_type = {'F', 'B'}
            lick_label = lick_type{1};
            section_start = find(contains(block, [lick_label ':']), 1, 'first');
            if isempty(section_start), continue; end

            lick_times = [];
            section_end = length(block);
            for i = section_start+1:length(block)
                line = strtrim(block{i});
                if contains(line, {'F:', 'B:', 'Subject:', 'Start Date:'})
                    section_end = i - 1;
                    break;
                end
            end
            for i = section_start+1:section_end
                line = strtrim(block{i});
                if contains(line, ':')
                    parts = strsplit(line, ':');
                    if length(parts) > 1
                        numbers = sscanf(parts{2}, '%f');
                        lick_times = [lick_times; numbers(:)];
                    end
                end
            end

            full_label = [char(subj_name) '_' lower(lick_label)];
            animal_ids{end+1} = subj_name;

            if isempty(lick_times)
                row = empty_row;            % 1 x nCols cell
                row{1} = full_label;        % Subject
                                
                % Fill each window's 9 metrics: Total=0, counts=0, sizes/durations=NaN
                for w = 1:nW
                    base = 2 + (w-1)*nM;  % start col for this window
                    row(base:base+nM-1) = {0, 0, NaN, NaN, NaN, 0, NaN, NaN, NaN};
                end

                row{end-3} = NaN;   % First_Meal_Size
                row{end-2} = NaN;   % First_Meal_Dur_s
                row{end-1} = NaN;   % Latency_to_Lick
                row{end}   = 'No';  % Flag_Early_Licks
                results = [results; row];

                all_cum_licks{end+1} = zeros(size(time_bins));
                subject_labels{end+1} = full_label;
                trial_timestamps(end+1,1) = {file_stamp};
                flag_list{end+1} = 'No';
                protected_zero_rows(end+1,1) = false;
                continue;
            end

            % === Filter pre-access licks ===
            had_early = any(lick_times < access_time);
            lick_times = lick_times(lick_times >= access_time);
            lick_times = lick_times - access_time; 

            if isempty(lick_times)
                row = empty_row;
                row{1} = full_label;
            
                for w = 1:nW
                    base = 2 + (w-1)*nM;
                    row(base:base+nM-1) = {0, 0, NaN, NaN, NaN, 0, NaN, NaN, NaN};
                end
            
                row{end-1} = NaN;                        % Latency_to_Lick
                row{end}   = tern(had_early,'Yes','No'); % Flag_Early_Licks
            
                results = [results; row];
            
                all_cum_licks{end+1} = zeros(size(time_bins));
                subject_labels{end+1} = full_label;
                trial_timestamps(end+1,1) = {file_stamp};
                flag_list{end+1} = tern(had_early,'Yes','No');
                protected_zero_rows(end+1,1) = false;
                continue;
            end

            % === Latency calculation ===
            latency_to_lick = lick_times(1); % first lick time after access

            Data = [lick_times, ones(size(lick_times))];
            cum_licks = sum(lick_times <= time_bins, 1);
            all_cum_licks{end+1} = cum_licks;
            subject_labels{end+1} = full_label;
            trial_timestamps(end+1,1) = {file_stamp};
            flag_list{end+1} = tern(had_early,'Yes','No');
            protected_zero_rows(end+1,1) = false;

            % === Compute metrics ===
            row_result = {full_label};
            for w = 1:size(time_windows, 1)
                t1 = time_windows(w, 1);
                t2 = time_windows(w, 2);
                Data_T = Data(Data(:,1) > t1 & Data(:,1) <= t2, :);
                total_licks = size(Data_T, 1);

                if isempty(Data_T)
                    row_result = [row_result, {0, 0, NaN, NaN, NaN, 0, NaN, NaN, NaN}];
                    continue;
                end

                time_diffs = [Inf; diff(Data_T(:,1))];
                bouts = [];
                bout_start = [];
                for i = 1:numel(time_diffs)
                    if isempty(bout_start)
                        bout_start = i;
                    elseif time_diffs(i) > ili
                        bout_size = i - bout_start;
                        if bout_size >= min_bout_length
                            bouts(end+1,:) = [bout_start, i-1];
                        end
                        bout_start = i;
                    end
                end
                if ~isempty(bout_start)
                    bout_size = numel(Data_T(:,1)) - bout_start + 1;
                    if bout_size >= min_bout_length
                        bouts(end+1,:) = [bout_start, numel(Data_T(:,1))];
                    end
                end

                num_bouts = size(bouts,1);
                bout_sizes = diff(bouts,1,2) + 1;
                mean_bout_size = mean(bout_sizes,'omitnan');

                if num_bouts > 0
                    bout_durations = Data_T(bouts(:,2),1) - Data_T(bouts(:,1),1);
                    mean_bout_duration   = mean(bout_durations,'omitnan');
                    median_bout_duration = median(bout_durations,'omitnan');
                else
                    mean_bout_duration   = NaN;
                    median_bout_duration = NaN;
                end

                     
                meals = [];
                meal_start = [];
                for i2 = 1:numel(time_diffs)
                    if isempty(meal_start)
                        meal_start = i2;
                    elseif time_diffs(i2) > imi
                       meal_size = i2 - meal_start;
                        if meal_size >= min_meal_length
                            meals(end+1,:) = [meal_start, i2-1];
                        end
                        meal_start = i2;
                    end
                end
                if ~isempty(meal_start)
                    meal_size = numel(Data_T(:,1)) - meal_start + 1;
                    if meal_size >= min_meal_length
                        meals(end+1,:) = [meal_start, numel(Data_T(:,1))];
                    end
                end

                num_meals = size(meals,1);
                meal_sizes = diff(meals,1,2) + 1;
                mean_meal_size = mean(meal_sizes,'omitnan');

                if num_meals > 0
                    meal_durations = Data_T(meals(:,2),1) - Data_T(meals(:,1),1);
                    mean_meal_duration   = mean(meal_durations,'omitnan');
                    median_meal_duration = median(meal_durations,'omitnan');
                else
                    mean_meal_duration   = NaN;
                    median_meal_duration = NaN;
                end
                                            
                
                row_result = [row_result, {total_licks, num_bouts, mean_bout_size, mean_bout_duration, median_bout_duration, ...
                num_meals, mean_meal_size, mean_meal_duration, median_meal_duration}];
            end

            time_diffs_all = [Inf; diff(Data(:,1))];
            meals_all = [];
            meal_start_all = [];
            for i3 = 1:numel(time_diffs_all)
                if isempty(meal_start_all)
                    meal_start_all = i3;
                elseif time_diffs_all(i3) > imi
                    ms = i3 - meal_start_all;
                    if ms >= min_meal_length
                        meals_all(end+1,:) = [meal_start_all, i3-1];
                    end
                    meal_start_all = i3;
                end
            end
            if ~isempty(meal_start_all)
                ms = numel(Data(:,1)) - meal_start_all + 1;
                if ms >= min_meal_length
                    meals_all(end+1,:) = [meal_start_all, numel(Data(:,1))];
                end
            end

            if ~isempty(meals_all)
                first_meal_size = meals_all(1,2) - meals_all(1,1) + 1;
                first_meal_dur  = Data(meals_all(1,2), 1) - Data(meals_all(1,1), 1);
            else
                first_meal_size = NaN;
                first_meal_dur  = NaN;
            end

            results = [results; [row_result, {first_meal_size, first_meal_dur, latency_to_lick, tern(had_early,'Yes','No')}]];
        end
    end
end

% === TABLE ===
results_table = cell2table(results, 'VariableNames', headers);

% === FILTERING RULES ===
base_ids = regexprep(results_table.Subject, '_(f|b)$', '');
unique_animals = unique(base_ids);
mask_keep = true(height(results_table), 1);
protected_zero_rows = false(height(results_table),1);

for a = 1:length(unique_animals)
    aid = unique_animals{a};
    idx = find(strcmp(base_ids, aid));
    if numel(idx) == 2
        iF = idx(contains(results_table.Subject(idx), '_f'));
        iB = idx(contains(results_table.Subject(idx), '_b'));
        TF = results_table.Total_Licks_All(iF);
        TB = results_table.Total_Licks_All(iB);

        % F=0 & B=0 → keep F, drop B, protect F
        if TF==0 && TB==0
            mask_keep(iB) = false;
            protected_zero_rows(iF) = true;
        % F=0 & B>0 → keep B only
        elseif TF==0 && TB>0
            mask_keep(iF) = false;
        % B=0 & F>0 → keep F only
        elseif TB==0 && TF>0
            mask_keep(iB) = false;
        end
    end
end

% === APPLY FILTERS ===
results_table = results_table(mask_keep,:);
subject_labels = subject_labels(mask_keep);
trial_timestamps = trial_timestamps(mask_keep);
all_cum_licks = all_cum_licks(mask_keep);
protected_zero_rows = protected_zero_rows(mask_keep);

% === REMOVE ZERO ROWS UNLESS PROTECTED ===
last_metric_col = 1 + nW*nM;   % Subject is col 1
data_numeric = results_table{:, 2:last_metric_col};
rows_all_zero = all(data_numeric==0 | isnan(data_numeric), 2);
rows_to_remove = rows_all_zero & ~protected_zero_rows;
results_table(rows_to_remove,:) = [];
subject_labels(rows_to_remove) = [];
trial_timestamps(rows_to_remove) = [];
all_cum_licks(rows_to_remove) = [];

% === SORTING ===
S = results_table.Subject;

subj_num = nan(height(results_table),1);
is_f = false(height(results_table),1);

for i = 1:height(results_table)
    % Capture the number after CGRP or CGPR
    tok = regexp(S{i}, 'CG(?:RP|PR)(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        subj_num(i) = str2double(tok{1});
    end

    is_f(i) = endsWith(S{i}, '_f');
end

[~, order] = sortrows([subj_num, ~is_f], [1 2]);   % number asc, f before b
results_table  = results_table(order,:);
subject_labels = subject_labels(order);
trial_timestamps = trial_timestamps(order);
all_cum_licks  = all_cum_licks(order);

% === SAVE EXCEL SUMMARY (CLEAN) ===
excel_name = ['lick_analysis_results_' folder_name '.xlsx'];
outname = fullfile(results_dir, excel_name);
if exist(outname,'file'), delete(outname); end

writetable(results_table, outname);

% === ADD PER-ANIMAL CUMULATIVE TABS ===
sheetCounts = containers.Map('KeyType','char','ValueType','double');

for i = 1:length(subject_labels)
    subj_name = upper(subject_labels{i});
    stamp = trial_timestamps{i};

    cum_trace = all_cum_licks{i}';
    output_table = table(time_bins', cum_trace, 'VariableNames', {'Time_s','Cumulative_Licks'});

    % Build base sheet name
    base = sprintf('%s_%s', subj_name, stamp);
    base = regexprep(base, '[^\w]', '_');   % Excel-safe

    % Truncate to Excel 31-char limit (leave room for suffix if needed)
    if length(base) > 31
        base = base(1:31);
    end

    % Ensure uniqueness (in case truncation causes collisions)
    if isKey(sheetCounts, base)
        sheetCounts(base) = sheetCounts(base) + 1;
    else
        sheetCounts(base) = 1;
    end
    k = sheetCounts(base);

    if k == 1
        sheet_name = base;
    else
        suffix = sprintf('_r%d', k);
        maxBase = 31 - length(suffix);
        sheet_name = [base(1:min(end,maxBase)) suffix];
    end

    writetable(output_table, outname, 'Sheet', sheet_name);
end

% === PLOTS ===
figure; hold on;
for i = 1:length(all_cum_licks)
    plot(time_bins, all_cum_licks{i}, 'LineWidth', 1.2);
end
xlabel('Time (s)'); ylabel('Cumulative Licks');
title('Cumulative Licks Over Time by Subject and Lickometer');
legend(subject_labels, 'Interpreter','none','Location','eastoutside');
xlim([0 7200]); grid on;
pdf_name = fullfile(results_dir, 'all_individual_cumulative_licks.pdf');
set(gcf,'PaperPositionMode','auto','PaperOrientation','landscape');
print(gcf,pdf_name,'-dpdf','-bestfit'); close;

% === POPULATION AVERAGE ===
cum_matrix = cell2mat(all_cum_licks');
mean_cum = mean(cum_matrix,1,'omitnan');
sem_cum = std(cum_matrix,0,1,'omitnan') ./ sqrt(size(cum_matrix,1));
figure; hold on;
fill([time_bins, fliplr(time_bins)], ...
     [mean_cum - sem_cum, fliplr(mean_cum + sem_cum)], ...
     [0.85 0.85 0.85], 'EdgeColor','none');
plot(time_bins, mean_cum, 'k', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Average Cumulative Licks');
title(sprintf('Population Average Cumulative Lick Distribution (N=%d)', size(cum_matrix,1)));
grid on; xlim([0 7200]); hold off;
pdf_pop = fullfile(results_dir,'population_average_cumulative_licks.pdf');
set(gcf,'PaperPositionMode','auto','PaperOrientation','landscape');
print(gcf,pdf_pop,'-dpdf','-bestfit'); close;

disp(['🎉 Analysis complete. Excel + plots saved in: ' results_dir]);

% === Helper ===
function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end
