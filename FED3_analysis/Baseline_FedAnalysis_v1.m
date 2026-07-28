%%FED3 Analysis for baseline intake during a user-defined set of hours at
%%the start of the dark cycle

% Reads all *.CSV files in DATA_DIR, matches sessions to the decoder
% spreadsheet, and computes within-window feeding metrics.
%
% OUTPUTS
%   FED_summary_<N>h.csv          – one row per decoder session
%   histograms/IPI_hist_<N>h.jpg  – pooled log-scale IPI histogram
%
% USAGE
%   Edit the USER PARAMETERS block, then run section-by-section (Ctrl+Enter)
%   or press Run.  All functions are at the bottom of this file.


clear; close all; clc;

data_dir      = 'Pathway to data';
decoder_path  = 'Pathway to decoder';

DARK_HOUR         = 17;    % dark-cycle onset (24-h clock, e.g. 17 = 5 PM)
WINDOW_HOURS      = 2;    % hours to analyze after dark onset (2–4 for DCZ)

MEAL_GAP_MIN      = 5;   % temporary guess; you'll tune using your histogram
MIN_MEAL_PELLETS = 3;   % 3 pellets = 0.06 g (adjustable)

% FED3 free-feeding baseline: ~1–3 Motor_Turns per pellet.
% >10 almost always indicates a blockage or repeated dispensing attempt.
JAM_TURNS_THRESH = 20;   % tune this

%% ── SECTION 1: Load & concatenate all CSVs ───────────────────────────────
fprintf('\n[1] Loading CSV files from: %s\n', data_dir);

csvFiles = dir(fullfile(data_dir, '*.CSV'));
FED  = struct();

allTables = cell(numel(csvFiles), 1);

% datetimeHeader = "MM:DD:YYYY hh:mm:ss";
% dtVar = matlab.lang.makeValidName(datetimeHeader);  % likely MM_DD_YYYYHh_mm_ss

for f = 1:numel(csvFiles)
    fpath = fullfile(data_dir, csvFiles(f).name);
    fprintf('    Reading: %s\n', csvFiles(f).name);

    try
        opts = detectImportOptions(fpath, 'PreserveVariableNames', true);

        % Timestamp column: read as string
        tsName = opts.VariableNames{1};
        opts   = setvartype(opts, tsName, 'char');

        % Retrieval_Time contains "Timed_out" strings – must be char/string,
        % not numeric
        if ismember('Retrieval_Time', opts.VariableNames)
            opts = setvartype(opts, 'Retrieval_Time', 'char');
        end

        T = readtable(fpath, opts);
        T.SourceFile = repmat({csvFiles(f).name}, height(T), 1);
        allTables{f} = T;

    catch ME
        warning('Could not read %s: %s', csvFiles(f).name, ME.message);
    end
end

% Remove empty cells (failed reads)
allTables = allTables(~cellfun(@isempty, allTables));
if isempty(allTables)
    error('No CSV files could be read.');
end

% Concatenate. synchronizeVars = false keeps all columns even if schemas differ.
ALL = vertcat(allTables{:});
fprintf('    Total rows loaded: %d\n', height(ALL));


% ── Parse timestamps ───────────────────────────────────────────────────────
% FED3 format: "M/D/YYYY HH:MM:SS"  (month and day may be 1 or 2 digits)
tsName  = ALL.Properties.VariableNames{1};
rawTS   = strtrim(string(ALL.(tsName)));

ALL.Timestamp = parseFedTimestamps(rawTS);

% Drop rows with unparseable timestamps
badTS = isnat(ALL.Timestamp);
if any(badTS)
    warning('%d rows had unparseable timestamps and were dropped.', sum(badTS));
    ALL = ALL(~badTS, :);
end

% Sort chronologically
[~, ord] = sort(ALL.Timestamp);
ALL = ALL(ord, :);

% Coerce key numeric columns (Retrieval_Time stays as char/string)
ALL = coerceNumeric(ALL, {'Device_Number','Motor_Turns','InterPelletInterval','Pellet_Count'});

devList = unique(ALL.Device_Number);
devList = devList(~isnan(devList));
fprintf('    Devices found: %s\n', num2str(devList', '%g '));


%% ── SECTION 2: Load decoder ───────────────────────────────────────────────
D = readtable(decoder_path, 'PreserveVariableNames', true);
D.Properties.VariableNames = strtrim(D.Properties.VariableNames);

D.Animal        = strtrim(string(D.Animal));
D.Condition     = upper(strtrim(string(D.Condition)));
D.Device_Number = double(D.Device_Number);

% NightDate: ensure it is a datetime with time zeroed
if ~isdatetime(D.NightDate)
    D.NightDate = datetime(D.NightDate);
end
D.NightDate = dateshift(D.NightDate, 'start', 'day');


% ── Build WinStart per row ─────────────────────────────────────────────────
% WindowStart in the decoder is OPTIONAL.  When absent (NaN / missing),
% default to DARK_HOUR on NightDate.
D.WinStart = NaT(height(D), 1);

for r = 1:height(D)
    nd = D.NightDate(r);    % midnight of the night

    if ~ismember('WindowStart', D.Properties.VariableNames)
        D.WinStart(r) = nd + hours(DARK_HOUR);
        continue
    end

    ws = D.WindowStart(r);  % may be duration, datetime, numeric, or string

    if isdatetime(ws) && ~isnat(ws)
        % Came in as a full datetime – keep only the time-of-day part
        D.WinStart(r) = nd + timeofday(ws);

    elseif isduration(ws) && ~isnan(hours(ws))
        D.WinStart(r) = nd + ws;

    elseif isnumeric(ws) && ~isnan(ws)
        % Excel fractional day (e.g. 0.708333 = 17:00)
        D.WinStart(r) = nd + days(ws);

    else
        % Missing / NaN / empty string → default to dark-cycle onset
        D.WinStart(r) = nd + hours(DARK_HOUR);
    end
end

D.WinEnd = D.WinStart + hours(WINDOW_HOURS);

disp(D(:, {'Animal','Device_Number','NightDate','Condition','WinStart','WinEnd'}));

%% ── SECTION 3: Per-session analysis ──────────────────────────────────────
outRows = {};

% IPI pools per condition (for histogram)
condList   = unique(D.Condition);
IPI_pool   = struct();
for c = 1:numel(condList)
    IPI_pool.(matlab.lang.makeValidName(condList(c))) = [];
end

for r = 1:height(D)

    animal = D.Animal(r);
    dev    = D.Device_Number(r);
    cond   = D.Condition(r);
    night  = D.NightDate(r);
    winS   = D.WinStart(r);
    winE   = D.WinEnd(r);


    % ── Subset to this device ─────────────────────────────────────────────
    devMask = ALL.Device_Number == dev;
    if ~any(devMask)
        fprintf('    WARNING: No data for device %g (%s) – skipping.\n', dev, animal);
        continue
    end
    Tdev = ALL(devMask, :);

    % ── Subset to analysis window ─────────────────────────────────────────
    winMask = (Tdev.Timestamp >= winS) & (Tdev.Timestamp < winE);
    Twin = Tdev(winMask, :);

    if isempty(Twin)
        fprintf('    WARNING: No data in window for device %g (%s) on %s – skipping.\n', ...
            dev, animal, datestr(night, 'yyyy-mm-dd'));
        continue
    end

    % ── Pellet rows only ──────────────────────────────────────────────────
    % IPI and Retrieval_Time are only recorded on 'Pellet' events.
    isPellet = strcmp(string(Twin.Event), 'Pellet');
    Tpel     = Twin(isPellet, :);
    totalPellets = height(Tpel);

    % ── Latency to first pellet (minutes) ────────────────────────
    % Time from window start (dark-cycle onset) to the first Pellet event.
    if isempty(Tpel)
        latencyFirstPellet_min = NaN;
    else
        latencyFirstPellet_min = seconds(Tpel.Timestamp(1) - winS) / 60;
    end

    % ── IPI (seconds) ─────────────────────────────────────────────────────
    ipi_raw = double(Tpel.InterPelletInterval);   % NaN where string/"Timed_out"
    if ~isempty(ipi_raw)
        ipi_raw(1) = NaN;   % zero out the cross-boundary IPI
    end
    ipi_valid = ipi_raw(isfinite(ipi_raw) & ipi_raw > 0);

    % Accumulate for histogram
    ck = matlab.lang.makeValidName(cond);
    IPI_pool.(ck) = [IPI_pool.(ck); ipi_valid(:)];


    % ── Meals ─────────────────────────────────────────────────────────────
    M = mealsFromIPI(ipi_valid, MEAL_GAP_MIN, MIN_MEAL_PELLETS);

    % ── Jams ──────────────────────────────────────────────────────────────
    J = detectJams(Tpel, JAM_TURNS_THRESH);

    % ── Build output row ──────────────────────────────────────────────────
    row = { ...
       animal, ...                             % 1  Animal
        dev, ...                                % 2  Device
        night, ...                              % 3  Date
        cond, ...                               % 4  Condition
        winS, ...                               % 5  WinStart
        winE, ...                               % 6  WinEnd
        totalPellets, ...                       % 8  TotalPellets
        latencyFirstPellet_min, ...             % 9  LatencyFirstPellet_min
        nanmedian_safe(ipi_valid), ...          % 11 Median_IPI_sec
        nanmean_safe(ipi_valid), ...            % 10 Mean_IPI_sec
        numel(ipi_valid), ...                   % 12 N_IPI
        M.nRawMeals, ...                        % 13 N_RawMeals
        nanmean_safe(M.rawSizes), ...           % 14 Mean_RawMealSize
        nanmean_safe(M.rawDurs)/60, ...         % 15 Mean_RawMealDur_min
        MIN_MEAL_PELLETS, ...                   % 16 MinMealPellets
        M.nQualMeals, ...                       % 17 N_QualMeals
        nanmean_safe(M.qualSizes), ...          % 18 Mean_QualMealSize
        nanmean_safe(M.qualDurs)/60, ...        % 19 Mean_QualMealDur_min
        M.firstSize, ...                        % 20 FirstMealSize
        M.firstDur/60, ...                      % 21 FirstMealDur_min
        numel(M.subseqSizes), ...               % 22 N_SubseqMeals
        nanmean_safe(M.subseqSizes), ...        % 23 Mean_SubseqMealSize
        nanmean_safe(M.subseqDurs)/60, ...      % 24 Mean_SubseqMealDur_min
        J.jamCount, ...                         % 25 JamCount
        J.jamResolved, ...                      % 26 JamResolved
        J.jamUnresolved, ...                    % 27 JamUnresolved
        J.totalJamTime_sec/60, ...              % 28 TotalJamTime_min
        J.jamTimestamps, ...                    % 29 JamTimestamps  (string)
        J.jamTurns, ...                         % 30 JamMotorTurns  (string)
        strjoin(arrayfun(@(x)sprintf('%.2f',x), ipi_valid(:)', 'UniformOutput',false), ';'), ...   % 31 IPI_sec
        strjoin(arrayfun(@(x)sprintf('%d',x),   M.qualSizes(:)', 'UniformOutput',false), ';'), ... % 32 QualMealSizes
        strjoin(arrayfun(@(x)sprintf('%.2f',x), M.qualDurs(:)'/60, 'UniformOutput',false), ';'), ...% 33 QualMealDurs_min
    };

    outRows(end+1, :) = row;

    fprintf('    %s | dev %g | %s | %s | %d pellets | %d qual-meals | %d jams\n', ...
        animal, dev, datestr(night,'yyyy-mm-dd'), cond, ...
        totalPellets, M.nQualMeals, J.jamCount);
end

if isempty(outRows)
    warning('No sessions produced output. Check that Device_Number values in the decoder match those in the CSV files.');
    return
end


%% ── SECTION 4: Build summary table & save ────────────────────────────────

varNames = { ...
    'Animal','Device','Date','Condition','WinStart','WinEnd', ...
    'TotalPellets','LatencyFirstPellet_min', ...
    'Median_IPI_sec','Mean_IPI_sec','N_IPI', ...
    'N_RawMeals','Mean_RawMealSize','Mean_RawMealDur_min', ...
    'MinMealPellets','N_QualMeals','Mean_QualMealSize','Mean_QualMealDur_min', ...
    'FirstMealSize','FirstMealDur_min', ...
    'N_SubseqMeals','Mean_SubseqMealSize','Mean_SubseqMealDur_min', ...
    'JamCount','JamResolved','JamUnresolved','TotalJamTime_min', ...
    'JamTimestamps','JamMotorTurns', ...
    'IPI_sec','QualMealSizes','QualMealDurs_min' ...
};


Summary = cell2table(outRows, 'VariableNames', varNames);

runStamp = datestr(now, 'yyyymmdd_HHMMSS');
outDir   = fullfile(data_dir, sprintf('FED_output_%dh_%s', WINDOW_HOURS, runStamp));
if ~exist(outDir, 'dir'), mkdir(outDir); end
 
outPath = fullfile(outDir, sprintf('FED_summary_%dh.csv', WINDOW_HOURS));
writetable(Summary, outPath);
fprintf('\n[4] Summary saved → %s\n', outPath);

%% ── SECTION 5: IPI histogram ──────────────────────────────────────────────

histDir = fullfile(data_dir, 'histograms');
if ~exist(histDir, 'dir'), mkdir(histDir); end

step     = 0.02;
logEdges = (-2 : step : 3) + step/2;
edges    = 10.^logEdges;
majTicks = [0.01 0.1 1 10 100 1000];

colors = lines(numel(condList));
fig = figure('Color','w','Position',[100 100 800 400]);
hold on;

for c = 1:numel(condList)
    ck   = matlab.lang.makeValidName(condList(c));
    if ~isfield(IPI_pool, ck), continue; end

    dMin = IPI_pool.(ck) / 60;          % seconds → minutes
    dMin = dMin(isfinite(dMin) & dMin > 0);
    if isempty(dMin), continue; end

    [N, E] = histcounts(dMin, edges);
    stairs(E(1:end-1), N, 'Color', colors(c,:), 'LineWidth', 2, ...
        'DisplayName', condList(c));
end

% Vertical line at meal-gap threshold
xline(MEAL_GAP_MIN, '--k', sprintf('Meal gap (%g min)', MEAL_GAP_MIN), ...
    'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom');

ax = gca;
ax.XScale = 'log';
ax.XTick  = majTicks;
xticklabels(ax, cellfun(@num2str, num2cell(majTicks), 'UniformOutput', false));
xlabel('Inter-pellet interval (min)');
ylabel('Pellet count');
title(sprintf('IPI distribution | %dh post-dark-cycle onset', WINDOW_HOURS), ...
    'Interpreter','none');
legend('Location','best');
box off;
ax.XMinorTick = 'off';

histPath = fullfile(histDir, sprintf('IPI_hist_%dh_by_condition.jpg', WINDOW_HOURS));
saveas(fig, histPath);
fprintf('[5] IPI histogram saved → %s\n', histPath);


%% ── SECTION 6: Quick condition-level summary ─────────────────────────────

fprintf('\n[6] Condition means:\n');
numericCols = {'TotalPellets','N_QualMeals','Mean_QualMealSize', ...
    'Mean_QualMealDur_min','FirstMealSize','FirstMealDur_min', ...
    'JamCount','TotalJamTime_min'};
numericCols = intersect(numericCols, Summary.Properties.VariableNames, 'stable');

for c = 1:numel(condList)
    idx = Summary.Condition == condList(c);
    fprintf('\n  %s (n=%d)\n', condList(c), sum(idx));
    for v = 1:numel(numericCols)
        vals = Summary.(numericCols{v})(idx);
        fprintf('    %-28s mean=%.2f  SEM=%.2f\n', ...
            numericCols{v}, mean(vals,'omitnan'), ...
            std(vals,'omitnan')/sqrt(sum(~isnan(vals))));
    end
end

fprintf('\nDone.  %d sessions in summary table.\n', height(Summary));


%%=========================================================================
%% LOCAL FUNCTIONS
%%=========================================================================

function t = parseFedTimestamps(s)
% Parse FED3 timestamp strings "M/D/YYYY HH:MM:SS" to datetime.
% Uses datetime() with InputFormat; falls back to NaT for bad rows.
t = NaT(size(s));
good = strlength(strtrim(s)) > 0 & ~ismissing(s);
if ~any(good), return; end

try
    % Vectorised parse (fast path – works when all rows share the same
    % width, e.g. two-digit month/day)
    t(good) = datetime(s(good), 'InputFormat', 'M/d/uuuu HH:mm:ss');
catch
    % Row-by-row fallback for mixed-width timestamps
    idx = find(good);
    for i = 1:numel(idx)
        try
            t(idx(i)) = datetime(s(idx(i)), 'InputFormat', 'M/d/uuuu HH:mm:ss');
        catch
            % leave as NaT
        end
    end
end
end


function T = coerceNumeric(T, colNames)
% For each column name in colNames, if it exists in T and is not already
% numeric, coerce it (non-numeric strings become NaN).
for k = 1:numel(colNames)
    cn = colNames{k};
    if ~ismember(cn, T.Properties.VariableNames), continue; end
    if isnumeric(T.(cn)), continue; end
    T.(cn) = double(string(T.(cn)));   % NaN for non-numeric strings
end
end


function M = mealsFromIPI(ipi_sec, mealGapMin, minMealPellets)
% MEALSFROIPI  Segment inter-pellet intervals into meals.
%
%   ipi_sec        – column vector, seconds (already finite & >0)
%   mealGapMin     – gap threshold in minutes
%   minMealPellets – minimum pellets to count as a "qualified" meal
%
%   Returns struct M with fields:
%     rawSizes, rawDurs         – all gap-defined meals
%     qualSizes, qualDurs       – meals meeting minMealPellets threshold
%     firstSize, firstDur       – first qualified meal
%     subseqSizes, subseqDurs   – qualified meals after the first
%     nRawMeals, nQualMeals

% Empty defaults
M = struct('rawSizes',[],'rawDurs',[],'qualSizes',[],'qualDurs',[], ...
    'firstSize',NaN,'firstDur',NaN, ...
    'subseqSizes',[],'subseqDurs',[], ...
    'nRawMeals',0,'nQualMeals',0);

ipi = ipi_sec(:);
ipi = ipi(isfinite(ipi) & ipi > 0);
if isempty(ipi), return; end

gapSec = mealGapMin * 60;
N      = numel(ipi);        % N intervals  →  N+1 pellets

% Pellet-level meal assignment (length N+1)
%   pellet 1 always starts meal 1
%   pellet k+1 starts a new meal when ipi(k) >= gapSec
newMeal  = [true; ipi >= gapSec];
mealID   = cumsum(newMeal);         % 1-indexed, length N+1
nMeals   = mealID(end);

% Raw meal sizes (pellet counts per meal)
rawSizes = accumarray(mealID, 1, [nMeals 1]);

% Meal durations = sum of within-meal IPIs
%   interval k connects pellet k to pellet k+1 → belongs to mealID(k+1)
intervalMeal = mealID(2:end);       % length N
within       = ipi < gapSec;
rawDurs      = accumarray(intervalMeal(within), ipi(within), [nMeals 1], @sum, 0);

M.rawSizes  = rawSizes;
M.rawDurs   = rawDurs;
M.nRawMeals = nMeals;

% Qualified meals
keep       = rawSizes >= minMealPellets;
qualSizes  = rawSizes(keep);
qualDurs   = rawDurs(keep);

M.qualSizes  = qualSizes;
M.qualDurs   = qualDurs;
M.nQualMeals = sum(keep);

if ~isempty(qualSizes)
    M.firstSize    = qualSizes(1);
    M.firstDur     = qualDurs(1);
    M.subseqSizes  = qualSizes(2:end);
    M.subseqDurs   = qualDurs(2:end);
end
end


function J = detectJams(Tpel, jamThresh)
% DETECTJAMS  Identify potential jam events from Motor_Turns.
%
%   Tpel       – table of Pellet-only rows within the analysis window
%   jamThresh  – Motor_Turns threshold (events >= thresh flagged)
%
%   A "jam" row is one where the motor had to spin many more times than
%   normal to dispense a pellet, indicating a blockage.  For each jam we
%   record the time until the *next* pellet (time-to-recovery, TTR).
%   Unresolved = no subsequent pellet in the window.
%
%   TotalJamTime = sum of all finite TTRs.  Use this to decide whether
%   to exclude the session (e.g. if >10% of window was jammed).

J = struct('jamCount',0,'jamResolved',0,'jamUnresolved',0, ...
    'totalJamTime_sec',0,'jamTimestamps',"", 'jamTurns',"");

if isempty(Tpel) || ~ismember('Motor_Turns', Tpel.Properties.VariableNames)
    return
end

mt      = double(Tpel.Motor_Turns);
ts      = Tpel.Timestamp;
jamIdx  = find(isfinite(mt) & mt >= jamThresh);

if isempty(jamIdx), return; end

J.jamCount = numel(jamIdx);
ttrList    = NaN(J.jamCount, 1);
tsStrs     = strings(J.jamCount, 1);
mtStrs     = strings(J.jamCount, 1);

for j = 1:J.jamCount
    i0 = jamIdx(j);
    t0 = ts(i0);

    % Next pellet (any motor-turn count) after this jam row
    future = find(ts > t0, 1, 'first');
    if ~isempty(future)
        ttrList(j) = seconds(ts(future) - t0);
    end

    tsStrs(j) = string(datestr(t0, 'mm/dd/yyyy HH:MM:SS'));
    mtStrs(j) = string(mt(i0));
end

J.jamResolved      = sum(isfinite(ttrList));
J.jamUnresolved    = sum(~isfinite(ttrList));
J.totalJamTime_sec = sum(ttrList, 'omitnan');
J.jamTimestamps    = strjoin(tsStrs, ';');
J.jamTurns         = strjoin(mtStrs, ';');
end


function v = nanmean_safe(x)
if isempty(x), v = NaN; else, v = mean(x, 'omitnan'); end
end

function v = nanmedian_safe(x)
if isempty(x), v = NaN; else, v = median(x, 'omitnan'); end
end