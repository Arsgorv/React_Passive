function trig = RP_sync_triggers_passive(datapath)
% RP_sync_triggers_passive
% Extracts fUS (F-*) and baphy (B-*) trigger times from a Baphy csv file
% and returns them on a common, zero-aligned timeline.
%
% INPUT
%   datapath : session folder, e.g.
%              '~/Arsenii/React_Passive/Processed_data/Edel/20220419_1_m_C'
%
% OUTPUT (struct 'trig')
%   trig.fus_time_s      : [nF x 1] fUS trigger times in seconds, t0 = first F
%   trig.baphy_time_s    : [nB x 1] baphy trigger times in seconds, t0 = first F
%   trig.fus_time_ts     : ts object, fUS times, units = 1e-4 s (tsd convention)
%   trig.baphy_time_ts   : ts object, baphy times, units = 1e-4 s
%   trig.csv_file        : full path to the csv used
%   trig.raw_table       : original table from readtable (for debugging)
%
% Expected csv format (no header):
%   col1: time in seconds (double)
%   col2: string code 'F-xxx', 'B-xxx', ...
%
% Arsenii Goriachenkov / React_Passive

%% locate csv file

% First try: .../baphy/*.csv
csv_list = dir(fullfile(datapath, 'baphy', '*.csv'));

% If nothing found, also try .../TrigFiles/*.csv (legacy location)
if isempty(csv_list)
    csv_list = dir(fullfile(datapath, 'TrigFiles', '*.csv'));
end

if isempty(csv_list)
    error('RP_sync_triggers_passive:NoCSV', ...
        'No *.csv trigger file found in %s/baphy or %s/TrigFiles', datapath, datapath);
end

if numel(csv_list) > 1
    % If there are many, you can add a smarter selector here based on name.
    % For now, just take the first and warn.
    warning('RP_sync_triggers_passive:MultipleCSV', ...
        'Multiple csv files found, using the first one: %s', csv_list(1).name);
end

csv_file = fullfile(csv_list(1).folder, csv_list(1).name);
disp(['[RP_sync_triggers_passive] Using trigger file: ' csv_file]);

%% read table

T = readtable(csv_file, 'ReadVariableNames', false);

if size(T,2) < 2
    error('RP_sync_triggers_passive:BadFormat', ...
        'Trigger csv must have at least 2 columns: time(s), code string.');
end

time_s = T{:,1};          % numeric, seconds
code   = T{:,2};          % cellstr or string

% Make sure 'code' is a cell array of char (for older MATLAB string tools)
if isstring(code)
    code = cellstr(code);
end

%% extract fUS and Baphy triggers

isF = strncmp(code, 'F-', 2);
isB = strncmp(code, 'B-', 2);

fus_time_s   = time_s(isF);
baphy_time_s = time_s(isB);

if isempty(fus_time_s)
    error('RP_sync_triggers_passive:NoFUS', ...
        'No F- entries found in %s', csv_file);
end

% Align both streams so that the first fUS trigger is time zero
t0 = fus_time_s(1);
fus_time_s   = fus_time_s   - t0;
baphy_time_s = baphy_time_s - t0;

%% convert to ts objects (tsd convention: units = 1e-4 s)

fus_time_ts   = ts(fus_time_s * 1e4);
baphy_time_ts = ts(baphy_time_s * 1e4);

%% pack output

trig = struct;
trig.fus_time_s      = fus_time_s(:);
trig.baphy_time_s    = baphy_time_s(:);
trig.fus_time_ts     = fus_time_ts;
trig.baphy_time_ts   = baphy_time_ts;
trig.csv_file        = csv_file;
trig.raw_table       = T;

out_file = fullfile(datapath, 'RP_sync.mat');
save(out_file, 'trig');

disp(['Saved trig to: ' out_file])
disp('----------------------------------')

end
