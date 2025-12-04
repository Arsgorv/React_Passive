function E = RP_build_epochs_passive_csv_only(datapath, trig_csv, B)
% RP_build_epochs_passive_csv_only
% Build intervalSet epochs for React Passive session when only
% Baphy/fUS csv triggers are available (no OpenEphys).
%
% Uses:
%   - trig_csv.fus_time_s    : fUS TTLs, seconds, t0 = first F
%   - trig_csv.baphy_time_s  : Baphy TTLs, seconds, t0 = first F
%   - B (from RP_parse_baphy_passive)
%
% OUTPUT E struct:
%   E.fus_blocks            : intervalSet of continuous fUS imaging blocks
%   E.trial_all             : intervalSet of all trials
%   E.trial_music/speech/ferret/silence : category-specific trial sets
%   E.stim_all              : intervalSet of all stimulus epochs
%   E.stim_music/speech/ferret/silence  : category-specific stim sets
%
% All times are in seconds on the csv time axis (t0 = first fUS TTL).

E = struct;

%% 1) fUS blocks from csv fUS TTLs (if present)

if isfield(trig_csv,'fus_time_s') && ~isempty(trig_csv.fus_time_s)
    t_fus = trig_csv.fus_time_s(:);
    if numel(t_fus) >= 2
        % Single continuous exposure block: from first to last TTL.
        exp_start = t_fus(1);
        exp_end   = t_fus(end);

        % In your tsd world you'll likely multiply by 1e4 when creating intervalSet.
        E.fus_exposure = intervalSet(exp_start, exp_end);
        E.fus_blocks   = E.fus_exposure;  % for backward compatibility

        fprintf('  [csv-only] fUS exposure: from %.2f to %.2f s\n', exp_start, exp_end);
    else
        warning('[RP_build_epochs_passive_csv_only] Not enough fUS TTLs to define exposure block.');
    end
else
    warning('[RP_build_epochs_passive_csv_only] No fus_time_s field; skipping fUS exposure.');
end

%% 2) Trial epochs from Baphy + csv Baphy TTLs

if ~isfield(trig_csv,'baphy_time_s') || isempty(trig_csv.baphy_time_s)
    error('[RP_build_epochs_passive_csv_only] trig_csv.baphy_time_s missing or empty.');
end

t_baphy = trig_csv.baphy_time_s(:);
n_b = B.n_trials;

[t_trial_onset, t_trial_offset] = detect_trial_onsets_from_baphy_ttl(t_baphy, n_b);

n_trials_detected = numel(t_trial_onset);
fprintf('  [csv-only] Trials: Baphy=%d, csv-detected=%d\n', n_b, n_trials_detected);

if n_trials_detected ~= n_b
    warning('[RP_build_epochs_passive_csv_only] Trial count mismatch: Baphy=%d, csv-detected=%d', ...
        n_b, n_trials_detected);
    n = min(n_b, n_trials_detected);
    t_trial_onset  = t_trial_onset(1:n);
    t_trial_offset = t_trial_offset(1:n);
else
    n = n_b;
end

trial_start_s = t_trial_onset(:);
trial_end_s   = t_trial_offset(:);

stim_start_s  = nan(n,1);
stim_end_s    = nan(n,1);

for it = 1:n
    t0 = t_trial_onset(it);  % csv time when trial starts

    % Purely stimulus epoch (same logic as in your RP_build_epochs_passive*)
    if ~isnan(B.t_stim_start(it)) && ~isnan(B.t_stim_stop(it))
        stim_start_s(it) = t0 + B.t_stim_start(it);
        stim_end_s(it)   = t0 + B.t_stim_stop(it);
    end
end

% drop any NaNs
bad = isnan(trial_start_s) | isnan(trial_end_s);
if any(bad)
    warning('[RP_build_epochs_passive_csv_only] Dropping %d trials with NaN start/stop.', sum(bad));
    trial_start_s(bad) = [];
    trial_end_s(bad)   = [];
    stim_start_s(bad)  = [];
    stim_end_s(bad)    = [];
    B.category(bad)    = [];
end

%% 3) Pack into intervalSets (same structure as RP_build_epochs_passive)

E.trial_all = intervalSet(trial_start_s*1e4, trial_end_s*1e4);

cat = B.category(:);
idx_music  = strcmp(cat,'music');
idx_speech = strcmp(cat,'speech');
idx_ferret = strcmp(cat,'ferret');

if ~isempty(trial_start_s(idx_music))
    E.trial_music  = intervalSet(trial_start_s(idx_music)*1e4,  trial_end_s(idx_music)*1e4);
end
if ~isempty(trial_start_s(idx_speech))
    E.trial_speech = intervalSet(trial_start_s(idx_speech)*1e4, trial_end_s(idx_speech)*1e4);
end
if ~isempty(trial_start_s(idx_ferret))
    E.trial_ferret = intervalSet(trial_start_s(idx_ferret)*1e4, trial_end_s(idx_ferret)*1e4);
end

% Stim epochs (only where defined)
good_stim = ~isnan(stim_start_s) & ~isnan(stim_end_s);
E.stim_all = intervalSet(stim_start_s(good_stim)*1e4, stim_end_s(good_stim)*1e4);

if ~isempty(trial_start_s(idx_music))
    E.stim_music  = intervalSet(stim_start_s(idx_music  & good_stim)*1e4, stim_end_s(idx_music  & good_stim)*1e4);
end
if ~isempty(trial_start_s(idx_speech))
    E.stim_speech = intervalSet(stim_start_s(idx_speech & good_stim)*1e4, stim_end_s(idx_speech & good_stim)*1e4);
end
if ~isempty(trial_start_s(idx_ferret))
    E.stim_ferret = intervalSet(stim_start_s(idx_ferret & good_stim)*1e4, stim_end_s(idx_ferret & good_stim)*1e4);
end

end
