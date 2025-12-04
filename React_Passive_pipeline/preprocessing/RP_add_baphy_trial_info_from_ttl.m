function trigOE = RP_add_baphy_trial_info_from_ttl(trigOE, trig_csv, B, datapath)
% RP_add_baphy_trial_info_from_ttl
% For Edel: use OE Baphy TTL bursts to define per-trial times:
%   - trigOE.baphy.t_s              : one onset time per trial (seconds, OE time)
%   - trigOE.baphy.trial_start_ts   : trial start in ts units (1e-4 s)
%   - trigOE.baphy.trial_stop_ts    : trial stop in ts units (1e-4 s)
% Also compares per-trial TTL counts OE vs csv (if csv is available).

if ~isfield(trigOE,'baphy') || ~isfield(trigOE.baphy,'t_raw_s') || isempty(trigOE.baphy.t_raw_s)
    warning('[RP_add_baphy_trial_info_from_ttl] No trigOE.baphy.t_raw_s for %s', datapath);
    return
end

t_oe = trigOE.baphy.t_raw_s(:);

% Segment OE Baphy TTLs into trials
[t_trial_onset, t_trial_offset, first_idx, last_idx] = ...
    detect_trial_onsets_from_baphy_ttl(datapath, t_oe, B.n_trials);

n_trials = numel(t_trial_onset);
if n_trials == 0
    warning('[RP_add_baphy_trial_info_from_ttl] No trials detected from OE Baphy TTLs in %s', datapath);
    return
end

n_per_trial_oe = last_idx - first_idx + 1;
n_mode_oe = mode(n_per_trial_oe);

fprintf('[RP_add_baphy_trial_info_from_ttl] OE Baphy TTLs: %d trials, mode TTL/trial = %d in %s\n', ...
    n_trials, n_mode_oe, datapath);

% Fill trigOE fields used by your RP_build_epochs_passive*
trigOE.baphy.t_s = t_trial_onset(:);                       % in seconds

trigOE.baphy.trial_start_ts = round(t_trial_onset(:)  * 1e4);  % ts units
trigOE.baphy.trial_stop_ts  = round(t_trial_offset(:)  * 1e4); % ts units

trigOE.baphy.n_ttl_per_trial_oe = n_per_trial_oe(:);
trigOE.baphy.n_ttl_mode_oe      = n_mode_oe;

% Optional: cross-check vs csv Baphy TTLs, if available
if nargin >= 2 && ~isempty(trig_csv) && isstruct(trig_csv) && ...
        isfield(trig_csv,'baphy_time_s') && ~isempty(trig_csv.baphy_time_s)

    t_csv = trig_csv.baphy_time_s(:);

    try
        [~, ~, first_idx_csv, last_idx_csv] = ...
            detect_trial_onsets_from_baphy_ttl(datapath, t_csv, B.n_trials);
    catch
        [~, ~, first_idx_csv, last_idx_csv] = ...
            detect_trial_onsets_from_baphy_ttl(datapath, t_csv, []);
    end

    if isempty(first_idx_csv)
        warning('[RP_add_baphy_trial_info_from_ttl] Could not segment csv Baphy TTLs into trials for %s', datapath);
        return
    end

    n_per_trial_csv = last_idx_csv - first_idx_csv + 1;
    n_mode_csv = mode(n_per_trial_csv);

    trigOE.baphy.n_ttl_per_trial_csv = n_per_trial_csv(:);
    trigOE.baphy.n_ttl_mode_csv      = n_mode_csv;

    bad_oe  = find(n_per_trial_oe  ~= n_mode_oe);
    bad_csv = find(n_per_trial_csv ~= n_mode_csv);

    if ~isempty(bad_oe)
        warning(['[RP_add_baphy_trial_info_from_ttl] OE Baphy TTLs: trials with non-modal ' ...
                 'TTL count (mode=%d) in %s: %s'], ...
            n_mode_oe, datapath, mat2str(bad_oe(:)'));
    end
    if ~isempty(bad_csv)
        warning(['[RP_add_baphy_trial_info_from_ttl] csv Baphy TTLs: trials with non-modal ' ...
                 'TTL count (mode=%d) in %s: %s'], ...
            n_mode_csv, datapath, mat2str(bad_csv(:)'));
    end

else
    % No csv for this session, nothing to validate against
end

end
