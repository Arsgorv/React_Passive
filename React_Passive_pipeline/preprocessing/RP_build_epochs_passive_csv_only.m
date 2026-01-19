function E = RP_build_epochs_passive_csv_only(datapath, trig_csv, B)
% RP_build_epochs_passive_csv_only
% Build intervalSet epochs for React Passive session when only
% Baphy/fUS csv triggers are available (no OpenEphys).
%
% Uses:
%   - trig_csv.fus_time_s    : concatenated fUS TTLs (Pre+Exp+Post), seconds
%   - trig_csv.baphy_time_s  : Baphy TTLs, seconds
%   - B (from RP_parse_baphy_passive)
%
% OUTPUT E struct:
%   E.fus_preexp / fus_exp / fus_postexp : fUS blocks (ts units)
%   E.fus_blocks                         : union of all 3 blocks
%   E.trial_all             : intervalSet of all trials
%   E.trial_music/speech/ferret/silence : category-specific trial sets
%   E.stim_all              : intervalSet of all stimulus epochs
%   E.stim_music/speech/ferret/silence  : category-specific stim sets

% RP_build_epochs_passive_csv_only
% Build intervalSet epochs for React Passive session when only
% Baphy/fUS csv triggers are available (no OpenEphys).
%
% Uses:
%   - fUS cat_tsd (concatenated PreExp+Exp+PostExp) to define fUS epochs
%   - exp_info to get number of frames per phase
%   - trig_csv.baphy_time_s  : Baphy TTLs, seconds (for trials)
%   - trig_csv.fus_time_s    : fUS TTLs for Exp only (sanity check)
%   - B (from RP_parse_baphy_passive)
%
% OUTPUT E struct:
%   E.fus_preexp / fus_exp / fus_postexp : fUS blocks (ts units)
%   E.fus_blocks                         : union of all 3 blocks
%   E.trial_all             : intervalSet of all trials
%   E.trial_music/speech/ferret/silence : category-specific trial sets
%   E.stim_all              : intervalSet of all stimulus epochs
%   E.stim_music/speech/ferret/silence  : category-specific stim sets

E = struct;

E.fus_preexp  = intervalSet([],[]);
E.fus_exp     = intervalSet([],[]);
E.fus_postexp = intervalSet([],[]);
E.fus_blocks  = intervalSet([],[]);

%% 1) fUS blocks from cat_tsd + exp_info sizes

n_pre = []; n_exp = []; n_post = [];
fus_dir  = fullfile(datapath,'fUS');
exp_file = fullfile(fus_dir,'exp_info.mat');
if exist(exp_file,'file')
    S = load(exp_file);
    if isfield(S,'exp_info')
        exp_info = S.exp_info;
        if isfield(exp_info,'PreExp')
            n_pre  = exp_info.PreExp.size(3);
            n_exp  = exp_info.Exp.size(3);
            n_post = exp_info.PostExp.size(3);
        elseif isfield(exp_info,'size')
            n_pre  = exp_info.size{1}(3);
            n_exp  = exp_info.size{2}(3);
            n_post = exp_info.size{3}(3);
        end
    end
end

cat_file = [];
if isfolder(fus_dir)
    d = dir(fullfile(fus_dir,'RP_data*_slice_*.mat'));
    if ~isempty(d)
        cat_file = fullfile(fus_dir, d(1).name);
    end
end

if ~isempty(cat_file) && ~isempty(n_pre) && ~isempty(n_exp) && ~isempty(n_post)
    try
        S = load(cat_file,'cat_tsd');
        if ~isfield(S,'cat_tsd') || ~isfield(S.cat_tsd,'data')
            warning('[RP_build_epochs_passive_csv_only] cat_tsd.data not found in %s', cat_file);
        else
            cat_tsd = S.cat_tsd;
            r = Range(cat_tsd.data);   % ts units, one per frame
            n_frames = numel(r);
            n_total_meta = n_pre + n_exp + n_post;

            if n_frames ~= n_total_meta
                warning('[csv-fUS] Frame count in cat_tsd (%d) != expected total (%d) in %s.', ...
                    n_frames, n_total_meta, datapath);
            end

            % frame indices per phase (clip to available frames)
            idx_pre  = 1 : min(n_pre, n_frames);
            idx_exp  = [];
            idx_post = [];
            if n_frames > n_pre
                start_exp = n_pre + 1;
                idx_exp = start_exp : min(start_exp + n_exp - 1, n_frames);
            end
            if n_frames > (n_pre + n_exp)
                start_post = n_pre + n_exp + 1;
                idx_post = start_post : min(start_post + n_post - 1, n_frames);
            end

            if ~isempty(idx_pre)
                t_pre_start_ts = r(idx_pre(1));
                t_pre_end_ts   = r(idx_pre(end));
                E.fus_preexp   = intervalSet(t_pre_start_ts, t_pre_end_ts);
            end
            if ~isempty(idx_exp)
                t_exp_start_ts = r(idx_exp(1));
                t_exp_end_ts   = r(idx_exp(end));
                E.fus_exp      = intervalSet(t_exp_start_ts, t_exp_end_ts);
            end
            if ~isempty(idx_post)
                t_post_start_ts = r(idx_post(1));
                t_post_end_ts   = r(idx_post(end));
                E.fus_postexp   = intervalSet(t_post_start_ts, t_post_end_ts);
            end

            % union of all non-empty blocks
            E.fus_blocks = union(E.fus_preexp, union(E.fus_exp, E.fus_postexp));

            fprintf('  [csv-only] fUS blocks from cat_tsd+exp_info: PreExp(%d frames), Exp(%d), PostExp(%d)\n', ...
                numel(idx_pre), numel(idx_exp), numel(idx_post));

            % optional: sanity check of Exp TTLs vs n_exp
            if isfield(trig_csv,'fus_time_s') && ~isempty(trig_csv.fus_time_s)
                n_ttl_exp = numel(trig_csv.fus_time_s);
                if n_ttl_exp ~= n_exp
                    warning('[csv-fUS] Exp TTL count (%d) != Exp frames in exp_info (%d) in %s.', ...
                        n_ttl_exp, n_exp, datapath);
                end
            end
        end
    catch ME
        warning('[RP_build_epochs_passive_csv_only] Failed to load cat_tsd from %s: %s', cat_file, ME.message);
    end
else
    warning('[RP_build_epochs_passive_csv_only] Missing cat_tsd or exp_info sizes; fUS epochs will not be phase-split.');
    % fallback: if we have fus_time_s, define a single Exp block from csv TTLs
    if isfield(trig_csv,'fus_time_s') && ~isempty(trig_csv.fus_time_s)
        t_fus = sort(trig_csv.fus_time_s(:));
        if numel(t_fus) >= 2
            t_start = t_fus(1);
            t_end   = t_fus(end);
            % here time is in seconds, convert to ts
            E.fus_exp    = intervalSet(t_start*1e4, t_end*1e4);
            E.fus_blocks = E.fus_exp;
        end
    end
end

%% 2) Trial epochs from Baphy + csv Baphy TTLs

if ~isfield(trig_csv,'baphy_time_s') || isempty(trig_csv.baphy_time_s)
    error('[RP_build_epochs_passive_csv_only] trig_csv.baphy_time_s missing or empty.');
end

t_baphy = trig_csv.baphy_time_s(:);
n_b = B.n_trials;

[t_trial_onset, t_trial_offset] = detect_trial_onsets_from_baphy_ttl(datapath,t_baphy, n_b);

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
