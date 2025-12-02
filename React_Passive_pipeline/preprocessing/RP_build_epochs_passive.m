function E = RP_build_epochs_passive(datapath, trigOE, B)
% RP_build_epochs_passive
% Build intervalSet epochs for React Passive session
% using:
%   - OE TTLs: trigOE.baphy (trial onsets), trigOE.fus (fUS TTLs)
%   - Baphy m-file parsed by RP_parse_baphy_passive (B)
%
% OUTPUT E struct:
%   E.fus_blocks            : intervalSet of continuous fUS imaging blocks
%
%   E.trial_all             : intervalSet of all trials
%   E.trial_music/speech/ferret/silence : category-specific trial sets
%
%   E.stim_all              : intervalSet of all stimulus epochs
%   E.stim_music/speech/ferret/silence  : category-specific stim sets
%
% All times are in seconds on the OE/ephys timeline.

[~, sess_name] = fileparts(datapath);
disp(['[RP_build_epochs_passive] Session: ' sess_name]);

E = struct;

%% 1) fUS blocks from fUS TTLs (if present)

if isfield(trigOE,'fus')
    if isfield(trigOE.fus,'t_reg_s') && ~isempty(trigOE.fus.t_reg_s)
        t_fus = trigOE.fus.t_reg_s(:);
    else
        t_fus = trigOE.fus.t_raw_s(:);
    end
    
    if numel(t_fus) >= 2
        iti = diff(t_fus);
        medITI = median(iti);
        gap_thr = 1.5 * medITI;   % gap > gap_thr => pause between imaging blocks
        
        gap_idx = find(iti > gap_thr);
        
        blk_start = t_fus(1);
        blk_starts = blk_start;
        blk_ends   = [];
        
        for g = 1:numel(gap_idx)
            blk_ends(end+1,1)   = t_fus(gap_idx(g));      % end at last pulse before gap
            blk_starts(end+1,1) = t_fus(gap_idx(g)+1);    % start at first pulse after gap
        end
        blk_ends(end+1,1) = t_fus(end);
        
        % small padding around blocks
        pad = 0.5 * medITI;
        blk_starts = blk_starts - pad;
        blk_ends   = blk_ends   + pad;
        
        E.fus_blocks = intervalSet(blk_starts, blk_ends);
        fprintf('  fUS blocks: %d blocks from %.2f to %.2f s\n', ...
            length(blk_starts), blk_starts(1), blk_ends(end));
    else
        warning('[RP_build_epochs_passive] Not enough fUS TTLs to define blocks.');
    end
else
    warning('[RP_build_epochs_passive] No trigOE.fus field; skipping fUS blocks.');
end

%% 2) Trial epochs from Baphy + Baphy TTLs

if ~isfield(trigOE,'baphy')
    error('[RP_build_epochs_passive] trigOE.baphy missing (no Baphy TTL channel).');
end

t_trial_onset = trigOE.baphy.t_raw_s(:);   % OE TTL times in seconds
n_oe = numel(t_trial_onset);
n_b = B.n_trials;

fprintf('  Trials: Baphy=%d, OE TTL=%d\n', n_b, n_oe);

if n_oe ~= n_b
    warning('[RP_build_epochs_passive] Trial count mismatch: Baphy=%d, OE=%d', n_b, n_oe);
    n = min(n_oe, n_b);
    t_trial_onset = t_trial_onset(1:n);
else
    n = n_b;
end

trial_start_s = nan(n,1);
trial_end_s   = nan(n,1);
stim_start_s  = nan(n,1);
stim_end_s    = nan(n,1);

for it = 1:n
    t0 = t_trial_onset(it);    % OE time when trial starts
    
    % durations from Baphy (relative to trial)
    t_pre_start  = B.t_pre_start(it);   % usually 0
    t_pre_stop   = B.t_pre_stop(it);
    t_stim_start = B.t_stim_start(it);
    t_stim_stop  = B.t_stim_stop(it);
    t_post_stop  = B.t_post_stop(it);
    
    % define trial start/end: from PreStimSilence start to PostStimSilence stop
    if ~isnan(t_pre_start)
        trial_start_s(it) = t0 + t_pre_start;
    else
        trial_start_s(it) = t0;  % fallback
    end
    
    if ~isnan(t_post_stop)
        trial_end_s(it) = t0 + t_post_stop;
    elseif ~isnan(t_stim_stop)
        trial_end_s(it) = t0 + t_stim_stop;
    else
        % fallback: 5 s trial if all else fails
        trial_end_s(it) = t0 + 5;
    end
    
    % stimulus epoch
    if ~isnan(t_stim_start) && ~isnan(t_stim_stop)
        stim_start_s(it) = t0 + t_stim_start;
        stim_end_s(it)   = t0 + t_stim_stop;
    end
end

% drop any NaNs
bad = isnan(trial_start_s) | isnan(trial_end_s);
if any(bad)
    warning('[RP_build_epochs_passive] Dropping %d trials with NaN start/stop.', sum(bad));
    trial_start_s(bad) = [];
    trial_end_s(bad)   = [];
    stim_start_s(bad)  = [];
    stim_end_s(bad)    = [];
    B.category(bad)    = [];
end

E.trial_all = intervalSet(trial_start_s, trial_end_s);

%% 3) Category-specific trial and stim epochs

cat = B.category(:);
n_trials = numel(cat);

idx_music  = strcmp(cat,'music');
idx_speech = strcmp(cat,'speech');
idx_ferret = strcmp(cat,'ferret');
idx_sil    = strcmp(cat,'silence');

E.trial_music  = intervalSet(trial_start_s(idx_music),  trial_end_s(idx_music));
E.trial_speech = intervalSet(trial_start_s(idx_speech), trial_end_s(idx_speech));
E.trial_ferret = intervalSet(trial_start_s(idx_ferret), trial_end_s(idx_ferret));
E.trial_silence= intervalSet(trial_start_s(idx_sil),    trial_end_s(idx_sil));

% Stim epochs (only where defined)
good_stim = ~isnan(stim_start_s) & ~isnan(stim_end_s);
E.stim_all = intervalSet(stim_start_s(good_stim), stim_end_s(good_stim));

E.stim_music  = intervalSet(stim_start_s(idx_music  & good_stim), stim_end_s(idx_music  & good_stim));
E.stim_speech = intervalSet(stim_start_s(idx_speech & good_stim), stim_end_s(idx_speech & good_stim));
E.stim_ferret = intervalSet(stim_start_s(idx_ferret & good_stim), stim_end_s(idx_ferret & good_stim));
E.stim_silence= intervalSet(stim_start_s(idx_sil    & good_stim), stim_end_s(idx_sil    & good_stim));

fprintf('  Trials kept: %d (music=%d, speech=%d, ferret=%d, silence=%d)\n', ...
    n_trials, sum(idx_music), sum(idx_speech), sum(idx_ferret), sum(idx_sil));

% quick sanity plot if you want
do_plot = 0;
if do_plot
    figure('Color','w');
    subplot(3,1,1);
    plot(trial_start_s, ones(size(trial_start_s)),'k.');
    hold on;
    plot(stim_start_s, 1.1*ones(size(stim_start_s)),'r.');
    xlabel('Time (s)'); ylim([0.9 1.2]);
    title('Trial and stim onsets');
    
    subplot(3,1,2);
    stairs(trial_start_s, strcmp(cat,'music'),'b'); hold on;
    stairs(trial_start_s, strcmp(cat,'speech'),'r');
    stairs(trial_start_s, strcmp(cat,'ferret'),'g');
    legend({'music','speech','ferret'});
    ylabel('category');
    
    if isfield(E,'fus_blocks')
        subplot(3,1,3);
        starts = Start(E.fus_blocks,'s');
        ends   = End(E.fus_blocks,'s');
        for i = 1:numel(starts)
            patch([starts(i) ends(i) ends(i) starts(i)], [0 0 1 1], 'c', 'FaceAlpha',0.3,'EdgeColor','none');
            hold on;
        end
        xlabel('Time (s)');
        ylabel('fUS blocks');
        ylim([0 1.1]);
    end
end

end
