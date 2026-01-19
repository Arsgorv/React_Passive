function E = RP_build_epochs_passive(datapath, trigOE, B)
% RP_build_epochs_passive
% Build intervalSet epochs for React Passive session
% using:
%   - OE TTLs: trigOE.baphy (trial onsets), trigOE.fus (fUS TTLs)
%   - Baphy m-file parsed by RP_parse_baphy_passive (B)
%
% OUTPUT E struct:
%   E.fus_preexp / fus_exp / fus_postexp : fUS imaging blocks (ts units)
%   E.fus_blocks                         : union of all fUS blocks
%
%   E.trial_all             : intervalSet of all trials
%   E.trial_music/speech/ferret/silence : category-specific trial sets
%
%   E.stim_all              : intervalSet of all stimulus epochs
%   E.stim_music/speech/ferret/silence  : category-specific stim sets
%
% All times are in seconds on the OE/ephys timeline, converted to ts (1e4).
        
[~, sess_name] = fileparts(datapath);
disp(['[RP_build_epochs_passive] Session: ' sess_name]);

E = struct;

E.fus_preexp  = intervalSet([],[]);
E.fus_exp     = intervalSet([],[]);
E.fus_postexp = intervalSet([],[]);
E.fus_blocks  = intervalSet([],[]);

%% 1) fUS blocks from OE TTLs (if present)

have_fus = isfield(trigOE,'fus') && ...
    ( (isfield(trigOE.fus,'t_s')     && ~isempty(trigOE.fus.t_s)) || ...
      (isfield(trigOE.fus,'t_raw_s') && ~isempty(trigOE.fus.t_raw_s)) );

if have_fus
    if isfield(trigOE.fus,'t_s') && ~isempty(trigOE.fus.t_s)
        t_fus = trigOE.fus.t_s(:);
    else
        t_fus = trigOE.fus.t_raw_s(:);
    end
    t_fus = sort(t_fus(:));

    if numel(t_fus) >= 2
        dt = diff(t_fus);
        medITI = median(dt);
        
        phase_gap_thr =  1.5 * medITI; % seconds
        gap_idx = find(dt > phase_gap_thr);
        seg_start = [1; gap_idx + 1];
        seg_end   = [gap_idx; numel(t_fus)];
        n_seg = numel(seg_start);

        % try to assign 3 segments: PreExp, Exp, PostExp
        if n_seg >= 3
            idx_pre  = seg_start(1):seg_end(1);
            idx_exp  = seg_start(2):seg_end(2);
            idx_post = seg_start(3):seg_end(3);

            t_pre_start  = t_fus(idx_pre(1));
            t_pre_end    = t_fus(idx_pre(end));
            t_exp_start  = t_fus(idx_exp(1));
            t_exp_end    = t_fus(idx_exp(end));
            t_post_start = t_fus(idx_post(1));
            t_post_end   = t_fus(idx_post(end));

            E.fus_preexp  = intervalSet(t_pre_start*1e4,  t_pre_end*1e4);
            E.fus_exp     = intervalSet(t_exp_start*1e4,  t_exp_end*1e4);
            E.fus_postexp = intervalSet(t_post_start*1e4, t_post_end*1e4);
            E.fus_blocks  = union(E.fus_preexp, union(E.fus_exp, E.fus_postexp));

            fprintf('  fUS blocks (OE): PreExp[%.2f–%.2f], Exp[%.2f–%.2f], PostExp[%.2f–%.2f] s\n', ...
                t_pre_start, t_pre_end, t_exp_start, t_exp_end, t_post_start, t_post_end);

            % optional: check counts against exp_info
            fus_dir  = fullfile(datapath,'fUS');
            exp_file = fullfile(fus_dir,'exp_info.mat');
            if exist(exp_file,'file')
                try
                    S = load(exp_file);
                    exp_info = S.exp_info;
                    if isfield(exp_info,'PreExp')
                        n_pre_exp  = exp_info.PreExp.size(3);
                        n_exp      = exp_info.Exp.size(3);
                        n_post_exp = exp_info.PostExp.size(3);
                    elseif isfield(exp_info,'size')
                        n_pre_exp  = exp_info.size{1}(3);
                        n_exp      = exp_info.size{2}(3);
                        n_post_exp = exp_info.size{3}(3);
                    else
                        n_pre_exp = []; n_exp = []; n_post_exp = [];
                    end

                    n_pre_ttl  = numel(idx_pre);
                    n_exp_ttl  = numel(idx_exp);
                    n_post_ttl = numel(idx_post);

                    if ~isempty(n_pre_exp) && n_pre_ttl ~= n_pre_exp
                        warning('[OE-fUS] PreExp frames mismatch in %s: TTL=%d, fUS=%d', ...
                            datapath, n_pre_ttl, n_pre_exp);
                    end
                    if ~isempty(n_exp) && n_exp_ttl ~= n_exp
                        warning('[OE-fUS] Exp frames mismatch in %s: TTL=%d, fUS=%d', ...
                            datapath, n_exp_ttl, n_exp);
                    end
                    if ~isempty(n_post_exp) && n_post_ttl ~= n_post_exp
                        warning('[OE-fUS] PostExp frames mismatch in %s: TTL=%d, fUS=%d', ...
                            datapath, n_post_ttl, n_post_exp);
                    end
                catch ME
                    warning('[OE-fUS] exp_info frame check failed in %s: %s', datapath, ME.message);
                end
            end

        else
            warning('[RP_build_epochs_passive] Could not detect 3 fUS blocks from OE TTLs (n_seg=%d).', n_seg);
            % We at least define a single fUS block spanning all TTLs
            t_start = t_fus(1);
            t_end   = t_fus(end);
            E.fus_exp    = intervalSet(t_start*1e4, t_end*1e4);
            E.fus_blocks = E.fus_exp;
        end

    else
        warning('[RP_build_epochs_passive] Not enough fUS TTLs to define blocks.');
    end
else
    warning('[RP_build_epochs_passive] No fUS TTLs in trigOE; fUS epochs must come from csv/data for this session.');
end

%% 2) Trial epochs from Baphy + Baphy TTLs

if ~isfield(trigOE,'baphy')
    error('[RP_build_epochs_passive] trigOE.baphy missing (no Baphy TTL channel).');
end

if isfield(trigOE.baphy,'t_s') && ~isempty(trigOE.baphy.t_s)
    t_trial_onset = trigOE.baphy.t_s(:);
else
    t_trial_onset = trigOE.baphy.t_raw_s(:); 
end

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

stim_start_s  = nan(n,1);
stim_end_s    = nan(n,1);

trial_start_s = trigOE.baphy.trial_start_ts/1e4;
trial_end_s = trigOE.baphy.trial_stop_ts/1e4;

for it = 1:n
    t0 = t_trial_onset(it);    % OE time when trial starts
    % stimulus epoch
    stim_start_s(it) = t0 + B.t_stim_start(it);
    stim_end_s(it)   = t0 + B.t_stim_stop(it);
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

E.trial_all = intervalSet(trial_start_s*1e4, trial_end_s*1e4);

%% 3) Category-specific trial and stim epochs

cat = B.category(:);
n_trials = numel(cat);

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

fprintf('  Trials kept: %d (music=%d, speech=%d, ferret=%d)\n', ...
    n_trials, sum(idx_music), sum(idx_speech), sum(idx_ferret));

% quick sanity plot if you want
do_plot = 1;
if do_plot
    figure('Color','w');
    subplot(3,1,1);
    plot(trial_start_s, ones(size(trial_start_s)),'k.');
    hold on;
    plot(stim_start_s, 1.1*ones(size(stim_start_s)),'r.');
    xlabel('Time (s)'); ylim([0 2]); legend({'trial onset'; 'stim onset'})
    title('Trial and stim onsets');
    
    subplot(3,1,2);
    if ~isempty(trial_start_s(idx_music))
        plot(trial_start_s(strcmp(cat,'music')), 1 ,'b*'); hold on;
    end
    if ~isempty(trial_start_s(idx_speech))
        
        plot(trial_start_s(strcmp(cat,'speech')),1, 'r*');
    end
    if ~isempty(trial_start_s(idx_ferret))
        plot(trial_start_s(strcmp(cat,'ferret')),1, 'g*');
    end
%     legend({'music','speech','ferret'});
    ylabel('category');
    
    if isfield(E,'fus_blocks')
        starts = Start(E.fus_blocks, 's');
        ends   = End(E.fus_blocks, 's');
        for i = 1:numel(starts)
            patch([starts(i) ends(i) ends(i) starts(i)], [0 0 1 1], 'k', 'FaceAlpha',0.1,'EdgeColor','none');
            hold on;
        end
        xlabel('Time (s)');
        ylabel('fUS blocks');
        ylim([0 1.1]);
    end
end

end
