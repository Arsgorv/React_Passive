function E = RP_build_epochs_tonotopy(datapath, trigOE, B)
% RP_build_epochs_passive
% Build intervalSet epochs for Tonotopy session
% using:
%   - OE TTLs: trigOE.baphy (trial onsets), trigOE.fus (fUS TTLs)
%   - Baphy m-file parsed by RP_parse_baphy_tonotopy (B)
%
% OUTPUT E struct:
%   E.fus_blocks            : intervalSet of continuous fUS imaging blocks
%
%   E.trial_all             : intervalSet of all trials
%   E.trial_200/400/800/1600/32000/6400/12800/etc : frequence-specific trial sets
%
%   E.stim_all              : intervalSet of all stimulus epochs
%   E.stim_200/400/800/1600/32000/6400/12800/etc/silence  : frequence-specific stim sets
%
% All times are in seconds on the OE/ephys timeline.

[~, sess_name] = fileparts(datapath);
disp(['[RP_build_epochs_passive] Session: ' sess_name]);

E = struct;

%% 1) fUS blocks from fUS TTLs (if present)

if isfield(trigOE,'fus')
    if isfield(trigOE.fus,'t_s') && ~isempty(trigOE.fus.t_s)
        t_fus = trigOE.fus.t_s(:);
    else
        t_fus = trigOE.fus.t_raw_s(:);
    end
    
    Slices={'A', 'B', 'C', 'D'};
    for iSlice = 1:length(Slices)
        fus_file = dir(strcat(datapath, '/fUS/RP_data_*slice_', Slices{iSlice}, '.mat'));
        load(fullfile(datapath, 'fUS', fus_file.name))
        cat_tsd_new = tsd(t_fus*10000, Data(cat_tsd.data));  
        cat_tsd.data = cat_tsd_new; 
        save(fullfile(datapath, 'fUS', fus_file.name), 'cat_tsd', '-append')
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
        
        E.fus_blocks = intervalSet(blk_starts*1e4, blk_ends*1e4);
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

if isfield(trigOE.baphy,'t_s') && ~isempty(trigOE.baphy.t_s)
    t_trial_onset = trigOE.baphy.t_s(:);
else
    t_trial_onset = trigOE.baphy.t_raw_s(:); 
end

n_oe = numel(t_trial_onset);
n_b = B.n_trials;

fprintf('  Trials: Baphy=%d, OE TTL=%d\n', n_b, n_oe);

if n_oe ~= n_b
    warning('[RP_build_epochs_tonotopy] Trial count mismatch: Baphy=%d, OE=%d', n_b, n_oe);
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
    warning('[RP_build_epochs_tonotopy] Dropping %d trials with NaN start/stop.', sum(bad));
    trial_start_s(bad) = [];
    trial_end_s(bad)   = [];
    stim_start_s(bad)  = [];
    stim_end_s(bad)    = [];
    B.frequence(bad)    = [];
end

E.trial_all = intervalSet(trial_start_s*1e4, trial_end_s*1e4);

%% 3) frequence-specific trial and stim epochs

cat = B.frequence(:);
n_trials = numel(cat);

cats = unique(cat);          % all unique categories (cell array of char)
nCats = numel(cats);


for iC = 1:nCats
    cname = cats{iC};
    idx   = strcmp(cat, cname);

    if any(idx)
        fieldname = ['trial_' cname];
        E.(fieldname) = intervalSet( ...
            trial_start_s(idx)*1e4, ...
            trial_end_s(idx)*1e4 );
    end
end

good_stim = ~isnan(stim_start_s) & ~isnan(stim_end_s);

% all stim epochs
E.stim_all = intervalSet( ...
    stim_start_s(good_stim)*1e4, ...
    stim_end_s(good_stim)*1e4 );

% per frequence category
for iC = 1:nCats
    cname = cats{iC};
    idx   = strcmp(cat, cname) & good_stim;

    if any(idx)
        fieldname = ['stim_' cname];
        E.(fieldname) = intervalSet( ...
            stim_start_s(idx)*1e4, ...
            stim_end_s(idx)*1e4 );
    end
end

fprintf('  Trials kept: %d\n', n_trials);
for iC = 1:nCats
    fprintf('    %s = %d\n', cats{iC}, sum(strcmp(cat, cats{iC})));
end


do_plot = 1;
if do_plot
    figure('Color','w');

    subplot(3,1,1);
    plot(trial_start_s, ones(size(trial_start_s)),'k.');
    hold on;
    plot(stim_start_s, 1.1*ones(size(stim_start_s)),'r.');
    xlabel('Time (s)');
    ylim([0 2]);
    legend({'trial onset','stim onset'});
    title('Trial and stim onsets');

    subplot(3,1,2); hold on;
    colors = lines(nCats);

    for iC = 1:nCats
        idx = strcmp(cat, cats{iC});
        plot(trial_start_s(idx), iC*ones(sum(idx),1), '*', ...
            'Color', colors(iC,:), 'DisplayName', cats{iC});
    end
    yticks(1:nCats);
    yticklabels(cats);
    ylabel('frequence');
    legend;

    if isfield(E,'fus_blocks')
        subplot(3,1,3); hold on;
        starts = Start(E.fus_blocks, 's');
        ends   = End(E.fus_blocks, 's');
        for i = 1:numel(starts)
            patch([starts(i) ends(i) ends(i) starts(i)], ...
                  [0 0 1 1], 'k', ...
                  'FaceAlpha',0.1,'EdgeColor','none');
        end
        xlabel('Time (s)');
        ylabel('fUS blocks');
        ylim([0 1.1]);
    end
end


end
