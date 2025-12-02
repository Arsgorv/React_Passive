function B = RP_parse_baphy_passive(datapath)
% RP_parse_baphy_passive
% Parse Baphy m-file for React Passive (Edel / Chabichou / Kosichka)
% and extract trial-wise stimulus identity and timing.
%
% OUTPUT struct B:
%   B.n_trials          : number of trials
%   B.trial_id          : [n_trials x 1] trial numbers (1..N)
%   B.stim_name         : cellstr, filename part from Stim note
%   B.cat_id            : numeric category id (e.g. 668, 650...)
%   B.category          : 'music' | 'speech' | 'ferret' | 'silence' | 'unknown'
%   B.is_silence        : 0/1
%   B.t_pre_start       : PreStimSilence start (s, rel to trial)
%   B.t_pre_stop        : PreStimSilence stop
%   B.t_stim_start      : Stim start
%   B.t_stim_stop       : Stim stop
%   B.t_post_start      : PostStimSilence start
%   B.t_post_stop       : PostStimSilence stop
%
% It assumes that running the .m file puts globalparams, exptparams,
% exptevents into the workspace (standard Baphy behaviour).

baphy_dir = fullfile(datapath,'baphy');
m_list = dir(fullfile(baphy_dir,'*.m'));
if isempty(m_list)
    error('RP_parse_baphy_passive:NoMfile','No *.m Baphy file found in %s',baphy_dir);
end
if numel(m_list) > 1
    warning('RP_parse_baphy_passive:MultipleMfile','Multiple m-files, using first: %s', m_list(1).name);
end

m_file = fullfile(m_list(1).folder, m_list(1).name);
disp(['[RP_parse_baphy_passive] Running m-file: ' m_file]);

run(m_file);   % defines globalparams, exptparams, exptevents

if ~exist('exptevents','var')
    error('RP_parse_baphy_passive:NoExptevents','exptevents not found after running %s',m_file);
end

trial_all = [exptevents.Trial];
trial_list = unique(trial_all);
trial_list(trial_list==0) = [];   % drop events not assigned to a trial, if any
n_trials = numel(trial_list);

B = struct;
B.n_trials     = n_trials;
B.trial_id     = trial_list(:);
B.stim_name    = cell(n_trials,1);
B.cat_id       = nan(n_trials,1);
B.category     = cell(n_trials,1);
B.is_silence   = zeros(n_trials,1);
B.t_pre_start  = nan(n_trials,1);
B.t_pre_stop   = nan(n_trials,1);
B.t_stim_start = nan(n_trials,1);
B.t_stim_stop  = nan(n_trials,1);
B.t_post_start = nan(n_trials,1);
B.t_post_stop  = nan(n_trials,1);

for it = 1:n_trials
    tr = trial_list(it);
    idx = find(trial_all == tr);
    
    notes = {exptevents(idx).Note};
    starts = [exptevents(idx).StartTime];
    stops  = [exptevents(idx).StopTime];
    
    % PreStimSilence
    p_idx = find(strncmp(notes,'PreStimSilence',14));
    if ~isempty(p_idx)
        B.t_pre_start(it) = starts(p_idx(1));
        B.t_pre_stop(it)  = stops(p_idx(1));
    end
    
    % Stim
    s_idx = find(strncmp(notes,'Stim ',5));
    if ~isempty(s_idx)
        note_stim = notes{s_idx(1)};
        B.t_stim_start(it) = starts(s_idx(1));
        B.t_stim_stop(it)  = stops(s_idx(1));
        
        % parse sound name: "Stim , NAME , ..."
        parts = strsplit(note_stim,' , ');
        if numel(parts) >= 2
            B.stim_name{it} = strtrim(parts{2});
        else
            B.stim_name{it} = note_stim;
        end
    else
        B.stim_name{it} = '';
    end
    
    % PostStimSilence
    q_idx = find(strncmp(notes,'PostStimSilence',15));
    if ~isempty(q_idx)
        B.t_post_start(it) = starts(q_idx(1));
        B.t_post_stop(it)  = stops(q_idx(1));
    end
    
    % category
    nm = B.stim_name{it};
    if isempty(nm)
        B.category{it} = 'unknown';
        continue;
    end
    
    if contains(nm,'Silence')
        B.is_silence(it) = 1;
        B.category{it}   = 'silence';
        continue;
    end
    
    % look for "catXXX" pattern
    tok = regexp(nm,'cat(\d+)','tokens');
    if ~isempty(tok)
        cid = str2double(tok{1}{1});
        B.cat_id(it) = cid;
        
        if cid == 668
            B.category{it} = 'ferret';
        elseif cid < 668
            B.category{it} = 'music';
        elseif cid > 668
            B.category{it} = 'speech';
        end
    else
        % fallback: use the snippet you mentioned if "catNNN" is not there
        B.cat_id(it)   = NaN;
        B.category{it} = 'unknown';
    end
end

% make category a column cellstr
B.category = B.category(:);

% simple summary
n_music  = sum(strcmp(B.category,'music'));
n_speech = sum(strcmp(B.category,'speech'));
n_ferret = sum(strcmp(B.category,'ferret'));
n_sil    = sum(B.is_silence);

disp(sprintf('[RP_parse_baphy_passive] Trials: %d (music=%d, speech=%d, ferret=%d, silence=%d)', ...
    n_trials, n_music, n_speech, n_ferret, n_sil));

end
