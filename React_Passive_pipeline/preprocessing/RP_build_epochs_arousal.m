function E = RP_build_epochs_arousal(datapath, trigOE)
% RP_build_epochs_arousal
% Build intervalSet epochs for React Passive session
% using:
%   - OE TTLs: trigOE.baphy (trial onsets), trigOE.fus (fUS TTLs)
%   - Baphy m-file parsed by RP_parse_baphy_passive (B)
%
% OUTPUT E struct:
%   E.fus_Arousal1, Arousal2, Arousal3, etc : fUS imaging blocks (ts units)
%   E.fus_blocks                         : union of all fUS blocks
%
% All times are in seconds on the OE/ephys timeline, converted to ts (1e4).
        
[~, sess_name] = fileparts(datapath);
disp(['[RP_build_epochs_passive] Session: ' sess_name]);

E = struct;
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
        phase_gap_thr = 1.5 * medITI;
        gap_idx = find(dt > phase_gap_thr);        
        seg_start = [1; gap_idx + 1];
        seg_end   = [gap_idx; numel(t_fus)];
        n_seg = numel(seg_start);
        
        % Getting N Arousal sessions
        fus_dir  = fullfile(datapath,'fUS');
        exp_file = fullfile(fus_dir,'exp_info.mat');        
        n_ArousalSessions = [];        
        if exist(exp_file,'file')
            S = load(exp_file);
            exp_info = S.exp_info;
            
            if isfield(exp_info,'size')
                nSessions = numel(exp_info.size);
                for iSess = 1:nSessions
                    n_ArousalSessions(iSess) = exp_info.size{iSess}(3);
                end
            else
                nSessions = 0;
            end
        else
            nSessions = 0;
        end
        
        % Check segmentation
        if n_seg >= 1            
            E.fus_blocks = intervalSet([],[]);  
            remove_idx = []; 
            for iSess = 1:n_seg                
                idx = seg_start(iSess):seg_end(iSess);                
                t_start = t_fus(idx(1));
                t_end   = t_fus(idx(end));
                fieldName = sprintf('Arousal_session%d', iSess);                
                E.(fieldName) = intervalSet(t_start*1e4, t_end*1e4);                
                % accumulate all blocks                
                if iSess == 1
                    E.fus_blocks = E.(fieldName);
                else
                    E.fus_blocks = union(E.fus_blocks, E.(fieldName));
                end                
                fprintf('  fUS block %d: [%.2f – %.2f] s (%d frames)\n', ...
                    iSess, t_start, t_end, numel(idx));
                
                % Compare with expected Nb of frames
                
                if iSess <= numel(n_ArousalSessions)                    
                    n_expected = n_ArousalSessions(iSess);
                    n_ttl = numel(idx);                    
                    if n_ttl ~= n_expected                        
                        warning('[OE-fUS] Arousal Session %d mismatch in %s: TTL=%d, fUS=%d', ...
                            iSess, datapath, n_ttl, n_expected);                        
                        % Remove the last ones
                        if n_expected < n_ttl
                            remove_idx = [remove_idx, idx(n_expected+1:end)];
                        end
                    end
                end
            end
            t_fus(remove_idx) = [];
            % Global consistency of nb of sessions detected warning
            if n_seg ~= nSessions
                warning('[OE-fUS] Segment count (%d) ≠ expected sessions (%d) in %s.', ...
                    n_seg, nSessions, datapath);
            end
            
        else
            warning('[RP_build_epochs_passive] Could not detect fUS blocks.');
            
            t_start = t_fus(1);
            t_end   = t_fus(end);
            
            E.fus_session1 = intervalSet(t_start*1e4, t_end*1e4);
            E.fus_blocks   = E.fus_session1;
        end
        
    else
        warning('[RP_build_epochs_passive] Not enough fUS TTLs to define blocks.');
    end
    
    Slices={'A', 'B', 'C', 'D'};
    for iSlice = 1:length(Slices)
        fus_file = dir(strcat(datapath, '/fUS/RP_data_*slice_', Slices{iSlice}, '.mat'));
        load(fullfile(datapath, 'fUS', fus_file.name), 'cat_tsd')
        if length(t_fus) == length(Range(cat_tsd.data))
            cat_tsd_new = tsd(t_fus*10000, Data(cat_tsd.data));  
            cat_tsd.data = cat_tsd_new; 
            save(fullfile(datapath, 'fUS', fus_file.name), 'cat_tsd', '-append')
        end
    end
    
else
    warning('[RP_build_epochs_passive] No fUS TTLs in trigOE; fUS epochs must come from csv/data for this session.');
end

end
