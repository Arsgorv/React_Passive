function Master_data_sync_preproc(sessions)
% Master_data_sync_preproc
% React_Passive_AG master wrapper: sync all data streams and build epochs.
%
% INPUT
%   sessions : cell array of session folders, e.g.
%       sessions = {
%         '/home/Arsenii/React_Passive/Processed_data/Edel/20220419_1_m_C'
%         '/home/Arsenii/React_Passive/Processed_data/Kosichka/20220309_1_m_A'
%         ...
%       };

if ~iscell(sessions)
    error('Input "sessions" must be a cell array of session paths');
end

for sess = 1:numel(sessions)
    datapath = sessions{sess};
    disp('------------------------------------------')
    disp(['Working on ' datapath])

    if ~isfolder(datapath)
        warning(['Folder not found, skipping: ' datapath])
        continue
    end

    try
        % 1) Passive-specific csv sync for Edel Exp phase
        trig_csv = [];
        if contains(datapath,'Edel')
            disp('  [1] Syncing Baphy–fUS from csv (RP_sync_triggers_passive)...')
            trig_csv = RP_sync_triggers_passive(datapath);
        else
            disp('  [1] No csv sync step for this animal (skip RP_sync_triggers_passive).')
        end

        % 2) Get channel configuration for this session
        disp('  [2] Loading trigger channel config (RP_get_trigger_cfg)...')
        cfg = get_trigger_config(datapath);

        % 3) Extract TTLs from OpenEphys (fUS / Baphy / video / respi / heart ...)
        disp('  [3] Extracting OE TTL triggers (extract_triggers_oe)...')
        trigOE = extract_triggers_oe(datapath, cfg);

        % 4) Parse Baphy m-file for React Passive (categories, timings, etc.)
        disp('  [4] Parsing Baphy file (RP_parse_baphy_passive)...')
        B = RP_parse_baphy_passive(datapath);

        % 5) Build epochs on the OE/common time axis
        disp('  [5] Building epochs (RP_build_epochs_passive)...')
        E = RP_build_epochs_passive(datapath, trigOE, B);

        % 6) Save everything in one place
        out_file = fullfile(datapath, 'RP_master_sync.mat');
        save(out_file, 'trig_csv', 'trigOE', 'B', 'E');
        disp(['  Saved master sync file: ' out_file])

    catch ME
        warning('Error in session %s:\n  %s', datapath, ME.message);
        disp(getReport(ME,'basic'))
        continue
    end
end

disp('Master_data_sync_preproc finished.')

end
