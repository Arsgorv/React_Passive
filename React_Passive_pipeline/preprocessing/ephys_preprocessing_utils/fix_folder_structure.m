function fix_folder_structure(datapath)
% 1) Fix OpenEphys folder structure inside [datapath filesep 'ephys']:
%    makes sure each *_PreExp / *_Exp / *_PostExp folder directly contains "recording1"
%    (and not nested in "Record Node*/experiment*").
%    Also removes now-useless "Record Node*/experiment*" folders.

if ~exist(datapath,'dir')
    error('datapath does not exist: %s', datapath);
end

ephysDir = fullfile(datapath, 'ephys');
if ~exist(ephysDir,'dir')
    error('No "ephys" folder found in datapath: %s', ephysDir);
end

fprintf('--- Step 1: Fixing OpenEphys folder structure ---\n');
fprintf('datapath: %s\n', datapath);
fprintf('ephysDir: %s\n\n', ephysDir);

D = dir(ephysDir);
recIdx = false(length(D),1);
for i = 1:length(D)
    if D(i).isdir && ~strcmp(D(i).name,'.') && ~strcmp(D(i).name,'..')
        recIdx(i) = true;
    end
end
recDirs = D(recIdx);

if isempty(recDirs)
    warning('No subfolders found in ephysDir = %s. Nothing to do.', ephysDir);
    return
end

for iRec = 1:length(recDirs)
    recName = recDirs(iRec).name;
    recPath = fullfile(ephysDir, recName);
    fprintf('  Recording folder: %s\n', recPath);

    rec1Dir = fullfile(recPath, 'recording1');
    if exist(rec1Dir,'dir')
        fprintf('    - recording1 already present at top level. OK.\n');
        continue
    end

    rnList = dir(fullfile(recPath, 'Record Node*'));
    if isempty(rnList)
        warning('    ! No "Record Node*" folder and no recording1 found in %s', recPath);
        continue
    end

    rnPath = fullfile(recPath, rnList(1).name);
    fprintf('    - Found %s\n', rnPath);

    expList = dir(fullfile(rnPath, 'experiment*'));
    if isempty(expList)
        warning('    ! No experiment* folder found inside %s. Skipping.', rnPath);
        continue
    end
    expPath = fullfile(rnPath, expList(1).name);
    fprintf('    - Found %s\n', expPath);

    recList2 = dir(fullfile(expPath, 'recording*'));
    if isempty(recList2)
        warning('    ! No recording* folder found inside %s. Skipping.', expPath);
        continue
    end

    srcRec = fullfile(expPath, recList2(1).name);
    dstRec = fullfile(recPath, recList2(1).name);

    if exist(dstRec,'dir')
        warning('    ! Destination %s already exists, not moving.', dstRec);
    else
        fprintf('    - Moving %s -> %s\n', srcRec, dstRec);
        movefile(srcRec, dstRec);
    end

    % Now clean up leftover experiment*/settings.xml and Record Node* if empty
    if exist(expPath,'dir')
        settingsFile = fullfile(expPath, 'settings.xml');
        if exist(settingsFile,'file')
            fprintf('    - Deleting %s\n', settingsFile);
            delete(settingsFile);
        end
        try
            rmdir(expPath);
            fprintf('    - Removed empty folder %s\n', expPath);
        catch
            fprintf('    - Could not remove %s (not empty?)\n', expPath);
        end
    end

    if exist(rnPath,'dir')
        try
            rmdir(rnPath);
            fprintf('    - Removed empty folder %s\n', rnPath);
        catch
            fprintf('    - Could not remove %s (not empty?)\n', rnPath);
        end
    end
end

fprintf('Step 1 completed.\n\n');
end