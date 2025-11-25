function copy_ExpeInfo(datapath)
% Copy ExpeInfo.mat into [datapath filesep 'ephys'] if it is not already there.

ephysDir = fullfile(datapath, 'ephys');
if ~exist(ephysDir,'dir')
    error('No "ephys" folder found in datapath: %s', ephysDir);
end

dest = fullfile(ephysDir, 'ExpeInfo.mat');
if exist(dest,'file')
    fprintf('  ExpeInfo.mat already exists in ephys. Skipping copy.\n\n');
    return
end

ExpeInfoSource = '';
if isempty(ExpeInfoSource)
    parentDir = fileparts(datapath);
    candidate = fullfile(parentDir, 'ExpeInfo.mat');
    if exist(candidate,'file')
        ExpeInfoSource = candidate;
        fprintf('  Using ExpeInfo source from parent folder: %s\n', ExpeInfoSource);
    else
        warning('  ! ExpeInfo.mat not found in parent folder: %s', candidate);
    end
end

if isempty(ExpeInfoSource)
    fprintf('  No ExpeInfo.mat copied (source not found).\n\n');
    return
end

fprintf('  Copying to: %s\n', dest);
copyfile(ExpeInfoSource, dest);
fprintf('Step 2 completed.\n\n');

end