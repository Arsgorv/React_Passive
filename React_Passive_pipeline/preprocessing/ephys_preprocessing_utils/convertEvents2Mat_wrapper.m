function convertEvents2Mat_wrapper(datapath, arg1, arg2)
% Run convertEvents2Mat.py on each recording1/continuous and recording1/events
% inside [datapath filesep 'ephys'].
%
% Windows call:
%   convertEvents2Mat_wrapper(datapath, pyExe, github_root)
%
% Linux call:
%   convertEvents2Mat_wrapper(datapath, github_root)

fprintf('--- Step 3: Running convertEvents2Mat.py ---\n');

if ~exist(datapath,'dir')
    error('datapath does not exist: %s', datapath);
end

ephysDir = fullfile(datapath, 'ephys');
if ~exist(ephysDir,'dir')
    error('No "ephys" folder found in datapath: %s', ephysDir);
end

% --------- parse inputs depending on OS ----------
if ispc
    if nargin < 3
        error('On Windows, call convertEvents2Mat_wrapper(datapath, pyExe, github_root)');
    end
    pyExe = arg1;
    github_root = arg2;
else
    if nargin < 2
        error('On Linux, call convertEvents2Mat_wrapper(datapath, github_root)');
    end
    github_root = arg1;
    % on Linux we assume "python" is on PATH
    pyExe = 'python';
end

if ~exist(github_root,'dir')
    warning('github_root folder does not exist: %s', github_root);
end

% --------- locate convertEvents2Mat.py ----------
defaultScript = fullfile(github_root, 'React_Passive', 'React_Passive_pipeline', ...
    'preprocessing', 'ephys_preprocessing_utils', 'convertEvents2Mat.py');

if exist(defaultScript,'file')
    pyScript = defaultScript;
    fprintf('  Found convertEvents2Mat.py at:\n    %s\n', pyScript);
else
    fprintf('  Could not find convertEvents2Mat.py at default location:\n    %s\n', defaultScript);
    fprintf('  Please select convertEvents2Mat.py manually.\n');
    [f,p] = uigetfile('*.py','Select convertEvents2Mat.py');
    if isequal(f,0)
        warning('  No script selected. Aborting Step 3.');
        return
    end
    pyScript = fullfile(p,f);
    fprintf('  Using script: %s\n', pyScript);
end

if ~exist(pyScript,'file')
    error('convertEvents2Mat.py not found at %s', pyScript);
end

% --------- check python executable (Windows sanity) ----------
if ispc
    if ~exist(pyExe,'file')
        error('Selected Python executable does not exist: %s', pyExe);
    end
    fprintf('  Using Python executable:\n    %s\n', pyExe);
else
    % Linux: trust PATH, but still print
    fprintf('  Using Python interpreter: %s (from PATH)\n', pyExe);
end

% --------- list recording folders ----------
D = dir(ephysDir);
recIdx = false(length(D),1);
for i = 1:length(D)
    if D(i).isdir && ~strcmp(D(i).name,'.') && ~strcmp(D(i).name,'..')
        recIdx(i) = true;
    end
end
recDirs = D(recIdx);

if isempty(recDirs)
    warning('No recording subfolders in ephysDir = %s. Nothing to do.', ephysDir);
    return
end

% --------- loop over recordings ----------
for iRec = 1:length(recDirs)
    recName = recDirs(iRec).name;
    recPath = fullfile(ephysDir, recName);
    rec1Dir = fullfile(recPath, 'recording1');

    fprintf('  Recording: %s\n', recPath);

    if ~exist(rec1Dir,'dir')
        fprintf('    - No recording1 folder, skipping.\n');
        continue
    end

    contDir = fullfile(rec1Dir, 'continuous');
    evtDir  = fullfile(rec1Dir, 'events');

    % continuous
    if exist(contDir,'dir')
        if ispc
            cmd1 = sprintf('"%s" "%s" -p "%s"', pyExe, pyScript, contDir);
        else
            cmd1 = sprintf('%s "%s" -p "%s"', pyExe, pyScript, contDir);
        end
        fprintf('    - Running (continuous): %s\n', cmd1);
        system(cmd1);
    else
        warning('    ! continuous folder not found: %s', contDir);
    end

    % events
    if exist(evtDir,'dir')
        if ispc
            cmd2 = sprintf('"%s" "%s" -p "%s"', pyExe, pyScript, evtDir);
        else
            cmd2 = sprintf('%s "%s" -p "%s"', pyExe, pyScript, evtDir);
        end
        fprintf('    - Running (events):     %s\n', cmd2);
        system(cmd2);
    else
        warning('    ! events folder not found: %s', evtDir);
    end
end

fprintf('Step 3 completed.\n\n');
end
