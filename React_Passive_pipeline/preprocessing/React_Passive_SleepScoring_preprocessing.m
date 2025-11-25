function React_Passive_SleepScoring_preprocessing(sessions)
% Full original preprocessing of the OB project can be found here: Ferret_ProcessData_BM.m

%% Prepare data
github_location = {'D:\Arsenii\GitHub\'; ''};
python_location = 'C:\Users\Arsenii Goriachenkov\.conda\envs\sleepscoring\python.exe';

for sess = 1:numel(sessions)
    disp(['Working on ' sessions{sess}])
    
    fix_folder_structure(sessions{sess})
    copy_ExpeInfo(sessions{sess})
    
    if ispc 
        % Windows LB1
        convertEvents2Mat_wrapper(datapath, python_location, github_location{1});
    else
        % Linux
        convertEvents2Mat_wrapper(datapath, github_location{2});
    end
end

%% PreProcess data
for sess = 1:numel(sessions)
    disp(['Working on ' sessions{sess}]) 
    cd(fullfile(sessions{sess}, 'ephys'))
    
    GUI_StepOne_ExperimentInfo
end

%% Calculate necessary spectrograms
for sess = 1:numel(sessions)
    disp(['Working on ' sessions{sess}])
    
    calculate_spectrograms(sessions{sess},'')
end

%% Do the SleepScoring
for sess = 1:numel(sessions)
    disp(['Working on ' sessions{sess}])
    cd(fullfile(sessions{sess}, 'ephys'))
    
    SleepScoring_Ferret_FV_BAMG('recompute', 1, 'full_ob', 1)
end

%% Calculate brain powers
sm_w = 0.1;
for sess = 1:numel(sessions)
    disp(['Working on ' sessions{sess}])
    
    calculate_brain_power(fullfile(sessions{sess}, 'ephys'), sm_w)
end

end