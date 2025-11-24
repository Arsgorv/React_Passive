function React_Passive_SleepScoring_preprocessing
% Don't forget to put ExpeInfo.mat in the target folder
Ferret_ProcessData_BM




%% Make LFP
load(fullfile(datapath, 'ephys', 'ExpeInfo.mat'))
BaseFileName = ['M' num2str(ExpeInfo.nmouse) '_' ExpeInfo.date '_' ExpeInfo.SessionType];
% Create the xml
WriteExpeInfoToXml(ExpeInfo)


% copyfile(fullfile(ExpeInfo.PreProcessingInfo.FolderForConcatenation_Ephys{f}, 'amplifier.xml'),...
%     fullfile(FinalFolder, [BaseFileName '-' sprintf('%02d',f) '.xml']));
% delete([ExpeInfo.PreProcessingInfo.FolderForConcatenation_Ephys{f} filesep 'amplifier.xml']);

SetCurrentSession([BaseFileName '.xml'])
MakeData_LFP_PluggedOnly([fullfile(datapath, 'ephys')],ExpeInfo)

%% Calculate smooth brain powers
for sess = 1:numel(session_dlc)
    disp('Calculating brain power...')
    disp([num2str(length(session_dlc)-sess + 1) '/' num2str(length(session_dlc)) ' left'])
    disp(['Running session: ' session_dlc{sess}])    
    calc_brain_gamma_powers(session_dlc{sess})
end








end