function React_Passive_SleepScoring_preprocessing
% Don't forget to put ExpeInfo.mat in the target folder
Ferret_ProcessData_BM
% quick version of sleep scoring, where I only calculate smooth brain powers
for sess = 1:numel(session_dlc)
    disp('Calculating brain power...')
    disp([num2str(length(session_dlc)-sess + 1) '/' num2str(length(session_dlc)) ' left'])
    disp(['Running session: ' session_dlc{sess}])    
    calc_brain_gamma_powers(session_dlc{sess})
end








end