function calculate_brain_power(datapath, sm_w)
% calculate_brain_power(datapath, sm_w)
% Computes gamma power (and optionally theta ratio + delta)
% for OB / HPC / PFC / ACx depending on the animal.
%
% INPUTS
%   datapath : folder containing LFPData/LFP*.mat
%   sm_w     : smoothing window in seconds (default 0.1)
%
% OUTPUT
%   Saves BrainPower into [datapath filesep 'SleepScoring_OBGamma.mat']
%
% BrainPower fields:
%   BrainPower.signal_names : cellstr, same order as BrainPower.Power
%   BrainPower.Power        : cell array of tsd objects

if nargin < 2 || isempty(sm_w)
    sm_w = 0.1;
end

if ~exist(datapath,'dir')
    error('calc_brain_gamma_powers: datapath does not exist: %s', datapath);
end

fprintf('\n--- calc_brain_gamma_powers ---\n');
fprintf('datapath: %s\n', datapath);

%% ------------------------------------------------------------------------
%  1) Animal-specific configuration
%     Define which channels and signals you want per animal
% -------------------------------------------------------------------------

cfg.gamma_ch   = [];   % vector of channels for gamma
cfg.gamma_names = {};  % corresponding names (same length as gamma_ch)
cfg.theta_ch   = [];   % single channel for theta/delta ratio
cfg.theta_name = '';   % name, e.g. 'HPC_theta'
cfg.delta_ch   = [];   % single channel for delta power
cfg.delta_name = '';   % name, e.g. 'OB_delta'

if contains(datapath, 'Labneh')
    % channels comment from your original code: [24 1] -> [OB HPC]
    cfg.gamma_ch    = [24 1];
    cfg.gamma_names = {'OB_gamma','HPC_gamma'};
    cfg.theta_ch    = 1;      % HPC channel
    cfg.theta_name  = 'HPC_theta';
    cfg.delta_ch    = 24;     % OB channel
    cfg.delta_name  = 'OB_delta';

elseif contains(datapath, 'Brynza')
    % [11 21 25 2] -> [OB HPC PFC ACx]
    cfg.gamma_ch    = [11 21 25 2];
    cfg.gamma_names = {'OB_gamma','HPC_gamma','PFC_gamma','ACx_gamma'};
    cfg.theta_ch    = 21;     % HPC
    cfg.theta_name  = 'HPC_theta';
    cfg.delta_ch    = 11;     % OB
    cfg.delta_name  = 'OB_delta';

elseif contains(datapath, 'Shropshire')
    % [21 18 12 65] -> [OB HPC PFC ACx]
    cfg.gamma_ch    = [21 18 12 65];
    cfg.gamma_names = {'OB_gamma','HPC_gamma','PFC_gamma','ACx_gamma'};
    cfg.theta_ch    = 18;     % HPC
    cfg.theta_name  = 'HPC_theta';
    cfg.delta_ch    = 21;     % OB
    cfg.delta_name  = 'OB_delta';

elseif contains(datapath, 'Tvorozhok')
    % only OB
    cfg.gamma_ch    = 12;
    cfg.gamma_names = {'OB_gamma'};
    cfg.delta_ch    = 12;
    cfg.delta_name  = 'OB_delta';

elseif contains(datapath, 'Kosichka')
    % only OB
    cfg.gamma_ch    = 4;
    cfg.gamma_names = {'OB_gamma'};
    cfg.delta_ch    = 4;
    cfg.delta_name  = 'OB_delta';

else
    error('Unknown animal in datapath. Add its config in calc_brain_gamma_powers.');
end

%% ------------------------------------------------------------------------
%  2) Init BrainPower struct
% -------------------------------------------------------------------------

BrainPower = struct();
BrainPower.signal_names = {};
BrainPower.Power        = {};

%% ------------------------------------------------------------------------
%  3) Gamma power for each configured channel
% -------------------------------------------------------------------------

fs = 1024; % hard-coded in original code

for i = 1:numel(cfg.gamma_ch)
    ch = cfg.gamma_ch(i);
    sigName = cfg.gamma_names{i};

    lfpFile = fullfile(datapath, 'LFPData', sprintf('LFP%d.mat', ch));
    if ~exist(lfpFile,'file')
        warning('LFP file not found for channel %d: %s (skipping %s)', ch, lfpFile, sigName);
        continue
    end

    L = load(lfpFile);
    LFP = L.LFP;

    FilGamma = FilterLFP(LFP,[40 60],fs);
    envGamma = abs(hilbert(Data(FilGamma)));

    dt = median(diff(Range(FilGamma,'s')));
    win = ceil(sm_w/dt);
    smoothed = runmean(envGamma, win);

    BrainPower.signal_names{end+1} = sigName;
    BrainPower.Power{end+1} = tsd(Range(FilGamma), smoothed);

    clear L LFP FilGamma envGamma smoothed
end

%% ------------------------------------------------------------------------
%  4) Delta power
% -------------------------------------------------------------------------

if ~isempty(cfg.delta_ch) && ~isempty(cfg.delta_name)
    lfpFile = fullfile(datapath, 'LFPData', sprintf('LFP%d.mat', cfg.delta_ch));
    if exist(lfpFile,'file')
        L = load(lfpFile);
        LFP = L.LFP;

        FilDelta = FilterLFP(LFP,[0.5 4],fs);
        envDelta = abs(hilbert(Data(FilDelta)));

        dt = median(diff(Range(FilDelta,'s')));
        win = ceil(sm_w/dt);
        smoothed = runmean(envDelta, win);

        BrainPower.signal_names{end+1} = cfg.delta_name;
        BrainPower.Power{end+1} = tsd(Range(FilDelta), smoothed);

        clear L LFP FilDelta envDelta smoothed
    else
        warning('Delta channel %d file not found: %s', cfg.delta_ch, lfpFile);
    end
end

%% ------------------------------------------------------------------------
%  5) Theta / delta ratio (theta power / delta power)
% -------------------------------------------------------------------------

if ~isempty(cfg.theta_ch) && ~isempty(cfg.theta_name)
    lfpFile = fullfile(datapath, 'LFPData', sprintf('LFP%d.mat', cfg.theta_ch));
    if exist(lfpFile,'file')
        L = load(lfpFile);
        LFP = L.LFP;

        Frequency{1} = [3 6];   % theta
        Frequency{2} = [0.2 3]; % delta-ish

        FilTheta = FilterLFP(LFP,Frequency{1},fs);
        FilDelta = FilterLFP(LFP,Frequency{2},fs);

        hilbert_theta = abs(hilbert(Data(FilTheta)));
        hilbert_delta = abs(hilbert(Data(FilDelta)));
        hilbert_delta(hilbert_delta < 10) = 10;  % avoid division by tiny numbers

        theta_ratio = hilbert_theta ./ hilbert_delta;
        ThetaRatioTSD = tsd(Range(FilTheta), theta_ratio);

        dt = median(diff(Range(ThetaRatioTSD,'s')));
        win = ceil(sm_w/dt);
        smoothed = runmean(Data(ThetaRatioTSD), win);

        BrainPower.signal_names{end+1} = cfg.theta_name;
        BrainPower.Power{end+1} = tsd(Range(ThetaRatioTSD), smoothed);

        clear L LFP FilTheta FilDelta hilbert_theta hilbert_delta theta_ratio ThetaRatioTSD smoothed Frequency
    else
        warning('Theta channel %d file not found: %s', cfg.theta_ch, lfpFile);
    end
end

%% ------------------------------------------------------------------------
%  6) Save
% -------------------------------------------------------------------------

outFile = fullfile(datapath, 'SleepScoring_OBGamma.mat');

if exist(outFile,'file')
    save(outFile, 'BrainPower', '-append', '-v7.3');
else
    save(outFile, 'BrainPower', '-v7.3');
end

fprintf('Saved BrainPower to %s\n', outFile);
disp('calc_brain_gamma_powers: done.');

end