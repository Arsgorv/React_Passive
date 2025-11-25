function calculate_spectrograms(datapath, animalName, channels)
% calculate_spectrograms(datapath, animalName, channels)
%
% INPUTS
%   datapath   : session folder where .lfp etc live
%   animalName : (optional) string, e.g. 'Shropshire','Edel','Brynza','Tvorozhok'
%   channels   : (optional) [AuCx OB HPC PFC] channel indices
%                1000 means "skip this structure"
%
% EXAMPLES
%   calculate_spectrograms('Z:\...\Kosichka\20251114_1_n_p05\ephys','Kosichka');
%   calculate_spectrograms(datapath, '', [65 21 18 12]);

if nargin < 1 || isempty(datapath)
    error('calculate_spectrograms: you must provide datapath');
end
if nargin < 2
    animalName = '';
end
if nargin < 3
    channels = [];
end

if ~exist(datapath,'dir')
    error('calculate_spectrograms: datapath does not exist: %s', datapath);
end

fprintf('\n--- calculate_spectrograms ---\n');
fprintf('datapath: %s\n', datapath);

%% 1) Default channels per animal
AnimalChannels = struct();
AnimalChannels.Shropshire = [65   21   18   12];    % [AuCx OB HPC PFC]
AnimalChannels.Edel       = [1000 26   24   1000];
AnimalChannels.Brynza     = [1    11   21   13];
AnimalChannels.Tvorozhok  = [1000 12   21   13];
AnimalChannels.Kosichka   = [1000  4 1000 1000];  % <- fill when known

%% 2) Determine channels

if isempty(channels)
    if ~isempty(animalName) && isfield(AnimalChannels, animalName)
        channels = AnimalChannels.(animalName);
        fprintf('Using channels for "%s": [AuCx OB HPC PFC] = [%d %d %d %d]\n', ...
            animalName, channels(1), channels(2), channels(3), channels(4));
    else
        % Try to guess from datapath
        fn = fieldnames(AnimalChannels);
        guessed = '';
        for k = 1:length(fn)
            if contains(datapath, fn{k}, 'IgnoreCase', true)
                guessed = fn{k};
                channels = AnimalChannels.(fn{k});
                fprintf('Guessed animal "%s" from path. Channels = [%d %d %d %d]\n', ...
                    guessed, channels(1), channels(2), channels(3), channels(4));
                break
            end
        end
        if isempty(channels)
            disp('Could not determine channels automatically.');
            disp('Enter channels as [AuCx OB HPC PFC], use 1000 to skip.');
            channels = input('channels = ');
        end
    end
else
    fprintf('Using user-provided channels [AuCx OB HPC PFC] = [%d %d %d %d]\n', ...
        channels(1), channels(2), channels(3), channels(4));
end

if numel(channels) ~= 4
    error('channels must be a 1x4 vector: [AuCx OB HPC PFC]');
end

AuCx_ch = channels(1);
OB_ch   = channels(2);
HPC_ch  = channels(3);
PFC_ch  = channels(4);

%% 3) Go to session folder

startDir = pwd;
cd(datapath);

fprintf('Running spectrograms in: %s\n', datapath);

%% 4) AuCx

if AuCx_ch ~= 1000
    if ~exist('AuCx_Low_Spectrum.mat','file')
        LowSpectrumSB(pwd, AuCx_ch, 'AuCx');
        disp('AuCx_Low done');
    else
        disp('AuCx_Low already exists');
    end

    if ~exist('AuCx_Middle_Spectrum.mat','file')
        MiddleSpectrum_BM(pwd, AuCx_ch, 'AuCx');
        disp('AuCx_Middle done');
    else
        disp('AuCx_Middle already exists');
    end
else
    disp('Skipping AuCx (channel = 1000)');
end

%% 5) OB (B)

if OB_ch ~= 1000
    if ~exist('B_UltraLow_Spectrum.mat','file')
        UltraLowSpectrumBM(pwd, OB_ch, 'B');
        disp('B_UltraLow done');
    else
        disp('B_UltraLow already exists');
    end

    if ~exist('B_Low_Spectrum.mat','file')
        LowSpectrumSB(pwd, OB_ch, 'B');
        disp('B_Low done');
    else
        disp('B_Low already exists');
    end

    if ~exist('B_Middle_Spectrum.mat','file')
        MiddleSpectrum_BM(pwd, OB_ch, 'B');
        disp('B_Middle done');
    else
        disp('B_Middle already exists');
    end

    if ~exist('B_High_Spectrum.mat','file')
        HighSpectrum(pwd, OB_ch, 'B');
        disp('B_High done');
    else
        disp('B_High already exists');
    end
else
    disp('Skipping OB (channel = 1000)');
end

%% 6) Hippocampus

if HPC_ch ~= 1000
    if ~exist('H_Low_Spectrum.mat','file')
        LowSpectrumSB(pwd, HPC_ch, 'H');
        disp('H_Low done');
    else
        disp('H_Low already exists');
    end

    if ~exist('H_Middle_Spectrum.mat','file')
        MiddleSpectrum_BM(pwd, HPC_ch, 'H');
        disp('H_Middle done');
    else
        disp('H_Middle already exists');
    end
else
    disp('Skipping HPC (channel = 1000)');
end

%% 7) PFC

if PFC_ch ~= 1000
    if ~exist('PFCx_Low_Spectrum.mat','file')
        LowSpectrumSB(pwd, PFC_ch, 'PFCx');
        disp('PFCx_Low done');
    else
        disp('PFCx_Low already exists');
    end

    if ~exist('PFCx_Middle_Spectrum.mat','file')
        MiddleSpectrum_BM(pwd, PFC_ch, 'PFCx');
        disp('PFCx_Middle done');
    else
        disp('PFCx_Middle already exists');
    end
else
    disp('Skipping PFC (channel = 1000)');
end

cd(startDir);
disp('calculate_spectrograms: done.');

end