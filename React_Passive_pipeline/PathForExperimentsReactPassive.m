function Dir = PathForExperimentsReactPassive(ferret_name, novelty, daytime, sound_category)

% PathForExperimentsReactPassive - Retrieves session paths and metadata for ReactPassive experiments.
% You can select several ferrets and several stimuli_protocols at the same time

% Inputs:
%   ferret_name     - {'all', 'Chabichou', 'Edel', 'Kosichka', ...} or cell array of ferret names
%   novelty         - {'all', 'first', 'second'} (default: 'all')
%   daytime         - {'all', 'morning', 'evening'} (default: 'all')
%   sound_category  - {'all', 'ferret-music', 'music-speech', 'ferret-speech'} (default: 'all')
%
                                %%%%%%% Examples of use %%%%%%% 
% Dir = PathForExperimentsReactPassive('Kosichka', 'all', 'all', {'ferret-music', 'ferret-speech'}); -- Will pull out all sessions with ferret-music and ferret-speech sounds for Kosichka
% Dir = PathForExperimentsReactPassive({'Chabichou', 'Edel'}, 'first', 'all', 'ferret-speech');      -- Will pull out first ferret-speech exposure sessions for  Edel and Chabichou
                                
%%                              
% Default parameter handling
if nargin < 1 || isempty(ferret_name)
    ferret_name = 'all';
end
if nargin < 2 || isempty(novelty)
    novelty = 'all';
end
if nargin < 3 || isempty(daytime)
    daytime = 'all';
end
if nargin < 4 || isempty(sound_category)
    sound_category = 'all';
end

% Ensure ferret_name and sound_category are cell arrays
if ~iscell(ferret_name)
    ferret_name = {ferret_name};
end
if ~iscell(sound_category)
    sound_category = {sound_category};
end

% Initialize Dir structure
Dir = struct();
Dir.path = {};
Dir.ExpeInfo = {};
Dir.name = {};

% Session data (add paths and session names)
data = {
    % Chabichou
    'Chabichou', 'first', 'morning', 'ferret-music', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {'20210329_1_m_C', '20210401_1_m_E', '20210408_1_m_H', '20210427_1_m_O', '20210429_1_m_Q', '20210505_1_m_S'};
    'Chabichou', 'first', 'morning', 'music-speech', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {'20210406_1_m_F', '20210413_1_m_I', '20210420_1_m_L', '20210422_1_m_N', '20210504_1_m_R', '20210512_1_m_U'};
    'Chabichou', 'first', 'morning', 'ferret-speech', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {'20210331_1_m_D', '20210407_1_m_G', '20210415_1_m_K', '20210421_1_m_M', '20210428_1_m_P', '20210511_1_m_T'};
    
    'Chabichou', 'first', 'evening', 'ferret-music', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {'20210414_1_n_J'};
    'Chabichou', 'first', 'evening', 'music-speech', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {};
    'Chabichou', 'first', 'evening', 'ferret-speech', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {};
    
    'Chabichou', 'second', 'morning', 'ferret-music', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {};
    'Chabichou', 'second', 'morning', 'music-speech', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {'20210414_2_m_I'};
    'Chabichou', 'second', 'morning', 'ferret-speech', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {};
    
    'Chabichou', 'second', 'evening', 'ferret-music', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {'20210330_2_n_C', '20210406_2_n_E', '20210409_2_n_H', '20210415_2_n_J', '20210427_2_n_O', '20210429_2_n_Q', '20210505_2_n_S'};
    'Chabichou', 'second', 'evening', 'music-speech', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {'20210407_2_n_F', '20210420_2_n_L', '20210422_2_n_N', '20210504_2_n_R', '20210512_2_n_U'};
    'Chabichou', 'second', 'evening', 'ferret-speech', 'Z:\Arsenii\React_Passive\Processed_data\Chabichou', {'20210401_2_n_D', '20210408_2_n_G', '20210416_2_n_K', '20210421_2_n_M', '20210428_2_n_P', '20210511_2_n_T'};
    
    % Edel
    'Edel', 'first', 'morning', 'ferret-music', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220419_1_m_C', '20220422_1_m_F', '20220519_1_m_U', '20220523_1_m_W'};
    'Edel', 'first', 'morning', 'music-speech', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220415_1_m_A', '20220420_1_m_D', '20220518_1_m_T', '20220520_1_m_V', '20220524_1_m_X'};
    'Edel', 'first', 'morning', 'ferret-speech', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220418_1_m_B', '20220421_1_m_E', '20220517_1_m_S'};
    
    'Edel', 'first', 'evening', 'ferret-music', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220429_1_n_J', '20220504_1_n_M', '20220509_1_n_O', '20220512_1_n_R'};
    'Edel', 'first', 'evening', 'music-speech', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220426_1_n_G', '20220503_1_n_L', '20220510_1_n_P'};
    'Edel', 'first', 'evening', 'ferret-speech', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220427_1_n_H', '20220428_1_n_I','20220502_1_n_K', '20220505_1_n_N', '20220511_1_n_Q'};
    
    'Edel', 'second', 'morning', 'ferret-music', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220430_2_m_J', '20220505_2_m_M', '20220510_2_m_O', '20220513_2_m_R'};
    'Edel', 'second', 'morning', 'music-speech', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220427_2_m_G', '20220504_2_m_L', '20220511_2_m_P'};
    'Edel', 'second', 'morning', 'ferret-speech', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220428_2_m_H', '20220429_2_m_I', '20220503_2_m_K', '20220506_2_m_N', '20220512_2_m_Q'};
    
    'Edel', 'second', 'evening', 'ferret-music', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220419_2_n_C', '20220422_2_n_F','20220519_2_n_U', '20220523_2_n_W'};
    'Edel', 'second', 'evening', 'music-speech', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220415_2_n_A', '20220420_2_n_D', '20220518_2_n_T', '20220520_2_n_V', '20220524_2_n_X'};
    'Edel', 'second', 'evening', 'ferret-speech', 'Z:\Arsenii\React_Passive\Processed_data\Edel', {'20220418_2_n_B', '20220421_2_n_E', '20220517_2_n_S'};

    % Kosichka
    'Kosichka', 'first', 'morning', 'ferret-music', basePaths.Kosichka, {};
    'Kosichka', 'first', 'morning', 'music-speech', basePaths.Kosichka, {'20251128_1_m_p14'};
    'Kosichka', 'first', 'morning', 'ferret-speech', basePaths.Kosichka, {'20251120_1_m_p06'};
    
    'Kosichka', 'first', 'evening', 'ferret-music', basePaths.Kosichka, {'20251117_1_n_p19', '20251124_1_n_p10', '20251215_1_n_p04'};
    'Kosichka', 'first', 'evening', 'music-speech', basePaths.Kosichka, {'20251114_1_n_p05', '20251208_1_n_p20'};
    'Kosichka', 'first', 'evening', 'ferret-speech', basePaths.Kosichka, {'20251201_1_n_p15'};
    
    'Kosichka', 'second', 'morning', 'ferret-music', basePaths.Kosichka, {};
    'Kosichka', 'second', 'morning', 'music-speech', basePaths.Kosichka, {};
    'Kosichka', 'second', 'morning', 'ferret-speech', basePaths.Kosichka, {'20251121_2_m_p06'};
    
    'Kosichka', 'second', 'evening', 'ferret-music', basePaths.Kosichka, {'20251118_2_n_p19', '20251125_2_n_p10'};
    'Kosichka', 'second', 'evening', 'music-speech', basePaths.Kosichka, {'20251208_1_n_p20'};
    'Kosichka', 'second', 'evening', 'ferret-speech', basePaths.Kosichka, {'20251202_2_n_p15'};
    
    % Ficello
    'Ficello', 'first', 'morning', 'ferret-music', basePaths.Ficello, {'20251212_1_m_p10', '20260105_1_m_p04'};
    'Ficello', 'first', 'morning', 'music-speech', basePaths.Ficello, {'20251204_1_m_p05'};
    'Ficello', 'first', 'morning', 'ferret-speech', basePaths.Ficello, {'20251215_1_m_p06', '20260128_1_m_p13'};
    
    'Ficello', 'first', 'evening', 'ferret-music', basePaths.Ficello, {'20260119_1_n_p19'};
    'Ficello', 'first', 'evening', 'music-speech', basePaths.Ficello, {'20260107_1_n_p20', '20260112_1_n_p14', '20260122_1_n_p8'};
    'Ficello', 'first', 'evening', 'ferret-speech', basePaths.Ficello, {'20260109_1_n_p15', '20260114_1_n_p28'};
    
    'Ficello', 'second', 'morning', 'ferret-music', basePaths.Ficello, {'20260121_2_m_p19'};
    'Ficello', 'second', 'morning', 'music-speech', basePaths.Ficello, {'20251205_2_m_p05', '20260113_2_m_p14'};
    'Ficello', 'second', 'morning', 'ferret-speech', basePaths.Ficello, {'20251216_2_m_p06'};
    
    'Ficello', 'second', 'evening', 'ferret-music', basePaths.Ficello, {'20260106_2_n_p04'};
    'Ficello', 'second', 'evening', 'music-speech', basePaths.Ficello, {'20260108_2_n_p20', '20260123_2_n_p8'};
    'Ficello', 'second', 'evening', 'ferret-speech', basePaths.Ficello, {'20260115_2_n_p28'};
    
    % Kiri
    };

% Filter sessions based on inputs
for i = 1:size(data, 1)
    if (any(strcmp('all', ferret_name)) || any(strcmp(data{i, 1}, ferret_name))) && ...
            (strcmp(novelty, 'all') || strcmp(data{i, 2}, novelty)) && ...
            (strcmp(daytime, 'all') || strcmp(data{i, 3}, daytime)) && ...
            (any(strcmp('all', sound_category)) || any(strcmp(data{i, 4}, sound_category)))
        
        % Add paths and names to Dir structure
        for j = 1:length(data{i, 6})
            session_name = data{i, 6}{j};
            session_path = fullfile(data{i, 5}, session_name);
            Dir.path{end+1} = session_path;
            Dir.name{end+1} = session_name;
            
            % Load ExpeInfo if available
            expe_info_path = fullfile(session_path, 'ExpeInfo.mat');
            if exist(expe_info_path, 'file')
                load(expe_info_path, 'ExpeInfo');
                Dir.ExpeInfo{end+1} = ExpeInfo;
            else
                Dir.ExpeInfo{end+1} = [];
                warning(['ExpeInfo.mat not found for session: ' session_name]);
            end
        end
    end
end

end
