function Dir = PathForExperimentsReactPassiveEphys(ferret_name, novelty, daytime, sound_category)

% PathForExperimentsReactPassiveEphys - Retrieves session paths and metadata for ReactPassive Ephys experiments.
% You can select several ferrets and several stimuli_protocols at the same time

% Inputs:
%   ferret_name     - {'all', 'Kiri', 'Ficello', ...} or cell array of ferret names
%   novelty         - {'all', 'first', 'second'} (default: 'all')
%   daytime         - {'all', 'morning', 'evening'} (default: 'all')
%   sound_category  - {'all', 'ferret-music', 'music-speech', 'ferret-speech'} (default: 'all')
%
                                %%%%%%% Examples of use %%%%%%% 
% Dir = PathForExperimentsReactPassive('Kiri', 'all', 'all', {'ferret-music', 'ferret-speech'}); -- Will pull out all sessions with ferret-music and ferret-speech sounds for Kosichka
% Dir = PathForExperimentsReactPassive({'Kiri', 'Ficello'}, 'first', 'all', 'ferret-speech');      -- Will pull out first ferret-speech exposure sessions for  Edel and Chabichou
                                
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
    % Ficello
    'Kiri', 'first', 'evening', 'ferret-music', 'Z:\Arsenii\React_Passive_ephys\Processed_data\Kiri', {'20251215'};

    % Kiri
    'Kiri', 'first', 'evening', 'ferret-music', 'Z:\Arsenii\React_Passive_ephys\Processed_data\Kiri', {'20251211', '20251218', '20251229', '20260113'};
    'Kiri', 'first', 'evening', 'music-speech', 'Z:\Arsenii\React_Passive_ephys\Processed_data\Kiri', {'20251216_m', '20251224', '20260109', '20260121'};
    'Kiri', 'first', 'evening', 'ferret-speech', 'Z:\Arsenii\React_Passive_ephys\Processed_data\Kiri', {'20251213', '20251222', '20251231', '20260116'};

    'Kiri', 'second', 'morning', 'ferret-music', 'Z:\Arsenii\React_Passive_ephys\Processed_data\Kiri', {'20251212', '20251219', '20251230', '20260114'};
    'Kiri', 'second', 'morning', 'music-speech', 'Z:\Arsenii\React_Passive_ephys\Processed_data\Kiri', {'20251216_n', '20251226', '20260112', '20260122'};
    'Kiri', 'second', 'morning', 'ferret-speech', 'Z:\Arsenii\React_Passive_ephys\Processed_data\Kiri', {'20251215', '20251223', '20260102', '20260116'};
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
