function Dir = PathForExperimentsArousal(ferret_name, type, daytime)

% PathForExperimentsReactPassive - Retrieves session paths and metadata for ReactPassive experiments.
% You can select several ferrets and several stimuli_protocols at the same time

% Inputs:
%   ferret_name     - {'all', 'Kosichka', 'Ficello' ...} or cell array of ferret names
%   type         - {'all', 'normal', 'atropine', 'saline'} (default: 'all')
%   daytime         - {'all', 'morning', 'evening'} (default: 'all')

                                %%%%%%% Examples of use %%%%%%% 
% Dir = PathForExperimentsArousal('Ficello', 'all', 'atropine'); -- Will
% pull out all sessions with atropine injection for Ficello
% Dir = PathForExperimentsArousal({'Kosichka', 'Ficello'}, 'normal',
% 'all');      -- Will pull out normal arousal sessions for Kosichka et
% Ficello
                                
%%                              
% Default parameter handling
if nargin < 1 || isempty(ferret_name)
    ferret_name = 'all';
end
if nargin < 2 || isempty(type)
    novelty = 'all';
end
if nargin < 3 || isempty(daytime)
    daytime = 'all';
end

% Ensure ferret_name and sound_category are cell arrays
if ~iscell(ferret_name)
    ferret_name = {ferret_name};
end
if ~iscell(type)
    type = {type};
end

% Initialize Dir structure
Dir = struct();
Dir.path = {};
Dir.ExpeInfo = {};
Dir.name = {};

% basePath according to OS (windows or Linux)
if ispc
    basePaths.Chabichou = 'Z:\Arsenii\OB_fUS_Arousal\Processed_data\Chabichou';
    basePaths.Edel      = 'Z:\Arsenii\OB_fUS_Arousal\Processed_data\Edel';
    basePaths.Kosichka  = 'Z:\Arsenii\OB_fUS_Arousal\Processed_data\Kosichka';
    basePaths.Ficello  = 'Z:\Arsenii\OB_fUS_Arousal\Processed_data\Kosichka';
elseif isunix
    basePaths.Chabichou = '/home/arsenii/data5/Arsenii/OB_fUS_Arousal/Processed_data/Chabichou';
    basePaths.Edel      = '/home/arsenii/data5/Arsenii/OB_fUS_Arousal/Processed_data/Edel';
    basePaths.Kosichka  = '/home/arsenii/data5/Arsenii/OB_fUS_Arousal/Processed_data/Kosichka';
    basePaths.Ficello  = '/home/arsenii/data5/Arsenii/OB_fUS_Arousal/Processed_data/Ficello';
else
    error('Adjust Data paths');
end

% Session data (add paths and session names)
data = {
    
    % Kosichka
    
    'Kosichka', 'normal', 'morning', basePaths.Kosichka, {};
    'Kosichka', 'atropine', 'morning' , basePaths.Kosichka, {};
    'Kosichka', 'saline', 'morning', basePaths.Kosichka, {};
    
    'Kosichka', 'normal', 'evening', basePaths.Kosichka, {'20251117', '20251118', '20251125', '20251208', '20251209'}; %corrupted sessions:'20251124', '20251201'
    'Kosichka', 'atropine', 'evening' , basePaths.Kosichka, {};
    'Kosichka', 'saline', 'evening', basePaths.Kosichka, {};
    
    % Ficello
    'Ficello', 'normal', 'morning', basePaths.Ficello, {'20251212','20251216', '20260105', '20260121', '20260128', '20260203', '20260204', '20260211', '20260217', '20260316', '20260317'};
    'Ficello', 'atropine', 'morning',  basePaths.Ficello, {};
    'Ficello', 'saline', 'morning',  basePaths.Ficello, {'20260408', '20260409'};
    
    'Ficello', 'normal', 'evening', basePaths.Ficello, {'20260106','20260107', '20260109', '20260112', '20260116', '20260122', '20260206', '20260209', '20260212', '20260213', '20260216','20260319', '20260320', '20260323'};
    'Ficello', 'atropine', 'evening',  basePaths.Ficello, {'20260327', '20260401', '20260403', '20260410', '20260417'};
    'Ficello', 'saline', 'evening',  basePaths.Ficello, {'20260402', '20260414', '20260416'};
    

    };

% Filter sessions based on inputs
for i = 1:size(data, 1)
    if (any(strcmp('all', ferret_name)) || any(strcmp(data{i, 1}, ferret_name))) && ...
            (strcmp(type, 'all') || strcmp(data{i, 2}, type)) && ...
            (strcmp(daytime, 'all') || strcmp(data{i, 3}, daytime))
        
        % Add paths and names to Dir structure
        for j = 1:length(data{i, 5})
            session_name = data{i, 5}{j};
            session_path = fullfile(data{i, 4}, session_name);
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
