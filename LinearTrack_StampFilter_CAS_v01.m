% Code written by Carrie Ann Stine
% Michael Bruchas Lab - UW

% Imports a linear track generated xlsx sheet of behavioral time stamps and 
% filters them down to a desired subset of time stamps (based on input 
% parameters).

% Writes a new PhAT compatible .csv file containing trimmed/desired time 
% stamps (BORIS format) as well as txt files of start times


%% Reset MatLab workspace - clears all variables and the command window

clear variables;  % clear all variables
close all;  % close all open graphs
clc   % clear command window

%% Select trials for extraction

% Enter the row numbers from your nameGenerator file of the trials you want to extract
trials = [20 21];

% Note: for this script, the row in your spreadsheet MUST have non-NA  
% values in the lintrack_vidname and lintrack_beheventname columns.

%%%%%%%%%%%%%%%%%%%%%%%%%%  EDIT THESE FIELDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set minimum session time to be reached before stamps are kept. Stamps 
% before this time will be excluded (good for cases where TTL triggers
% when arduino program starts for example, this will trim out the start
% TTL)
trim.minstart = []; % to ignore this operation, set to [];

% Set minimum time in seconds that needs to pass between kept stamps
trim.interval = []; % to ignore this operation, set to [];

% Set minimum distance from end of session (in seconds) that last time
% stamp must be.
trim.laststamp = []; % to ignore this operation, set to [];

% Create an identifier for the adjusted time stamps to be appended to 
% the filename 
trim.name = '_filtered';

% Set to 1 if you want to save the timestamp files
save_data = 1; % change to 0 if you just want to look at data in Matlab 

%% Set paths to time stamps
% point to timestamp and function locations (enter file path)
%enter name of your excel sheet containing trial info
nameGenerator = 'nameGenerator_sample_v02.csv'; 

% point to spreadsheet (enter file path)
path_to_data.spreadsheet = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/CodeToSHARE/20231102_PhAT_extract_CAS/name_spreadsheets';
    addpath(path_to_data.spreadsheet);
    
% point to folder where you want to store data
path_to_data.store = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/CodeToSHARE/20231102_PhAT_extract_CAS/data';
    addpath(genpath(path_to_data.store));
% identify the subfolders within your main extract data folder where the 
% different saved files should be written to. These names will be
% appended to the path_to_data.store filepath
    path_to_data.path_figure = 'figures_rawtrace';
    path_to_data.path_casmat = 'matfiles_extracted';
    %path_to_data.path_bruchasmat = 'Bruchas_matfiles_extracted';
    path_to_data.path_PhATmat = 'PhAT_matfiles';
    path_to_data.path_PhATstream = 'PhAT_stream_csv';
    path_to_data.path_PhATstamp = 'PhAT_stamp_csv';
    path_to_data.path_txtstamps = 'timestamps';
    path_to_data.path_linbehevents = 'lineartrack_behevents';

    
%enter name of your excel sheet containing naming info
nameGen1 = readtable(nameGenerator); 
nameGen = table2cell(nameGen1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for a = trials
   %% Set current block/trial identifiers
    curr_trial = a-1;
  
    % pull tank name from nameGenerator spreadsheet
    tankname = nameGen{curr_trial,1};
    
    % pull block name from nameGenerator spreadsheet
    blockname = nameGen{curr_trial,2};
    
    % create variable containing current trial's row information
    tempID = string(nameGen(curr_trial,:));
    
    % pull the linear track behavioral event video name
    lintrack_vidname = char(tempID(1,20));
    lintrack_eventname = char(tempID(1,21));
    
    % create variable containing current trial's session length
    sessionlength = str2double(tempID(1,8))*60;
    
    
    % create file name to open associated mat file
    filename1 = strjoin(tempID(3:5),'');
    filename = char(filename1);
    
    matname = append(filename,'.','mat');
    
    PhATmatname = append('PhAT_',filename,'.','mat');

    % set file name for new trimmed PhAT stamps (BORIS format)
    PhAT_BORIStrimname = append('PhAT_BORIS',filename,trim.name,'_stamps','.','csv');

    %% Open the linear track event spreadsheet
    % generate the path name to the behavioral events spreadsheet, format:
    % each sheet named for the behavior, col 1 = vidname, col 2 = start
    % time, col 3 = stop time.
    xlsx_filepath = append(path_to_data.store,'/', path_to_data.path_linbehevents,'/',lintrack_eventname);
    
    % get the sheet names for the xlsx file
    [~, sheetnames] = xlsfinfo(xlsx_filepath);
    
    % initialize a struct to store the behavioral event sheets
    allbehevents = cell2struct(cell(size(sheetnames)), sheetnames, 2);
    
    % read and save each sheet into the behevents struct
    for i = 1:numel(sheetnames)
        allbehevents.(sheetnames{i}) = readtable(xlsx_filepath, 'Sheet', sheetnames{i});
    end
    
    %% Separate out current trial's behavioral events
    % multiple different trials (videos) have events stored in the
    % behavioral events spreadsheet. Save a separate struct containing
    % only the current trial's events
    for i = 1:numel(sheetnames)
        clear temp_table;
        temp_table = table2cell(allbehevents.(sheetnames{i}));
        vidname_idx = zeros(size(temp_table,1),1);
        
        % find indices for rows matching the current trial's vidname
        for k = 1:numel(vidname_idx)
            if strcmp(temp_table{k,1},lintrack_vidname)
                vidname_idx(k) = 1;
            end
        end

        % save current trial's behavioral events in separate struct
        behevents.(sheetnames{i}) = cell2mat(temp_table((vidname_idx==1),2:3));
    end
    
    %% Perform trim operation on the raw time stamps
    filt_start = {};
    filt_interval = {};
    filt_end = {};
    
    % exclude time stamps happening prior to the min start cutoff
    for i = 1:numel(sheetnames)
        % set a temporary array containing the set of time stamps currently
        % undergoing evaluation
        clear temp_ts;
        temp_ts = behevents.(sheetnames{i});
        if ~isempty(trim.minstart)
           % find the index of the first stamp above the min threshold
           temp_newstart_idx = find(temp_ts(:,1) > trim.minstart,1,'first'); 
           % save all stamps above the min threshold to a new array
           filt_start.(sheetnames{i}) = temp_ts(temp_newstart_idx:end,:);
        
        % fill the first filtered array with all time stamps if no min
        % start filter was specified.   
        elseif isempty(trim.minstart)
           filt_start.(sheetnames{i}) = behevents.(sheetnames{i});
        end
    end
    
    % exclude time stamps where the preceding stamp is within the interval
    % window
    for i = 1:numel(sheetnames)
        % set a temporary array containing the set of time stamps currently
        % undergoing evaluation
        clear temp_ts;
        temp_ts = filt_start.(sheetnames{i});
        if ~isempty(trim.interval) && ~isempty(temp_ts)
            % initialize variables to store bout information
            combined_bouts = zeros(0,2);
            current_bout = temp_ts(1,:);
            
            % loop through each time stamp
            for k = 2:size(temp_ts,1)
                int_currentrow = temp_ts(k,:);
                int_starttime = int_currentrow(1);
                int_endtime = int_currentrow(2);
                
                % check if the current stamp is within the interval for the
                % current bout
                if (int_starttime - current_bout(2)) <= trim.interval
                    % extend the current bout end time
                    current_bout(2) = int_endtime;
                else
                    % save the current bout to the combined bouts array
                    combined_bouts(end+1,:) = current_bout;
                    % start a new bout with the current time stamp
                    current_bout = int_currentrow;
                end
            end
            % add the last bout to the combined_bouts array
            combined_bouts(end+1,:) = current_bout;
            % save the combined bouts to the filtered interval struct
            filt_interval.(sheetnames{i}) = combined_bouts;
            
        else
            % fill the second filtered array with all time stamps if no
            % interval filter was specified.
            filt_interval.(sheetnames{i}) = filt_start.(sheetnames{i});
        end
    end
    
    
    
    % exclude time stamps happening after the max cutoff window
    for i = 1:numel(sheetnames)
        % set a temporary array containing the set of time stamps currently
        % undergoing evaluation
        clear temp_ts;
        temp_ts = behevents.(sheetnames{i});
        if ~isempty(trim.laststamp)
           trim.maxend = sessionlength - trim.laststamp;
           % find the index of the last stamp below the max threshold
           temp_newend_idx = find(temp_ts(:,1) < trim.maxend,1,'last'); 
           % Save all stamps below the max threshold to a new array
           filt_end.(sheetnames{i}) = temp_ts(1:temp_newend_idx,:);
        
        % fill the first filtered array with all time stamps if no last
        % stamo filter was specified.   
        elseif isempty(trim.laststamp)
           filt_end.(sheetnames{i}) = filt_interval.(sheetnames{i});
        end
    end
    
    %% Store filtered time stamps as a text file for each behavior
    
    % Initialize empty variables to store the start times for each type of
    % behavior
    for i = 1:numel(sheetnames)
        if ~isempty(filt_end.(sheetnames{i}))
            ts.(sheetnames{i}) = filt_end.(sheetnames{i})(:,1);
            ts_name.(sheetnames{i}) = append(filename,'_',sheetnames{i},trim.name,'_stamps','.','txt');
        end
    end
    
    %% Store filtered time stamps in PhAT BORIS format
    boris = {};
    boris.id = {}; % name of behavior, from spreadsheet
    boris.type = {}; % BORIS event type (state or point) store name
   
    % save behavior name and type (all 'state' from linear track) to struct
    for i = 1:numel(sheetnames)
        if ~isempty(filt_end.(sheetnames{i}))
            boris.id{end+1} = sheetnames{i};
            boris.type{end+1} = 'state';
        end
    end
    
    % initialize empty BORIS table with 'Time', 'Behavior', and 'Status'
    % columns
    boris_table = table([], [], [], 'VariableNames',{'Time', 'Behavior', 'Status'});
    
    % run through all behaviors and save as a BORIS table
    if ~isempty(boris.id)
        for i = 1:numel(boris.id)
            filt_identifier = boris.id{i};

            % Get the raw onset and offset timestamps for the current boris.id
            boris.table_time{i} = [filt_end.(filt_identifier)(:,1); filt_end.(filt_identifier)(:,2)];

            % Store the number of onset and offset timestamps as a variable
            numel_start = size(filt_end.(filt_identifier)(:,1),1);
            numel_stop = size(filt_end.(filt_identifier)(:,2),1);
            
            % Create a cell array of 'Status' values containing 'START'
            % for all onset markers and 'STOP' for all offset
            boris.table_status{i} = [cellstr(repmat('START',numel_start,1)); cellstr(repmat('STOP',numel_stop,1))];
            
            % Create a cell array of 'Behavior' values the same size as
            % boris.table_time
            boris.table_behavior{i} = cellstr(repmat(boris.id{i}, size(boris.table_time{i})));
            
            % Create a table for the current behavior
            boris.table_ind{i} = table(boris.table_time{i}, boris.table_behavior{i}, boris.table_status{i}, 'VariableNames', {'Time', 'Behavior', 'Status'});
            
            % Concatenate the temporary table with the main data table
            boris_table = [boris_table; boris.table_ind{i}];
            
        end
        
        % Sort BORIS table by 'Time' column in ascending order
        boris_table = sortrows(boris_table,'Time');
        
        % Save BORIS table to block_stamps struct
        block_stamps.boris_filt = boris_table;
    end
    
    %% Save filtered time stamps as a csv in PhAT BORIS format
    % define the PhAT stamp csv file name and location
    boris_path = append(path_to_data.store,'/',path_to_data.path_PhATstamp);
        boris_name = append(boris_path,'/',PhAT_BORIStrimname);

    % save the downsampled time stamp data as a csv file
    if save_data == 1
        writetable(boris_table, boris_name);
    end
    
    %% Write new txt to include new filtered stamps
    % define the CAS stamp txt file location
    stamp_textpath = append(path_to_data.store,'/',path_to_data.path_txtstamps);
    
    % get the field names of the ts struct (txt file formatted stamps)
    ts_id = fieldnames(ts);
    
    % run through each field of the ts struct and save as a txt file
    for i = 1:numel(ts_id)
        ts_path.(ts_id{i}) = append(stamp_textpath,'/',ts_name.(ts_id{i}));
        if save_data == 1
            writematrix(ts.(ts_id{i}), ts_path.(ts_id{i}))
        end
    end
end


    
    
    
    
    
    
    
    
    