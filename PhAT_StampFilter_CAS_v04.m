% Code written by Carrie Ann Stine
% Michael Bruchas Lab - UW

% Imports a PhAT compatible .csv file of behavioral time stamps and filters
% them down to a desired subset of time stamps (based on input parameters).

% Writes a new PhAT compatible .csv file containing trimmed/desired time 
% stamps.


% UPDATES
% 11/1/23: v03 Updated to filter BORIS formatted .csv files

% 11/6/23: v04 Updated to identify point (discrete) vs state (duration 
% period) TTL time stamps

%% Reset MatLab workspace - clears all variables and the command window

clear variables;  % clear all variables
close all;  % close all open graphs
clc   % clear command window

%% Select trials for extraction

% Enter the row numbers from your nameGenerator file of the trials you want to extract
trials = [17 18];

%%%%%%%%%%%%%%%%%%%%%%%%%%  EDIT THESE FIELDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set source parameters
% IF INPUT IS FROM A PHAT MAT FILE WITH TTLS
% Specify whether incoming data is from a PhAT formatted mat file or a
% manually entered txt file
input.phatmat = 1; % 1 = stamps from PhAT mat file, 0 = other format 
% Define which time stamp series trim should be applied to.
input.series = {'TTL_1', 'TTL_2'}; 

% IF INPUT IS A SINGLE TEXT FILE
input.txt = 0; % 1 = stamps from a manually generated txt file, 0 = other 
    % Set the name of the behavior marked by the txt file stamps
    input.beh = []; % if input.txt = 0, set to [];
    
    % Specify whether txt file to filter is txt file named in
    % spreadsheet (auto = 1) or if you will manually input txt file
    % name (manual = 1)
    input.txt_auto = 0;
    input.txt_manual = 0;
        % Set the name of the source txt file if  input.txt_manual = 1
        input.txtfile = []; 
        % if input.txt_manual = 0, set to [];

%% Set trimming parameters
% Set minimum session time to be reached before stamps are kept. Stamps 
% before this time will be excluded (good for cases where TTL triggers
% when arduino program starts, this will trim out the start TTL)
trim.minstart = 50; % to ignore this operation, set to [];

% Set minimum time in seconds that needs to pass between kept stamps
trim.interval = 10; % to ignore this operation, set to [];

% Set minimum distance from end of session (in seconds) that last time
% stamp must be.
trim.laststamp = 30; % to ignore this operation, set to [];

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
    
    % create variable containing current trial's session length
    sessionlength = str2double(tempID(1,8))*60;
    
    % create file names to save mat, figure, and TTL files
    filename1 = strjoin(tempID(3:5),'');
    filename = char(filename1);
    
    matname = append(filename,'.','mat');
    
    PhATmatname = append('PhAT_',filename,'.','mat');
    PhATcsvname = append('PhAT_',filename,'.','csv');

    Cam1textname = append(filename,'_Cam1_stamps','.','txt');
    % set file name for new trimmed cam1 stamps
    Cam1trimname = append(filename,trim.name, '_Cam1_stamps','.','txt');
    
    
        TTL1_beh = char(tempID(16));
    TTL1textname = append(filename, '_',TTL1_beh,'_stamps','.','txt');
    % set file name for new trimmed TTL1 stamps
    TTL1trimname = append(filename,trim.name, '_',TTL1_beh,'_stamps','.','txt');
    
        TTL2_beh = char(tempID(17));
    TTL2textname = append(filename, '_',TTL2_beh,'_stamps','.','txt');
    % set file name for new trimmed TTL1 stamps
    TTL2trimname = append(filename,trim.name, '_',TTL2_beh,'_stamps','.','txt');
    
    % set file name for new trimmed PhAT stamps (alternative format)
    PhAT_alttrimname = append('PhAT_',filename,trim.name,'_stamps','.','csv');
    
    % set file name for new trimmed PhAT stamps (BORIS format)
    PhAT_BORIStrimname = append('PhAT_BORIS',filename,trim.name,'_stamps','.','csv');
    
    % set file name for new trimmed stamps from an input text file
    if input.txt == 1
        behtextname = append(filename,trim.name,'_',input.beh,'_stamps','.','txt');
    end
        
    %% Open the untrimmed time stamp source file
    % Load the PhAT mat file containing the block_stamps struct
    if input.phatmat == 1
        load(PhATmatname)
        trim.series = input.series;
    % Load the txt file containing the single time series
    elseif input.txt == 1
        % Load the text file, automatically from spreadsheet or manually
        if input.txt_auto == 1
            auto_txtfilename = char(tempID(1,6));
            rawts.TTL_1 = dlmread(auto_txtfilename);
        elseif input.txt_manual == 1
            rawts.TTL_1 = dlmread(input.txtfile);
        end
        % Specify number of time stamp series to be evaluated, set second
        % time stamp id to empty.
        rawts.TTL_2 = [];
        trim.series = {'TTL_1'};
    end
    
    %% Separate time series to be trimmed from PhAT mat file
    if input.phatmat == 1
        % Initialize variables
        rawts.Cam1 = []; rawts.TTL_1 = []; rawts.TTL_2 = [];

        % Check if 'Cam1', 'TTL1', and 'TTL2' are found in block_stamps.id
        Cam1_idx = find(strcmp(block_stamps.id, 'Cam1'));
        TTL1_idx = find(strcmp(block_stamps.id, 'TTL_1'));
        TTL2_idx = find(strcmp(block_stamps.id, 'TTL_2'));

        % Check if 'Cam1', 'TTL1', and 'TTL2' are found in trim.series
        Cam1_trim = find(strcmp(trim.series, 'Cam1'));
        TTL1_trim = find(strcmp(trim.series, 'TTL_1'));
        TTL2_trim = find(strcmp(trim.series, 'TTL_2'));

        % Extract time stamps for Cam1
        if ~isempty(Cam1_idx) && ~isempty(Cam1_trim)
            rawts.Cam1 = block_stamps.rawdata{Cam1_idx};
            beh_type.Cam1 = block_stamps.type{Cam1_idx};
        end

        % Extract time stamps for TTL1
        if ~isempty(TTL1_idx) && ~isempty(TTL1_trim)
            rawts.TTL_1 = block_stamps.rawdata{TTL1_idx};
            beh_type.TTL_1 = block_stamps.type{TTL1_idx};
        end

        % Extract time stamps for TTL2
        if ~isempty(TTL2_idx) && ~isempty(TTL2_trim)
            rawts.TTL_2 = block_stamps.rawdata{TTL2_idx};
            beh_type.TTL_2 = block_stamps.type{TTL2_idx};
        end
    end
   
    %% Perform trim operation on raw time stamps
    filt_start = {};
    filt_interval = {};
    filt_end = {};
    
    % Exclude time stamps happening prior to the min start cutoff
    for i = 1:numel(trim.series)
        % Set a temporary array containing the set of time stamps currently
        % undergoing evaluation
        temp_ts = rawts.(trim.series{i})(:,1);
        if ~isempty(trim.minstart)
           % Find the index of the first stamp above the min threshold
           temp_newstart_idx = find(temp_ts > trim.minstart,1,'first'); 
           % Save all stamps above the min threshold to a new array
           filt_start.(trim.series{i})(:,1) = temp_ts(temp_newstart_idx:end);
           
           % Save stop times for state (duration) epoch types
           if strcmp(beh_type.(trim.series{i}), 'state')
               filt_start.(trim.series{i})(:,2) = rawts.(trim.series{i})(temp_newstart_idx:end,2);
           end
           
       elseif isempty(trim.minstart)
           % Fill the first filtered array with all time stamps if no
           % min start filter was specified.
           filt_start.(trim.series{i}) = rawts.(trim.series{i});
       end
    end
    
    clear temp_ts
    
    % Exclude time stamps where the preceding stamp is within the interval
    % window
    for i = 1:numel(trim.series)
        % Set a temporary array containing the set of time stamps currently
        % undergoing evaluation
        temp_ts = filt_start.(trim.series{i})(:,1);
       if ~isempty(trim.interval) && ~isempty(temp_ts)
           % Initialize array of zeros the size of the time stamp series
           temp_keepidx = zeros(size(temp_ts));
           % Keep the first time stamp
           temp_keepidx(1) = 1;
           % Iterate through the timestamps and set temp_keepidx to 1 if 
           % the time difference is greater than trim.interval
           for k = 2:numel(temp_ts)
               % subtract the previous time stamp from the current time
               % stamp to determine size of their interval
               time_diff = temp_ts(k) - temp_ts(k-1);
               % mark time series index with 1 if the time_diff interval is
               % at least as large as the desired trim interval
               if time_diff > trim.interval
                   temp_keepidx(k) = 1;
               end
           end
           
           % Save stamps with difference greater than interval threshold
           filt_interval.(trim.series{i})(:,1) = temp_ts(temp_keepidx==1);
           
           % find the bout stop times for state (duration) behaviors
           if strcmp(beh_type.(trim.series{i}), 'state')
               % define a temporary array containing the stop times
               % undergoing evaluation
               temp_stop_ts = filt_start.(trim.series{i})(:,2);
               % make an array containing the indices of kept time stamps
               ts_kept = find(temp_keepidx==1);
               
               % initiate an array to contain the number of events in
               % between each bout
               kept_diff = zeros(numel(ts_kept),1);
               for s = 2:numel(ts_kept)
                   kept_diff(s-1) = ts_kept(s) - ts_kept(s-1);
               end
               % for last kept time stamp, subtract total number of time
               % stamps from last kept to identify where next actual time
               % stamp would theoretically start (one space after last
               % zero)
               kept_diff(end) = numel(temp_ts) - ts_kept(end) + 1;
               
               % add the identifier for the stop index (kept_diff) to the
               % identifier for the start index (ts_kept) and subtract 1 to
               % find the indices of the stop times
               temp_stopidx = ts_kept + (kept_diff-1);
               temp_stop_val = temp_stop_ts(temp_stopidx);
               % save stop times for each bout
               filt_interval.(trim.series{i})(:,2) = temp_stop_val;
           end
       else
           % Fill the second filtered array with all time stamps if no
           % interval filter was specified.
           filt_interval.(trim.series{i}) = filt_start.(trim.series{i});
       end
    end
    
    clear temp_ts       
    
    % Exclude time stamps happening after the max cutoff window
    for i = 1:numel(trim.series)
        % Set a temporary array containing the set of time stamps currently
        % undergoing evaluation
        temp_ts = filt_interval.(trim.series{i})(:,1);
        if ~isempty(trim.laststamp)
            trim.maxend = sessionlength - trim.laststamp;
            % Find the index of the last stamp below the max threshold
           temp_newend_idx = find(temp_ts < trim.maxend,1,'last'); 
           % Save all stamps below the max threshold to a new array
           filt_end.(trim.series{i}) = temp_ts(1:temp_newend_idx);
            
           % Save stop times for state (duration) epoch types
           if strcmp(beh_type.(trim.series{i}), 'state')
               filt_end.(trim.series{i})(:,2) = filt_interval.(trim.series{i})(1:temp_newend_idx,2);
           end
           
        elseif isempty(trim.laststamp)
           % Fill the last filtered array with all time stamps if no
           % last stamp filter was specified.
           filt_end.(trim.series{i}) = filt_interval.(trim.series{i});
        end
    end
    
    %% Update block_stamps to include new filtered stamps
    if input.phatmat == 1
        % Initialize new field to contain filtered time stamps
        block_stamps.filtdata = cell(size(block_stamps.id));
            % Populate first cell with time series data (doesn't change)
            block_stamps.filtdata{1} = block_stamps.rawdata{1};
            
        % Initialize new field to contain PhAT format time stamp series
        block_stamps.filt_traces_downsamp = zeros(size(block_stamps.traces_downsamp));
            % Populate first column with time series data (doesn't change)
            block_stamps.filt_traces_downsamp(:,1) = block_stamps.traces_downsamp(:,1);
            
        % Iterate through possible time stamp sources 
        for i = 2:numel(block_stamps.id)
            filt_identifier = block_stamps.id{i};
            % Check if the identifier exists in trim.series
            if any(strcmp(filt_identifier, trim.series))
                % If stamp id exists in trim.series and was filtered, 
                % populate filtdata with the new filtered stamps.
                block_stamps.filtdata{i} = filt_end.(filt_identifier);
            else
                % If it doesn't exist in trim.series, use the data from 
                % block_stamps.rawdata
                block_stamps.filtdata{i} = block_stamps.rawdata{i}; 
            end
        end
        
        % Generate the phAT format time stamp series using filtered
        % stamps
        Dts_downsampled = block_stamps.filt_traces_downsamp(:,1); 
        for i = 2:numel(block_stamps.id)
            % if epoch type is 'state' (has a duration), define time stamp
            % start and end blocks
            if strcmp(block_stamps.type{i}, 'state')
                % pull out start and stop times for each behavior
                ts_current_start = block_stamps.filtdata{i}(:,1);
                ts_current_stop = block_stamps.filtdata{i}(:,2);
                
                % run through each event stamp for the current epoch and find
                % closest corresponding start and stop index in Dts
                for t = 1:length(ts_current_start)
                    clear start_index; clear end_index;
                    [~,start_index] = min(abs(Dts_downsampled - ts_current_start(t)));
                    [~,stop_index] = min(abs(Dts_downsampled - ts_current_stop(t)));
               
                    if start_index == stop_index
                        stop_index = stop_index + 1;
                    end
                
                    % set the corresponding time points in trace array to 1
                    block_stamps.filt_traces_downsamp(start_index:stop_index,i) = 1;
                end
                
            else
           % define the current time stamp array based on which iteration 
           % the loop is on
            ts_current = block_stamps.filtdata{i};
        
                % run through each event stamp for the current epoch and find 
                % closest corresponding index in the time series
                for t = 1:length(ts_current)
                    % find index of closest value in Dts
                    [~,index] = min(abs(Dts_downsampled - ts_current(t)));

                    % set the corresponding element in trace array to 1
                    block_stamps.filt_traces_downsamp(index,i) = 1;
                end
            end
        end
    end
    
    
    %% Store filtered camera, epoch/time stamp data in CAS standard format
    if input.phatmat == 1
        % Initialize variables
        ts_Cam1 = []; ts_TTL1 = []; ts_TTL2 = [];
        
        % Extract time stamps for Cam1
        if ~isempty(Cam1_trim)
            ts_Cam1(:,1) = filt_end.Cam1(:,1);
        end
        
        % Extract time stamps for TTL1
        if ~isempty(TTL1_trim)
            ts_TTL1(:,1) = filt_end.TTL_1(:,1);
        end
        
        % Extract time stamps for TTL2
        if ~isempty(TTL2_trim)
            ts_TTL2(:,1) = filt_end.TTL_2(:,1);
        end
        
    elseif input.txt == 1
        % Extract time stamps for TTL1
        ts_TTL1(:,1) = filt_end.TTL_1(:,1);
    end
   
    %% Store filtered camera, epoch/time stamp data in PhAT BORIS form
    boris = {};
    boris.id = {}; % TTL numeric ID 
    boris.val = {}; % TTL store name
    boris.type = {}; % BORIS event type (state or point) store name
    boris.beh = {}; % name of behavior, from spreadsheet
    
    % pull epoch names from name Generator spreadsheet
    TTL_1_id = tempID{1,14};
    TTL_2_id = tempID{1,15};
    
    % pull epoch types from spreadsheet (state = event with duration, 
    % point = event with no duration)
    TTL_1_type = tempID{1,18};
    TTL_2_type = tempID{1,19};
    
    % save potential time stamp attributes in matched cells
    boris_id = {'TTL_1', 'TTL_2'};    
    boris_val = {TTL_1_id, TTL_2_id};
    boris_type = {TTL_1_type, TTL_2_type};
    boris_beh = {TTL1_beh, TTL2_beh};
   
    
    % iterate through potential boris_vals
    for i = 1:numel(boris_id)
        if input.phatmat == 1
            % check if the current boris id is in the filtered time stamp set.
            % Also check that any time stamps remain after filtering
            if any(strcmp(trim.series,boris_id{i})) && ~isempty(filt_end.(boris_id{i}))
                % append the non-'NA' value to the block_stream structs
                boris.id{end+1} = boris_id{i};
                boris.val{end+1} = boris_val{i};
                boris.type{end+1} = boris_type{i};
                boris.beh{end+1} = boris_beh{i};
            end
        elseif input.txt == 1 && ~isempty(filt_end.TTL_1)
            boris.id{1} = 'TTL_1';
            boris.val{1} = TTL_1_id;
            boris.type{1} = 'point';
            if isempty(input.beh)
                boris.beh{1} = TTL1_beh;
            else
                boris.beh{1} = input.beh;
            end
        end 
    end
 
    %% Create the BORIS table
    % Initialize empty boris table with 'Time', 'Behavior', and 'Status'
    % columns
    boris_table = table([], [], [], 'VariableNames',{'Time', 'Behavior', 'Status'});
    
    % populate boris table with raw time stamp data
    if ~isempty(boris.id)
        for i = 1:numel(boris.val)
            filt_identifier = boris.id{i};
            
            % if the behavior has start and stop times, evaluate both in
            % this loop
            if strcmp(boris.type{i}, 'state')
                % Get the raw onset and offset timestamps for the current boris.val
                boris.table_time{i} = [filt_end.(filt_identifier)(:,1); filt_end.(filt_identifier)(:,2)];
                
                % Store the number of onset and offset timestamps as a variable
                numel_start = size(filt_end.(filt_identifier)(:,1),1);
                numel_stop = size(filt_end.(filt_identifier)(:,2),1);
                
                % Create a cell array of 'Status' values containing 'START'
                % for all onset markers and 'STOP' for all offset
                boris.table_status{i} = [cellstr(repmat('START',numel_start,1)); cellstr(repmat('STOP',numel_stop,1))];
                
                % Create a cell array of 'Behavior' values the same size as
                % boris.table_time
                boris.table_behavior{i} = cellstr(repmat(boris.beh{i}, size(boris.table_time{i})));
            
            % if the behavior is a 'point behavior (no stop time), evaluate
            % using this loop
            else
                
                % Get the raw timestamps for the current boris.val
                boris.table_time{i} = filt_end.(filt_identifier)(:,1);
                
                % Create a cell array of 'Status' values with the same size as boris.rawdata{i}
                boris.table_status{i} = cellstr(repmat('POINT', size(boris.table_time{i})));
                
                % Create a cell array of 'Behavior' values with the same size as boris.rawdata{i}
                boris.table_behavior{i} = cellstr(repmat(boris.beh{i}, size(boris.table_time{i})));
                
            end
            % Create a table for the current boris.val and behavior
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
   
    
    %% Overwrite PhAT mat file to include new filtered stamps
    if input.phatmat == 1
        % define the PhAT data mat file name and location
        PhATmat_path = append(path_to_data.store,'/',path_to_data.path_PhATmat);
        PhATmat_name = append(PhATmat_path,'/',PhATmatname);
    
        % save streaming and time stamp data as a mat file
        if save_data == 1
            save(PhATmat_name, 'block_streams', 'block_stamps');
        end
    end
    
    %% Write new txt file to include new filtered stamps
    % define the CAS stamp txt file location
    stamp_textpath = append(path_to_data.store,'/',path_to_data.path_txtstamps);
    
    if input.phatmat == 1
        % for Cam1, TTL1, and TTL2, save non-empty arrays as a txt file
    
        if ~isempty(Cam1_trim)
            Cam1_name = append(stamp_textpath,'/',Cam1trimname);
            if save_data == 1
                writematrix(ts_Cam1,Cam1_name)
            end
        end
        
        if ~isempty(TTL1_trim)
            TTL1_name = append(stamp_textpath,'/',TTL1trimname);
            if save_data == 1
            	writematrix(ts_TTL1,TTL1_name)
            end
        end
        
        if ~isempty(TTL2_trim)
            TTL2_name = append(stamp_textpath,'/',TTL2trimname);
            if save_data == 1 
                writematrix(ts_TTL2,TTL2_name)
            end
        end
        
    elseif input.txt == 1
        TTL1_name = append(stamp_textpath,'/',behtextname);
        if save_data == 1
            writematrix(ts_TTL1,TTL1_name)
        end
    end
    
    %% Reset variables for next block in the loop
    % clear loop generated variables to free up memory space
    if a ~= trials(end)
        clearvars -except a trials nameGen path_to_data save_data input trim
    end
    
end
    