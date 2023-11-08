% Code written by Carrie Ann Stine
% Michael Bruchas Lab - UW

% Extracts photometry data and up to 2 TTLs from TDT Data Tank.
% Outputs photometry data as a .mat file and TTL time series as individual
% .txt files, in format compatible with my custom photometry analysis.

% Also saves data and TTLs in PhAT compatible .csv files, for use with the  
% Donaldson Lab photometry analysis GUI.

% REQUIRES TDTbin2mat.m file to be in the path!


% UPDATES
% 11/1/23: v04 Updated to save BORIS formatted TTLs instead of Alternative
% format

% 11/2/23: v05 Updated to save variables in Bruchas standard format

% 11/6/23: v06 Updated to identify point (discrete) vs state (duration 
% period) TTL time stamps


%% Reset MatLab workspace - clears all variables and the command window

clear variables;  % clear all variables
close all;  % close all open graphs
clc   % clear command window

%%%%%%%%%%%%%%%%%%%%%%%%%  EDIT THESE FIELDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select trials for extraction
% Enter the row numbers from your nameGenerator file of the trials you want to extract
trials = [2:5 7:10];

%% Setup locations of files and spreadsheets
%enter name of your excel sheet containing trial info
nameGenerator = 'nameGenerator_sample_v02.csv'; 

%change to factor you want to downsample data at
resampf = 100; 

 %change to 1 if you want to export files, change to 0 if you just want 
 % to look at data without saving mat and csv files/figures
save_data = 1;

% point to tank (enter file path)
path_to_data.tank = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/CodeToSHARE/20231102_PhAT_extract_CAS/sample_data/Tanks';

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
    path_to_data.path_casmat = 'CAS_matfiles_extracted';
    path_to_data.path_bruchasmat = 'Bruchas_matfiles_extracted';
    path_to_data.path_PhATmat = 'PhAT_matfiles';
    path_to_data.path_PhATstream = 'PhAT_stream_csv';
    path_to_data.path_PhATstamp = 'PhAT_stamp_csv';
    path_to_data.path_txtstamps = 'timestamps';
    
    
%add any other paths to other locations accessed by this script
% path_1 = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/CodeToSHARE/20231102_PhAT_extract_CS/functions';
    % addpath(path_1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% set up Tank name, variables to extract
tankdir = [path_to_data.tank];

% import csv file's trial info into a table
nameGen1 = readtable(nameGenerator);
nameGen = table2cell(nameGen1);
   
%% Run through loop to extract signal for each trial    

for a = trials  
    %% Set current block/trial identifiers
    curr_trial = a-1;
    
    % pull tank name from nameGenerator spreadsheet
    tankname = nameGen{curr_trial,1};
    
    % pull block name from nameGenerator spreadsheet
    blockname = nameGen{curr_trial,2};
    
    % define path to block (block = location of Synapse folder for a given trial)
    block_path = append(tankdir,'/',tankname,'/',blockname);
    
    
    % create variable containing current trial's row information
    tempID = string(nameGen(curr_trial,:));
    
    % create file names to save mat, figure, and TTL files
    filename1 = strjoin(tempID(3:5),'');
    filename = char(filename1);
    
    % create filename for CAS compatible analysis
    cas_matname = append(filename,'.','mat');
    
    % create filename for Bruchas lab compatible analysis
    bruchas_matname = append(filename,'_Bruchas','.','mat');
    
    % create filename for PhAT formatted mat + csv files (data and timestamps)
    PhATmatname = append('PhAT_',filename,'.','mat');
    PhATcsvname = append('PhAT_',filename,'.','csv');
    PhAT_boris_stampname = append('PhAT_BORIS',filename,'_stamps','.','csv');
    
    % create filename for raw data figure
    figname = append(filename,'.','svg');
    
    % extract name for behavior 1 from spreadsheet
    TTL1_beh = char(tempID(16));
    % create filename for behavior 1's time stamps as a text file
    TTL1textname = append(filename, '_',TTL1_beh,'_stamps','.','txt');
    
    % extract name for behavior 2 from spreadsheet
    TTL2_beh = char(tempID(17));
    % create filename for behavior 2's time stamps as a text fil
    TTL2textname = append(filename, '_',TTL2_beh,'_stamps','.','txt');
   
    
    
    %% Extract tdt streams/data
    tdt_data = TDTbin2mat(block_path);
        %Alldata = TDTbin2mat(block_path, 'TYPE' , {'streams'}); %alternate extraction that pulls only streams, not epocs/other data
    
    %% Define expected storenames 
    % initialize block_streams to contain stream ID names and corresponding values
    block_streams.id = {};
    block_streams.val = {};
    block_streams.color = {};
    
    % append each channel name with x, to match TDTbin2mat extracted format
    isos_id = append('x', tempID{1,11});
    green_id = append('x', tempID{1,12});
    red_id = append('x', tempID{1,13});
    
    % save potential channel information in matched cells (ie isosbestic, isos_id, and k all stored in idx 1)
    channel_id = {'Isosbestic', 'Green', 'Red'};
    channel_val = {isos_id, green_id, red_id};
    channel_color = {'k', 'b', 'r'};
    
    % iterate through potential channel_vals
    for i = 1:3
        % check if the value is not 'xNA'
        if ~strcmp(channel_val{i},'xNA')
            % append the non-'xNA' value to the block_stream structs
            block_streams.id{end+1} = channel_id{i};
            block_streams.val{end+1} = channel_val{i};
            block_streams.color{end+1} = channel_color{i};
        end
    end
    
    %% Store streaming channel data in a formatted struct
    % initialize block_streams.rawdata to store streaming data for corresponding channel
    block_streams.rawdata = {};
    block_streams.traces = [];
    
    block_length = zeros(1,size(block_streams.val,2));
    % find streaming length of each channel
    for i = 1:size(block_streams.val,2)
        block_length(i) = length(tdt_data.streams.(block_streams.val{i}).data);
    end
    
    % define length of streaming data to trim all streaming data to same length
    stream_length = min(block_length);
    
    % extract channel data
    % store each channel's data in the block_streams struct (transposed, save as column)
    for i = 1:size(block_streams.val,2)
        % save entire raw data set from Synapse, as cell
        block_streams.rawdata{i} = transpose(tdt_data.streams.(block_streams.val{i}).data);
        % save trimmed to equal length data trace from Synapse as array
        block_streams.traces(:,i) = transpose(tdt_data.streams.(block_streams.val{i}).data(1:stream_length));
    end
    
    % extract sampling rate using first channel (same across all channels so shouldn't matter which channel is used)
    Fs = tdt_data.streams.(block_streams.val{1}).fs;
    % generate total length of recording, in seconds. Rounded for Dts calc
    total_time = (stream_length - 1)/Fs;
    % generate Dts variable (transposed, save as column)
    Dts = (0:1/Fs:total_time)';
    
    % add Dts as a variable in the block_streams format
    block_streams.id{end+1} = 'Timestamp';
    block_streams.val{end+1} = 'Dts';
    block_streams.color{end+1} = 'NA';
    block_streams.rawdata{end+1} = Dts;
    block_streams.traces(:,end+1) = Dts;
    
    
    
    %% Store streaming data in PhAT form
    % define array formatted for PhAT analysis 
    PhAT_data_full = block_streams.traces;
    
    % downsample PhAT array so csv size is readable
        % initialize PhAT_array as an empty cell array
        PhAT_data_downsampled = cell(1, size(PhAT_data_full, 2));

    % Downsample each column of PhAT_array_full
    for i = 1:size(PhAT_data_full, 2)
        PhAT_data_downsampled{i} = downsample(PhAT_data_full(:, i), resampf);
    end
    
    % define associated column titles for PhAT array
    PhAT_stream_coltitles = block_streams.id;
    
    % store column titles and data in single array
    % convert downsampled cell data to an array
    PhAT_stream_array = cell2mat(PhAT_data_downsampled);
        % save downsampled streaming data to block_streams struct
        block_streams.traces_downsamp = PhAT_stream_array;
    
    % store data array as a table with appropriate column titles
    PhAT_stream_table = array2table(PhAT_stream_array,'VariableNames', PhAT_stream_coltitles);

    % clear PhAT_stream_array to save memory space
    clear PhAT_stream_array
    clear PhAT_data_full
    
    %% Define struct to store camera, epoch/time stamp data
    % initialize stamps to contain time stamp ID names and corresponding values
    block_stamps.id = {};
        block_stamps.id{1} = 'Time';
    block_stamps.val = {};
        block_stamps.val{1} = 'Dts';
    block_stamps.type = {};
        block_stamps.type{1} = 'NA';
    block_stamps.beh = {};
        block_stamps.beh{1} = 'NA';
    
    % pull epoch names from name Generator spreadsheet
    TTL_1_id = tempID{1,14};
    TTL_2_id = tempID{1,15};
    
    % pull epoch types from spreadsheet (state = event with duration, 
    % point = event with no duration)
    TTL_1_type = tempID{1,18};
    TTL_2_type = tempID{1,19};
    
    
    % pull name of behavior associated with epoch type from spreadsheet
    TTL_1_beh = tempID{1,16};
    TTL_2_beh = tempID{1,17};
    
    % save potential epoch information in matched cells
    ts_id = {'Cam1', 'TTL_1', 'TTL_2'};    
    ts_val = {'Cam1', TTL_1_id, TTL_2_id};
    ts_type = {'point', TTL_1_type, TTL_2_type};
    ts_beh = {'NA', TTL_1_beh, TTL_2_beh};
    
    % iterate through potential ts_vals, only store if a valid TTL is found
    for i = 1:3
        % check if the value is not 'NA'
        if ~strcmp(ts_val{i},'NA')
            % append the non-'NA' value to the block_stream structs
            block_stamps.id{end+1} = ts_id{i};
            block_stamps.val{end+1} = ts_val{i};
            block_stamps.type{end+1} = ts_type{i};
            block_stamps.beh{end+1} = ts_beh{i};
        end
    end
    
    %% Store camera, epoch/time stamp data in a formatted struct
    % initialize block_stamps.rawdata to store time stamps for corresponding channel as cell and array
    block_stamps.rawdata = {};
        % populate entire trace time series data
        block_stamps.rawdata{1} = Dts;
        
    % Intialize block_stamps.traces_downsamp to store time stamps as series
    % of 0s (null) and 1s (event occurred) across entire time array
    block_stamps.traces_downsamp = [];
        % populate entire trace time series data
        %block_stamps.traces(:,1) = ts;
        Dts_downsampled = downsample(Dts,resampf);
        block_stamps.traces_downsamp(:,1) = Dts_downsampled;
    
    % extract time stamp data
    % store camera and each TTL data in the block_stamps struct
    for i = 2:size(block_stamps.val,2)
        % save entire raw time series from Synapse, as a cell
        block_stamps.rawdata{i}(:,1) = tdt_data.epocs.(block_stamps.val{i}).onset;
        if strcmp(block_stamps.type{i}, 'state')
            block_stamps.rawdata{i}(:,2) = tdt_data.epocs.(block_stamps.val{i}).offset;
        end
    end

    % create an array of 0s for each actual time stamp series
    block_stamps.traces_downsamp(:,2:size(block_stamps.val,2)) = 0;
    
    % mark times where events occur with 1s for each time stamp type
    for i = 2:size(block_stamps.val,2)
        % if epoch type is 'state' (has a duration), define time stamp
        % start and end blocks
        if strcmp(block_stamps.type{i}, 'state')
            % pull out start and stop times for each behavior
            ts_current_start = block_stamps.rawdata{i}(:,1);
            ts_current_stop = block_stamps.rawdata{i}(:,2);
            
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
                block_stamps.traces_downsamp(start_index:stop_index,i) = 1;
                
            end
            
            
            
        else
            % define the current time stamp array based on which iteration
            % the loop is on
            ts_current = block_stamps.rawdata{i};
        
            % run through each event stamp for the current epoch and find 
            % closest corresponding index in Dts
            for t = 1:length(ts_current)
                % find index of closest value in Dts
                [~,index] = min(abs(Dts_downsampled - ts_current(t)));

                % set the corresponding element in trace array to 1
                block_stamps.traces_downsamp(index,i) = 1;
            end
        end
    end
  
    %% Store camera, epoch/time stamp data in BORIS format
    % Initialize boris struct to save time stamps in boris format
    % (Time, Behavior, and Status columns where Time only includes times
    % where event did occur)
    boris = {};
    boris.id = {}; % TTL numeric ID 
    boris.val = {}; % TTL store name
    boris.type = {}; % BORIS event type (state or point) store name
    boris.beh = {}; % name of behavior, from spreadsheet
    
    % save potential time stamp attributes in matched cells
    boris_id = {'TTL_1', 'TTL_2'};    
    boris_val = {TTL_1_id, TTL_2_id};
    boris_type = {TTL_1_type, TTL_2_type};
    boris_beh = {TTL1_beh, TTL2_beh};
   
 %%   
    % iterate through potential boris_vals
    for i = 1:numel(boris_id)
        % check if the value is not 'NA'
        if ~strcmp(boris_val{i},'NA')
            % append the non-'NA' value to the block_stream structs
            boris.id{end+1} = boris_id{i};
            boris.val{end+1} = boris_val{i};
            boris.type{end+1} = boris_type{i};
            boris.beh{end+1} = boris_beh{i};
        end
    end
    %% Create the BORIS table
    % Initialize empty boris table with 'Time', 'Behavior', and 'Status'
    % columns
    boris_table = table([], [], [], 'VariableNames',{'Time', 'Behavior', 'Status'});
        
        
    % populate boris.rawdata with raw time stamp data
    if ~isempty(boris.id)
        for i = 1:numel(boris.val)
            if strcmp(boris.type{i}, 'state')
                % Get the raw onset and offset timestamps for the current boris.val
                boris.table_time{i} = [tdt_data.epocs.(boris.val{i}).onset; tdt_data.epocs.(boris.val{i}).offset];
                
                % Store the number of onset and offset timestamps as a variable
                numel_start = size(tdt_data.epocs.(boris.val{i}).onset,1);
                numel_stop = size(tdt_data.epocs.(boris.val{i}).offset,1);
                
                % Create a cell array of 'Status' values containing 'START'
                % for all onset markers and 'STOP' for all offset
                boris.table_status{i} = [cellstr(repmat('START',numel_start,1)); cellstr(repmat('STOP',numel_stop,1))];
            
                % Create a cell array of 'Behavior' values the same size as
                % boris.table_time
                boris.table_behavior{i} = cellstr(repmat(boris.beh{i}, size(boris.table_time{i})));
            else
                % Get the raw timestamps for the current boris.val
                boris.table_time{i} = tdt_data.epocs.(boris.val{i}).onset;

                % Create a cell array of 'Status' values with the same data as boris.rawdata{i}
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
        block_stamps.boris = boris_table;
    end
          
    %% clear tdt_data variable to save memory space
    clear tdt_data
    
    %% Extract mat variables into CAS standard format
    % Initialize variables
    raw405 = []; raw470 = []; raw565 = []; rawsub405 = []; rawsub565 = [];
    
    % Check if 'Isosbestic', 'Green', and 'Red' are found in block_streams.id
    isosbestic_idx = find(strcmp(block_streams.id, 'Isosbestic'));
    green_idx = find(strcmp(block_streams.id, 'Green'));
    red_idx = find(strcmp(block_streams.id, 'Red'));

    % Extract data for Isosbestic
    if ~isempty(isosbestic_idx)
        raw405(:,1) = block_streams.traces(:, isosbestic_idx);
    end
    % Extract data for Green
    if ~isempty(green_idx)
        raw470(:,1) = block_streams.traces(:, green_idx);
    end
    % Extract data for Red
    if ~isempty(red_idx)
        raw565(:,1) = block_streams.traces(:, red_idx);
    end
    % Calculate Green - Isosbestic
    if ~isempty(isosbestic_idx) && ~isempty(green_idx)
        rawsub405(:,1) = raw470 - raw405;
    end
    % Calculate Green - Red
    if ~isempty(red_idx) && ~isempty(green_idx)
        rawsub565(:,1) = raw470 - raw565;
    end
    

    %% Store mat variables in CAS standard format
    % Store variables
    FPvarnames1 = {'Dts' 'raw470' 'raw405' 'raw565' 'rawsub405' 'rawsub565'};
    FPvariables1 = {Dts raw470 raw405 raw565 rawsub405 rawsub565};
    nonEmptyIndices = false(1, numel(FPvarnames1));
    %nonEmptyIndices = zeros(1,size(FPvarnames1,2));
    
    % Loop through each variable and check if it's a non-empty array
    for i = 1:numel(FPvarnames1)
        varName = FPvarnames1{i};
        varValue = FPvariables1{i};
        
        % Check if the current variable is empty
        if ~isempty(varValue)
            nonEmptyIndices(i) = true;
        end
    end
        
    % Trim FPvarnames and FPvariables to only include non-empty variables
    FPvarnames = FPvarnames1(nonEmptyIndices);
    
    FPvariables2 = FPvariables1(nonEmptyIndices);
        FPvariables = cell2mat(FPvariables2);
    
    % Store the non-empty variables and their ids in a struct 'FPvalues'
    for i = 1:numel(FPvarnames)
        FPvalues.(string(FPvarnames(i))) = FPvariables(:,i);
    end
    
    %% Store mat variables in Bruchas standard format
    % generate subdat variable (same as rawsub405)
    subdat = raw470-raw405;
    
    % find a polynomial fit for the isosbestic subtracted or raw signal 
    p1 = polyfit(Dts,subdat,4); %Fit to dF
    p2 = polyfit(Dts,raw470,4); %Fit to raw 470 For sensors
    % %Evaluate the polynomial on a finer grid and plot the results.
    
    % fit the polynomial curve to the time series
    fitcurve = polyval(p1,Dts);
    fitcurve2 = polyval(p2,Dts);
    
    % subtract the fitted curve to correct for photobleaching
    dataFilt = subdat - fitcurve; % Fit to isosbestic corrected
    dataFilt2 = raw470 - fitcurve2; % Fit to raw 470 (for sensors)
    
    % calculate dF/F
    % this gives deltaF/F with isobestic correction
    normDat1 = (dataFilt - median(dataFilt))./abs(median(raw405)); 
    normDat = normDat1.*100; % make a percentage
    data1 = normDat;
    
    % this gives deltaF/F for only using 470 data
    normDat2 = (dataFilt2 - median(dataFilt2))./abs(median(raw470)); 
    normDat2 = normDat2.*100;
    data2 = normDat2;
    
    %% Store epoch/time stamp data in CAS standard format
    % Initialize variables
    ts_TTL1 = []; ts_TTL2 = [];
    
    % Check if 'TTL1' and 'TTL2' are found in block_stamps.id
    TTL1_idx = find(strcmp(block_stamps.id, 'TTL_1'));
    TTL2_idx = find(strcmp(block_stamps.id, 'TTL_2'));

    % Extract time stamps for TTL1
    if ~isempty(TTL1_idx)
        ts_TTL1(:,1) = block_stamps.rawdata{TTL1_idx}(:,1);
    end
    
    % Extract time stamps for TTL2
    if ~isempty(TTL2_idx)
        ts_TTL2(:,1) = block_stamps.rawdata{TTL2_idx}(:,1);
    end
    
     %% Plot raw signal
     figure(1)
     hold on
     for i = 1:(size(block_streams.val,2) - 1) %exclude Dts
         plot(Dts_downsampled,block_streams.traces_downsamp(:,i),'Color',block_streams.color{i});
     end
     % plot signal minus isosbestic
     if ~isempty(isosbestic_idx) && ~isempty(green_idx)
         plot(Dts,rawsub405,'Color','g')
     end
     
     title('all channels (raw)')
     
     xlabel('time (s)')
     ylabel('amplitude')
     hold off
     
     % save figure to specified location
     fig_path = append(path_to_data.store,'/',path_to_data.path_figure);    
     full_file_fig = fullfile(fig_path,figname);
     
     if save_data == 1
        saveas(gcf,full_file_fig,'svg');
     end
     
    
    %% Save streaming data as a csv in PhAT format
    % define the PhAT data csv file name and location
    stream_path = append(path_to_data.store,'/',path_to_data.path_PhATstream);
        stream_name = append(stream_path,'/',PhATcsvname);
    
    % save the downsampled data as a csv file
    if save_data == 1
        writetable(PhAT_stream_table, stream_name);
    end
    
 
    %% Save time stamps as a csv in PhAT BORIS format
    % define the PhAT BORIS stamp csv file name and location
    boris_path = append(path_to_data.store,'/',path_to_data.path_PhATstamp);
        boris_name = append(boris_path,'/',PhAT_boris_stampname);
        
    % save the BORIS time stamp data as a csv file
    if save_data == 1 && ~isempty(boris_table)
        writetable(boris_table, boris_name);
    end
    
    %% Save streaming data and time stamps as a mat file in PhAT format
    % define the PhAT data mat file name and location
    PhATmat_path = append(path_to_data.store,'/',path_to_data.path_PhATmat);
        PhATmat_name = append(PhATmat_path,'/',PhATmatname);
    
    % save streaming and time stamp data as a mat file
    if save_data == 1
        save(PhATmat_name, 'block_streams', 'block_stamps');
    end
     
    
    %% Save streaming data as a mat file in CAS standard format
    % define the CAS data mat file name and location
    mat_path = append(path_to_data.store,'/',path_to_data.path_casmat);
        full_file_mat = fullfile(mat_path,cas_matname);
    
    % save streaming data as a mat file
    if save_data == 1
        save(full_file_mat,'Fs','FPvalues','FPvarnames');
    end
    
    
    %% Save streaming data as a mat file in Bruchas lab standard format
    % define the Bruchas data mat file name and location
    bruchas_mat_path = append(path_to_data.store,'/',path_to_data.path_bruchasmat);
        full_file_bruchas_mat = fullfile(bruchas_mat_path,bruchas_matname);
    
    % save streaming data as a mat file
    if save_data == 1
        save(full_file_bruchas_mat,'subdat','Dts', 'data1','raw470','data2','raw405', 'raw565'); %data1 is df/f fitted, data2 is only 470 fitted​
    end
    
    %% Save time stamps as a txt file
    % define the stamp txt file location
    stamp_textpath = append(path_to_data.store,'/',path_to_data.path_txtstamps);
    
    % for TTL1 and TTL2 save non-empty arrays as a txt file
    if ~isempty(TTL1_idx)
        TTL1_name = append(stamp_textpath,'/',TTL1textname);
        if save_data == 1
            writematrix(ts_TTL1,TTL1_name)
        end
    end
    
    if ~isempty(TTL2_idx)
        TTL2_name = append(stamp_textpath,'/',TTL2textname);
        if save_data == 1
            writematrix(ts_TTL2,TTL2_name)
        end
    end
    
    
    
    %% Reset variables for next block in the loop
    % clear loop generated variables to free up memory space
    if a ~= trials(end)
        close all
        clearvars -except a trials tankdir nameGen path_to_data resampf save_data
    end
end
   
