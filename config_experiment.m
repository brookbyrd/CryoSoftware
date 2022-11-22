function[config, flatfield] = config_experiment(config)
%% CONFIG_EXPERIMENT Function dedicated to setting up all parameters for running all Bioslice data processing sequences
%%
% * INPUTS-*
% * CONFIG: returns a config variable which has subfields of experimental
%  parameters necessary to process macrotome data in the main loop. 
% * FLATFIELD: returns the flatfield variable which has stored flatfield
% images for each channel at every depth. 
%%
% *Author: Brook Byrd $ $Date: 6-11-19 $*
% *Copyright: Dartmouth College*

%% SET UP FILEPATH CONFIGURATIONS
% Set up filepath dependencies to use during central main loop processing
% code
currentPath = pwd;
config.currentPath = currentPath; 
addpath(config.currentPath);
addpath('color');
addpath('SegmentationFit');

 % load in the camera barrel-distortion parameters, last calibrated 7-3-19
config.params = load('distortion_param.mat','params');

% Select raw data for a new configuration
disp('Select all folder path(s) for this data');
config.filePath = string(uigetdir2()); % Select the filepaths each time bc it can change

% Name the specimen being processed (autofills as folder selection name) 
if ~isfield(config, 'studyName')
    fileSplit = split(config.filePath,filesep);
    studyName = inputdlg({'Name this study: '},'Channel selection',[ 1 35],string(fileSplit(end)));
    config.studyName = string(studyName);
end
 

%% SELECT CHANNELS TO PROCESS
if ~isfield(config, 'channels')
    cd(config.filePath(1));
    direct = dir('20*');
    cd(direct(1).name)
    
    % SEARCH FOR ALL CHANNEL FOLDERS IN THE MAIN DATA DIRECTORY
    channelDir = dir(); keepconfig.channels = string;
    
        for ch = 1:length(channelDir)
            match = string(intersect(channelDir(ch).name, {'405','470', '530', '532', '635', '760','WL'}));
            if (~isempty(match))
                keepconfig.channels = [keepconfig.channels, match];
            end
        end
    % Combine possible channels into one string
    allOneString = sprintf('%s,' , keepconfig.channels);
    allOneString = string(allOneString(2:end-1)); % get rid of the comma
    
    % Provide the string of channels to to user to allow them to chose in a
    % text box
    input = inputdlg({'Enter channels to process'},'Channel selection',[ 1 35],string(allOneString));
    channelArray = cellfun(@(x)regexp(x,',','split'),input,'UniformOutput',0);
    config.channels = string(channelArray{1,1});
    
end

%% DISCOVERY ALL WAVELENGTHS USED
% Find lambda for each channel and save into configurations for later
config.numberOfSlices = zeros(length(config.channels),1);

% Iterate across all channels
for c = 1:length(config.channels)
    
    % Iterate across all data folders previously selected by the user
    for f = 1:length(config.filePath)
        
        % Direct into the raw data files 
        cd(char(strcat(config.filePath(f))));
        date = dir('*20*');
        for d = 1:length(date)
            
            try
                % cd to raw data files under a given date, and a given
                % channel name
                cd(char(strcat(config.filePath(f)))); cd(date(d).name); cd(config.channels(c));
                files =  dir('*.tif');
                fileInfo = {};
                
                %  Iterate through the files, and populate fileInfo using
                %  the parsed data from each file's name.
                sliceCount = 1;
                for slice = 1:length(files)
                    % There should be 8 pieces of data
                    if(length(split(files(slice).name,'_')) == 8)
                        fileInfo(sliceCount,:) = split(files(slice).name,'_');
                        sliceCount = sliceCount + 1;
                    end
                end
                
                % Pull out the important bits of detail from all files in
                % the fileInfo string matrix.
                if ~isempty(fileInfo)
                    
                    %Convert file info into strings that can be called on
                    fileInfo = string(fileInfo);
                    
                    % pull out the unique waves collected in file
                    wave = unique(fileInfo(:,3));
                    
                    % pull out the unique SpatialCalib waves
                    config.numberOfSlices(c,1) = config.numberOfSlices(c,1) + length(unique(fileInfo(:,2)));
                    
                    % pull out starting slice
                    if(d == 1)
                        config.startingSlice = cellfun(@(x)sscanf(x,'Slice%d'),fileInfo(1,2));
                    end
                    
                    %record the wavelengths as a channel-specific field
                    field = strcat('lambda_',char(config.channels(c)));
                    wavelengths.(field) = cellfun(@(x)sscanf(x,'LCTF%dnm'),wave);
                    
                end
                
                % clear file info before going into a new raw data folder 
                clear wave; clear fileInfo; clear files;
            catch
            end
            
        end
    end
    
    % set up the wavelengths 
    config.wavelengths = wavelengths;
    cd(currentPath);
    
end

%% FEED IN ALL EXPERIMENTAL PARAMETERS
% Give the slice thickness
if ~isfield(config, 'sliceThickness')
    input = inputdlg({'Enter slice thickness'},'Slice thickness (um)',[ 1 35],{'100'});
    config.sliceThickness = str2num(char(input));
end

% Give the starting depth of the data to do a proper FF correction
if ~isfield(config, 'startDepth')
    input = inputdlg({'Enter slice starting depth'},'Start depth (um)',[ 1 35],{'9000'});
    config.startDepth = str2num(char(input)) - 1; % Accounts for an indexing error in LABview
end

% Give the starting slice of the data
input = inputdlg({'Enter starting slice number'},'Slice #:',[ 1 35],{char(num2str(config.startingSlice))});
config.startingSlice = str2num(char(input));

% Set the range for integration
newRange = setRange_app(config.channels);
waitfor(newRange,'running',0);

% Gather the agents and ranges from the setRange app
config.agents = newRange.agents;
config.range = newRange.range;
delete(newRange);

% Set which rep to use per channel
config.rep = 2*ones([1,length(config.channels)]);

% Set up for possible slice-specific transforms that need to be added due to physical shifts that occured in the OCT block mounted on the stage during H&E-slide recovery. (10/20/22)
% config.sliceTransform(:,:,range) = repmat([1 0 0; 0 1 0; tx ty 1],1,1,length(range));
% Positive ty is down and positive tx moves it right
config.sliceTransform = repmat(eye(3,3),1,1,config.numberOfSlices(1));
 
close all;
%%  SET UP THE SAVEPATH FOR PROCESSED DATA

% Identify the location to save this processed data 
if ~isfield(config, 'savePath')
    % Selecting where to save the data
    disp('Select where to save processed data');
    % Collect folder info for reading and writing
    config.savePath =  string(uigetdir2());
    
end

% Make save path directory and folder
cd(strcat(config.savePath)); mkdir(char(strcat(config.studyName,'_processed')));
cd(char(strcat(config.studyName,'_processed')));
config.newSavePath = pwd;

% Create a folder structure to deposite processed image stacks
mkdir('RGB'); cd('RGB'); mkdir('RGB_stacks');


%% OPTIONAL CROPPING
currentPath = pwd;

% Option to select a raw file for cropping 
config.cropQ = questdlg('Crop image?');
if(strcmp(config.cropQ,'No'));
    config.channel_rect = [0,0,2048,2048];
    config.channel_roi = ones(2048);
elseif ~isfield(config, 'channel_rect')
    cd(config.filePath(1));
    disp('Select an image to crop');
    [imFile,impath] = uigetfile('.tif');
    im = imread(strcat(impath,imFile));
    im = undistortImage(im,params,'OutputView','same');
    % Save the cropped image parameters to be applied to every processed
    % image.
    [ config.channel_roi,  config.channel_rect] = imcrop(mat2gray(im));
    close all;
end

cd(currentPath);


%% SELECT NIR-I CAMERA TRANSFORM
% Optional NIR-I Channel transform selection to set-up the transform for registering the NIR-I camera to the visible camera's FOV. 
cd(currentPath);
if (~isempty(find(config.channels == "760")))
    
    % Offer the user the ability to hand select an NIR-I registration
    config.tform_Q = questdlg('Select a different NIR-I camera registration transform?');
    if(strcmp(config.tform_Q ,'No'))
        config.transformFile = ('tform_760.mat');  %Use the default transform 
    else
        % Allow the user to navigate to find a different registration transform file
        disp('Select an IR camera registration transform .mat file');
        [tFile,tpath] = uigetfile('.mat');
        
        % Save the transform file as part of the config file info. 
        config.transformFile = (strcat(tpath,tFile));  
    end
end
cd(currentPath);

%% SETUP FLATFIELD PATH
config.flatfieldQ = questdlg('Perform flatfield corrections?');
if(strcmp(config.flatfieldQ ,'Yes'))
    
    % Give the option to select a different folder
    config.FF_Q = questdlg('Select a different folder for flatfield calibration data?');
    if(strcmp(config.FF_Q ,'No'))
        config.FFfilePath = config.filePath(1);
    else
        disp('Select Flatfield folder directory');
        cd(config.filePath(1));
        config.FFfilePath = string(uigetdir2());
    end
end

%% ANALYZE SPATIAL CALIBRATION
config.spatialCalibQ = questdlg('Perform Color Card calibrations?');
if(strcmp(config.spatialCalibQ,'Yes'))
    
    % Allow user to select a specific color card folder
    config.Spatial_Q = questdlg('Is the Color Card data in a different location as FF calibration data?');
    if(strcmp(config.Spatial_Q ,'No'))
        config.SQfilePath = config.FFfilePath;
    else
        cd(config.filePath);
        disp('Select Spatial Calibration folder directory');
        config.SQfilePath = string(uigetdir2());
    end
    
    % run spatial calibration program
    [config] = process_SpatialCalib(config);
    
end
drawnow;
disp('User input complete');

%% SAVE ALL CONFIGURATIONS
disp('Saving configurations');
save(char(strcat(config.newSavePath,filesep, config.studyName,'_config.mat')), 'config', '-v7.3');

%% PROCESS FLATFIELD
currentPath = pwd;
addpath(currentPath);
if ~exist('flatfield')
    
    % %%Interpolate Flatfield
    if(strcmp(config.flatfieldQ,'Yes'))
        
        % Generate the flatfield for each channel using
        % interpolateFlatfield_correctDepth.m
        for c = 1:length(config.channels)
            disp(strcat('Interpolating Channel ',config.channels(c),' Flatfield'));
            field = char(strcat('Channel_',config.channels(c),'_mat'));
            [ff_mat] = interpolateFlatfield_correctDepth(config, config.channels(c), currentPath);
            flatfield.(field) = (ff_mat);
        end
    else
        % If a flatfield folder is not selected, just generate a flatfield
        % which consists of all ones.
        for c = 1:length(config.channels)
            field = char(strcat('Channel_',config.channels(c),'_mat'));
            flatfield.(field) = ones(size(config.channel_roi,1),size(config.channel_roi,2),config.numberOfSlices(c));
        end
        
    end
    
    % Write out the processed flatfield to keep record
    writeProcessedFlatfield(config,flatfield);
    
    % Save flatfield variable as a .mat file
    save(char(strcat(config.newSavePath,filesep, config.studyName,'_flatfield.mat')), 'flatfield', '-v7.3');
    
end

end