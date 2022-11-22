function main_loop(config,flatfield)
% SET UP FILE PATHS
% Setting up function pathscd
currentPath = pwd;
addpath(currentPath);
addpath('color');
addpath('SegmentationFit');

%% Load/process configurations
% Choose data files or load in pre-existing figures
if ~exist('config')
    config = struct;
    
    % SET-UP configurations for processing
    [config, flatfield] = config_experiment(config);
end

%% MAIN PROCESSING LOOP
% Run through all the channels, make integrated and hyperspectral stacks

% FOR each channel
for c = 1:length(config.channels)        
    processAllStacks(config, flatfield, config.channels(c), currentPath,2,config.startingSlice,0); % rep2
end

%% THRESHELD PROCESSING LOOP
% Run through all the channels, make integrated and hyperspectral stacks

% FOR each channel
for c = 1:length(config.channels)        
    processAllStacks_threshold(config, flatfield, config.channels(c), currentPath,2,config.startingSlice,0); % rep2
end

%% Run this loop For unsorted images only
for c = 1:length(config.channels)
   processAllStacks_unsorted(config, flatfield, config.channels(c), config.range(c,:), currentPath,2,config.startingSlice,0); % rep2
end

end
