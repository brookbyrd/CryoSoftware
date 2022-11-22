function  writeProcessedFlatfield(config, flatfield)

channels = config.channels;
filePath = config.filePath;
channel_rect = config.channel_rect;
params = config.params;
savePath = config.newSavePath;



% SET UP FILE PATHS
% Setting up function paths

currentPath = pwd;
addpath(currentPath);
addpath(strcat(currentPath,filesep, 'color'));


for c = 1:length(config.channels)
    channel = config.channels(c);
    
    mkdir(char(strcat(savePath,filesep, channel, filesep, channel, '_processed_flatfield_interp3')));
    
    disp(strcat('Writing Channel ',config.channels(c),' Flatfield'));
    field = char(strcat('Channel_',config.channels(c),'_mat'));
    flatfield_mat = flatfield.(field);
    
    for z = 1:size(flatfield_mat,3)
        
        imwrite(uint16(flatfield_mat(:,:,z)), char(strcat(savePath,'\',channel,'\', channel,'_processed_flatfield_interp3\Channel', channel, '_Slice_',num2str(z), '_','flatfield.tif')));
        slice = z
    end
end


end