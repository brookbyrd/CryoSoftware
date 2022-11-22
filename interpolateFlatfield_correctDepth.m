function[interp_Frames] =  interpolateFlatfield_correctDepth(config, channel, currentPath)
%% interpolateFlatfield_correctDepth Function generates a makima interpolated 3D flatfield based on 3-plane sampled flatfield from the experimental data files.

%%
% * INPUTS-
% * CONFIG: experiment-specific config variable which has subfields of experimental
%  parameters necessary to generate flatfields.
% * CHANNEL: Channel which is being processed. Only one channel's flatfield
% is generated at a time.
% * CURRENTPATH: PATH WHICH HAS ALL NECESSARY FILEPATH DEPENDENCIES
%%
% * OUTPUTS-
% * INTERP_FRAMES: returns the flatfield variable which has stored flatfield
% images for each channel at every depth.


% *Author: Brook Byrd $ $Date: 4-18-20 $*
% *Copyright: Dartmouth College*

%% LOAD IN SAMPLED FLATFIELD IMAGES 
% Storage for interpolated flatfield frames
interp_Frames = [];

% Direct to file path
fNum = 1
cd(char(config.FFfilePath(fNum)));

% Collect all the date folders
date = dir('*20*');
for day = 1:length(date)
    dateNum(day) = datenum(date(day,1).name,'yyyymmmdd');
end
[dateNum,I] = sort(dateNum);

% Provide an update
disp(strcat('Processing the flatfield of: ', channel));

% Load in the NIR-I transform if it is the 760 channel
cd(currentPath);
if(strcmp(channel,"760"))
    load(strcat('tform_760.mat'));
    disp(strcat('loading: ','tform_',channel,'_default.mat'));
else
    load('tform_identity.mat');
end


% Assumes the FF was all taken on the first day of experimentation.
cd(char(strcat(config.FFfilePath(fNum), filesep , date(I(1)).name, filesep, 'Flatfield', filesep, channel)))
disp(strcat('Loading flatfields of channel: ',channel));

% Work around to make sure all files are in the correct format
raw_files =  dir('*.tif');
counter = 1; %
for f = 1:length(raw_files)
    %check to make sure the formatting is correct
    if(length(split(raw_files(f).name,'_')) == 8)
        files(counter) = raw_files(f); % remove any files that do not fit this format
        fileInfo(counter,:) = split(files(counter).name,'_');
        counter = counter + 1;
    end
end

%Convert file info into strings that can be called on
fileInfo = string(fileInfo);

% pull out the unique repitions
repitions = unique(fileInfo(:,8));

% pull out the last 3 slices of the acquired FF images
allSlices = [allSlices; unique(fileInfo(:,2))];
selectedSlices = unique(fileInfo(:,2));
if (length(selectedSlices) > 2)
    selectedSlices = selectedSlices(end - 2:end);
else
    selectedSlices = selectedSlices(end);
end

% pull out the unique LCTF waves for this particular channel
wave = unique(fileInfo(:,3));

%% PROCESS THREE SAMPLED FLATFIELD IMAGES 
% Iterate over each image, using the first wavelength as the basis for
% the flatfield image.

% Slice counting index
index = 1;

% Storage for sampled flatfield slices
allSlices = [];

for slice = 1:length(selectedSlices)
    
    % pull out the file indices for this slice
    sliceIndexArray = find(fileInfo(:,2) == selectedSlices(slice));
    
    % pull out the historical slice number
    sliceLabel = selectedSlices(slice);
    N = cellfun(@(x)sscanf(x,'Slice%d'),{sliceLabel});
    
    %Figure out the slice number and appropriate affine transforms
    z_microns =  5000 + 15000*(slice - 1);
    tform_z = affine2d([1.3e-6*z_microns + 0.99,0,0;0,1.3e-6*z_microns + 0.99,0; -0.0015*z_microns + 7.3,-0.0014*z_microns + 6.8,1]);
    
    % pull out the dark image from the last LCTF
    % wavelength used
    OffIndex = intersect(intersect(find(strcmp(fileInfo(:,7),'OFF')), find(strcmp(fileInfo(:,8),'rep2.tif'))),find(strcmp(fileInfo(:,2), selectedSlices(slice))));
    
    % just in case the file is corrupted
    imOff = imread(files(OffIndex(1)).name);
    
    % Gather the exposure time in milliseconds of the
    % specific instance and repititions
    D = regexp(fileInfo(OffIndex,6),'\d*','Match');
    
    %Exposure time expressed in ms so adjusting to for cps
    exptime  = str2num(char(D{1,1}))./1e3;
    
    % pull out the slice within the time set that has the
    % correct wavelength for the given LCTF wavelength
    OnIndex = intersect(intersect(intersect(find(strcmp(fileInfo(:,7),'ON')), find(strcmp(fileInfo(:,8),'rep2.tif'))),find(strcmp(fileInfo(:,3), wave(1)))),find(strcmp(fileInfo(:,2), selectedSlices(slice))));
    imOn = imread(files(OnIndex(1)).name);
    
    %perform basic image subtract
    temp_image =  ((double(imOn - imOff)./double(exptime)));%./double(flatfield_norm));
    
    % perform distortion correction based on barrel distortion
    temp_image = undistortImage(temp_image,config.params,'OutputView','same');
    
    %register all channels together
    temp_image = imwarp((temp_image), mytform, 'OutputView', imref2d(size(temp_image)));%, 'OutputView',imref2d(size(temp_image)));
    
    % adjust for scale and transform based on stage
    temp_image = imwarp((temp_image), tform_z, 'OutputView', imref2d(size(temp_image)));%, 'OutputView',imref2d(size(temp_image)));
    
    % apply smoothing function to take out the paper artifacts
    temp_image = imgaussfilt(temp_image,8);
    
    % crop down the flatfield to save on time
    temp_image_rect = imcrop(temp_image, config.channel_rect);
    
    % save the processed flatfield image in a stack
    flatfield(:,:,index) = (temp_image_rect);
    
    % progress through all selectedSlices
    index = index + 1;
    
end

%% INTERPOLATE VOLUME FROM PROCESSEd FLATFIELD IMAGES

% If there are three or more slices present, assume the standard three
% flatfields at 35,000, 20,000 and 5,000 have been acquired and proceed
% with Makima 3D interpolation
if(length(allSlices) > 2)
    interpSlices = [config.startDepth:config.sliceThickness:config.sliceThickness*config.numberOfSlices(1) - config.sliceThickness + config.startDepth];
    [Xq,Yq,Zq] = meshgrid(1:size(flatfield,2),1:size(flatfield,1),config.startDepth:config.sliceThickness:config.sliceThickness*config.numberOfSlices(1) - config.sliceThickness + config.startDepth);
    
    % Concatenate the 3 planes to create a volume for interpolation
    X = cat(3, Xq(:,:,1), Xq(:,:,1), Xq(:,:,1));
    Y = cat(3, Yq(:,:,1), Yq(:,:,1),  Yq(:,:,1));
    
    % Z is the height of the acquired flatfield
    Z = cat(3, 5000*ones([size(flatfield,1), size(flatfield,2)]),20000*ones([size(flatfield,1), size(flatfield,2)]), 35000*ones([size(flatfield,1),size(flatfield,2)]));
    
    % Feed in the last image and the two proceeding (assuming last is 5000)
    V =  cat(3,squeeze(flatfield(:,:,size(flatfield,3))),squeeze(flatfield(:,:,size(flatfield,3)-1)),squeeze(flatfield(:,:,size(flatfield,3)-2)));
    
    % Perform the makima interpolation
    interp_Frames = interp3(X,Y,Z,V,Xq,Yq,Zq,'makima');
    disp('Multi-slice FF detected');
    
% In the condition of having only one FF slice, just use that one FF as the whole volume   
else 
    disp('Single-slice FF detected');
    
    % no need for interpolation, just return the single flatfield image at
    % every depeth
    for n = 1:config.numberOfSlices(1)
        interp_Frames(:,:,n) = double(flatfield(:,:,size(flatfield,3)));
    end
    
end

% Provides an update for which channel is processed
disp(strcat('Completed channel: ', channel));
cd(currentPath);

end
