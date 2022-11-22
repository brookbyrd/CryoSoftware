function processAllStacks(config,flatfield, channel, currentPath, rep, beginSlice, countFlag, ME)
%% PROCESSALLSTACKS Core function for basic processing of Bioslice Macrotome Data
%
% PROCESSALLSTACKS(CONFIG, FLATFIELD, CHANNEL, CURRENTPATH, REP,
% BEGINSLICE, COUNTFLAG, ME) processes the full stack of images from a


% *Author: Brook Byrd $ $Date: 6-11-19 $*
% *Copyright: Dartmouth College*

%% SET UP DEPENDENT FILEPATHS an EXPERIMENTAL PARAMETERS FOR PROCESSING
% Current path is the folder with all processing functions
cd(currentPath);
addpath(currentPath);
%Color weights to use for converting light to D65 illuminant
load('color_weights_D65.mat');
% Transform for merging two channels, default is the identity transform
load('tform_identity.mat');

% Use the second repitition unless specified
if ~exist('rep')
    rep = 2;
end

% Use the first slice unless specified
if ~exist('beginSlice')
    beginSlice = config.startingSlice
end

% Allow the user to report signal to background images
if ~exist('countFlag')
    countFlag = 0;
end

% Set the channelRange to the range specified for the given channel
channelRange = config.range(find(config.channels == channel),:);

% Allow the user to preselect ME value for converting processed image intensities to fit within a uint-16 range by using a Base 10 scale
if ~exist('ME')
    switch channel
        
        % modify the uint16 conversion to suite each channel
        case "760"
            ME = 1;
        case "470"
            ME = 1;
        case "635"
            ME = 1;
        case "WL"
            ME = 1;
        case "530"
            ME = 2;
        case "532"
            ME = 1;
        otherwise
    end
    
end
%
% load in the flatfield array matrix
field = char(strcat('Channel_',channel,'_mat'));
flatfield_array = getfield(flatfield,field);
%% CREATE FOLDER STRUCTURE FOR PROCESSED DATA
disp(strcat('Processing Channel ',char(channel)));

% Make the necessary folders for each channel
cd(config.newSavePath);
mkdir(char(channel));
cd(char(channel));
mkdir(char(strcat(channel,'_processed_LCTF_stacks'))); % Hyperspectral image file directory
mkdir(char(strcat(channel,'_nextImage_corrected_stacks'))); % Next-image corrected file directory

% NO LCTF exists in NIR-I channel so don't need to make an integrated stack
% for this channel
if~(strcmp(channel, "760"))
    mkdir(char(strcat(channel,'_integratedStack'))); % Integrated stack file directors
else
    % Load in the 760 transform file if this is channel 760
    load(config.transformFile);
end

% Flag to allow the images to be recorded in units of absolute pixel counts
% at each voxel.
if (countFlag)
    mkdir(strcat(config.newSavePath, filesep, channel, filesep, channel, '_processed_Counts_stacks'));
end

%% LOAD FILES IN EACH RAW DATA FOLDER
% z is the cumulative slice counter
for fNum = 1:length(config.filePath)
    
    % Examine the first raw data folder directory  provided by the user
    cd(char(config.filePath(fNum)));
    date = dir('*20*'); %dates when data was acquired
    for day = 1:length(date)
        dateNum(day) = datenum(date(day,1).name,'yyyymmmdd');
    end
    % sort dates correctly by year, month, and date.
    [dateNum,I] = sort(dateNum);
    
    % Loop through all dates in the sorted order
    for d = 1:length(date)
        
        % change directory to date folder of a given file path
        cd(strcat(config.filePath(fNum), filesep, date(I(d)).name));
        
        if(isdir(channel))
            % Select the correct channel folder director
            cd(strcat(config.filePath(fNum), filesep, date(I(d)).name,filesep,channel));
            
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
            slices = [];
            
            if exist('fileInfo','var')
                %Convert file info into strings that can be called on
                fileInfo = string(fileInfo);
                
                % pull out the unique repitions
                repitions = unique(fileInfo(:,8));
                
                % pull out the unique slices
                slices = unique(fileInfo(:,2));
            end
            
            currentFolder = pwd;
            
            %% MAIN PROCESSING LOOP
            % Iterate over all slices in the folder
            for slice = 1:(length(slices))
                cd(currentFolder);
                tic
                
                % pull out the file indices for this slice
                sliceIndexArray = find(fileInfo(:,2) == slices(slice));
                
                % pull out the historical slice number
                sliceLabel = slices(slice);
                N = cellfun(@(x)sscanf(x,'Slice%d'),{sliceLabel});
                
                % Assign the startingSlice to the first slice number N, if
                % it hasn't already been assigned in the configExperiment.m
                % program.
                if ~isfield(config, 'startingSlice')
                    config.startingSlice = N;
                end
                
                %Figure out the slice number and appropriate affine transforms
                z =  N - config.startingSlice + 1; % This counts from 1 to the number of slices
                z_microns = config.startDepth + config.sliceThickness*z;
                
                % Linear function for affine transform to correctly scale
                % and shift field as it moves away from the camera.
                tform_z = affine2d([1.3e-6*z_microns + 0.99,0,0;0,1.3e-6*z_microns + 0.99,0; -0.0015*z_microns + 7.3,-0.0014*z_microns + 6.8,1]);
                
                if(N < beginSlice)
                    % skip over the processing because the slice number is
                    % lower than the assigned starting slice.
                    disp(strcat('Skipping Slice:', {' '}, num2str(N)));
                else
                    
                    % select the flatfield file per slice
                    flatfield_temp = (squeeze(flatfield_array(:,:,z)));
                    
                    % pull out the unique LCTF waves for this particular slice
                    wave = unique(fileInfo(sliceIndexArray,3));
                    lambda = cellfun(@(x)sscanf(x,'LCTF%dnm'),wave);
                    
                    % select the range to integrate over
                    bandpass = [channelRange(1):channelRange(end)];
                    
                    % pull out which wavelengths to integrate
                    available_waves = intersect(lambda(:),bandpass(:));
                    available_indexes = find(ismember(lambda(:),bandpass(:)));
                    
                    % pull out the dark image from the last LCTF
                    % wavelength used
                    OffIndex = intersect(intersect(find(strcmp(fileInfo(:,7),'OFF')), find(strcmp(fileInfo(:,8),strcat('rep',num2str(rep),'.tif')))),find(strcmp(fileInfo(:,2), slices(slice))));
                    
                    if(~isempty(OffIndex))
                        
                        % Read in off image using the correct OffIndex
                        try
                            imOff = imread(files(OffIndex(1)).name);
                        catch
                            warning(strcat('Slice', {' '}, sliceLabel, {' '}, 'does not exist'));
                        end
                        
                        % Gather the exposure time in seconds of the
                        % specific instance
                        D = cellfun(@(x)sscanf(x,'exp%dms'),(fileInfo(OffIndex,6)));
                        exptime  = D./1e3; %Exposure time expressed in cps assuming it is the same for laser off and on
                        
                        % Iterate over all wavelengths
                        for wavecount = 1:length(wave)
                            
                            % pull out the specific image that is the
                            % correct wavelength for a given LCTF wavelength
                            OnIndex = intersect(intersect(intersect(find(strcmp(fileInfo(:,7),'ON')), find(strcmp(fileInfo(:,8),strcat('rep',num2str(rep),'.tif')))),find(strcmp(fileInfo(:,3), wave(wavecount)))),find(strcmp(fileInfo(:,2), slices(slice))));
                            
                            % skip over wavelength if it does not exist
                            if(~isempty(OnIndex))
                                imOn = loadtiff(files(OnIndex(1)).name);
                            end
                            
                            %Perform basic image corrections for signal background and exposure time.
                            temp_image =  ((double(imOn - imOff))/double(exptime(1)));
                            
                            
                            % Perform radial distortion correction
                            temp_image = undistortImage(temp_image,config.params,'OutputView','same');
                            
                            % Register all channels together
                            temp_image = imwarp((temp_image), mytform, 'OutputView', imref2d(size(temp_image)));%,
                            
                            % Adjust for scale and transform based on stage
                            temp_image = imwarp((temp_image), tform_z, 'OutputView', imref2d(size(temp_image)));%,
                            
                            
                            % crop the image down to save on memory
                            temp_image_rect = imcrop(temp_image, config.channel_rect);
                            
                            % Apply flatfield correction and then scale up by a
                            % Multiply the ratio of image intensity divided by the flatfield by a factor of x1000
                            try
                                temp_image_rect = 1000*(double(temp_image_rect)./double(flatfield_temp));
                            catch
                            end
                            
                            % CIE D65 Scaling for RGB reconstruction to
                            % make the WL hyperspectral image look like it
                            % came from a D65 illuminant
                            if(strcmp(cellstr(channel),'WL'))
                                scalingFactor = weightedLambda(find(weightedLambda(:,1) ==  sscanf(wave(wavecount),'LCTF%dnm')),2);
                            else
                                scalingFactor = 1;
                            end
                            
                            % Scale the matrix by the appropriate scalingFactor
                            stack(:,:,wavecount) =  (temp_image_rect).*double(scalingFactor);
                            
                        end
                        
                        % Apply a slice-specific transform to the entire stack if
                        % config.sliceTransform exists, which is a
                        % transform necessary to adjust for physical shifts
                        % in the OCT block location.
                        if isfield(config,'sliceTransform')
                            offset = affine2d;
                            offset.T = config.sliceTransform(:,:,z);
                            stack = imwarp(stack, offset, 'OutputView', imref2d(size(temp_image_rect)));%,
                        end
                        
                        % Calculate raw stack counts based on the exposure
                        % time and flatfield correction.
                        flatfield_temp = (squeeze(flatfield_array(:,:,z)));
                        stackCounts = (double(stack)/1000).*double(flatfield_temp)*exptime(1) + 100;
                        maxCounts = max(stackCounts,[],3);
                        
                        % Count map is a logical map which shows when
                        % voxels achieve a signal greater than the lower
                        % detection limit of 300 raw counts per voxel.
                        countMap = zeros(size(maxCounts));
                        countMap(find(maxCounts > 300)) = 1;
                        
                        % CONVERT THE PROCESSED STACK ON A BASE 10
                        % using ME as the power coefficient.
                        cd(config.newSavePath);
                        if(ME < 0)
                            fullStack = convertTo16bit(stack.*10^abs(ME), 65535);
                        else
                            fullStack = convertTo16bit(stack, 65535*10^ME);
                        end
                        
                        %% SAVE HYPERSPECTRAL IMAGE FILES
                        % Code to save the 3D uint16 hyperspectral images in the
                        % LCTF directory for each slice.
                        saveLCTFPath =  char(strcat(config.newSavePath, filesep, channel, filesep, channel, '_processed_LCTF_stacks'));
                        saveLCTFFile = char(strcat('Channel', channel,'_', sliceLabel, '_','exp',num2str(D(1)),'ms_rep',num2str(rep),'_','cpsE',num2str(ME),'_cropped_stack.tif'));
                        cd(saveLCTFPath);
                        options.overwrite = true;
                        options.message = 1;
                        saveastiff(uint16(fullStack),char(saveLCTFFile),options);
                        cd(currentFolder);
                        
                        % write the count stack if countFlag = 1
                        if(countFlag)
                            saveCountPath =  char(strcat(config.newSavePath, filesep, channel, filesep, channel, '_processed_Counts_stacks'));
                            saveCountFile = char(strcat('Channel', channel,'_', sliceLabel, '_','exp',num2str(D(1)),'ms_rep',num2str(rep),'_','cpsE',num2str(ME),'_Count_stack.tif'));
                            cd(saveCountPath);
                            saveastiff(uint16(stackCounts),char(saveCountFile),options);
                            cd(currentFolder);
                        end
                        
                        
                        %%  INTEGRATE HYPERSPECTRAL IMAGES OVER SIMULATED BANDPASS
                        if(length(available_indexes) > 1)
                            
                            %Integrate over necessary wavelengths
                            % We divide the range by 10 nm to keep dx at 10 nm instead
                            % of 1 nm.
                            integralStack = trapz(available_waves, double(stack(:,:,available_indexes)),3)./(10);
                            
                            % convert integrated staack to 16 bit accounting for +/- ME
                            if(ME < 0)
                                dataIntegrated = convertTo16bit(integralStack.*10^abs(ME), 65535);
                            else
                                dataIntegrated = convertTo16bit(integralStack, 65535*10^ME);
                            end
                            
                            % Code to save the integrated image
                            saveIntegratedPath =  char(strcat(config.newSavePath, filesep, channel, filesep, channel, '_integratedStack')); options.message = 0;
                            saveIntegratedFile =  char(strcat(sliceLabel,'_integrated_exp',num2str(D),'ms_rep',num2str(rep),'_','cpsE',num2str(ME),'_', num2str(channelRange(1)),'_', num2str(channelRange(2)),'nm.tif')); options.message = 0;
                            cd(saveIntegratedPath);
                            options.overwrite = true;
                            options.message = 0;
                            saveastiff(uint16(dataIntegrated) ,char(saveIntegratedFile),options);
                            cd(currentFolder);
                        end
                        
                        
                        
                        %% NEXT-IMAGE CORRECTION
                        % Function using calculated muscle attenuation
                        % (post-correction for cylindrical volume)
                        
                        if exist('topIm','var')
                            switch channel
                                case "470"
                                    % Using Microsphere phantom data
                                    bottomIm = uint16(dataIntegrated);
                                    mu = 0.0037; %based off microspheres
                                    s1 = 0.36;  %based off microspheres
                                    s2 = 91.15;  %based off microspheres
                                    
                                case "530"
                                    bottomIm = uint16(dataIntegrated);
                                    
                                    mu = 0.0032; %based off microspheres
                                    s1 = 0.46;  %based off microspheres
                                    s2 = 91.42; %based off microspheres
                                case "532"
                                    bottomIm = uint16(dataIntegrated);
                                    
                                    mu = 0.0032; %based off microspheres
                                    s1 = 0.46;  %based off microspheres
                                    s2 = 91.42;  %based off microspheres
                                    
                                case "635"
                                    % Measured mut (post-correction for cylindrical volume)
                                    bottomIm = uint16(dataIntegrated);
                                    
                                    mu = 0.000847; % Measured in vivo
                                    s1 = 0.37;  %based off microspheres
                                    s2 = 98.89;  %based off microspheres
                                    
                                case "760"
                                    % Measured mut (post-correction for cylindrical volume)
                                    bottomIm = uint16(fullStack);
                                    mu = 0.000682;  % Measured in vivo
                                    s1 = 1.506;  %Based off
                                    s2 = 4047;  %based off microspheres
                                    
                                case "WL"
                                    % Using Alexandrakis, 2005 value for muscle
                                    % attentuation at 520 nm
                                    bottomIm = uint16(dataIntegrated);
                                    mu = 0.0028;  %based off microspheres (530)
                                    s1 = 0.36;  %based off microspheres (530)
                                    s2 = 91.15;  %based off microspheres (530)
                                otherwise
                                    bottomIm = [];
                            end
                            
                            % generate voxel-wise attenuation coefficient
                            % Kernel radius
                            n = 15;
                            
                            % generate the Gaussian 2D spread
                            h = generateGausPsf(config, n, s1, s2);
                            %imshow(h,[]);
                            
                            % Attenuate from source and fluorophore
                            mu_a = exp(-2*mu*config.sliceThickness);
                            
                            % filter the top and bottom images
                            % topIm was already filtered
                            %bottomIm = medfilt2(double(bottomIm),[3 3]);
                            
                            %%convolve the image below
                            gausIm = conv2(bottomIm,h/(sum(h(:))),'same');
                            
                            % subtract off top image and rescale to maintain
                            % global intensities
                            % rescaleFactor = double(1/(1- mu_a));
                            rescaleFactor = 1;
                            Rs = rescaleFactor*(double(topIm) - mu_a.*gausIm);
                            
                            
                            % Save the file in the NextImage subfolder
                            saveNextImagePath =  char(strcat(config.newSavePath, filesep, channel, filesep, channel, '_nextImage_corrected_stacks'));
                            saveNextImageFile = char(strcat('Channel', channel,'_', sliceLabel, '_','exp',num2str(D(1)),'ms_rep',num2str(rep),'_','cpsE',num2str(ME),'_', num2str(channelRange(1)),'_', num2str(channelRange(2)),'nm_','mut_',num2str(mu,2),'_nextImage_corrected.tif'));
                            cd(saveNextImagePath);
                            options.overwrite = true;
                            options.message = 1;
                            saveastiff(uint16(imboxfilt(medfilt2(Rs,[3 3]),[3 3])),char(saveNextImageFile),options);
                            
                            cd(currentFolder);
                            
                            % reset the topIm for the next image
                            topIm = bottomIm;
                            
                            
                        else % If no top image exists, it is the first image so just save the image as the top image for next-image corrections on the second image.
                            switch channel
                                case "470"
                                    topIm = uint16(dataIntegrated);
                                case "530"
                                    topIm = uint16(dataIntegrated);
                                case "532"
                                    topIm = uint16(dataIntegrated);
                                case "635"
                                    topIm = uint16(dataIntegrated);
                                case "760"
                                    topIm = uint16(fullStack);
                                case "WL"
                                    topIm = uint16(dataIntegrated);
                                otherwise
                                    topIm = [];
                            end
                            
                        end
                        
                        %% RGB remixing
                        if(strcmp(channel,'WL'))
                            cd(config.newSavePath);
                            
                            % If no maximum signal is predetermined in this
                            % channel, use a pre-set maximum for the
                            % spectre2rgb.m function
                            if(~isfield(config,'maximum'))
                                disp('no max set'); % used for setting global maxima
                                config.maximum = [77968.1503850000];
                                [rgbStack] = spectrum2rgb(lambda, fullStack,config.maximum);
                                save(char(strcat(config.newSavePath,filesep, config.studyName,'_config.mat')), 'config', '-v7.3');
                            else
                                rgbStack = spectrum2rgb(lambda, fullStack, config.maximum);
                            end
                            
                            % Use pre-calibrated color card to perform
                            % black and white balancing, by applying the
                            % predetermined offset and slope in each R, G,
                            % and B channel.
                            disp('Using pre-calibrated color card');
                            m = config.colorConsts(1,[1:3]);
                            b = config.colorConsts(2,[1:3]);
                            
                            rgbStack(:,:,1) = m(1)*rgbStack(:,:,1) + b(1);
                            rgbStack(:,:,2) = m(2)*rgbStack(:,:,2) + b(2);
                            rgbStack(:,:,3) = m(3)*rgbStack(:,:,3) + b(3);
                            imshow(rgbStack);
                            RGBfilename = char(strcat(config.newSavePath,filesep, 'RGB', filesep, 'RGB_stacks', filesep, sliceLabel,'_rgb_twopoint_corrected.tif'));
                            
                            imwrite(rgbStack, RGBfilename);
                            cd(currentFolder);
                        end
                        
                    else
                        disp(strcat('Could not find Slice:', {' '}, num2str(N)));
                    end
                    %%
                    
                    toc
                    disp(strcat('Processed ',{' '}, channel,{' '}, 'Slice', num2str(N)));
                    disp(strcat('Slice ',{' '}, num2str(z),{' '}, 'out of', {' '}, num2str(config.numberOfSlices(1))));
                    
                end
            end
            
            % Clear slices and file info before moving onto the next data
            % folder.
            clear('files')
            clear('fileInfo');
            clear('slices');
            clear('raw_files');
            
        else
            % Safety check to make sure data folder exists.
            disp(strcat('Folder: ',channel,' does not exist'));
        end
        
    end
end

