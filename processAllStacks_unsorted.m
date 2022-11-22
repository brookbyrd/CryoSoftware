function processAllStacks(config,flatfield, channel, channelRange, currentPath, rep, beginSlice, sbrFlag)
% Function for complete processing of macrotome LCTF images
% Created 6-11-19
% Modified 2-11-20
cd(currentPath);
addpath(currentPath);
load('color_weights_D65.mat');
load('tform_identity.mat');
ME = 0;

% Must set this every time
expTime = 1500;
lambda = [510:10:720];
%%
% Load in weighting factors for illumination source

% load in all configurations
newSavePath = config.newSavePath;
wavelengths = config.wavelengths;
numberOfSlices = config.numberOfSlices;
channel_rect = config.channel_rect;
OCT_rect = config.OCT_rect;
studyName = config.studyName;
channels = config.channels;
range = config.range;
params = config.params;
filePath = config.filePath;


if ~exist('rep')
    rep = 2;
end

if ~exist('beginSlice')
    beginSlice = config.startingSlice
end


if ~exist('sbrFlag')
    sbrFlag = 0;
end

%% Create folders for saving
disp(strcat('Processing Channel ',char(channel)));

cd(config.newSavePath);
mkdir(char(channel));
cd(char(channel));
mkdir(char(strcat(channel,'_processed_LCTF_stacks')));

if~(strcmp(channel, "760"))  % NO LCTF exists in IR channel
    mkdir(char(strcat(channel,'_integratedStack')));
end

if (sbrFlag)
    mkdir(char(strcat(channel,'_processed_SBR_stacks')));
end


% load in the necessary channel registration transforms, IR Camera has
% to be flipped and transfered
cd(currentPath);
switch channel
    
    % modify the uint16 conversion to suite each channel
    case "760"
        
        % we need to load in a different transform
        if isfield(config, 'transformFile')
            disp(strcat('loading: ',config.transformFile));
            load(config.transformFile);
        else
            disp(strcat('loading: ','tform_760.mat'));
            load(char(strcat('tform_760.mat')));
        end
        
        ME = 0;
    case "470"
        ME = 0;
        
    case "635"
        ME = 0;
        
    case "WL"
        ME = 1;
        
    case "530"
        ME = 0;
        
    case "532"
        ME = 1;
    otherwise
end

%%
% load in the flatfield array matrix
field = char(strcat('Channel_',channel,'_mat'));
flatfield_array = getfield(flatfield,field);

%%
% Direct to file path
%%z = 1; % z is the cumulative slice counter
for fNum = 1
    
    cd(char(filePath(fNum)));
    
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
    
    cd(currentPath);
    
    %%
    if exist('fileInfo','var')
        %Convert file info into strings that can be called on
        fileInfo = string(fileInfo);
        
        % pull out the unique repitions
        repitions = unique(fileInfo(:,8));
        
        % pull out the slices
        slices = unique(fileInfo(:,2));
        
        exposures = unique(fileInfo(:,6));
        
        for e = 1:length(exposures)
            expTimes(e) = cellfun(@(x)sscanf(x,'exp%dms'),{exposures(e)});
        end
        expCorrect = exposures(find(expTimes == expTime));
    end
    currentFolder = pwd;
    %%
    % Iterate over all slices in the folder
    for slice = 1:(length(slices))
        cd(currentFolder);
        tic
        
        % pull out the file indices for this slice
        sliceIndexArray = find(fileInfo(:,2) == slices(slice));
        % pull out the historical slice number
        sliceLabel = slices(slice)
        N = cellfun(@(x)sscanf(x,'Slice%d'),{sliceLabel});
        
        if ~isfield(config, 'startingSlice')
            config.startingSlice = N;
        end
        
        z =  N - config.startingSlice + 1
        
        if(N < beginSlice)
            % skip over the processing
            disp(strcat('Skipping Slice:', {' '}, num2str(N)));
        else
            
            try
                % select the flatfield file per slice
                flatfield_temp = (squeeze(flatfield_array(:,:,z)));
            catch
            end
            % pull out the unique LCTF waves for this particular slice
            wave = unique(fileInfo(sliceIndexArray,3));
            %lambda = cellfun(@(x)sscanf(x,'LCTF%dnm'),wave);
            
            
            % select the range to integrate over
            bandpass = [channelRange(1):channelRange(end)];
            
            % pull out which wavelengths to integrate
            available_waves = intersect(lambda(:),bandpass(:));
            available_indexes = find(ismember(lambda(:),bandpass(:)));
            
            
            % pull out the dark image from the last LCTF
            % wavelength used
            OffIndex = intersect(intersect(find(strcmp(fileInfo(:,7),'OFF')), find(strcmp(fileInfo(:,6),expCorrect))),find(strcmp(fileInfo(:,2), slices(slice))));
            
            
            if(~isempty(OffIndex))
                cd(char(filePath(fNum)));
                
                try
                    imOff = imread(files(OffIndex(1)).name);
                catch
                    warning(strcat('Slice', {' '}, sliceLabel, {' '}, 'does not exist'));
                end
                
                % Gather the exposure time in seconds of the
                % specific instance
                D = cellfun(@(x)sscanf(x,'exp%dms'),(fileInfo(OffIndex,6)));
                exptime  = D./1e6; %Exposure time expressed in cps assuming it is the same for laser off and on
                
                
                % Iterate over all wavelengths
                w = 1;
                for wavecount = 1:length(wave)
                    
                    
                    % pull out the slice within the time set that has the
                    % correct wavelength for the given LCTF wavelength
                    OnIndex = intersect(intersect(intersect(find(strcmp(fileInfo(:,7),'ON')), find(strcmp(fileInfo(:,6),expCorrect))),find(strcmp(fileInfo(:,3), wave(wavecount)))),find(strcmp(fileInfo(:,2), slices(slice))));
                    
                    % skip over wavelength if it does not exist
                    if(~isempty(OnIndex))
                        try
                            imOn = loadtiff(files(OnIndex(1)).name);
                        catch
                        end
                        
                        
                        
                        % try
                        %perform basic image corrections
                        temp_image =  ((double(imOn - imOff))/double(exptime(1)));
                        
                        %                         catch % in case imOff or imOn do not exist;
                        %                             warning(strcat('Slice', {' '}, slices(slice), {' '}, 'does not exist'));
                        %                             temp_image = ones(size(flatfield_temp));
                        %                         end
                        
                        % perform radial distortion correction
                        temp_image = undistortImage(temp_image,params,'OutputView','same');
                        
                        %register all channels together
                        temp_image = imwarp((temp_image), mytform, 'OutputView', imref2d(size(temp_image)));%, 'OutputView',imref2d(size(temp_image)));
                        
                        % crop the image down to save on memory
                        temp_image_rect = imcrop(temp_image, channel_rect);
                        
                        % apply flatfield correction
                        temp_image_rect = (double(temp_image_rect)./double(flatfield_temp));
                        
                        
                        % CIE D65 Scaling for RGB reconstruction
                        if(strcmp(cellstr(channel),'WL'))
                            scalingFactor = weightedLambda(find(weightedLambda(:,1) ==  sscanf(wave(wavecount),'LCTF%dnm')),2);
                        else
                            scalingFactor = 1.0;
                        end
                        
                        % save the slice data as a matrix
                        stack(:,:,w) =  (temp_image_rect).*double(scalingFactor);
                        
                        % write the SBR stack if necessary
                        if(sbrFlag)
                            sbrStack(:,:,w) = imcrop(double(imOn)./double(imOff), channel_rect);
                        end
                        
                        w = w + 1;
                    end
                end
                
                %%
                cd(newSavePath);
                % Hyperspectral stack saving
                if(ME < 0)
                    fullStack = convertTo16bit(stack.*10^abs(ME), 65535);
                else
                    fullStack = convertTo16bit(stack, 65535*10^ME);
                end
                
                saveLCTFPath =  char(strcat(newSavePath, filesep, channel, filesep, channel, '_processed_LCTF_stacks'));
                %saveLCTFFile = char(strcat('Channel', channel, '_Slice_', num2str(z), '_','exp',num2str(D(1)),'ms_rep',num2str(rep),'_','cpsE',num2str(ME),'_cropped_stack.tif'));
                saveLCTFFile = char(strcat('Channel', channel,'_', sliceLabel, '_','exp',num2str(D(1)),'ms_rep',num2str(rep),'_','cpsE',num2str(ME),'_cropped_stack.tif'));
                
                cd(saveLCTFPath);
                
                options.overwrite = true;
                options.message = 1;
                saveastiff(uint16(fullStack),char(saveLCTFFile),options);
                cd(currentFolder);
                
                % write the SBR stack if necessary
                if(sbrFlag)
                    saveSBRPath =  char(strcat(newSavePath, filesep, channel, filesep, channel, '_processed_SBR_stacks'));
                    saveSBRFile = char(strcat('Channel', channel,'_', sliceLabel, '_','exp',num2str(D(1)),'ms_rep',num2str(rep),'_','cpsE',num2str(ME),'_SBR_stack.tif'));
                    cd(saveSBRPath);
                    saveastiff(uint16(10.*sbrStack),char(saveSBRFile),options);
                    cd(currentFolder);
                end
                %% Integrate the bandpass
                if(length(available_indexes) > 1)
                    
                    
                    %Integrate over necessary wavelengths
                    % We divide the range by 10 nm to keep dx at 10 nm instead
                    % of 1 nm.
                    %
                    
                    integralStack = trapz(available_waves, stack(:,:,available_indexes),3)./(10);
                    
                    % convert integrate to 16 bit
                    if(ME < 0)
                        dataIntegrated = convertTo16bit(integralStack.*10^abs(ME), 65535);
                    else
                        dataIntegrated = convertTo16bit(integralStack, 65535*10^ME);
                    end
                    
                    % save the integrated image
                    saveIntegratedPath =  char(strcat(newSavePath, filesep, channel, filesep, channel, '_integratedStack')); options.message = 0;
                    %saveIntegratedFile =  char(strcat('Slice_', num2str(z),'_integrated_exp',num2str(D),'ms_rep',num2str(rep),'_','cpsE',num2str(ME),'_', num2str(channelRange(1)),'_', num2str(channelRange(2)),'nm.tif')); options.message = 0;
                    saveIntegratedFile =  char(strcat(sliceLabel,'_integrated_exp',num2str(D),'ms_rep',num2str(rep),'_','cpsE',num2str(ME),'_', num2str(channelRange(1)),'_', num2str(channelRange(2)),'nm.tif')); options.message = 0;
                    
                    cd(saveIntegratedPath);
                    options.overwrite = true;
                    options.message = 0;
                    saveastiff(uint16(dataIntegrated) ,char(saveIntegratedFile),options);
                    cd(currentFolder);
                end
                
                                    %% Next image correction
                    % Function using calculated muscle attenuation
                    % (post-correction for cylindrical volume)
                    
                    if exist('topIm')
                        switch channel
                            case "470"
                                % Using Microsphere phantom data
                                bottomIm = uint16(dataIntegrated);
                                mu = 0.0037;
                                s1 = 0.36;
                                s2 = 91.15;
                                
                            case "530"
                                bottomIm = uint16(dataIntegrated);
                                
                                mu = 0.0032;
                                s1 = 0.46;
                                s2 = 91.42;
                            case "532"
                                bottomIm = uint16(dataIntegrated);
                                
                                mu = 0.0032;
                                s1 = 0.46;
                                s2 = 91.42;
                            case "635"
                                % Measured mut (post-correction for cylindrical volume)
                                bottomIm = uint16(dataIntegrated);
                                
                                mu = 0.0028;
                                s1 = 0.37;
                                s2 = 98.89;
                                
                            case "760"
                                % Measured mut (post-correction for cylindrical volume)
                                bottomIm = uint16(fullStack);
                                
                                mu = 0.0012;
                                s1 = 0.37;
                                s2 = 98.89;
                            case "WL"
                                % Using Alexandrakis, 2005 value for muscle
                                % attentuation at 520 nm
                                bottomIm = uint16(dataIntegrated);
                                mu = 0.0037;
                                s1 = 0.36;
                                s2 = 91.15;
                            otherwise
                                bottomIm = [];
                        end
                        
                        % generate voxel-wise attenuation coefficient
                        %%
                        n = 15;
                        
                        % generate the Gaussian 2D spread
                        h = generateGausPsf(config, n,  s1, s2);
                        imshow(h,[]);
                        
                        % Attenuate from source and fluorophore
                        mu_a = exp(-2*mu*config.sliceThickness);
                        
                        % filter the top and bottom images
                        %topIm = medfilt2(double(topIm),[3 3]); topIm was
                        % already filtered
                        bottomIm = medfilt2(double(bottomIm),[3 3]);
                        
                        %%convolve the image below
                        gausIm = conv2(bottomIm,h/(sum(h(:))),'same');
                        
                        % subtract off top image
                        rescaleFactor = double(max(bottomIm(:))./(max(gausIm(:)))); % rescaling factor to ensure no loss from convolution.                        
                         Rs = (double(topIm) - mu_a.*rescaleFactor*gausIm);
                        
                        % Save the file in the NextImage subfolder
                        saveNextImagePath =  char(strcat(newSavePath, filesep, channel, filesep, channel, '_nextImage_corrected_stacks'));
                        saveNextImageFile = char(strcat('Channel', channel,'_', sliceLabel, '_','exp',num2str(D(1)),'ms_rep',num2str(rep),'_','cpsE',num2str(ME),'mut_',num2str(mu,2),'_nextImage_corrected.tif'));
                        
                        cd(saveNextImagePath);
                        options.overwrite = true;
                        options.message = 1;
                        saveastiff(uint16(Rs),char(saveNextImageFile),options);
                        
                        cd(currentFolder);
                        
                        % reset the topIm for the next image
                        topIm = bottomIm;
                        
                    else
                        switch channel
                            case "470"
                                topIm = uint16(dataIntegrated);
                                topIm = medfilt2(double(topIm),[3 3]);
                            case "530"
                                topIm = uint16(dataIntegrated);
                                topIm = medfilt2(double(topIm),[3 3]);
                            case "532"
                                topIm = uint16(dataIntegrated);
                                topIm = medfilt2(double(topIm),[3 3]);
                            case "635"
                                topIm = uint16(dataIntegrated);
                                topIm = medfilt2(double(topIm),[3 3]);
                            case "760"
                                topIm = uint16(fullStack);
                                topIm = medfilt2(double(topIm),[3 3]);
                            case "WL"
                                topIm = uint16(dataIntegrated);
                                topIm = medfilt2(double(topIm),[3 3]);
                            otherwise
                                topIm = [];
                        end
                        
                    end
                
                %% RGB remixing
                if(strcmp(channel,'WL'))
                    
                    if(~isfield(config,'maximum'))
                        disp('no max set'); % used for setting global maxima
                        [rgbStack,config.maximum] = spectrum2rgb(lambda, fullStack);
                    else
                        rgbStack = spectrum2rgb(lambda, fullStack, config.maximum);
                    end
                    
                    
                    if isfield(config, 'colorConsts')
                        disp('Using pre-calibrated color card');
                        m = config.colorConsts(1,[1:3]);
                        b = config.colorConsts(2,[1:3]);
                        
                        rgbStack(:,:,1) = m(1)*rgbStack(:,:,1) + b(1);
                        rgbStack(:,:,2) = m(2)*rgbStack(:,:,2) + b(2);
                        rgbStack(:,:,3) = m(3)*rgbStack(:,:,3) + b(3);
                        imshow(rgbStack);
                        RGBfilename = char(strcat(newSavePath,filesep, 'RGB', filesep, 'RGB_stacks', filesep, sliceLabel,'_rgb_twopoint_corrected.tif'));
                        
                    else
                        % Perform white balancing
                        if(~exist('Rmax'))
                            if(config.OCT_rect(1) ~= [0])
                                
                                disp('Using white balancing');
                                Rmax = mean(imcrop(rgbStack(:,:,1),OCT_rect),[1 2]);
                                Gmax = mean(imcrop(rgbStack(:,:,2),OCT_rect),[1 2]);
                                Bmax = mean(imcrop(rgbStack(:,:,3),OCT_rect),[1 2]);
                            else
                                Rmax = max(max(rgbStack(:,:,1)));
                                Gmax = max(max(rgbStack(:,:,2)));
                                Bmax = max(max(rgbStack(:,:,3)));
                            end
                        end
                        %
                        rgbStack(:,:,1) =  rgbStack(:,:,1)./Rmax;
                        rgbStack(:,:,2) = rgbStack(:,:,2)./Gmax;
                        rgbStack(:,:,3) = rgbStack(:,:,3)./Bmax;
                        imshow(rgbStack);
                        RGBfilename = char(strcat(newSavePath,filesep, 'RGB', filesep, 'RGB_stacks', filesep, sliceLabel,'_rgb_Rmax_corrected.tif'));
                        
                    end
                    
                    
                    imwrite(rgbStack, RGBfilename);
                    cd(currentFolder);
                end
                
            else
                disp(strcat('Could not find Slice:', {' '}, num2str(N)));
            end
            %%
            %                 % remove stacks from local memory
            %                 stack = [];
            %                 OffIndex = [];
            %                 fullStack = [];
            %                 dataIntegrated = [];
            %                 integralStack = [];
            %                 rgbStack = [];
            %                 temp_image_stack = [];
            %                 flatfield_temp = [];
            %
            
            toc
            % move onto the next slice
            %disp(strcat('Processed ',{' '}, channel,{' '}, 'Slice ', {' '}, num2str(z)));
            
            %disp(strcat('Processed ',{' '}, channel,{' '}, sliceLabel));
            disp(strcat('Processed ',{' '}, channel,{' '}, 'Slice', num2str(N)));
            disp(strcat('Slice ',{' '}, num2str(z),{' '}, 'out of', {' '}, num2str(config.numberOfSlices(1))));
            
        end
        
        
    end
    
    
    clear('files')
    clear('fileInfo');
    clear('slices');
    clear('raw_files');
    %catch
    %end
end




end











