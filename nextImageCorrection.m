function nextImageCorrection(config,channelAgent, bodyRegion)
% Function for doing nextImage corrections following spectral fitting
addpath('H:\Shared drives\MacroImage_software\Final Analysis Package (All Platforms)');

% add cases to load
switch channelAgent
    
    case "PEG550"
        channel = "470";
        
        
    case "PEG5000"
        channel = "470";
        
        
    case "FITC"
        channel = "470";
        
    case "GFP"
        channel = "470";
        
    case "YFP"
        channel = "470";
        
    case "AF488"
        channel = "470";
        
    case "Bodipy"
        channel = "470";
        
    case "tdTomato"
        channel = "530";
        
    case "OF550"
        channel = "530";
        
    case "Rhd550"
        channel = "530";
        
    case "Rhd1k"
        channel = "530";
        
    case "Rhd5k"
        channel = "530";
        
    case "PPIX"
        channel = "635";
        
    case "AF647"
        channel = "635";
        
    case "AF680"
        channel = "635";
        
    case "IR680"
        channel = "635";
        
    case "OF650"
        channel = "635";
        
    otherwise
end

%%
mkdir(char(strcat(config.newSavePath, filesep, channel, filesep, channel,'_', channelAgent, '_unmixed_nextImage_corrected_stacks')));
cd(strcat(config.newSavePath,filesep,channel,filesep,channel,'_SpecFit_',bodyRegion));
currentFolder = pwd;

files = dir(strcat('*',channelAgent,'_LLS_unmixed_stack.tif'));

for i = 1:length(files)
    filenames(i) = string(files(i).name);
end

filenum = cellfun(@(x)sscanf(x,strcat('Channel',channel,'_Slice%d_exp1500s_cropped_stack.tif')),filenames,'UniformOutput',false);
[~,Sidx] = sort(cell2mat(filenum));
sortedFilenames = filenames(Sidx);

%%

for f = 1:length(sortedFilenames)
    
    cd(strcat(config.newSavePath,filesep,channel,filesep,channel,'_SpecFit_',bodyRegion));
    
    dataUnmixed = imread(sortedFilenames(f));
    
    if exist('topIm')
        
        switch channel
            case "470"
                % Using Microsphere phantom data
                bottomIm = uint16(dataUnmixed);
                
                mu = 0.0037; %based off microspheres
                s1 = 0.36;  %based off microspheres
                s2 = 91.15;  %based off microspheres
                
            case "530"
                bottomIm = uint16(dataUnmixed);
                
                mu = 0.0032; %based off microspheres
                s1 = 0.46;  %based off microspheres
                s2 = 91.42; %based off microspheres
            case "532"
                bottomIm = uint16(dataUnmixed);
                
                mu = 0.0032; %based off microspheres
                s1 = 0.46;  %based off microspheres
                s2 = 91.42;  %based off microspheres
                
            case "635"
                % Measured mut (post-correction for cylindrical volume)
                bottomIm = uint16(dataUnmixed);
                
                mu = 0.000847; % Measured in vivo
                s1 = 0.37;  %based off microspheres
                s2 = 98.89;  %based off microspheres
                
            case "760"
                % Measured mut (post-correction for cylindrical volume)
                bottomIm = uint16(fullStack);
                
                mu = 0.000682;  % Measured in vivo
                s1 = 0.37;  %based off microspheres
                s2 = 98.89;  %based off microspheres
                
            case "WL"
                % Using Alexandrakis, 2005 value for muscle
                % attentuation at 520 nm
                bottomIm = uint16(dataUnmixed);
                mu = 0.0037;  %based off microspheres (530)
                s1 = 0.36;  %based off microspheres (530)
                s2 = 91.15;  %based off microspheres (530)
            otherwise
                bottomIm = [];
        end
        
        % generate voxel-wise attenuation coefficient
        %%
        n = 15;
        
        % generate the Gaussian 2D spread
        h = generateGausPsf(config, n, s1, s2);
        
        
        % Attenuate from source and fluorophore
        mu_a = exp(-2*mu*config.sliceThickness);
        
        % filter the top and bottom images
        %topIm = medfilt2(double(topIm),[3 3]); topIm was
        % already filtered
        bottomIm = medfilt2(double(bottomIm),[3 3]);
        
        %%convolve the image below
        gausIm = conv2(bottomIm,h/(sum(h(:))),'same');
        
        % subtract off top image and rescale to maintain
        % global intensities
        % rescaleFactor = double(1/(1- mu_a));
        rescaleFactor = 1; % We are no longer rescaling
        Rs = rescaleFactor*(double(topIm) - mu_a.*gausIm);
        
        % Save the file in the NextImage subfolder
        saveNextImagePath =  char(strcat(config.newSavePath, filesep, channel, filesep, channel,'_',channelAgent, '_unmixed_nextImage_corrected_stacks'));
        saveNextImageFile = char(strcat('Channel',channel,'_Slice',num2str(f),'_mut_',num2str(mu,2),'_nextImage_corrected.tif'));
        
        cd(saveNextImagePath);
        options.overwrite = true;
        options.message = 1;
        imwrite(uint16(Rs),char(saveNextImageFile));
        
        cd(currentFolder);
        
        % reset the topIm for the next image
        topIm = bottomIm;
        
        
    else
        switch channel
            case "470"
                topIm = uint16(dataUnmixed);
                topIm = medfilt2(double(topIm),[3 3]);
            case "530"
                topIm = uint16(dataUnmixed);
                topIm = medfilt2(double(topIm),[3 3]);
            case "532"
                topIm = uint16(dataUnmixed);
                topIm = medfilt2(double(topIm),[3 3]);
            case "635"
                topIm = uint16(dataUnmixed);
                topIm = medfilt2(double(topIm),[3 3]);
            case "760"
                topIm = uint16(fullStack);
                topIm = medfilt2(double(topIm),[3 3]);
            case "WL"
                topIm = uint16(dataUnmixed);
                topIm = medfilt2(double(topIm),[3 3]);
            otherwise
                topIm = [];
        end
        
    end
    
    
end

end