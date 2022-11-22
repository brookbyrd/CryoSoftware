function colorFeatureExtraction(config)

% Direct user to color image example
cd(strcat(config.newSavePath,filesep,'RGB',filesep,'RGB_stacks'));
mkdir(strcat(config.newSavePath,filesep,'RGB',filesep,'RGB_color_filtered'));

% %% Select the color feature of interest
% disp('Select an image to identify color feature');
% [imFile,impath] = uigetfile('.tif');
% im = imread(strcat(impath,imFile));
% disp('Select feature color');
% BW = roipoly(im);
% 
% rIm = double(im(:,:,1));
% gIm = double(im(:,:,2));
% bIm = double(im(:,:,3));
% sumRGB = double(rIm + gIm + bIm);
% 
% Rv = mean(rIm(BW));
% Gv = mean(gIm(BW));
% Bv = mean(bIm(BW));

%%

rgbImageStack = dir('*twopoint_corrected.tif');

for f = 1:length(rgbImageStack)
    
    % read in the image files
    im = imread(rgbImageStack(f).name);
    
    % Break out image into channels
    rIm = double(im(:,:,1));
    gIm = double(im(:,:,2));
    bIm = double(im(:,:,3));
    sumRGB = double(rIm + gIm + bIm);
    
    % Color of Ventricles
    Rv = 130;
    Gv = 108;
    Bv = 95;
    
    % Color of the brain
    Rb = 153;
    Gb = 135;
    Bb = 129;
    
    % Weighting ratios of ventricles
    cR = Rv/(Rv + Gv + Bv);
    cG = Gv/(Rv + Gv + Bv);
    cB = Bv/(Rv + Gv + Bv);
    
    % Weighted combination the R,G, and B feature filters
    % intensityImage = SUM of Contribution Ratio * (Background RGB channel - image RGB channel)
    intensityIm = 255*(cR*(Rb - rIm)+ cG*(Gb - gIm)+ cB*(Bb - bIm))./sumRGB;
    %intensityIm = 255*(cR*(1 - abs(Rv - rIm)./Rv)+ cG*(1- abs(Gv - gIm)./Gv) + cB*(1 - abs(Bv - bIm)./Bv));%./sumRGB;
    [Gmag,Gdir] = imgradient(intensityIm);
    
    imwrite(uint8(intensityIm),strcat(config.newSavePath,filesep,'RGB',filesep,'RGB_color_filtered',filesep,rgbImageStack(f).name,'_ventricle_filtered.tif'));
    
    
end

end