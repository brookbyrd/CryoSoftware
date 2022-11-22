function nextImageCorrect_exponential(config, channel, filename)

addpath(pwd);
%%
% Function using median Alexandrakis values for attenuation
%%
switch channel
    case "470"
        %filename = "470_integratedStack.tif";
        %mu = 0.0050204;
        %mu = 0.0072;
         mu = 0.0034;
    case "635"
        %filename = "635_integratedStack.tif";
         mu = 0.0008377;
         %mu = 0.0016;
    case "760"
        %filename = "760_processed_LCTF_stacks.tif";
        %mu = 0.00058614;
        mu = 0.0009;
    otherwise
end

%%
mu_x = mu*42;
mu_z = mu*config.sliceThickness;
mu_a = exp(-2*mu_z);

n = 25;
psf = generatePsf(n, mu_x);

% normalize psf
psf = psf./sum(psf(:));

%%
cd(config.newSavePath);
delete(strcat(filename, '_nextImage_exponential_m_',num2str(mu),'.tiff'));

%%
for s = 1:length(imfinfo(filename))-1
    
    %%
    topIm = imread(filename,s);
    bottomIm = imread(filename,s+1);
    
    topIm = medfilt2(double(topIm),[3 3]);
    bottomIm = medfilt2(double(bottomIm),[3 3]);
    
    %convolve the image below
    gausIm = convn(bottomIm,psf,'valid');
     
    % Subtract off the bottom image from the top to get residual (Rs)
    Rs = double(topIm) - mu_a.*double(bottomIm);

    
    imwrite(uint16(Rs), strcat(filename, '_nextImage_exponential_m_',num2str(mu),'.tiff'),'Writemode', 'append');
    %clear topIm; clear bottomIm; clear Rs; clear  gausIm; clear scale; clear offset;
    s
end
end

%%