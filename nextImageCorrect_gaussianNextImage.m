function nextImageCorrect_gaussianNextImage(config, channel, filename)

addpath(pwd);
%%
% Function using median Alexandrakis values for attenuation
%%
switch channel
    case "470"
         mu = 0.0037; %Measured by MS
         s1 = 0.36;
         s2 = 91.15;
    case "530"    
         mu = 0.0032; %Measured by MS
         s1 = 0.46;
         s2 = 91.42;
    case "635"    
         %mu = 0.0028; % MS
         %mu = 0.0012; % Alexandrakis adipose
         mu = 0.000847; % Measured in vivo
      
         s1 = 0.37;
         s2 = 98.89;
    case "760"
         %mu = 0.0009;  % Alexandrakis adipose
         mu = 0.000682;  % Measured in vivo
        
         s1 = 0.37;
         s2 = 98.89;
    otherwise
end

%%
n = 15;

% generate the Gaussian 2D spread
h = generateGausPsf(config, n, s1, s2);

% Attenuate from source and fluorophore
mu_a = exp(-2*mu*config.sliceThickness);


%%
cd(config.newSavePath);
delete(strcat(filename, '_nextImage_gaussian_m_',num2str(10000*mu),'.tiff'));

%%
for s = 1:length(imfinfo(filename))-1
    
    %%
    topIm = imread(filename,s);
    bottomIm = imread(filename,s+1);
    
    
    topIm = medfilt2(double(topIm),[3 3]);
    bottomIm = medfilt2(double(bottomIm),[3 3]);
   
    %%
    %%convolve the image below
    gausIm = conv2(bottomIm,h/(sum(h(:))),'same');
    
    % subtract off top image
    Rs = double(topIm) - mu_a.*gausIm;
    
    imwrite(uint16(Rs), strcat(filename, '_nextImage_gaussian_m_',num2str(mu),'s1_',num2str(s1),'_s2_',num2str(s2),'.tiff'),'Writemode', 'append');
    %clear topIm; clear bottomIm; clear Rs; clear  gausIm; clear scale; clear offset;
    s
end
end

%%