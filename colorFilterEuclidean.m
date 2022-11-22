function colorFilterEuclidean(config, flag)

% Brain weighted = 0.5R 0.3G 0.2B
%Lung weighted = 0.508R 0.24G 0.2B

if(flag)
    % Self select ROI
    cd(strcat(config.newSavePath, filesep,'RGB',filesep,'RGB_stacks'));

    disp('Select an image for color ROI');
    [imFile,impath] = uigetfile('.tif');
    im = imread(strcat(impath,imFile));
    [BW] = roipoly(mat2gray(im));
    close all;

    %
    R = im(:,:,1);
    G = im(:,:,2);
    B = im(:,:,3);

    % Find mean value of R, G, and B values for selection
    cR = mean(double(R(BW)));
    cG = mean(double(G(BW)));
    cB = mean(double(B(BW)));
else
    % Ventricles
    segmentName = "ventricles";
    cR = 130;
    cG = 108;
    cB = 95;
end

if ~exist('segmentName')
    segmentName = inputdlg({'Name this segment ROI '},'Segment selection:',[ 1 35]);
    segmentName = string(segmentName);
end

savePath = config.newSavePath;
mkdir(char(strcat(savePath,'/RGB/RGB_filtered_',segmentName)));
cd(char(strcat(savePath,'/RGB/RGB_stacks')));

%%
files = dir('*.tif');

for i = 1:length(files)
    filenames(i) = string(files(i).name);
end

filenum = cellfun(@(x)sscanf(x,'Slice%d_rgb_corrected.tif'),filenames,'UniformOutput',false);
[~,Sidx] = sort(cell2mat(filenum));
sortedFilenames = filenames(Sidx);

for f = 1:length(sortedFilenames)
    
    rgbStack = imread(sortedFilenames(f));
    
    % Compute Euclidean Distance between RGB pixels and ideal color
    deltaR2 = (double(rgbStack(:,:,1)) - cR).^2;
    deltaG2 = (double(rgbStack(:,:,2)) - cG).^2;
    deltaB2 = (double(rgbStack(:,:,3)) - cB).^2;
    
    % Equation pulled from: https://www.tdx.cat/bitstream/handle/10803/6189/07Jvl07de11.pdf;jsessionid=11F0A58F28C195CE50BD3760C98C7321?sequence=7
    deltaC = sqrt(deltaR2 + deltaG2 + deltaB2);
    
    % Subtract distance from 255 to invert the distance image
    EucInvertedStack = 100*(255-deltaC);
    
    imwrite(uint16(EucInvertedStack) , char(strcat(savePath,'/RGB/RGB_filtered_',segmentName,filesep,sortedFilenames(f),'_',segmentName,'_filtered.tif')));
    
    clear rgbStack; clear deltaC; clear EucInvertedStack;
    
    slice = f;
end


end

