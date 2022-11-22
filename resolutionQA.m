
%% ANALYZE SPATIAL CALIBRATION
function resolutionQA(config)
% Process written to quickly run through spatial calibration measures
% Created 7-20-20
addpath(config.currentPath);
cd(config.currentPath);
addpath('color');
load('color_weights_D65.mat');
load('tform_default.mat');
ME = 1;
load('distortion_param.mat');
mkdir(strcat(config.newSavePath, filesep, 'RGB',filesep,'RGB_stacks_spatialCalib',filesep,'ResolutionQA'));
 
%% COLLECT resolution target stack
try
    cd(config.newSavePath);
    cd(char(strcat('RGB', filesep,'RGB_stacks_spatialCalib')));
    files = dir('Slice*');
    for i = 1:length(files)
        filenames(i) = string(files(i).name);
    end
    
    filenum = cellfun(@(x)sscanf(x,'Slice_%d_rgb.tif'),filenames,'UniformOutput',false);
    [~,Sidx] = sort(cell2mat(filenum));
    sortedFilenames = filenames(Sidx);
    
    disp('Select the top image');
    [imFile,impath] = uigetfile('.tif');
    startingInd = find(strcmp(sortedFilenames,imFile));
    
    clear im;
    ind = 1;
    for s = startingInd:length(sortedFilenames)
        temp_image = imread(sortedFilenames(s));
        im(:,:,ind) = rgb2gray(temp_image);
        ind = ind + 1;
    end
    
    %% Register images
    [optimizer, metric] = imregconfig('multimodal');
    optimizer.InitialRadius = 0.009;
    optimizer.Epsilon = 1.5e-4;
    optimizer.GrowthFactor = 1.01;
    optimizer.MaximumIterations = 1000;
    
    tform500.T = [1 0 0; 0 1 0; 0 0 1];
    tform2000 = imregtform(im(:,:,2), im(:,:,1), 'affine', optimizer,metric);
    tform3500 = imregtform(im(:,:,3), im(:,:,1), 'affine', optimizer,metric);
    
    %register all images together
    registered500 = im(:,:,1);
    registered2000 = imwarp(im(:,:,2), tform2000, 'OutputView', imref2d(size(im(:,:,2))));%, 'OutputView',imref2d(size(temp_image)));
    registered3500 = imwarp(im(:,:,3), tform3500, 'OutputView', imref2d(size(im(:,:,2))));%, 'OutputView',imref2d(size(temp_image)));
    
    %%
    subplot(2,2,1);
    imshowpair(im(:,:,1), im(:,:,2)); title('2000 registration');
    subplot(2,2,2);
    imshowpair(im(:,:,1),  registered2000);
    subplot(2,2,3);
    imshowpair(im(:,:,1), im(:,:,3)); title('3500 registration');
    subplot(2,2,4);
    imshowpair(im(:,:,1), registered3500);
    
    %% Interpolate Registration
    interpTransforms = [config.startDepth:config.sliceThickness:config.sliceThickness*config.numberOfSlices(1) - config.sliceThickness + config.startDepth];
    [Xq,Yq,Zq] = meshgrid(1:3,1:3,config.startDepth:config.sliceThickness:config.sliceThickness*config.numberOfSlices(1) - config.sliceThickness + config.startDepth);
    
    X = cat(3, Xq(:,:,1), Xq(:,:,1), Xq(:,:,1));
    Y = cat(3, Yq(:,:,1), Yq(:,:,1),  Yq(:,:,1));
    Z = cat(3, 5000*ones(3,3),20000*ones(3,3), 35000*ones(3,3));
    
    V =  cat(3,tform500.T,tform2000.T,tform3500.T);
    
    config.interpTransformsScaled = interp3(X,Y,Z,V,Xq,Yq,Zq,'makima');
    save(char(strcat(config.newSavePath,filesep, config.studyName,'_config.mat')), 'config', '-v7.3');
    disp('Set to go transform interpolation');
    
    %% Profile Analysis
    
    %%
    disp('Crop Image left to right, top to bottom');
    
    for count = 1:9
        
        figure(4);
        [square, rect(count,:)] =  imcrop(rgb2gray(im));
        count
    end
    
    
    
    %%
    clear element_minima; clear limit;
    disp('Draw profile left to right, top to bottom');
    for count = 1:9
        
        
        figure(1);
        imshow(registered500,[]);
        clear cx, clear cy; clear c;
        [cx,cy,c] = improfile();
        x_point(1,count) = cx(1);    x_point(2,count) = cx(end);
        y_point(1,count) = cy(1);    y_point(2,count) = cy(end);
        
        
        count
        
        figure(2);
        subplot(3,3,count);
        plot(c); hold on;
        
        x = 1:length(c);
        TF = islocalmin(c,'MinSeparation',3,'MinProminence',1);
        TF_space = islocalmin(c,'MinSeparation',26);
        plot(x,c,x(TF),c(TF),'b*'); hold on;
        plot(x,c,x(TF_space),c(TF_space),'r*');
        spaces = find(TF_space == 1);
        minima = find(TF == 1);
        hold on;
        
        for element = 1:min(length(spaces),6)
            element_minima500(element,count) = length(intersect(find(TF == 1),[spaces(element)-15:spaces(element)+30]));
        end
        
        limitingElement500 = (min(find(element_minima500(:,count)<3)))
        if(isempty(limitingElement500))
            title(strcat('All resolutions met in Group 1'));
        else
            title(strcat('Resolution Limit: Group 1 Element',{' '}, num2str(limitingElement500)));
        end
        
    end
    saveas(gcf,strcat(config.newSavePath, filesep, 'RGB',filesep,'RGB_stacks_spatialCalib',filesep,'ResolutionQA', filesep,'500_curves.png'));

    %%
    registered500Copy = registered500;
    for count = 1:9
        registered500Copy = insertShape( registered500Copy ,'Line',[x_point(1,count) y_point(1,count) x_point(2,count) y_point(2,count)],'LineWidth',2,'Color','blue');
        
        figure(3);
        subplot(3,3,count);
        square =  imcrop(rgb2gray(registered500Copy),rect(count,:));
        imshow(square,[]);
        limitingElement500 = (min(find(element_minima500(:,count)<3)))
        if(isempty(limitingElement500))
            title(strcat('All resolutions met in Group 1'));
            group1ElementLimit(count,1) = 0;
        else
            title(strcat('Resolution Limit: Group 1 Element',{' '}, num2str(limitingElement500)));
            group1ElementLimit(count,1) = limitingElement500;
        end
        hold on;
        
    end
     saveas(gcf,strcat(config.newSavePath, filesep, 'RGB',filesep,'RGB_stacks_spatialCalib',filesep,'ResolutionQA', filesep,'500_ROI_squares.png'));
 
    figure(5);
    imshow(registered500Copy,[]); hold on; title('5000 ROI selections');
    saveas(gcf,strcat(config.newSavePath, filesep, 'RGB',filesep,'RGB_stacks_spatialCalib',filesep,'ResolutionQA', filesep,'500_ROI.png'));
 
    %% Compute Resolution in 2nd and 3rd image 
    
      %%
    clear element_minima; clear limit;
    for count = 1:9
        
       
        [cx,cy,c] = improfile(registered2000,[x_point(1,count),x_point(2,count)],[y_point(1,count),y_point(2,count)]);
        
        count
        
        figure(4);
        subplot(3,3,count);
        plot(c); hold on;
        
        x = 1:length(c);
        TF = islocalmin(c,'MinSeparation',3,'MinProminence',1);
        TF_space = islocalmin(c,'MinSeparation',26);
        plot(x,c,x(TF),c(TF),'b*'); hold on;
        plot(x,c,x(TF_space),c(TF_space),'r*');
        spaces = find(TF_space == 1);
        minima = find(TF == 1);
        hold on;
        
        for element = 1:min(length(spaces),6)
            element_minima2000(element,count) = length(intersect(find(TF == 1),[spaces(element)-15:spaces(element)+30]));
        end
        
        limitingElement2000 = (min(find(element_minima2000(:,count)<3)))
        if(isempty(limitingElement2000))
            title(strcat('All resolutions met in Group 1'));
        else
            title(strcat('Group 1 Element',{' '}, num2str(limitingElement2000)));
        end
        
    end
    saveas(gcf,strcat(config.newSavePath, filesep, 'RGB',filesep,'RGB_stacks_spatialCalib',filesep,'ResolutionQA', filesep,'2000_curves.png'));

     registered2000Copy = registered2000;
    for count = 1:9
         registered2000Copy = insertShape(  registered2000Copy ,'Line',[x_point(1,count) y_point(1,count) x_point(2,count) y_point(2,count)],'LineWidth',2,'Color','blue');
        
        figure(5);
        subplot(3,3,count);
        square =  imcrop( rgb2gray(registered2000Copy),rect(count,:));
        imshow(square,[]);
        limitingElement2000 = (min(find(element_minima2000(:,count)<3)))
        if(isempty(limitingElement2000))
            title(strcat('All resolutions met in Group 1'));
            group1ElementLimit(count,2) = 0;
        else
            title(strcat('Group 1 Element',{' '}, num2str(limitingElement2000)));
            group1ElementLimit(count,2) = limitingElement2000;
        end
        hold on;
        
    end
    saveas(gcf,strcat(config.newSavePath, filesep, 'RGB',filesep,'RGB_stacks_spatialCalib',filesep,'ResolutionQA', filesep,'2000_ROI_squares.png'));
 
    figure(6);
    imshow(registered2000Copy,[]); hold on; title('2000 ROI selections');
    saveas(gcf,strcat(config.newSavePath, filesep, 'RGB',filesep,'RGB_stacks_spatialCalib',filesep,'ResolutionQA', filesep,'2000_ROI.png'));
 
    %%
    
       %%
    clear element_minima; clear limit;
    for count = 1:9
        
       
        [cx,cy,c] = improfile(registered3500,[x_point(1,count),x_point(2,count)],[y_point(1,count),y_point(2,count)]);
        
        count
        
        figure(7);
        subplot(3,3,count);
        plot(c); hold on;
        
        x = 1:length(c);
        TF = islocalmin(c,'MinSeparation',3,'MinProminence',1);
        TF_space = islocalmin(c,'MinSeparation',26);
        plot(x,c,x(TF),c(TF),'b*'); hold on;
        plot(x,c,x(TF_space),c(TF_space),'r*');
        spaces = find(TF_space == 1);
        minima = find(TF == 1);
        hold on;
        
        for element = 1:min(length(spaces),6)
            element_minima3500(element,count) = length(intersect(find(TF == 1),[spaces(element)-15:spaces(element)+30]));
        end
        
        limitingElement3500 = (min(find(element_minima3500(:,count)<3)))
        if(isempty(limitingElement3500))
            title(strcat('All resolutions met in Group 1'));
            group1ElementLimit(count,3) = 0;
        else
            title(strcat('Group 1 Element',{' '}, num2str(limitingElement3500)));
            group1ElementLimit(count,3) = limitingElement3500;
        end
        
    end
    saveas(gcf,strcat(config.newSavePath, filesep, 'RGB',filesep,'RGB_stacks_spatialCalib',filesep,'ResolutionQA', filesep,'3500_curves.png'));
 
    
     registered3500Copy = registered3500;
    for count = 1:9
         registered3500Copy = insertShape(  registered3500Copy ,'Line',[x_point(1,count) y_point(1,count) x_point(2,count) y_point(2,count)],'LineWidth',2,'Color','blue');
        
        figure(8);
        subplot(3,3,count);
        square =  imcrop( rgb2gray(registered3500Copy),rect(count,:));
        imshow(square,[]);
        limitingElement3500 = (min(find(element_minima3500(:,count)<3)))
        if(isempty(limitingElement3500))
            title(strcat('All resolutions met in Group 1'));
        else
            title(strcat('Group 1 Element',{' '}, num2str(limitingElement3500)));
        end
        hold on;
        
    end
    saveas(gcf,strcat(config.newSavePath, filesep, 'RGB',filesep,'RGB_stacks_spatialCalib',filesep,'ResolutionQA', filesep,'3500_ROI_squares.png'));
  
    figure(9);
    imshow(registered3500Copy,[]); hold on; title('3500 ROI selections'); 
    saveas(gcf,strcat(config.newSavePath, filesep, 'RGB',filesep,'RGB_stacks_spatialCalib',filesep,'ResolutionQA', filesep,'3500_ROI.png'));
    
    save(strcat(config.newSavePath, filesep, 'RGB',filesep,'RGB_stacks_spatialCalib',filesep,'ResolutionQA', filesep,'group1ElementLimit.mat','group1ElementLimit'));
    disp('QA complete');
catch
    disp('Run process_spatialCalib first');
end

