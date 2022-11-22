% METHOD for registering together the same image as it goes up and down the
% stage

% Written by BKB 7-9-21

info = imfinfo('WL_integratedStack.tif');

for i = 1:length(info)
    stack_res(:,:,i) = imread('WL_integratedStack.tif',i);
end


%%

% Register all images to the fixed image
fixed = stack_res(:,:,1);

[optimizer, metric] = imregconfig('monomodal');
optimizer.MinimumStepLength = 5e-4;
optimizer.MaximumIterations = 30;

tform_start = affine2d([1,0,0;0,1,0; 0,0,1]);

% use the SURF method to detect like features to automatically match
for i = 2:length(info)
    
    % register the next slice
    moving = stack_res(:,:,i);
    
    % using the SURF method
    boxImage = moving;
    sceneImage= fixed; % Scene is fixed
    
    % gather the points and features of a SURF registration
    boxPoints = detectSURFFeatures(boxImage);
    scenePoints = detectSURFFeatures(sceneImage);
    [boxFeatures, boxPoints] = extractFeatures(boxImage, boxPoints,'FeatureSize',128,'Method','KAZE','Upright',true);
    [sceneFeatures, scenePoints] = extractFeatures(sceneImage, scenePoints,'FeatureSize',128,'Method','KAZE','Upright',true);
    
    % Collect the matching metric in cupMetric, 0 represents a perfect match
    [boxPairs, boxMetric] = matchFeatures(boxFeatures, sceneFeatures,'MaxRatio',1,'MatchThreshold',90);
    
    % Organize the pairs into a table with one column representing the match
    % metric
    totalPairs = [single(boxPairs),boxMetric];
    % Sort the rows by the match metric
    totalPairs = sortrows(totalPairs,3);
    
    % Pick out the best 50 matchings indices
    n = 200;
    best50 = totalPairs(1:n,1:2);
    matchedBoxPoints = boxPoints(best50(:, 1));
    matchedScenePoints = scenePoints(best50(:, 2));
    
    figure;
    showMatchedFeatures(boxImage, sceneImage, matchedBoxPoints, ...
        matchedScenePoints, 'montage');
    title('Putatively Matched Points (Including Outliers)');
       
    
    % estimate the transform (can be scaling and offset)
    [tform, inlierBoxPoints, inlierScenePoints] = ...
        estimateGeometricTransform(matchedBoxPoints, matchedScenePoints, 'similarity');
    
    %tform = imregtform(moving,fixed,'similarity', optimizer, metric);
    
    % store the entire stack of transforms
    transform_stack(:,:,i) = tform.T;
    
    % reregister the image
    movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
    stack_res_reregistered(:,:,i) = movingRegistered;
    imshowpair(fixed, movingRegistered);
    
    % Write the reregistered image
    %imwrite(uint16(movingRegistered),'WL_integratedStack_reregistered_SURF.tiff','WriteMode','append');
    
    i
    
end


%% Plot results
transform_stack(:,:,1) = [1,0,0;0,1,0; 0,0,1];
z = 5000:500:5000 + 500*50;

scatter3(39*transform_stack(3,1,:),39*transform_stack(3,2,:),z,squeeze(single(1500*(transform_stack(1,1,:)-0.99))),squeeze(single(1500*(transform_stack(1,1,:)-0.99))),'filled');
set(gca, 'Zdir', 'reverse');
xlabel('X-direction (um)');
ylabel('Y-direction (um)');
zlabel('Z-direction (um)');

ylim([-1024*39, 1024*39]);
xlim([-1024*39, 1024*39]);
pbaspect([1 1 1]);


%%
z_range = 5000:1:35000;
xtrans = squeeze(transform_stack(3,1,:));
ytrans = squeeze(transform_stack(3,2,:));
xscale = squeeze(transform_stack(1,1,:));
yscale = squeeze(transform_stack(2,2,:));
xshear = squeeze(transform_stack(2,1,:));
yshear = squeeze(transform_stack(1,2,:));

figure;
subplot(2,1,1);
mdl1 = fitlm(z,xtrans); 
m1 = mdl1.Coefficients.Estimate(2);
b1 = mdl1.Coefficients.Estimate(1);
mdl2 = fitlm(z,ytrans); 
m2 = mdl2.Coefficients.Estimate(2);
b2 = mdl2.Coefficients.Estimate(1);

scatter(z,xtrans,'s','LineWidth',3); hold on; 
scatter(z,ytrans,'o','LineWidth',3); hold on;
legend({'X transform','Y transform'});
plot(z_range, z_range*m1 + b1,'r'); hold on;
plot(z_range, z_range*m2 + b2,'r');
title('Translation');
ylabel('pixels');
xlabel('Distance from top (um)');
title(strcat('X transform: y =', num2str(m1,2), 'x + ',num2str(b1,2),  ', Y transform y =', num2str(m2,2), 'x + ',num2str(b2,2)));


subplot(2,1,2);
mdl3 = fitlm(z,xscale); 
m3 = mdl3.Coefficients.Estimate(2);
b3 = mdl3.Coefficients.Estimate(1);
scatter(z,xscale,'md','LineWidth',3); hold on; 
legend({'Scale Factor'});
plot(z_range,z_range*m3 + b3,'r');
title(strcat('Scale: y =', num2str(m3,2), 'x + ',num2str(b3,2)));
xlabel('Distance from top (um)');
ylabel('Scale');
