function[mytform] = register_760()
%%
  if(~exist('imMoving'))
         disp('Select a moving raw image');
        [imMoving,imMovingpath] = uigetfile('.tif');
        
        disp('Select a fixed raw image');     
        [imFixed,imFixedpath] = uigetfile('.tif');
  end
  
  moving = imread(strcat(imMovingpath, imMoving));
 
  fixed = imread(strcat(imFixedpath, imFixed));
    
  %%
  subplot(1,2,1);
  imshow(moving,[0 5000]);  
  subplot(1,2,2);
  imshow(fixed,[0 10000]);
  
  movingMat = mat2gray(moving,[0 5000]);
  fixedMat = mat2gray(fixed,[0 10000]);
  %%
  cpselect(movingMat, fixedMat);
  
%cpselect(mat2gray(movingMat,[0 1000]), mat2gray(fixed));
%%
%Variables have been cdreated in the base workspace.
mytform = fitgeotrans(movingPoints, fixedPoints,'similarity'); 


%%
registeredMat = imwarp((movingMat), mytform,'OutputView',imref2d(size(fixed)));
imshowpair(fixedMat, registeredMat);

end