function nrrdReformater(config, filename_nrrd, filename_tiff)
addpath(pwd);
cd(config.newSavePath);

temp = split(filename_tiff,filesep);
filename_save = string(temp(end));
[img] = nrrdread(filename_nrrd);
 
%%
for z = 1:length(imfinfo(filename_tiff))
    oldVolume(:,:,z,:) = imread(filename_tiff,z);
end

rotVolume = fliplr(rot90(rot90(rot90(oldVolume))));

img.pixelData = rotVolume ;
%%
nrrdwrite(strcat(filename_save,'.nrrd'), img);

end       