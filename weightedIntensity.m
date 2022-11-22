function weightedIntensity(savePath, channel, lambda, weights, z_range)
%filePath_WL_LCTF = 'S:\lazaret\Members\Byrd\Data_Analysis\Macrotome\Data Analysis\BrookWork\CurrentDataAnalysis\M630_processed\WL_processed_LCTF_stacks';


mkdir(char(strcat(savePath,'/',channel,'/',channel,'_weightedStack')));
cd(char(strcat(savePath,'/',channel,'/',channel,'_processed_LCTF_stacks')));

files = dir('*.tif');

for i = 1:length(files)
    filenames(i) = string(files(i).name);
end

filenum = cellfun(@(x)sscanf(x,'Channel470_Slice_%d_exp500s_cropped_stack.tif'),filenames,'UniformOutput',false);
[~,Sidx] = sort(cell2mat(filenum));
sortedFilenames = filenames(Sidx);

weightedStack = [];
available_waves = intersect(lambda(:),weights(:,1));
available_indexes = find(ismember(lambda(:),weights(:,1)));


for s = z_range(1):1:z_range(end)
  
        filename = sortedFilenames(s)
        weightedStack = zeros(size(imread(filename)));
        
        % Load up the slices used in a stack
        index = 1;
        for i = [available_indexes(1):available_indexes(end)]
            % weight of wavelength we are looking at
            weight = weights(find(available_waves(i) == weights(:,1)),2);
            tempStack = imread(filename,i);
            tempStack = double(tempStack).*double(weight);
           
            % keep summing up the weighted image
            weightedStack = weightedStack + tempStack;
            index = index + 1;
            
        end
        
            
    
                          data = uint32(( weightedStack));
                            t =  Tiff(char(strcat(strcat(savePath,'/',channel,'/',channel,'_weightedStack/', 'Slice_', num2str(s),'_weighted.tif'))),'w');
                                
                            setTag(t,'Photometric',Tiff.Photometric.MinIsBlack);
                            setTag(t,'Compression',Tiff.Compression.None);
                            setTag(t,'BitsPerSample',32);
                            setTag(t,'SamplesPerPixel',1);
                            setTag(t,'SampleFormat',Tiff.SampleFormat.UInt);
                            setTag(t,'ImageLength',size(weightedStack,1));
                            setTag(t,'ImageWidth',size(weightedStack,2));
                            setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
                            
                            write(t,data);
                            close(t);
                            
                            
    %   imwrite(rgbsubStack, char(strcat(savePath,'/RGB/RGB_stacks_prime/', 'Slice_', num2str(s),'_rgb_corrected.tif')));
    %imwrite(uint16(weightedStack), char(strcat(savePath,'/',channel,'/',channel,'_weightedStack/', 'Slice_', num2str(s),'_weighted.tif')));
    
    clear weightedStack;
    
    slice = z_range(s)
        

end

    
end

