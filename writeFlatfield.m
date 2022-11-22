function  writeFlatfield(config)

channels = config.channels;
filePath = config.filePath;
channel_rect = config.channel_rect;
params = config.params;
savePath = config.newSavePath;



% SET UP FILE PATHS
% Setting up function paths

currentPath = pwd;
addpath(currentPath);
addpath(strcat(currentPath,filesep, 'color'));
 

% Direct to file path
for fNum = 1:length(filePath)
    
    cd(char(filePath(fNum)));
    date = dir('*20*');
    for day = 1:length(date)
        dateNum(day) = datenum(date(day,1).name,'yyyymmmdd');
    end
    [dateNum,I] = sort(dateNum);

    % iterate over each passed in channel
    for c = 1:length(channels)
        %try
            z = 1;
            channel = channels(c);
            allSlices = [];
            
            
             switch channel
            % modify the uint16 conversion to suite each channel
            
            case "470"
                centeredWave = 'LCTF510nm';
            case "530"
                centeredWave = 'LCTF570nm';
            case "635"
                centeredWave = 'LCTF660nm';
            case "760"
                centeredWave = 'LCTF0nm';
            case "WL"
                centeredWave = 'LCTF600nm';
            otherwise
             end
        
        
             
             
            disp(strcat('Processing the flatfield of: ', channel));
            
            if(strcmp(channel,"760"))
                load(strcat('tform_',channel,'_M677_v2.mat'));
                disp(strcat('loading: ','tform_',channel,'_M677_v2.mat'));
            else
                load('tform_default.mat');
            end
            
             mkdir(char(strcat(savePath,filesep, channel, filesep, channel, '_raw_flatfield')));
  
            
            % Slice counting index
            index = 1;
            for d = 2:length(date)
                
                
                cd(char(strcat(filePath(fNum),filesep, date(I(d)).name, filesep,'Flatfield', filesep,  channel)))
                
                
                % Work around to make sure all files are in the correct format
                raw_files =  dir('*.tif');
                counter = 1; %
                for f = 1:length(raw_files)
                    %check to make sure the formatting is correct
                    if(length(split(raw_files(f).name,'_')) == 8)
                        files(counter) = raw_files(f); % remove any files that do not fit this format
                        fileInfo(counter,:) = split(files(counter).name,'_');
                        counter = counter + 1;
                    end
                end
                
                %Convert file info into strings that can be called on
                fileInfo = string(fileInfo);
                
                % pull out the unique repitions
                repitions = unique(fileInfo(:,8));
                
                % pull out the slices
                allSlices = [allSlices; unique(fileInfo(:,2))];
                selectedSlices = unique(fileInfo(:,2));
                
                % pull out the unique LCTF waves for this particular slice
                wave = unique(fileInfo(:,3));
                
                %%
                disp(strcat('Loading flatfields of channel: ',channel));
                
                % Iterate over each centered wavelength (selecting the centered for LCTF)
                for wavecount = 1;
                    
                    
                    % length of slices is only 2 every time
                    for slice = 1:length(selectedSlices)
                        
                        
                        % pull out the file indices for this slice
                        sliceIndexArray = find(fileInfo(:,2) == selectedSlices(slice));
                        
                        
                        % pull out the dark image from the last LCTF
                        % wavelength used
                        OffIndex = intersect(intersect(find(strcmp(fileInfo(:,7),'OFF')), find(strcmp(fileInfo(:,8),'rep2.tif'))),find(strcmp(fileInfo(:,2), selectedSlices(slice))));
                        OffLowerIndex = intersect(intersect(find(strcmp(fileInfo(:,7),'OFF')), find(strcmp(fileInfo(:,8),'rep1.tif'))),find(strcmp(fileInfo(:,2), selectedSlices(slice))));
                        
                        
%                         try % just in case the file is corrupted
                            
                            imOff = imread(files(OffIndex(1)).name);
                            
                            % Gather the exposure time in seconds of the
                            % specific instance
                            D = regexp(fileInfo(OffIndex,6),'\d*','Match');
                            exptime  = str2num(char(D{1,1}))./1000; %Exposure time expressed in seconds assuming it is the same for laser off and on
                            
                            
                            
                            
                            % pull out the slice within the time set that has the
                            % correct wavelength for the given LCTF wavelength
                            OnIndex = intersect(intersect(intersect(find(strcmp(fileInfo(:,7),'ON')), find(strcmp(fileInfo(:,8),'rep2.tif'))),find(strcmp(fileInfo(:,3), wave(wavecount)))),find(strcmp(fileInfo(:,2), selectedSlices(slice))));
                            OnLowerIndex = intersect(intersect(intersect(find(strcmp(fileInfo(:,7),'ON')), find(strcmp(fileInfo(:,8),'rep1.tif'))),find(strcmp(fileInfo(:,3), wave(wavecount)))),find(strcmp(fileInfo(:,2), selectedSlices(slice))));
                            
                            try
                                imOn = imread(files(OnIndex(1)).name);
                            catch
                                %warning(strcat((files(OnLowerIndex(1)).name),' is an invalid file'));
                            end
                            
                            
                            
                            %perform basic image corrections
                            temp_image =  ((double(imOn - imOff)./double(exptime)));%./double(flatfield_norm));
                            
                            %
                            % perform distortion correction
                            temp_image = undistortImage(temp_image,params,'OutputView','same');
                            
                            %                     %register all channels together
                            temp_image = imwarp((temp_image), mytform);
                            
                            % apply smoothing function
                            %temp_image = imgaussfilt(temp_image,8);
                            
                            % crop down the flatfield to save on time
                            temp_image_rect = imcrop(temp_image, channel_rect);
                            
                            
                            % crop the image down to save on memory
                            flatfield(:,:,index) = (temp_image_rect);%, channel_rect);
                            index = index + 1;
                            
                            
                            data = uint32(10000*(temp_image_rect));
                            t =  Tiff(char(strcat(savePath,'\',channel,'\', channel,'_raw_flatfield\Channel', channel, '_Slice_', num2str(wave(wavecount)), '_', num2str(z), '_','exp',strcat(D{1}),'ms_cropped_stack.tif')),'w');
                                
                            setTag(t,'Photometric',Tiff.Photometric.MinIsBlack);
                            setTag(t,'Compression',Tiff.Compression.None);
                            setTag(t,'BitsPerSample',32);
                            setTag(t,'SamplesPerPixel',1);
                            setTag(t,'SampleFormat',Tiff.SampleFormat.UInt);
                            setTag(t,'ImageLength',size(temp_image_rect,1));
                            setTag(t,'ImageWidth',size(temp_image_rect,2));
                            setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
                            
                            write(t,data);
                            close(t);
                            
                            
                        z = z+1;
                        slice
                        end
                    end
                end
            end
            
            
          
%         catch
%         end
        
        disp(strcat('Completed channel: ', channel));
        
    end
end
    
    
