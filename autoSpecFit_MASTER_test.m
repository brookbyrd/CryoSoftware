function autoSpecFit_MASTER(config, channelAgent, bodyRegion, flag)
%% Function for spectral fitting with a library for each agent
% @author Brook Byrd (10-21-21)

% Inputs: config = experimental configuration
% channelAgent (optional) = agent of interest
% flag (optional) = 1 (plot) or 0 (process all slices)

instance = struct;

fittingFolder = pwd;
addpath(fittingFolder);
addpath('G:\Shared drives\MacroImage_software\Final Analysis Package (All Platforms)');
load('basis.mat','biofit_basis_470_table_M778','tdTomato_530','biofit_basis_AF680_table_DDSI');

if(~exist('flag','var'))
    %flag = 1; % PLOT
    flag = 0; % Paralellyzed run
end
%% GUI for selecting  agents and ROI

app = spectralSelection_app(config);
waitfor(app,'running',0);
close all;

% Collect the selected data
bodyRegion = app.bodyRegion;
channelAgent = app.channelAgent;

% Load in the wavelengths
filePath = config.newSavePath;
wavelengths = app.wavelengths;
lambda = app.lambda;

% Set up the appropriate indexing and integral ranges
available_indexes = app.available_indexes;
available_stack_indexes = app.available_stack_indexes;
available_waves = app.available_waves;

% select the range to integrate over
bandpass = app.bandpass;
integral_range = app.integral_range
integral_waves= app.integral_waves;


% Plot the results of the spectral bases in use.
A = app.A;
A_names = app.A_names;

figure;
plot(lambda(available_indexes),A./max(A),'-o','LineWidth',4);
legend( A_names);
ylim([0 1.2]);
ylabel('Normalized Spectra'); xlabel('Wavelength (nm)');


%% Call in file info and pull out relevant data
mkdir(strcat(config.newSavePath,filesep,channel,filesep,channel,'_SpecFit_',bodyRegion));
mkdir(strcat(config.newSavePath,filesep,channel,filesep,channel,'_SpecFit_',bodyRegion,'_nM'));
cd(strcat(config.newSavePath,filesep,channel,filesep,channel,'_processed_LCTF_stacks'));
files = dir('*E0_cropped_stack.tif');

for i = 1:length(files)
    filenames(i) = string(files(i).name);
end

filenum = cellfun(@(x)sscanf(x,strcat('Channel',channel,'_Slice%d_exp1500s_cropped_stack.tif')),filenames,'UniformOutput',false);
[~,Sidx] = sort(cell2mat(filenum));
sortedFilenames = filenames(Sidx);

temp_var = cellfun(@(x)sscanf(x,strcat('Channel',channel,'_Slice%d_exp%dms_rep%d_cpsE%d_cropped_stack.tif')),sortedFilenames(1),'UniformOutput',false);
temp_var = cell2mat(temp_var);
ME = temp_var(4);


%% Create a region specific ROI
cd(strcat(config.newSavePath));

try
    load(strcat(bodyRegion,'ROI.mat'));
catch
    cd(strcat(config.newSavePath,filesep, channel , filesep,channel,'_processed_LCTF_stacks'));
    slice = imread(sortedFilenames(round(length(sortedFilenames)/2)));
    imshow(slice,[]);
    disp(strcat('Select ',{' '}, bodyRegion,{' '}, 'from image'));
    ROI = drawrectangle();
    bw = createMask(ROI);
    cd(strcat(config.newSavePath));
    save(strcat(bodyRegion,'ROI.mat'),"ROI","bw");
    close all;
end

cd(strcat(filePath,filesep, channel , filesep, channel,'_processed_LCTF_stacks'));

%% PLOTTING CODE
if(flag) % plot the spectral fit of a selected point in a loop
    
    % Set this number to the slice of interest
    %f = round(length(sortedFilenames)/2);
    f = 106;
    running = true;
    
    % Load up the slices used in a stack
    cd(strcat(filePath,filesep, channel , filesep, channel,'_processed_LCTF_stacks'));
    index = 1;
    for i = 1:length(available_stack_indexes)
        temp = imread(sortedFilenames(f),available_stack_indexes(i));
        temp = medfilt2(temp,[3 3]);
        stack(:,:,index) = temp;
        index = index + 1;
    end
    
    %%
    % Create a figure
    figure1 = figure('Color',[1 1 1],'Units','normalized','Position',[0.25,0.25,0.7,0.5]);
    
    while(running)
        subplot1 = subplot(1,2,1,'Parent',figure1);
        imshow(squeeze(stack(:,:,3)),[0 1000]); hold on;
        
        % plot the previous ROI
        try
            plot(x2, y2, 'r.-', 'MarkerSize', 15);
        catch
        end
        
        % Select an ROI
        [BW, x2, y2] = roipoly();
        
        % collect the spectra
        if(length(x2) == 1)
            data = double(squeeze(stack(round(y2),round(x2),:)));
        else
            for j = 1:size(stack,3)
                slice470 = squeeze(stack(:,:,j));
                data(j) = median(slice470(find(BW == 1)));
            end
            data = data';
        end
        
        % Collect indices that are above p = 0.01
        % F_spec_2: uses only top contibuting base + agent spectra if it is
        % statistically shown to be present ( p > 0.01)
        [Fpd, indices, spec_filtered, u_filtered] = F_spec_2(A,  data);

        % Plot the final spectral fit
        subplot2 = subplot(1,2,2,'Parent',figure1);
        
        x_axis = [1:size(A,1)]';
        Legend=cell(length(indices) + 2,1);
        Legend{1}='Data';
        Legend{2} = 'Fit';
        fit = sum(spec_filtered(1,:,:),3);
        
        plot(x_axis,data, '-k',x_axis,fit,':','LineWidth',3);
        hold on
        for j = 1:length(indices)
            plot(x_axis,spec_filtered(1,:,j),'-*','LineWidth',3);
            Legend{j+2}=strcat(A_names(indices(j)));
        end
        legend(Legend);
        xlabel('Wavelength [nm]');
        ylabel('Intensity [AU]');
        
        set(gca,'linew',2,'FontSize',20);
        xticks(available_indexes(1:2:end));
        xticklabels(available_waves(1:2:end));
        hold off;
        clear data;
    end
else
    %%
    %% Begin paralellyzed group
    
    % % Configure paralizable work
    delete(gcp('nocreate'));
    parpool('local',8); %current work station is an 8 core
    poolobj = gcp('nocreate');
    
    % send each slice to a different processor in an 8 core processor
    parfor f = 1:length(sortedFilenames)
        
        % update the progress bar
        %waitbar(double(f-1)./length(sortedFilenames),'Spectral fitting progress:'); drawnow;
        
        %%
        % Create a storage place
        stack = [];
        filename = sortedFilenames(f);
        data = [];
        topIm = [];
        
        %%
        cd(strcat(filePath,filesep, channel , filesep, channel,'_processed_LCTF_stacks'));
        % Load up the slices used in a stack
        index = 1;
        for i = 1:length(available_stack_indexes)
            temp = imread(filename,available_stack_indexes(i));
            temp = medfilt2(temp,[3 3]);
            stack(:,:,index) = temp;
            index = index + 1;
        end
        
        % create a map to save the integrated spectral values
        u_map = zeros(size(temp,1),size(temp,2),size(A,2));
        
        tic
        for r = 1:size(stack,1)
            for c = 1:size(stack,2)
                % Only proceed if the region is within the ROI
                
                if(bw(r,c) == 1)
                    
                    % Collect data spectra
                    data = double(squeeze(stack(r,c,:)));
                    indices = [1:size(A,2)];
                    
                    % Check F-stat, if p>0.01,
                    % then remove the lower contributor
                    % Spec filtered is the weighted base contributions
                    % u_filtered is the weighting scheme
                    [Fpd, indices, spec_filtered, u_filtered] = F_spec_2(A,  data);
                    
                    % save the unmixed signals
                    u_trap = trapz(integral_waves, spec_filtered(1,integral_range,:))./(10);
                    
                    % store the integrated values
                    u_map(r,c,indices) = u_trap;
                end
            end
        end
        toc
        
        
        % Write out the unmixed spectral stacks in integrated RFUs
        for ind = 1:size(u_map,3)
            imwrite(uint16(imboxfilt(medfilt2(u_map(:,:,ind),[4 4]),3)),strcat(filePath,filesep,channel, filesep, channel,'_SpecFit_', bodyRegion, filesep,filename,'_Slice_',num2str(f), '_',num2str(bandpass(1)),'_',num2str(bandpass(end)),'nm_ME',num2str(ME),'_', A_names(ind),'_F_spec_2_unmixed_stack.tif'));
        end
        
        
        %%Provide nM conversion coefficients
        % added cases to convert RFU to nM concentrations
        switch channelAgent
            case "PEG550"
                m = 0.2033;
                b = 0;
            case "PEG5000"
            case "FITC"
                m = 0.0462;
                b = 0.485;
            case "AF488"
                m = 0.0783;
                b = -12.016;
            case "Bodipy"
                m = 0.2741;
                b = -14.41;
            case "AF647"
                m = 0.0231;
                b = 0;
            case "AF680"
            case "IR680"
                m = 0.0291;
                b = 0;
            otherwise
        end
        % Adjust the slope based off of ME scaling
        m = m.*10^ME;
        
        %  Write image in nM concentration units
        
        % Save the file in the NextImage subfolder
        Rs_nM = (u_map(:,:,1)*m + b);
        
        % clean up the image using a gaussian and sharpening filter
        Rs_nM = imboxfilt(medfilt2(Rs_nM,[4 4]),3);
        
        % Write out the nM converted next-image corrected slice
        imwrite(uint16(Rs_nM), strcat(filePath,filesep,channel, filesep, channel,'_SpecFit_', bodyRegion,'_nM', filesep,filename,'_Slice_',num2str(f), '_',num2str(bandpass(1)),'_',num2str(bandpass(end)),'nm_ME',num2str(ME),'_', A_names(1),'_nM_F_spec_2_stack_new.tif'));
        
        % Slice number
        slice = f
        
    end
    
end
% % selection
%         function Test(app, event)
%             disp('test worked');
%             
% 
%         end
end

% Create ValueChangedFcn callback:
function selection(instance,dd,ax)

val = dd.Value

switch dd.Value
        case 'Head'
            instance.bodyRegion = "head"
        case 'Hind Tumor'
            instance.bodyRegion = "hind"
        case 'GI System'
            instance.bodyRegion = "gi"
        case 'Phantom'
           instance.bodyRegion = "phantom"
        otherwise
end
   

updatePlot(instance,dd,ax);

end

function updatePlot(instance,dd,ax)

val = dd.Value

switch dd.Value
        case 'Head'
            instance.bodyRegion = "head"
        case 'Hind Tumor'
            instance.bodyRegion = "hind"
        case 'GI System'
            instance.bodyRegion = "gi"
        case 'Phantom'
           instance.bodyRegion = "phantom"
        otherwise
end
   

end
