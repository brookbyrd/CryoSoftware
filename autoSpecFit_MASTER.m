function autoSpecFit_MASTER(config, channelAgent, bodyRegion, flag, ME)
%% Function for spectral fitting with a library for each agent
% @author Brook Byrd (10-21-21)

% Inputs: config = experimental configuration
% channelAgent (optional) = agent of interest
% bodyRegion = region of interest (head,gi,hind, phantom)
% flag (optional) = 1 (plot) or 0 (process all slices)
% ME (optional) = 0 (default) (ME you want to process from the LCTFs)

fittingFolder = pwd;
addpath(fittingFolder);
load('basis.mat');
cd ..
addpath(pwd);

if(~exist('flag','var'))
    flag = 1; % PLOT
    %flag = 0; % Paralellyzed run
end
%% GUI for selecting  agents and ROI

if(exist('channelAgent') && exist('bodyRegion'))
    
    % if channel agent and bodyRegion already exist, surpass the GUI
    disp(strcat('Experiment name:',{' '},config.studyName));
    disp(strcat('Processing agent:',{' '},channelAgent,{' '}, 'using the',{' '},bodyRegion,{' '},'spectral bases'));
else
    fig = uifigure('Name','Spectral Fitting Selection','Position',[400 200 300 150]);
    txt = uilabel(fig,'Position',[30 120 250 20],'Text','Select experiment type for spectral fitting:');
    
    txtAgent = uilabel(fig,'Position',[30 90 100 20],'Text','Agent:');
    txtROI = uilabel(fig,'Position',[30 50 100 20],'Text','ROI:');
    
    
    if (~exist('channelAgent'))
        cc = uidropdown(fig,'Position',[70 90 100 20],'Items',{'PEG550','PEG5000','FITC','GFP','YFP','AF488','Bodipy','tdTomato','OF550','Rhd550','Rhd1k','Rhd5k','PPIX','AF647','AF680','IR680','OF650'},...
            'Value','IR680');
    else
        cc = uidropdown(fig,'Position',[70 90 100 20],'Items',{'PEG550','PEG5000','FITC','GFP','YFP','AF488','Bodipy','tdTomato','OF550','Rhd550','Rhd1k','Rhd5k','PPIX','AF647','AF680','IR680','OF650'},...
            'Value',channelAgent);
    end
    
    dd = uidropdown(fig,'Position',[70 50 100 20],'Items',{'Head','Hind Tumor','GI System','Phantom','GFP Head'},...
        'Value',"Head");
    
    % Accept button closes everything down
    acceptButton = uibutton(fig,'state','Position',[70 10 100 20],'Text','Accept set-up');
    waitfor(acceptButton,'Value');

    % Collect user inputs
    channelAgent = string(cc.Value)
    
    switch dd.Value
        case 'Head'
            bodyRegion = "head"
        case 'Hind Tumor'
            bodyRegion = "hind"
        case 'GI System'
            bodyRegion = "gi"
        case 'Phantom'
            bodyRegion = "phantom"
        case 'GFP Head'
            bodyRegion = "gfp"
        otherwise
    end
    
    delete(fig);
end

%%
% Load in the wavelengths
wavelengths = config.wavelengths;
filePath = config.newSavePath;

clear autofl;clear basis;

% add cases to load
switch channelAgent
    
    case "PEG550"
        channel = "470";
        lambda = config.wavelengths.lambda_470;
        table = biofit_basis_470_table_M778;
        basis(:,1) = (biofit_basis_470_table_M778{:,2});
        
    case "PEG5000"
        channel = "470";
        lambda = config.wavelengths.lambda_470;
        table = biofit_basis_470_table_M778;
        basis(:,1) = (biofit_basis_470_table_M778{:,3});
        
        
    case "FITC"
        channel = "470";
        lambda = config.wavelengths.lambda_470;
        table = biofit_basis_470_table_M778;
        basis(:,1) = (biofit_basis_470_table_M778{:,4});
        
    case "GFP"
        channel = "470";
        lambda = config.wavelengths.lambda_470;
        table = biofit_basis_470_table_M778;
        basis(:,1) = (biofit_basis_470_table_M778{:,6});
        
    case "YFP"
        channel = "470";
        lambda = config.wavelengths.lambda_470;
        table = biofit_basis_470_table_M778;
        basis(:,1) = (biofit_basis_470_table_M778{:,7});
        
    case "AF488"
        channel = "470";
        lambda = config.wavelengths.lambda_470;
        table = biofit_basis_470_table_M778;
        basis(:,1) = (biofit_basis_470_table_M778{:,8});
        
    case "Bodipy"
        channel = "470";
        lambda = config.wavelengths.lambda_470;
        table = biofit_basis_470_table_M778;
        basis(:,1) = (biofit_basis_470_table_M778{:,5});
        
    case "tdTomato"
        channel = "530";
        lambda = config.wavelengths.lambda_530;
        table = tdTomato_530;
        basis(:,1) = (tdTomato_530{:,2});
        
    case "OF550"
        channel = "530";
        lambda = config.wavelengths.lambda_530;
        table = tdTomato_530;
        basis(:,1) = (tdTomato_530{1:14,7});
        
    case "Rhd550"
        channel = "530";
        lambda = config.wavelengths.lambda_530;
        table = pegRhd_530_5nm;
        basis(:,1) = (pegRhd_530_5nm{:,2});
     
    case "Rhd1k"
        channel = "530";
        lambda = config.wavelengths.lambda_530;
        table = pegRhd_530_5nm;
        basis(:,1) = (pegRhd_530_5nm{:,3});
    
    case "Rhd5k"
        channel = "530";
        lambda = config.wavelengths.lambda_530;
        table = pegRhd_530_5nm;
        basis(:,1) = (pegRhd_530_5nm{:,4});
        
    case "PPIX"
        channel = "635";
        lambda = config.wavelengths.lambda_635;
        table = test_phantom635;
        basis(:,1) = (test_phantom635{:,5});
        
    case "AF647"
        channel = "635";
        lambda = config.wavelengths.lambda_635;
        table = biofit_basis_IRDye680_table_M688;
        basis(:,1) = (biofit_basis_IRDye680_table_M688{:,16});
        
    case "AF680"
        channel = "635";
        lambda = config.wavelengths.lambda_635;
        table = biofit_basis_AF680_table_DDSI;
        basis(:,1) = (biofit_basis_AF680_table_DDSI{:,2});
        
    case "IR680"
        channel = "635";
        lambda = config.wavelengths.lambda_635;
        table = biofit_basis_IRDye680_table_M688;
        basis(:,1) = (biofit_basis_IRDye680_table_M688{:,2});
        
    case "OF650"
        channel = "635";
        lambda = config.wavelengths.lambda_635;
        table = biofit_basis_IRDye680_table_M688;
        basis(:,1) = (biofit_basis_IRDye680_table_M688{:,17});
        
        
    otherwise
end

% Set up the appropriate indexing and integral ranges
available_indexes = find(ismember(table{:,1},lambda(:)));
available_stack_indexes = find(ismember(lambda(:),table{:,1}));
available_waves = intersect(lambda(:),table{:,1});

% select the range to integrate over
bandpass = [config.range(find(config.channels == channel),1):config.range(find(config.channels == channel),2)];
integral_range= find(ismember(available_waves(:),bandpass(:)));
integral_waves= intersect(available_waves(:),bandpass(:));

% assign autofl bases to each channel
switch channel
    case "470"
        switch bodyRegion
            case "head"
                autofl(:,1) = (biofit_basis_470_table_M778{available_indexes,17}); % Brain
                autofl(:,2) = (biofit_basis_470_table_M778{available_indexes,16}); % Necrosis
                autofl(:,3) = (biofit_basis_470_table_M778{available_indexes,15}); % Muscle
                autofl_names = ["brain","necrosis","muscle"];
            case "hind"
                autofl(:,1) = (biofit_basis_470_table_M778{available_indexes,11}); % Hind Tumor
                autofl_names = ["hind"];
            case  "gi"
                autofl(:,1) = (biofit_basis_470_table_M778{available_indexes,10}); % Stomach
                autofl(:,2) = (biofit_basis_470_table_M778{available_indexes,12}); % Food
                autofl(:,3) = (biofit_basis_470_table_M778{available_indexes,13}); % GI
                autofl(:,4) = (biofit_basis_470_table_M778{available_indexes,15}); % Small Intestine
                
                autofl_names = ["Stomach","Colon","Upper GI","Muscle"];
            case  "phantom"
                autofl(:,1) = (biofit_basis_FITC_table_M688{available_indexes,10}); % Muscle
                autofl_names = ["well"];
            otherwise
        end
    case "530"
        switch bodyRegion
            case "head"
                autofl(:,1) = (pegRhd_530_5nm{available_indexes,6}); % Brain
                autofl(:,2) = (pegRhd_530_5nm{available_indexes,7}); % Nec
                autofl(:,3) = (pegRhd_530_5nm{available_indexes,9}); % Sinus
                autofl_names = ["brain","necrosis","sinus"];
           case "gfp"
                autofl(:,1) = (pegRhd_530_5nm{available_indexes,6}); % Brain
                autofl(:,2) = (pegRhd_530_5nm{available_indexes,7}); % Nec
                autofl(:,3) = (pegRhd_530_5nm{available_indexes,9}); % Sinus
                autofl(:,4) = (pegRhd_530_5nm{available_indexes,11}); % GFP
                autofl_names = ["brain","necrosis","sinus","gfp"];
            case "hind"
                autofl(:,1) = (pegRhd_530_5nm{available_indexes,8}); % Muscle
                autofl_names = ["muscle"];
            case  "gi"
                autofl(:,1) = (pegRhd_530_5nm{available_indexes,8}); % Muscle
                autofl(:,2) = (pegRhd_530_5nm{available_indexes,12}); % Stomach
                autofl_names = ["muscle","stomach"];
            case  "phantom"
                autofl(:,1) = (pegRhd_530_5nm{available_indexes,5}); % Well
                autofl_names = ["well"];
            otherwise
        end
    case "635"
        switch bodyRegion
            case "head"
                
 %               autofl(:,1) = (test_phantom635{available_indexes,4}); % Brain
               autofl(:,1) = (biofit_basis_IRDye680_table_M688{available_indexes,10}); % Brain
                autofl(:,2) = (biofit_basis_IRDye680_table_M688{available_indexes,8}); % Necrosis
                autofl(:,3) = (biofit_basis_IRDye680_table_M688{available_indexes,11}); % Muscle
                 autofl_names = ["brain","necrosis","body"];
            case "hind"
                autofl(:,1) = (biofit_basis_AF680_table_DDSI{available_indexes,3}); % Brain
                autofl_names = ["hind"];
            case  "gi"
                autofl(:,1) = (biofit_basis_IRDye680_table_M688{available_indexes,3}); % Brain
                autofl(:,2) = (biofit_basis_IRDye680_table_M688{available_indexes,7}); % Stomach
                autofl(:,3) = (biofit_basis_IRDye680_table_M688{available_indexes,12}); % Bowels
                autofl_names = ["brain","stomach","bowels"];
            case  "phantom"
                autofl(:,1) = (biofit_basis_IRDye680_table_M688{available_indexes,9}); % Brain
                autofl_names = ["well"];
            otherwise
        end
    otherwise
end

% Plot the results of the spectral bases in use.
figure('Name','Spectral Bases');
A = [basis(available_indexes),autofl];
A_names = [channelAgent,autofl_names];


plot(lambda(available_stack_indexes),A./max(A),'-o','LineWidth',4);
legend(A_names);
ylim([0 1.2]);
ylabel('Normalized Spectra'); xlabel('Wavelength (nm)');
set(gca,'FontSize',24);


%% Call in file info and pull out relevant data
mkdir(strcat(config.newSavePath,filesep,channel,filesep,channel,'_SpecFit_',bodyRegion));
mkdir(strcat(config.newSavePath,filesep,channel,filesep,channel,'_SpecFit_',bodyRegion,'_nM'));
cd(strcat(config.newSavePath,filesep,channel,filesep,channel,'_processed_LCTF_stacks'));

%%
if ~exist('ME') 
    ME = 1;
end

files = dir(strcat('*E',num2str(ME,1),'_cropped_stack.tif'));

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
    imshow(slice,[0 100]);
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
    f = 3;
    running = true;
    
    % Load up the slices used in a stack
    cd(strcat(filePath,filesep, channel , filesep, channel,'_processed_LCTF_stacks'));
    index = 1;
    for i = 1:length(available_stack_indexes)
        temp = imread(sortedFilenames(f),available_stack_indexes(i));
        temp = medfilt2(temp,[5 5]);
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
            data = double(data)';
        end
        
        data = reshape(data,length(data),1);
        
        % Collect indices that are above p = 0.01
        % F_spec_2: uses only top contibuting base + agent spectra if it is
        % statistically shown to be present ( p > 0.01)
        [spec_filtered, u_filtered, chi2, r2] = LLS_specfit(A,  data,0);
        
        
        % Plot the final spectral fit
        subplot2 = subplot(1,2,2,'Parent',figure1);
        
        x_axis = [1:size(A,1)]';
        Legend=cell(size(A,2)+ 2,1);
        Legend{1}='Data';
        Legend{2} = 'Fit';
        fit = sum(spec_filtered,2);
        
        plot(x_axis,data, '-k',x_axis,fit,':','LineWidth',3);
        hold on
        for j = 1:(size(A,2))
            plot(x_axis,spec_filtered(:,j),'-*','LineWidth',3);
            Legend{j+2}=strcat(A_names(j));
        end
        legend(Legend);
        xlabel('Wavelength [nm]');
        ylabel('Intensity [AU]');

        title('R^2 GOF = ', num2str(r2,2));
        set(gca,'linew',2,'FontSize',20);
        xticks(available_indexes(1:2:end));
        xticklabels(available_waves(1:2:end));
        hold off;
        
     %   clear data;
    end
else
    %%
    %% Begin paralellyzed group
    
%     % % Configure paralizable work
    delete(gcp('nocreate'));
    parpool('local',8); %current work station is an 8 core
    poolobj = gcp('nocreate');
%     
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
        chi2_map = zeros(size(temp,1),size(temp,2));
        r2_map = zeros(size(temp,1),size(temp,2));
        unique_map = zeros(size(temp,1),size(temp,2));
         
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
                    %[Fpd, indices, spec_filtered, u_filtered] = F_spec_2(A,  data);
                    %[Fpdf, indices] = F_spec(A, data);
                    [spec_filtered, u_filtered, chi2, r2, uniqueFlag] = LLS_specfit(A,  data,0);
                    
                    % save the unmixed signals
                    %u_trap = trapz(integral_waves, spec_filtered(1,integral_range,:))./(10);
                    u_trap = trapz(integral_waves, spec_filtered(integral_range,:))./(10);
                    
                    % store the integrated values
                    u_map(r,c,indices) = u_trap;
                    chi2_map(r,c) = chi2;
                    r2_map(r,c) = r2;
                    unique_map(r,c) = uniqueFlag;
                end
            end
        end
        toc
        
       percentage =  sum(unique_map(:))/  sum(bw(:))
      
       imwrite(uint16(chi2_map*1000),strcat(filePath,filesep,channel, filesep, channel,'_SpecFit_', bodyRegion, filesep,filename,'_Slice_',num2str(f), '_',num2str(bandpass(1)),'_',num2str(bandpass(end)),'nm_ME',num2str(ME),'_chi2_dofx1000_LLS_unmixed_stack.tif'));
       imwrite(uint16(r2_map*1000),strcat(filePath,filesep,channel, filesep, channel,'_SpecFit_', bodyRegion, filesep,filename,'_Slice_',num2str(f), '_',num2str(bandpass(1)),'_',num2str(bandpass(end)),'nm_ME',num2str(ME),'_r2x1000_LLS_unmixed_stack.tif'));
       imwrite(uint16(unique_map*1000),strcat(filePath,filesep,channel, filesep, channel,'_SpecFit_', bodyRegion, filesep,filename,'_Slice_',num2str(f), '_',num2str(bandpass(1)),'_',num2str(bandpass(end)),'nm_ME',num2str(ME),'_unique_LLS_unmixed_stack.tif'));
      
        % Write out the unmixed spectral stacks in integrated RFUs
        for ind = 1:size(u_map,3)
            imwrite(uint16(imboxfilt(medfilt2(u_map(:,:,ind),[4 4]),3)),strcat(filePath,filesep,channel, filesep, channel,'_SpecFit_', bodyRegion, filesep,filename,'_Slice_',num2str(f), '_',num2str(bandpass(1)),'_',num2str(bandpass(end)),'nm_ME',num2str(ME),'_', A_names(ind),'_LLS_unmixed_stack.tif'));
        end
        
%         if ismember(channelAgent,["PEG550","PEG5000","FITC","Bodipy","OF550","Rhd550","OF650","AF647","AF680","IR680"])
%             %%Provide nM conversion coefficients
%             % added cases to convert RFU to nM concentrations
%             switch channelAgent
%                 case "PEG550"
%                     m = 0.2033;
%                     b = 0;
%                 case "PEG5000"
%                     m = 0.2033;
%                     b = 0;
%                 case "FITC"
%                     m = 0.0462;
%                     b = 0.485;
%                 case "AF488"
%                     m = 0.0783;
%                     b = -12.016;
%                 case "Bodipy"
%                     m = 0.2741;
%                     b = -14.41;
%                 case "OF550"
%                     m = 0.0127;
%                     b = 0.2429;
%                 case "Rhd550"
%                     m = 0.0133;
%                     b = -0.902
%                 case "OF650"
%                     m = 0.0473;
%                     b = 0.53;
%                 case "AF647"
%                     m = 0.0231;
%                     b = 0;
%                 case "AF680"
%                 case "IR680"
%                     m = 0.0291;
%                     b = 0;
%                 otherwise
%             end
%             % Adjust the slope based off of ME scaling
%             m = m.*10^ME;
%             
%             %  Write image in nM concentration units
%             
%             % Save the file in the NextImage subfolder
%             Rs_nM = (u_map(:,:,1)*m + b);
%             
%             % clean up the image using a gaussian and sharpening filter
%             Rs_nM = imboxfilt(medfilt2(Rs_nM,[4 4]),3);
%             
%             % Write out the nM converted next-image corrected slice
%             imwrite(uint16(Rs_nM), strcat(filePath,filesep,channel, filesep, channel,'_SpecFit_', bodyRegion,'_nM', filesep,filename,'_Slice_',num2str(f), '_',num2str(bandpass(1)),'_',num2str(bandpass(end)),'nm_ME',num2str(ME),'_', A_names(1),'_nM_LLS_stack.tif'));
%         end
        
        % Slice number
        slice = f
        
    end
    
    %%
    % Apply Next-Image correction to newly processed images
    nextImageCorrection(config,channelAgent, bodyRegion)
    
    
end

end


