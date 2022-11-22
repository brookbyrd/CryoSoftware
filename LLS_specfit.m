function [spec,u, reduced_chi_squared, rsq, uniqueFlag, mRPE] = LLS_specfit(A,data,flag,basis_names)

% Performs least squares fitting between measured spectrum (data), and
% measured basis spectra (in matrix A).

% input basis spectra in matrix A with columns denoting spectra.
% input data as column vector
% set flag to 1 to plot results when the fitting is done
% output data is matrix A with each column scaled by its minimization
% coefficient
% S Davis 2013
% Optimized to do spectral fit for any number of basis spectrums and
% generate plot for the fit (Brook Byrd, 2019)

for column = 1:size(A,2)
    A(:,column) = A(:,column)./sum(A(:,column));
end

[npix,m] = size(A);

H = A'*A;
b = A'*data;
u = H\b;

% Constrain negative fits

for i = 1:m
    if u(i) < 0
        act(i) = 0;
    else
        act(i) = 1;
    end
end

if sum(act) < m
%     disp('Negative fit detected');
    clear u
    ind = find(act == 0);
    A2 = A;
    A2(:,ind) = [];
    H = A2'*A2;
    b = A2'*data;
    u = H\b;
    k = 1;
    for i = 1:m
        if act(i) == 1
            u2(i) = u(k);
            k = k+1;
        else
            u2(i) = 0;
        end
    end
    u = u2;
end



% Collect unmixed spectra
spec=[];
for i = 1:m
    spec = [spec, u(i)*A(:,i)];
end

% Constrain to positive values of unmixed spectra
spec(spec<0) = 0;

fit = sum(spec,2);

%% Solution Uniqueness test
% Check for solution uniqueness
[Xsub,idx]=licols(spec,1e-2);

uniqueBases = length(idx);
presentBases = length(find(spec(1,[1:m]) ~= 0));

if(length(idx) ~= length(find(spec(1,[1:m]) ~= 0)))
    uniqueFlag = 0;
    %warning('Uniqueness is not maintained, choosen bases are not independent');
else 
    uniqueFlag = 1;
end


%% Calculate the Least Square Error (Residual sum of squares)

LSE_array  = (data - fit);
LSE = sum((data - fit).^2);

%Normalize fit
fit_norm = (fit - min(fit(:)))./(max(fit(:)) - min(fit(:)));
data_norm = (data - min(fit(:)))./(max(fit(:)) - min(fit(:)));

% Also defined as std^2 at https://en.wikipedia.org/wiki/Ordinary_least_squares
LSE_norm = sum((data_norm - fit_norm).^2)./(npix);

% Chi squared statistic https://towardsdatascience.com/chi-square-test-for-feature-selection-in-machine-learning-206b1f0b8223
chi_squared = sum((data_norm - fit_norm).^2)./var(data_norm);

% Reduced chi statistic (supposed to be less biased)
reduced_chi_squared = chi_squared./(npix - m);

%% Evaluating different error metrics
% Chad Kanick shown here: https://www.nature.com/articles/s41598-017-09727-8
mRPE_array = (abs(data - fit)./fit);%/npix;
mRPE = sum(abs(data - fit)./fit)/npix;

%Normalize fit over DOF
fit_max = fit./max(fit(:));
data_max = data./max(fit(:));

LSE_max_array = (data_max - fit_max);
LSE_max = sum((data_max - fit_max).^2)./(npix - m);

fit_auc = fit./sum(fit(:));
data_auc = data./sum(fit(:));

LSE_auc_array = (data_auc - fit_auc);
LSE_auc = sum((data_auc - fit_auc).^2)./(npix - m);


%% Calculate RMSE
% rmse_norm = sqrt(mean((fit_norm(:)-data_norm(:)).^2));
% rmse = sqrt(mean((fit(:)-data(:)).^2));
% 
y_hat = (1/length(data_norm(:))).*sum(data_norm(:));
ss_t = sum((data_norm(:) - y_hat).^2);
ss_res = sum((data_norm(:) - fit_norm(:)).^2);

rsq = 1 - ss_res./ss_t;
%%
% Plot flag
if flag == 1
    fit = sum(spec,2);
    
    x_axis = [1:npix]';
    Legend=cell(m+2,1);
    Legend{1}='Data';
    Legend{2}='Fit';
    plot(x_axis,data, '-k',x_axis,fit,':','LineWidth',3);
    hold on
    for j = 1:m
        subplot(1,m,j)
        plot(x_axis,spec(:,j),'-*','LineWidth',3);
        Legend{j+2}=strcat('Basis', num2str(j));
    end
   if(exist('basis_names'))
       for i = 1:length(basis_names)
           Legend(i+2) = (basis_names(i));
       end
   end
    ylabel('counts/s'); xlabel('pixel number');
    legend(Legend)
    set(gca,'linew',2)
    %set(gca,'fontsize',30)
    set(gcf,'color','w');
    %title(strcat('Normalized Error: ',{' '},num2str(LSE_auc,2)));
    hold off;
end

