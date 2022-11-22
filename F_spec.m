function [Fpdf, indices, LSE_initial, errorFlag] = F_spec(A, data)
% Function written by Brook to return the indices which are shown to
% significantly contribute to the observed variance.
%%

[spec,u, LSE_initial, LSE_max, LSE_auc, mRPE] = LLS_specfit( A, data ,0);  
 u_trap = trapz(spec(:,:))./(10);  
 
 %% Collect tree based on smallest to largest contributions

 % index ranks the contributions
 clear tree
[~,index] = sort(u_trap,'descend');
indices = [1:size(A,2)];


i = 1;
for count = size(A,2):-1:1
    tree(i,1:count) = indices(index(1:end - i + 1));
    i = i + 1;
end


%% Calculate the RSS for all combos
j = size(A,2);
%figure('Name','Spectral Fits');
for t = 1:size(tree,1)
    %subplot(1, size(tree,1),t);
    clear matrix
    matrix = A(:,tree(t,1:j));
    %tree_names = cellstr(A_names(tree(t,1:j)));
    plotflag = 0;
    [spec,u, LSE(t), ~ , LSE_auc(t)] = LLS_specfit( matrix, data,plotflag); 
    j = j - 1;
    %pause;
end
%%
LSE = LSE';
LSE_auc  = LSE_auc';
%%
% Calculate the F stat for these given combinations
k = 1;
for r = size(A,2):-1:2   
    
    num_DOF = 1;
    num = (LSE(r) - LSE(r-1))./(1);
    
    denom_DOF = size(A,1) - k - 1;
    denom = (LSE(r-1))./(denom_DOF);

    F(r,1) = num./denom;
    Fpdf(r,1) = fpdf(F(r,1),(num_DOF),(denom_DOF));
%     subplot(1, size(tree,1),r);
%     text(0,0, strcat('F stat:', {' '},num2str(Fpdf(r,1),2)),'FontSize',24); 
    k = k + 1;
end


%% Return the indices necessary
filtered_tree = tree;
Fpdf;
for r = 2:size(tree,1)
    % if P is not significant then the variance between the two fits is not significant     
    if(Fpdf(r,1) > 1E-8 || Fpdf(r,1) == 0 || isnan(Fpdf(r,1))) 
        
        % The lower the p, the more values are considered non-significant
         act(r-1,1) = 1;
    else
        act(r-1,1) = 0;
    end
end

%ind =  find(act == 1);
ind = 1:max(find(act == 1));
filtered_tree(ind,:) = [];
indices = intersect(1:size(A,2),unique(filtered_tree));

end