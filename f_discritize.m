function [ T_discrete_temp_data ] = f_discritize( csv, n_levels )
% temp_data: n by m matrix 
% n_levels : # for gap
% discrete_temp_data : n by m matrix, discritized.

% if n_levels = 2 -> val(discrete_temp_data) = { 0, 1 }
% if n_levels = 3 -> val(discrete_temp_data) = { 0, 1, 2 }
T_table = readtable(csv,'ReadRowNames',true,'ReadVariableNames',true);

temp_data = table2array(T_table);
agis = T_table.Properties.RowNames;
var = T_table.Properties.VariableNames;

[n_genes, T] = size(temp_data);
discrete_temp_data = zeros(n_genes,T);
%% Discretize original expression data in n_levels levels

   if n_levels == 2
        m = repmat( mean(temp_data,2),1,T);
        discrete_temp_data = temp_data >= m;
   elseif n_levels > 2 
    q = quantile(temp_data,n_levels-1,2); % quantile(X,N,dim)
    for i = 1 : T
        for j = 1:n_genes
            if temp_data(j,i) < q(j,1)
                discrete_temp_data(j,i) = 0;
            end
            for k = 2:n_levels-1
                if temp_data(j,i) > q(j,k-1) && temp_data(j,i) < q(j,k)
                    discrete_temp_data(j,i) = k-1;
                end
            end
            if temp_data(j,i) > q(j,k)
                discrete_temp_data(j,i) = k;
            end
        end
    end
   end

    T_discrete_temp_data = array2table(discrete_temp_data,'RowNames',agis,'VariableNames',var);
    [pathstr,name,ext] = fileparts(csv); 

    mkdir(sprintf('%s/Dscrtz/',pathstr));
    writetable(T_discrete_temp_data,sprintf('%s/Dscrtz/%s-Dscrtz%s',pathstr,name,ext)...
        ,'WriteRowNames',true,'WriteVariableNames',true);
end

