function [ Table_interp ] = f_interp( csv )
global outputDir
    % 0 0.25 0.5,1, 4, 12 ,24
    T = [ 0 0.25 0.5,1, 4, 12 ,24 ];
    T_new = 0 : 0.25 : 24;
    
    T_table = readtable(csv,'ReadRowNames',true,'ReadVariableNames',true);
    
    Data = table2array(T_table);
    agis = T_table.Properties.RowNames;
    
    [n_genes,~] = size(Data); 
    
    expr = [];
    for i = 1:3:21%7 time points; 3 replicates;
         expr = [expr sum(Data(:,i:i+2),2)];
    end
    expr = 1/3*expr;
    
    
    Matrix_interp = zeros(n_genes,length(T_new));
    for i = 1 : n_genes
        Matrix_interp(i,:) = interp1(T,expr(i,:),T_new,'linear');
    end
    
%     Table_interp = array2table(Matrix_interp,'RowNames',agis,'VariableNames',...
%         strcat('T',num2strcell(T_new,'%2.2f')));
%     
    Table_interp = array2table(Matrix_interp,'RowNames',agis,'VariableNames',...
        strcat('IT',num2strcell(1:length(T_new),'%d')));
    
    [pathstr,name,ext] = fileparts(csv); 

    mkdir(sprintf('%s/Interp/',outputDir));
    writetable(Table_interp,sprintf('%s/Interp/%s-Interp%s',outputDir,name,ext)...
        ,'WriteRowNames',true,'WriteVariableNames',true);
end

