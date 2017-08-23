csv = '140269-4.xlsx';
T = readtable(csv,...
     'ReadVariableNames',true,'ReadRowNames',true);
T_mapping = [];
for i = 1 : size(T,1)
      if ~strcmp(T{i,1},'')
          T_mapping=[T_mapping ; T(i,:)];
      end
end

T_mapping = T_mapping(:,1);
T_mapping.Properties.VariableNames = {'ntf'}; % other notification for each gene

writetable(T_mapping,'table_mapping.csv','WriteRowNames',true,'WriteVariableNames',true);

%% test f_tranlate
[agis_new,T_map] = f_tranlate({'AT1G68640','AT4G33950','aetejjjjj'},'table_mapping.csv');