myDir = './'; %gets directory
myFiles = dir(fullfile(myDir,'TFs-DEGs-time-Activation*.csv')); %gets all wav files in struct

geneID = [];
at = [];
for i = 1 : length(myFiles)
    T = readtable(sprintf('%s%s',myDir,myFiles(i).name),...
        'ReadRowNames',true,'ReadVariableNames',true);
    tmp = T.Properties.RowNames;
    geneID = [geneID ; tmp];
    at = [ at ; i * ones(length(tmp),1) ];
end

T_g2at = table(geneID,at);
writetable(T_g2at,'Table_TF2at.csv','WriteRowNames',true,'WriteVariableNames',true);