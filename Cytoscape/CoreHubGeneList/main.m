%%
myFiles = dir(fullfile('.','CoreHubList*.txt')); 
for i = 1 : length(myFiles)
    unix(sprintf('awk -F "." ''{print $1}'' %s | sort -u > GO_%s',myFiles(i).name, myFiles(i).name));
end

%%
myFiles = dir(fullfile('.','GO_CoreHubList*.txt'));
for i = 1 : length(myFile) 
    unix(sprintf('grep -f %s ATH_GO_GOSLIM.txt | sort -u > BPD_%s',myFiles(i).name,myFiles(i).name));
end

%%
myFiles = dir(fullfile('.','BPD_GO_CoreHubList*.txt'));
for i = 1 : length(myFile) 
    unix(sprintf('grep ''0009873'' %s | sort -u > ETHYLENE_%s',myFiles(i).name,myFiles(i).name));
end

%%
myFiles = dir(fullfile('.','ETHYLENE_BPD_GO_CoreHubList*.txt'));
for i = 1 : length(myFile)
    T = readtable(myFiles(i).name,'ReadRowNames',false,'ReadVariableNames',false);
    tmp = unique(T(:,1));
    writetable(tmp,sprintf('GO_GENE_%s',myFiles(i).name),'WriteRowNames',false,'WriteVariableNames',false);
end

%%
myFiles = dir(fullfile('.','GO*.txt'));
for i = 1 : length(myFiles)
    unix(sprintf('wc -l %s >> GO_Summary.txt',myFiles(i).name));
end

%%
FilesGL = dir(fullfile('.','GO_GENE_ETHYLENE_BPD_GO_CoreHubList*.txt'));
FilesHUBLIST = dir(fullfile('.','CoreHubList*.txt'));
for i = 1 : length(FilesGL)
    gl = FilesGL(i).name;
    hbl = FilesHUBLIST(i).name;
    unix(sprintf('grep -f %s %s | sort -u > CYTOSCAPE_%s',gl,hbl,hbl));
end



