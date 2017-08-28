myFile = dir(fullfile('./','A_Stage*/hubgobjs*.txt'));
for i = 1 : length(myFile)
    display(myFile(i).folder);
    display(myFile(i).name);
end

for i = 1 : length(myFile)
    file = sprintf('%s/%s',myFile(i).folder,myFile(i).name);
    file_go = strcat(myFile(i).folder,'/hubgenes_go');
    unix(sprintf('awk -F "." ''{print $1}'' %s | sort -u > %s.txt',file, file_go));
    unix(sprintf('awk ''END{print NR}'' %s.txt > %s.txt',file_go,'hubgenes_num'));
end