%% Example1
A.matrix = 1:9;
A.at = 4;
B.matrix = A.matrix + 10*ones(1,9);
C = cell(1,9);
for i = 1 : 9
    B.at = i;
    C{i} = f_alignment(A,B,0);
end

%% Example2
A.matrix = [1:9;11:19;21:29];
A.at = 4;
B.matrix = [31:39;41:49];
C = cell(1,9);
for i = 1 : 9
    B.at = i;
    C{i} = f_alignment(A,B,0);
end

for i = 1 : 9
    display(C{i}.matrix);
end


