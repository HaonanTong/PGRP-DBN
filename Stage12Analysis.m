%% Interpolation on three time course data
myDir = './A_Stage12/'; %gets directory
myFiles = dir(fullfile(myDir,'*.csv')); %gets all wav files in struct

for i = 1 : length(myFiles)
    f_interp(sprintf('%s%s',myDir,myFiles(i).name));
end

%% Discretize
myDir2 = './A_Stage12/Interp/';
myFiles2 = dir(fullfile(myDir2,'*.csv'));
n_levels = 3;

for i = 1 : length(myFiles2)
    csv = sprintf('%s%s',myDir2,myFiles2(i).name);
    f_discritize( csv, n_levels );
end

%% 
myDir2 = './A_Stage12/Interp/Dscrtz/';
% T.5 TFs-DEGs-time-Activation2-up-Interp-Dscrtz.csv    
csv = sprintf('%s%s',myDir2,'TFs-DEGs-time-Activation2-up-Interp-Dscrtz.csv');
T_table2 = readtable(csv,'ReadRowNames',true,'ReadVariableNames',true);

TF_at2.tp = 0.5; %T0.25
TF_at2.it = 3; %IT3;
TF_at2.stage = 1;
TF_at2.dmtrx = table2array(T_table2); 
TF_at2.table = T_table2; 
TF_at2.ngenes = size( TF_at2.dmtrx , 1 );
TF_at2.glist = T_table2.Properties.RowNames;


% T1 TFs-DEGs-time-Activation3-up-Interp-Dscrtz.csv
csv = sprintf('%s%s',myDir2,'TFs-DEGs-time-Activation3-up-Interp-Dscrtz.csv');
T_table3 = readtable(csv,'ReadRowNames',true,'ReadVariableNames',true);

TF_at3.tp = 1; %T1
TF_at3.it = 5; %IT5;
TF_at3.stage = 1;
TF_at3.dmtrx = table2array(T_table3); 
TF_at3.table = T_table3; 
TF_at3.ngenes = size( TF_at3.dmtrx , 1 );
TF_at3.glist = T_table3.Properties.RowNames;


% T4 TFs-DEGs-time-Activation4-up-Interp-Dscrtz.csv
csv = sprintf('%s%s',myDir2,'TFs-DEGs-time-Activation4-up-Interp-Dscrtz.csv');
T_table4 = readtable(csv,'ReadRowNames',true,'ReadVariableNames',true);

TF_at4.tp = 4; %T4
TF_at4.it = 17; %IT17;
TF_at4.stage = 2;
TF_at4.dmtrx = table2array(T_table4); 
TF_at4.table = T_table4; 
TF_at4.ngenes = size( TF_at4.dmtrx , 1 );
TF_at4.glist = T_table4.Properties.RowNames;


PoPaGList = [ TF_at2.glist ; TF_at3.glist];
target = cell(1,TF_at4.ngenes);

A.matrix = TF_at2.dmtrx;
A.at = TF_at2.it;
B.matrix = TF_at3.dmtrx;
B.at = TF_at3.it;
[ PPaMtrx ] = f_alignment( A, B, 0 );

rmax = 3;
ESS = 1E-10;
for i = 1 : TF_at4.ngenes % for each potential target
    % alignment potential regulator
    Atarget = TF_at4.dmtrx(i,:);
    B.matrix = Atarget;
    B.at = TF_at4.it;
    ADM = f_alignment(PPaMtrx, B, 0);
    % get topology
    [ topologyi, ~ ] = f_getTopology( ADM.matrix, n_levels, rmax, ESS ); 
    target{i}.name = TF_at4.glist{i};
    target{i}.indx = topologyi.indx;
    target{i}.BDeu = topologyi.BDeu;
    target{i}.Pa = PoPaGList(target{i}.indx);
    target{i}.PoPaGList = PoPaGList;
end

% save('Analysis_Stage12_TFs_UP.mat');
%% Analysis

% target{27} and target{29} has one parent others have 2 parents.

BDeu = zeros(TF_at4.ngenes,1);
for i = 1 : TF_at4.ngenes % for each potential target
    BDeu(i) = target{i}.BDeu;
end

BDeu_max = max(BDeu);
indx = BDeu == BDeu_max; % corresponds to strongest causal relationship;

DBN_target = target(indx);
fprintf('Totally have %d sub-structure\n',length(DBN));

for i = 1 : length(DBN_target)
   fprintf('-----------------------\n');
   fprintf('structure: %d\n',i);
   fprintf('Regulators:\n');
   disp(DBN_target{i}.Pa);
   fprintf('Target: %s\n',DBN_target{i}.name);
   fprintf('BDeu score: %G\n',DBN_target{i}.BDeu);  
end

% grep 3 structure.
structure = cell(1:length(DBN_target));
for i = 1 : length(DBN_target)
    structure{i} = [ DBN_target{i}.Pa ; DBN_target{i}.name ]; % last one is target;
end

% output expression profiles.
T = readtable('kat-rpkm-expression.csv','ReadVariableNames',true,'ReadRowNames',true);
for i = 1 : length(DBN_target)
    sub_T = T(structure{i},:);
    
    writetable(sub_T,sprintf('./A_Stage12/Sub_Table%d.csv',i)...
        ,'WriteRowNames',true,'WriteVariableNames',true);
end

% Visulization
myDir = './A_Stage12/'; %gets directory
myFiles = dir(fullfile(myDir,'Sub_Table*.csv')); %gets all wav files in struct
for i = 1 : length(myFiles)
    csv = sprintf('%s%s',myDir,myFiles(i).name);
    f_plotTable2( csv, [], 'Normalized legend' );
end
