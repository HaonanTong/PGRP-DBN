% T = readtable('timecourse.csv','ReadRowNames',true,'ReadRowNames',true);
% The examples is trying to test f_getTopology function. 
%% Example 1
ADM = [ zeros(1,3) ones(1,3) zeros(1,3) ones(1,3) ];
ADM = [ADM; ADM];

n_levelS = 2;
rmax = 1;

topology = f_getTopology( ADM, n_levels, rmax );
%-2.7726;


% ESS = 1E-15;
% j = 1;
% BDeu_config(j) = log(gamma(ESS/2)/gamma(6+ESS/2)) ...
%    + log(gamma(6+ESS/2)/gamma(ESS/4)) +...
%    log(gamma(ESS/2)/gamma(ESS/4));
    
%% Example2
ADM = [ ones(1,3) zeros(1,3) ones(1,3) zeros(1,3) ];
ADM = [ ADM ; zeros(1,3) ones(1,3) zeros(1,3) ones(1,3) ];

n_levelS = 2;
rmax = 1;

topology = f_getTopology( ADM, n_levels, rmax );

%% Example 3
% G1->G3
G1 = [ 0 0 1 1 0 0 1 1 ];
G2 = [ 0 1 0 1 0 1 0 1 ];
G3 = [ 1 1 1 0 0 0 1 0 ];
ADM = [G1;G3];
n_levelS = 2;
rmax = 1;
ESS = 1E-15;

topology = f_getTopology( ADM, n_levels, rmax );
topology.BDeu

% G2 -> G3
ADM = [G2;G3];
n_levelS = 2;
rmax = 1;

topology = f_getTopology( ADM, n_levels, rmax );
topology.BDeu

% G1 & G2 -> G3
% ESS = 1E-4;
BDeu_j1 = log(gamma(ESS/4)/gamma(2+ESS/4)) + ...
    log(gamma(1+ESS/4)/gamma(ESS/8)) *2 ;
BDeu_j2 = BDeu_j1;
BDeu_j3 = log(gamma(ESS/4)/gamma(2+ESS/4)) +...
    log(gamma(2+ESS/4)/gamma(ESS/8)) + ...
    log(gamma(ESS/4)/gamma(ESS/8));
BDeu_j4 = BDeu_j3;
BDeu = BDeu_j1 + BDeu_j2 + BDeu_j3 + BDeu_j4

%% Example 3

G1 = [ 0 0 1 1 0 0 1 1 ];
G2 = [ 0 1 0 1 0 1 0 1 ];
G3 = [ 1 1 1 0 0 0 1 0 ];

n_levels = 2;
ESS = 1E-15;
%ESS = 1E-10;
%ESS = 1E-7;
%ESS = 1E-4;
%ESS = 1E-1;
BDeu13 = f_calculateScore( [G1;G3], n_levels, ESS );
BDeu23 = f_calculateScore( [G2;G3], n_levels, ESS );
BDeu123 = f_calculateScore( [G1;G2;G3], n_levels, ESS );

ADM = [G1;G2;G3];
n_levels = 2;
rmax = 2;

[N1k] = f_getNijk([0 0]',ADM,n_levels);
[N2k] = f_getNijk([0 1]',ADM,n_levels);
[N3k] = f_getNijk([1 0]',ADM,n_levels);
[N4k] = f_getNijk([1 1]',ADM,n_levels);


G1 = [ 0 0 1 1 0 0 1 1 ];
G2 = [ 0 1 0 1 0 1 0 1 ];
G3 = [ 1 1 1 0 0 0 1 0 ];
[ BDeu ] = f_calculateScore( [G1;G3], 2 );
[ BDeu ] = f_calculateScore( [G1;G2;G3], 2 );
[topology,BDeu_memory] = f_getTopology( ADM, n_levels, rmax );

%% Example 4
G1 = [ 0 0 0 1 1 1 0 0 0 1 1 1 1 1 1 ];
G2 = [ 1 1 1 0 0 0 1 1 1 0 0 0 1 0 1 ];
G3 = [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ];
G4 = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];
G5 = [ 0 0 0 1 1 1 0 0 0 1 1 1 1 1 1 ];

% my Prediction G1->G5
% G2 -> G5
% G1,G2 -> G5
%BDeu15 = f_calculateScore( [G1;G5], n_levels, ESS );
%BDeu25 = f_calculateScore( [G2;G5], n_levels, ESS );
%BDeu125 = f_calculateScore( [G1;G2;G5], n_levels, ESS );

ADM = [G1;G2;G3;G4;G5];
n_levelS = 2;
rmax = 3;
ESS = 1E-15;

[topology,BDeu_Memory] = f_getTopology( ADM, n_levels, rmax, ESS );

for i = 1 : rmax
    for j = 1 : nchoosek(4,i)
        BDeu_Memory{i,j}.BDeu
    end
end