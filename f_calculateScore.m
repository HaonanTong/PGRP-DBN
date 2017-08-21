function [ BDeu ] = f_calculateScore( ScoreMatrix, n_levels, ESS )
% topology is determined by ScoreMatrix
[n,~]=size(ScoreMatrix);
nPa = n-1;% # of parent in the topology

config_Pa = f_getConfiguration(n_levels, nPa);

[~, nConfig] = size(config_Pa);
C = ESS/nConfig;

Nijk_config = cell(1,nConfig);
BDeu_config = zeros(1,nConfig);
for j = 1 : nConfig
    Nijk_config{j} = f_getNijk(config_Pa(:,j),ScoreMatrix,n_levels);
    BDeu_config(j) = log(gamma(C)/gamma(sum(Nijk_config{j})+C));
    for k = 1 : n_levels
        BDeu_config(j) = BDeu_config(j) + ...
            log(gamma(Nijk_config{j}(k) + C)/gamma(C/n_levels));
    end
end

BDeu = sum(BDeu_config);
end