function [ Nij ] = f_getNijk( config, ScoreMatrix, n_levels )
% Nij(k) = # of x_i = k & Parents of x_i is in jth config
% config : jth configuration fo pa of x_i, a vector
[ngenes, tp] = size(ScoreMatrix);
nReg = ngenes - 1 ;

if length(config)~=nReg
    error('Confict with function design!');
end

% find time point indx that has jth indx
tmp = repmat(config,1,tp)==ScoreMatrix(1:nReg,:);
indx_confj = ones(1,tp);
for i = 1 : nReg
    indx_confj = indx_confj .* tmp(i,:);
end

% find time point indx x has kth value
indx_x = cell(1, n_levels); % indx_x{k} contains index info of x=k
for k = 1 : n_levels
    indx_x{k} = ScoreMatrix(end,:) == (k-1);
end

% find time point indx x has kth value $ parents of x has jth config
Nij = zeros(1,n_levels);
indx_confj_x = cell(1,n_levels);
for k = 1 : n_levels
    indx_confj_x{k} = indx_x{k} .* indx_confj;
    Nij(k) = sum(indx_confj_x{k});
end



end

