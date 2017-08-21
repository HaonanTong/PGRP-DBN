function [ topology, BDeu_Memory ] = f_getTopology( ADM, n_levels, rmax, ESS )
% topology is a r by tp matrix, where r-1 regulators regulates rth gene
fprintf('Estimating topology based on BDeu score with equivalent sample size %G\n', ESS);
[n,~] = size(ADM);
nPReg = n - 1;

BDeu = -Inf;
BDeu_Memory = cell(1,rmax);
for nPa = 1 : rmax % # of parents
    Cmbntn_Pa = nchoosek(1:nPReg,nPa);
    [nCmb,~] = size(Cmbntn_Pa);
    for cmb = 1 : nCmb
        % get ScoreMatrix
        indx = Cmbntn_Pa(cmb,:);
        ScorePaMatrix = ADM(indx,:);
        ScoreMatrix = [ ScorePaMatrix ; ADM(end,:) ];
        BDeu_new = f_calculateScore( ScoreMatrix, n_levels, ESS );
        if BDeu_new > BDeu
            BDeu = BDeu_new;
            topology.indx = indx;
            topology.BDeu = BDeu;
        end
        BDeu_Memory{nPa,cmb}.ScoreMatrix = ScoreMatrix;
        BDeu_Memory{nPa,cmb}.Indx = indx;
        BDeu_Memory{nPa,cmb}.BDeu = BDeu_new;
    end
end


end

