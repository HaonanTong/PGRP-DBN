function [agis_new, T_map] = f_tranlate(agis,transfile)
% Translation file has to have Variable Names.
% agis : cell string 1 by n
% transfile : path + .csv file
%       Example of Translation file opened in MS Excell might look as follows:
%       |    ORF    |    ntf     |
%       | AT1G66340 |    ETR1    |
%       |    ...    |    ...     |
%       | AT3G20770 |    EIN3    |
%       |    ...    |    ...     |
% agis_new : cell string contains correponding new name, cell string
    T_ORF2CNames = readtable(transfile,...
     'ReadVariableNames',true,'ReadRowNames',true);
 
    agis_org = agis;
 
    T_ORF2CNames.Properties.VariableNames = {'ntf'};
    agis_new = cell(1,length(agis));
    for i = 1 : length(agis)
        agis{i} = agis{i}(1:9);    % filter the dot
        if any(ismember(agis{i},T_ORF2CNames.Properties.RowNames))
            agis_new(i) = T_ORF2CNames(agis{i},:).ntf;
        else
            agis_new{i} = agis{i};
        end
    end

    agis_new = agis_new';
    T_map = table(agis_org',agis',agis_new,'VariableNames',{'obj','orf','nrf'});
    return

end

