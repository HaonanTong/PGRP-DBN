function [ fig, ngene, expr,plotData, agis, agis_new ] = ...
    f_plotTable2( csv, transfile, varargin )
warning('off')
% #####################################
% @PGRP
% Haonan Tong
% #####################################

% [ fig, ngene, expr, agis ] = plotTable( csv, isVariableNames )
% fig - figure of profiles in csv file
% ngene - # of genes in the file
% expr - profiles of genes in RPKM at each time point
% agis - gene list
% csv - exel file with 3 replicates for each time point, totally 6 time
% point
%       Example of .csv file opened in MS Excell might look as follows:
%       |    AIG    |  R1T1  |  R2T1  |  R3T1  | ... 
%       |AT1G01010.1| 1.2691 | 1.7789 | 3.0794 | ...
%       |    ...    |  ...   |  ...   |   ...  | ...
%       |AT5G67360.1| 150.93 | 243.19 |   ...  | ...
%       |    ...    |  ...   |  ...   |   ...  | ... 
% transfile - Translation file could be empty but if not has to have Variable Names.
%       Example of Translation file opened in MS Excell might look as follows:
%       |    ORF    | OtherNames |
%       | AT1G66340 |    ETR1    |
%       |    ...    |    ...     |
%       | AT3G20770 |    EIN3    |
%       |    ...    |    ...     |
% Plot Pattern Selection
%   'Mean Plot' - log2 ratio plot of profile for each gene and mean of the
%           gene list
%   'Normalized' - Normalized each profile of gene corresponds to the mean.
%   'Normalized legend' - Normalized each profile of gene corresponds to the mean.

%   'Network Analysis log' - plot profile of each genes on a single
%       figure and take log ratio onto the initial time point ; 
%   'Plot Individually log' - plot profile of each genes individually and take log ratio onto the initial time
%   'Plot Individually' - plot profile of each genes individually in RPKM unit 
%   'No Variable Name' - Expression Profile Table not contains Variable
%                   Names
if nargin < 2
    fprintf('#####################################\n');

    fprintf('\ntwo parameter required for the function\n');
    fprintf(' path of .csv file as the first parameter\n');
    fprintf(' if the file contains VariableNames as second parameter\n');
    fprintf(' example:  [ fig, ngene, expr, agis ] = plotTable(''kat-rpkm-expression.csv'', 1)\n\n');
    fprintf('#####################################\n');

    return;
end




fprintf('Preparing ploting\n');

% csv_token = strtok(csv,'.');
% segments = string(0);
% remain = csv_token;
% while (remain ~= "")
%    [token,remain] = strtok(remain, '/');
%    segments = [segments ; token];
% end
% csv_token = segments(end);

[pathstr,csv_token,~] = fileparts(csv);
if isempty(pathstr)
    dir_Figures = 'Figures';
else
    dir_Figures = strcat(pathstr,'/Figures');
end
mkdir(dir_Figures);
 
% fig = cell(nargin-2);
counter = 0;

%% Read File
% =========  Read File ===========
fprintf(sprintf(' Reading file %s...\n',csv))
if any(strcmp(varargin,'No Variable Name'))
    isVariableNames = 0;
else
    isVariableNames = 1;
end

if isVariableNames == 1
    T = readtable(csv,...
     'ReadVariableNames',true);
else
    T = readtable(csv,...
     'ReadVariableNames',false);
end

if isempty(T)
    fprintf(' File %s is empty!\n',csv);
    ngene = 0;
    fig = [];
    return
end


 % summary(T);
Data = table2array(T(:,2:end));
agis = table2array(T(:,1));
    
if ~isempty(transfile)
   agis_new = f_tranlate(agis,transfile);
else
    agis_new = agis;
end

[ngene,~] = size(Data);

fprintf(' Calucating Expression Profiles at each time point...\n')
expr = [];
for i = 1:3:21%7 time points; 3 replicates;
    expr = [expr sum(Data(:,i:i+2),2)];
end
expr = 1/3*expr;

expr_bar = mean(expr,2);
expr_n_m = expr - repmat(expr_bar,1,7);
expr_n_mv = expr_n_m./repmat(sqrt(var(expr_n_m,0,2)),1,7);


fprintf(' Calucating log ratio expression pattern at each time point...\n')
tmp = [];
for i = 2:7
    tmp = [tmp log2( expr(:,i)./expr(:,1) )];
end
plotData = tmp;
plotData = [ zeros(size(plotData,1),1) plotData ];

if nargin == 2
    fig = [];
    return
end

%% Plot
fprintf(' Preparing for plotting...\n')
% =========  'Mean Plot' ===========
if any(strcmp(varargin,'Mean Plot'))
    counter = counter + 1;
    fprintf(' Calculating mean and visualizing...\n')

    %% Plot
    fig{counter} = figure;x = 0 : 1 : 6;
    hold on;axis([0 6 -2 5])
    plot(x, plotData','Color','[.4,.4,.4]');
    plot(x,mean( plotData),'Color','r','LineWidth',4);
    xticks(0:6)
    xticklabels({'0','0.25','0.5','1','4','12','24'})
    title(sprintf...
        ( 'plot of expression file\n "%s" with %d genes"',...
        csv, ngene))
    xlabel('Ethylene treatment(hrs)');
    ylabel('Expression-log2ratio(reference at 0 hrs)');
    set(gca,'fontsize',14);
    

    print(fig{counter},sprintf('%s/%s-Mean-Plot',dir_Figures,csv_token),'-dpng');
    fprintf(' Success!\n')
    fprintf('#####################################\n');

end


% =========  'Normalized' ===========
% Zero mean standard deviation
if any(strcmp(varargin,'Normalized'))
    counter = counter + 1;
    fprintf(' Calculating mean and visualizing...\n')

    %% Plot
    fig{counter} = figure;x = 0 : 1 : 6;
    hold on;axis([0 6 -3 3])
    plot(x, expr_n_mv');
    %plot(x,mean( expr_n_mv),'Color','r','LineWidth',4);
    xticks(0:6)
    xticklabels({'0','0.25','0.5','1','4','12','24'})
    title(sprintf...
        ( 'plot of expression file\n "%s" with %d genes"',...
        csv, ngene))
    xlabel('Ethylene treatment(hrs)');
    ylabel('Normalized');
    set(gca,'fontsize',14);
    grid on

    print(fig{counter},sprintf('%s/%s-Normalized',dir_Figures,csv_token),'-dpng');
    fprintf(' Success!\n')
    fprintf('#####################################\n');

end


% =========  'Normalized legend' ===========
% Zero mean standard deviation
if any(strcmp(varargin,'Normalized legend'))
    counter = counter + 1;
    fprintf(' Calculating mean and visualizing...\n')

    %% Plot
    fig{counter} = figure;x = 0 : 1 : 6;
    hold on;axis([0 6 -3 3])
    plot(x, expr_n_mv');
    %plot(x,mean( expr_n_mv),'Color','r','LineWidth',4);
    xticks(0:6)
    xticklabels({'0','0.25','0.5','1','4','12','24'})
    title(sprintf...
        ( 'plot of expression file\n "%s" with %d genes"',...
        csv, ngene))
    xlabel('Ethylene treatment(hrs)');
    ylabel('Normalized');
    set(gca,'fontsize',14);
    grid on

    legend(agis,'Location','best');
    
    print(fig{counter},sprintf('%s/%s-Normalized-legend',dir_Figures,csv_token),'-dpng');
    fprintf(' Success!\n')
    fprintf('#####################################\n');

end
% =========  'Network Analysis log' ===========
if any(strcmp(varargin,'Network Analysis log'))
    counter = counter + 1;
    fprintf(' Analyzing toy network and visualizing...\n')

    fig{counter} = figure; x = 0:6;
    hold on;axis([0 6 -2 5]);grid on;
    for i = 1 : ngene
       plot(x, plotData(i,:),'LineWidth',2);
    end
    legend(agis_new,'Location','best')
    xticks(0:6)
    xticklabels({'0','0.25','0.5','1','4','12','24'})
    title(sprintf( 'Network Analysis of expression file\n "%s"', csv))
    xlabel('Ethylene treatment(hrs)');
    ylabel('Expression-log2ratio(reference at 0 hrs)');
    set(gca,'fontsize',14);
    
    print(fig{counter},'./Figures/Network-Analysis-log','-dpng');

end

% =========  'Plot Individually log' ===========
if any(strcmp(varargin,'Plot Individually log'))
    counter = counter + 1;
    fprintf(' Individually plot log ratio figures...\n')
    
    fig_ind = cell(ngene);
    
    for i = 1 : ngene
        fig_ind{i} = figure; x = 0:6;
        hold on;axis([0 6 -2 5])
        plot(x, plotData(i,:),'LineWidth',2);
        legend(agis_new{i},'Location','best')
        xticks(0:6)
        xticklabels({'0','0.25','0.5','1','4','12','24'})
        title(sprintf( 'Individually analysis of gene\n "%s"', agis{i}))
        xlabel('Ethylene treatment(hrs)');
        ylabel('Expression-log2ratio(reference at 0 hrs)');
        set(gca,'fontsize',14);
        
        print(fig_ind{i},sprintf('./Figures/Network-Analysis-log%d',i),'-dpng');
    end
    fig{counter} = fig_ind;
end

% =========  'Plot Individually' ===========
if any(strcmp(varargin,'Plot Individually'))
    counter = counter + 1;
    fprintf(' Individually plot figures...\n')
    
    fig_ind = cell(ngene);
    for i = 1 : ngene
        fig_ind{i} = figure; x = 0:6;
        hold on;grid on
        plot(x, expr(i,:),'LineWidth',2);
        for j = x
            plot(j*ones(1,3), Data(i,3*j+1:3*j+3),'k*')
        end
        legend(agis_new{i},'Location','best')
        xticks(0:6)
        xticklabels({'0','0.25','0.5','1','4','12','24'})
        title(sprintf( 'Individually analysis of gene\n "%s"', agis{i}))
        xlabel('Ethylene treatment(hrs)');
        ylabel('Expression level(PKRM)');
        set(gca,'fontsize',14);
        print(fig_ind{i},sprintf('./Figures/Network-Analysis%d',i),'-dpng');
    end
    fig{counter} = fig_ind;
    fprintf('DONE!\n')
end

return
end


function [agis_new] = f_tranlate(agis,transfile)
% Translation file has to have Variable Names.
%       Example of Translation file opened in MS Excell might look as follows:
%       |    ORF    | OtherNames |
%       | AT1G66340 |    ETR1    |
%       |    ...    |    ...     |
%       | AT3G20770 |    EIN3    |
%       |    ...    |    ...     |
    T_ORF2CNames = readtable(transfile,...
     'ReadVariableNames',true,'ReadRowNames',true);

    for i = 1 : length(agis)
        agis{i} = agis{i}(1:9);
    end
 
    agis_new = T_ORF2CNames(agis,:).OtherNames;

    return
end

