% This script compares IC measurements before and after the low
% concentration curve is implemented
% Written by Haihui Zhu, Oct, 2023

close all; clear; clc
addpath('UtilityFunctions')
addpath('/storage1/fs1/rvmartin/Active/haihuizhu/1.code') % for statistc function

% Setup directories 
direc = find_root_dir(1);

direc_sampling = strcat(direc,'/Site_Sampling'); 
site_details = readtable(strcat(direc_sampling,'/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));

% IC ions
species = {'Na','K','NH4','Ca','Mg','Li','Br','Cl','F','NO2','NO3','PO4','SO4'};
fontsize = 10;


% path to the OLD master files 
direc_master1 = strcat(direc,'/Analysis_Data/Master_files/Backups/Backup_2023-10-09');
% path to the NEW master files 
direc_master2 = strcat(direc,'/Analysis_Data/Master_files');
% path to save figures
direc_Figure = strcat(direc,'/Analysis_Data/Master_files/Figure_IC_v4-vs-v5');

for loc = 1:length(Site_codes)
    % read the OLD file
    master_file = sprintf('%s/%s_master.csv',direc_master1,Site_codes{loc});
    [Titles,Filter_ID, ~, ~, ~, ~,~, ~, ...
        ~, ~, Master_IC, ~, ~,...
        ~, ~, ~] = ReadMaster(master_file,Site_codes{loc});

    % read the NEW file
    master_file = sprintf('%s/%s_master.csv',direc_master2,Site_codes{loc});
    [~,~, ~, ~, ~, ~,~, ~, ...
        ~, ~, Master_IC2, ~, ~,...
        ~, ~, ~] = ReadMaster(master_file,Site_codes{loc});

    IC_titles = Titles(find(contains(Titles,'IC_')));
    x = 1:size(Master_IC2,1);

    figure('Position',[1 1 1000 1000])
    for sp = 1:length(species)

        ind = find(contains(IC_titles,species{sp}));
        x = Master_IC(:,ind);
        y = Master_IC2(:,ind);
        if ~isempty(x)
            rng = [-0.1*max(round(max(x)*1.3),0.3) max(round(max(x)*1.3),0.3) ];

            subplot(4,4,sp)
            scatter(x,y,10,'filled')
            xlabel('original')
            ylabel('updated')
            title (sprintf('%s (%s)',species{sp},'ug'))
            hold on

            % formatting
            set(gca, ...
                'Box'         , 'on'     , ...
                'TickDir'     , 'in'     , ...
                'TickLength'  , [.02 .02] , ...
                'XMinorTick'  , 'on'      , ...
                'YMinorTick'  , 'on'      , ...
                'YGrid'       , 'off'     , ...
                'XColor'      , 'k', ...
                'YColor'      , 'k', ...
                'Xlim'        , rng, ...
                'Ylim'        , rng, ...
                'YTick'       , [rng(1), 0:(rng(2)/4):rng(2)], ...
                'XTick'       , [rng(1), 0:(rng(2)/4):rng(2)], ...
                'LineWidth'   , 1         );
            pbaspect([1 1 1])

            % adding statistics
            % calc statistics
            Ind = find(~isnan(x)& ~isnan(y));
            Xt=x(Ind);
            Yt=y(Ind);
            [ta1,m,b] = getStatics(Xt,Yt);

            xl = -0.3:0.01*rng(2):rng(2);
            yl = m*xl + b;
            plot(xl,yl,'-r','Linewidth',0.8)
            hold on
            plot(xl,zeros(size(yl)),'--k','Linewidth',0.8)
            hold on

            text(rng(1)+0.05*(rng(2)-rng(1)),rng(1)+0.97*(rng(2)-rng(1)),ta1, ...
                'horizontalalignment','left','verticalalignment','top','color','k','fontsize',fontsize+2);
        end
    end
    figname = sprintf('%s/%s_revison5_vs_revision4.png',direc_Figure,Site_codes{loc});
    saveas(gcf,figname)
   

    % reprocessing time series
    figure('Position',[1 1 1000 1000])
    for sp = 1:length(species)

        ind = find(contains(IC_titles,species{sp}));
        x = Master_IC(:,ind);
        y = Master_IC2(:,ind);
        if ~isempty(x)
            subplot(7,2,sp)
            plot(1:length(x),x, '.','linewidth',5,'displayname','original')
            hold on
            plot(1:length(y),y, '.','linewidth',5,'displayname','updated')
            xticks(1:10:length(x))
            xlabel('filiters')
            ylabel(sprintf('%s mass (%s)',species{sp},'ug'))
            legend
            hold on
        end
    end
    figname = sprintf('%s/%s_revison5_vs_revision4_timeseries.png',direc_Figure,Site_codes{loc});
    saveas(gcf,figname)
   
end

