% This script probe if there is possibly any leakage in SPARTAN cartridges.
% We want to look at the most recent two year's PM2.5 mass and Blank
% filter mass and see if there is any correlation at any site. 

% ---- Set directories ----
addpath('/storage1/fs1/rvmartin/Active/haihuizhu/6.SPARTAN')
direc = '/storage1/fs1/rvmartin/Active/SPARTAN-shared'; % For Haihui 
% addpath('/Volumes/rvmartin/Active/haihuizhu/6.SPARTAN')
% direc= '/Volumes/rvmartin/Active/SPARTAN-shared'; % For Haihui 

direc_master = strcat(direc,'/Analysis_Data/Master_files');
direc_blank = strcat(direc,'/Analysis_Data/blank');

% ----- user input ------
NoBadFilter = 1; 
for tyear = 999 % set to 999 and mannualy change line 61 to makes figures for all data after 2020
%-------------   SITE DETAILS   --------------
site_details = readtable(strcat(direc,'/Site_Sampling/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));

%-------------   Sampling Methods & Flag   --------------
Sampling_Parameters_Methods = strcat(direc,'/SOPs/Public SOPs/Sampling_Parameters_Methods_2.3.xlsx');
Flags = readtable(Sampling_Parameters_Methods,'Sheet','Flags');
Flags = Flags.Flag;

r1 = nan(length(Site_cities),1);
r2 = nan(length(Site_cities),1);
n1 = nan(length(Site_cities),1);
n2 = nan(length(Site_cities),1);
siteid = [];
site = cell(0,0);
cartridge = cell(0,0);
filter1 = [];
filter2 = [];
filter3 = [];
filter4 = [];
filter5 = [];
filter6 = [];
blank_mass = [];
blank2total = [];
odd2even = [];
avg2to6 = [];

sheetname = 'all_filters'; 
fnamenote = sprintf('_%d', tyear);

for loc = 1:length(Site_codes)

    % Find if there is a master file
    master_file = sprintf('%s/%s_master.csv',direc_master,Site_codes{loc});
    [Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, projectIDs_master,Master_hours, Master_masstype, ...
        Master_dates, Master_mass, Master_IC, Master_ICP, Master_XRF,...
        Master_carbon, Master_Method, Master_flags] = ReadMaster(master_file,Site_codes{loc});

    uniquecart = unique(CartridgeIDs_master);

    if isempty(uniquecart)
        continue
    else
        yearid = find(Master_dates(:,5)> 2020); % change here
        if isempty(yearid)
            fprintf('%s has no data for %d\n',Site_cities{loc},tyear)
        else
            uniquecart = unique(CartridgeIDs_master(yearid));

            % 0 = blank
            % 1 = PM2.5
            % 2 = PM10
            % 3 = PMcoarse (nuclepore)
            % 4 = unknown/void, nuclepore filter saturated
            % 5 = negative mass
            % 6 = invalid flow rates
            
            pm25mass = nan(length(uniquecart),6);
            blankmass = nan(size(uniquecart));
            for cid = 1:length(uniquecart)

                tcarind = find(ismember(CartridgeIDs_master,uniquecart{cid}));
                if length(tcarind) == 8
                    tmasstype = Master_masstype(tcarind);
                    tmass = Master_mass(tcarind,:);
                    tflag = Master_flags(tcarind);
                    
                    PM25ind = find(tmasstype==1);
                    pm25mass(cid,PM25ind) = tmass(PM25ind)';
%                     flags = contains(tflag,Flags); % flag it if filter has any of the conditions
                    flags = zeros(size(tflag)); % flag it if filter has any of the conditions
                    flags(isnan(tmass(:,2))) = 1; % flag it if flow rate is nan
                    
                    if NoBadFilter == 1
                        pm25mass(cid,flags==1) = nan;
                        sheetname = 'no_bad_filter';
                        fnamenote = sprintf('_%d_NoBadFilter',tyear);
                    end

                    % find blank mass
                    ind = find(tmasstype==0 & flags==0);
                    if ~isempty(ind) &&  sum(isnan(pm25mass(cid,:))) < 6
                        blankmass(cid) = tmass(ind,1);
                    end

                else
                    warning('%s has %d filters\n',uniquecart{cid},length(tcarind) == 8)
                end
            end
 
            [r1(loc),n1(loc)] = makescatter(mean(pm25mass,2,'omitnan'),blankmass,Site_codes{loc});
            % save figure
            sfname = sprintf('%s/%s_blank_vs_PM2.5_Scatter%s.png',direc_blank,Site_codes{loc},fnamenote);
            saveas(gcf,sfname)
            fprintf('%s saved\n',sfname)
            close

            % add data point to table
            ind = find(~isnan(blankmass));
            if ~isempty(ind)
                siteid(end+1:end+length(ind),1) = loc;
                site(end+1:end+length(ind),1) = Site_cities(loc);
                cartridge(end+1:end+length(ind),1) = uniquecart(ind);
                filter1(end+1:end+length(ind),1) = pm25mass(ind,1);
                filter2(end+1:end+length(ind),1) = pm25mass(ind,2);
                filter3(end+1:end+length(ind),1) = pm25mass(ind,3);
                filter4(end+1:end+length(ind),1) = pm25mass(ind,4);
                filter5(end+1:end+length(ind),1) = pm25mass(ind,5);
                filter6(end+1:end+length(ind),1) = pm25mass(ind,6);
                blank_mass(end+1:end+length(ind),1) = blankmass(ind);

                % ratio between blank and total
                blank2total(end+1:end+length(ind),1) = blankmass(ind)./sum(pm25mass(ind,:),2,'omitnan');
                odd2even(end+1:end+length(ind),1) = sum(pm25mass(ind,[1 3 5]),2,'omitnan') ./ sum(pm25mass(ind,[2 4 6]),2,'omitnan') ;
                avg2to6(end+1:end+length(ind),1) = mean(pm25mass(ind,2:6),2,'omitnan');
            end
        end
    end
end

% bar_chart (r1,r2,n1,n2,Site_cities)
% % save figure
% sfname = sprintf('%s/Bar_Correlations_blank_vs_PM2.5%s.png',direc_blank,fnamenote);
% saveas(gcf,sfname)
% fprintf('%s saved\n',sfname)
% close 

% print to spreadsheet
t = table(site,cartridge,blank_mass,filter1,filter2,filter3,filter4,filter5,filter6,blank2total,odd2even,avg2to6);
outfname = sprintf('%s/Sheet_blank_vs_PM2.5_Mass_%s.xlsx',direc_blank,num2str(tyear));
writetable(t,outfname,'Sheet',sheetname)

% making bar + std of those matrice
% blank to toal mass
thetitle = 'Blank to Other Filters';
bar_std(siteid,blank2total,Site_cities,thetitle)
sfname = sprintf('%s/Bar_Blank2TotMass%s.png',direc_blank,fnamenote);
saveas(gcf,sfname)
fprintf('%s saved\n',sfname)
close 
% odd to even fiilters ratio
thetitle = 'Odd to Even Filters';
bar_std(siteid,odd2even,Site_cities,thetitle)
sfname = sprintf('%s/Bar_Odd2Even%s.png',direc_blank,fnamenote);
saveas(gcf,sfname)
fprintf('%s saved\n',sfname)
close
% filter 1 vs others
bar_std_double(siteid,filter1,avg2to6,Site_cities,direc_blank,fnamenote)
end

%% Function
function bar_std(siteid,blank2total,Site_cities,thetitle)
    numsites = length(Site_cities);
    avgs = nan(numsites,1);
    medians = nan(numsites,1);
    stds = nan(numsites,1);
    N = nan(numsites,1);
    blank2total(blank2total== Inf) = nan;
    for id = 1:numsites
        N(id,1) = length(find(siteid == id));
        avgs(id,1)  = mean(blank2total(siteid == id),'omitnan');
        stds(id,1)  = std(blank2total(siteid == id),'omitnan')./sqrt(N(id,1));
        medians(id,1)  = median(blank2total(siteid == id),'omitnan');
    end

    % screen out nan
    ind = find(isnan(avgs));
    avgs(ind) = [];
    stds(ind) = [];
    medians(ind) = [];
    N(ind) = [];
    Site_cities(ind) = [];

    figure('Position',[10 10 800 400])
    x = 1:length(avgs);
    bar(x,avgs)                
    hold on

    er = errorbar(x,avgs,stds',stds');
%     er.Color = [0 0 0];
    er.LineStyle = 'none';
    hold on
    
    % adding median
    plot(x,medians,'.','MarkerSize',20,'Color','#EDB120')
    hold on

    xl = 0:length(avgs)+1;
    avgall = mean(avgs);
    medall = median(medians);
    plot(xl,avgall.*ones(size(xl)),'--r','linewidth',1)
    str = sprintf('SPARTAN Avg = %6.4f\nSPARTAN Median = %6.4f\n',avgall,medall);
    text(0.05,0.95,str,'units','normalized','horiz','left','verti','top','fontsize',16)
    hold on

    % plot the R = 1 line
    if avgall > 1
        plot(xl, ones(size(xl)),':k','linewidth',1.5)
        hold on
    end

    % lable N
    ylim2 = ylim;
    for id = x
        text(id,0,num2str(N(id)),'horiz','center','verti','top','fontsize',10)
        hold on
    end

    xticks(x)
    xticklabels(Site_cities)
    set(gca,'fontsize',16)
%     title(thetitle)
    ylabel('Ratio (unitless)')
    legend({'mean','std error','median','SPARTAN mean','1:1'},'Location','northeast')
    legend('boxoff')
    hold off
end

function bar_std_double(siteid,filter1,avg2to6,Site_cities,direc_blank,fnamenote)
    numsites = length(Site_cities);
    
    avg1 = nan(numsites,1);
    std1 = nan(numsites,1);
    med1 = nan(numsites,1);
    avg2 = nan(numsites,1);
    std2 = nan(numsites,1);
    med2 = nan(numsites,1);
    ratios = nan(numsites,1);
    stdr = nan(numsites,1);
    medr = nan(numsites,1);
    N = nan(numsites,1);

    filter1(filter1== Inf) = nan;
    avg2to6(avg2to6== Inf) = nan;

    for id = 1:numsites
        N(id,1) = length(find(siteid == id));

        avg1(id,1)  = mean(filter1(siteid == id),'omitnan');
        std1(id,1)  = std(filter1(siteid == id),'omitnan')./ sqrt(N(id,1));
        med1(id,1) = median(filter1(siteid == id),'omitnan');
        avg2(id,1)  = mean(avg2to6(siteid == id),'omitnan');
        std2(id,1)  = std(avg2to6(siteid == id),'omitnan')./ sqrt(N(id,1));
        med2(id,1) = median(filter1(siteid == id),'omitnan');

        ratios(id,1) = mean(filter1(siteid == id)./avg2to6(siteid == id),'omitnan');
        stdr(id,1) = std(filter1(siteid == id)./avg2to6(siteid == id),'omitnan') ./ sqrt(N(id,1));
        medr(id,1) = median(filter1(siteid == id)./avg2to6(siteid == id),'omitnan');

    end

    % screen out nan
    ind = find(isnan(avg1) | isnan(avg2));
    avgs  = [avg1,avg2];
    avgs(ind,:) = [];
    std1(ind) = [];
    std2(ind) = [];
    med1(ind) = [];
    med2(ind) = [];
    N(ind) = [];
    Site_cities(ind) = [];
    ratios(ind) = [];
    stdr(ind) = [];
    medr(ind) = [];

    % ---- plotting filter 1 and other filters ----
    thetitle = 'Filter 1 vs. Others';
    figure('Position',[10 10 800 400])
    x = 1:length(avgs);
    bar(x,avgs)                
    hold on

    xe = x-0.125;
    er = errorbar(xe,avgs(:,1),std1',std1');
%     er.Color = [0 0 0];
    er.LineStyle = 'none';
    hold on
    % adding median
    plot(xe,med1,'.','MarkerSize',14,'Color','#EDB120')
    hold on

    hold on
    xe = x+0.125;
    er = errorbar(xe,avgs(:,2),std2',std2');
%     er.Color = [0 0 0];
    er.LineStyle = 'none';
    hold on
    % adding median
    plot(xe,med2,'.','MarkerSize',14,'Color','#EDB120')
    hold on
   
    % lable N
    for id = x
        text(id,-0.05,num2str(N(id)),'horiz','center','verti','top','fontsize',10)
        hold on
    end

    xticks(x)
    xticklabels(Site_cities)
    ylabel(sprintf('Mass (%s)','\mug'))
    set(gca,'fontsize',16)
%     title(thetitle)
    legend({'Filter 1','Filter 2 to 6','std error','median'},'location','northeast')
    legend('boxoff') 

    hold off

    sfname = sprintf('%s/Bar_Filter1_and_Others%s.png',direc_blank,fnamenote);
    saveas(gcf,sfname)
    fprintf('%s saved\n',sfname)
    close

    % ---- plotting the ratio ----
    figure('Position',[10 10 800 400])
    thetitle = 'Filter 1 - Others Ratio';
    bar(x,ratios)
    hold on

    er = errorbar(x,ratios,stdr',stdr');
%     er.Color = [0 0 0];
    er.LineStyle = 'none';
    hold on
    % adding median
    plot(x,medr,'.','MarkerSize',20,'Color','#EDB120')
    hold on

    % plot the network avg
    avgall = mean(ratios);
    medall = median(medr);
    plot(x,avgall.*ones(size(x)),'--r','linewidth',1.5)
    str = sprintf('SPARTAN Avg = %4.2f\nSPARTAN Median = %4.2f\n',avgall,medall);
    text(0.05,0.95,str,'units','normalized','horiz','left','verti','top','fontsize',16)
    hold on

    % plot the R = 1 line
    plot(x, ones(size(x)),':k')
    hold on 

    % lable N
    for id = x
        text(id,-0.01,num2str(N(id)),'horiz','center','verti','top','fontsize',10)
        hold on
    end

    xticks(x)
    xticklabels(Site_cities)
    set(gca,'fontsize',16)
%     title(thetitle)
    ylabel('Ratio (unitless)')
    legend({'mean','std error','median','SPARTAN mean','1:1'},'Location','northeast')
    legend('boxoff')
    hold off

    sfname = sprintf('%s/Bar_Filter1_TO_Others%s.png',direc_blank,fnamenote);
    saveas(gcf,sfname)
    fprintf('%s saved\n',sfname)
    close

end

function [R1,  N1] = makescatter(pm25mass,blankmass,city)

    fz = 12;
    dotsize = 10;
    figure('Position',[10 10 800 400])

% subplot 1
    subplot(1,2,1)
    scatter(pm25mass,blankmass,dotsize,"black","filled")
    hold on
    % adjust x and/or y limite if needed
    xlim2 = xlim;
    ylim2 = ylim;

    xlabel(sprintf('%s mass (%s)', 'PM_{2.5}','\mug/m^3'), 'fontsize', fz, 'fontweight', 'bold');
    ylabel(sprintf('%s mass (%s)', 'Blank filter','\mug/m^3'), 'fontsize', fz, 'fontweight', 'bold');
    hold on

    % adding trend and p value
    [str1, m, b, R1, p, N1] = getStatics(pm25mass, blankmass);

    xl = 0:1/10 * (xlim2(2) - xlim2(1)):xlim2(2);
    yl = m * xl + b;
    plot(xl, yl, '--k', 'Linewidth', 1)
    hold on
    text(xlim2(2) - 0.05 * (xlim2(2) - xlim2(1)), ylim2(1) + 0.97 * (ylim2(2) - ylim2(1)), str1, ...
        'horizontalalignment', 'right', 'verticalalignment', 'top', 'color', 'k', 'fontsize', fz);
    text(0.05 * (xlim2(2) - xlim2(1)), ylim2(1) + 0.97 * (ylim2(2) - ylim2(1)), city, ...
        'horizontalalignment', 'left', 'verticalalignment', 'top', 'color', 'k', 'fontsize', fz+2);
  

end


function [ta1, m, b, R, p, N] = getStatics(X, Y)

    nonnan = find(~isnan(X) & ~isnan(Y));
    X = X(nonnan);
    Y = Y(nonnan);

    if length(Y) < 3
        fprintf('No enough input data for scatter\n')
        %     tN = sprintf('N = %d',numel(X));
        ta1 = 'NaN';
        m = NaN;
        b = NaN;
        R = NaN;
        p = NaN;
        N = NaN;
    else
        % calc statistics
        [R, P] = corrcoef(X, Y, 'rows', 'complete');
        R = R(1, 2); % R2
        p = P(1, 2); % p-value
        b2 = regress(Y, [X ones(length(Y), 1)]);
        m = b2(1);
        b = b2(2);
        pm = '+-';

        N = sum(length(X));

        tslope = sprintf('y = %.4fx%s%.4f', m, pm((b < 0) + 1), abs(b));
        tr = ['r = ' sprintf('%.2f', R)];
        tp = ['p = ' sprintf('%.3f', p)];
        tN = ['N = ' sprintf('%d', N)];

        ta1 = sprintf('%s\n%s\n%s\n%s', tslope, tr, tp, tN);
    end

end

function bar_chart (R1, R2,n1,n2, Site_cities)
    fz = 8;
    load('comp_color.mat')
    Colors = {sulfate, ss};

    % remove sites with no data
    ind = find(~isnan(R1));
    R1 = R1(ind);
    R2 = R2(ind);
    n1 = n1(ind);
    n2 = n2(ind);
    Site_cities = Site_cities(ind);

    [R1,I] = sort(R1);
    R2 = R2(I);
    Site_cities = Site_cities(I);
    n1 = n1(I);
    n2 = n2(I);

    allR = [R1, R2];

    figure('Position',[10 10 900 400])
    b = bar(allR);
    for ss = 1:size(allR,2)
        b(ss).FaceColor = Colors{ss};
    end
    ylim([-1 1.2])

    xticks(1:length(Site_cities))
    xticklabels(Site_cities)
    ylabel('Correlation')

    legend ({'All','No bad filters'}, 'location', 'best')
    set(gca,'fontsize',fz+6)

    % adding # of data points for each site
    for i = 1:length(R2)
        yloc = max(allR(i,:));
        str = sprintf('N1: %d\nN2: %d',n1(i),n2(i));
        if yloc <= 0 
            text(i,0.08,str,'fontsize',fz,'Color','k','HorizontalAlignment','center') 
        else
            text(i,yloc+0.08,str,'fontsize',fz,'Color','k','HorizontalAlignment','center') 
        end
        hold on
    end


end