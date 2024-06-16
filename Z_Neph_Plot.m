% This script making plots comparing total scatter under the shared dir and under my own dir
% Haihui 
% 2023-Nov
clear 
close all

S_AllData = 0;
S_StatisFig = 1;

% MyDir  = '/storage1/fs1/rvmartin/Active/haihuizhu/6.SPARTAN'; 
ShareDir = '/storage1/fs1/rvmartin/Active/SPARTAN-shared';

TrackerFileName = sprintf('%s/Public_Data/Neph_Processed/File_condition_tracker.xlsx', ShareDir);

site_details = readtable(sprintf('%s/Site_Sampling/Site_details.xlsx', ShareDir), 'PreserveVariableNames', true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));

if S_AllData == 1

for loc =  1:length(Site_codes)
    fprintf('Plotting for %s\n', Site_cities{loc})
    SiteCode = Site_codes{loc};

    % Scatter Data
    % read data
    fname1 = sprintf('%s/Public_Data/Neph_Processed_20231218/%s/%s_%s_NephPM25_daily.csv', ShareDir, SiteCode, SiteCode, Site_cities{loc});
    fname2 = sprintf('%s/Public_Data/Neph_Processed/%s/%s_%s_NephPM25_daily.csv', ShareDir, SiteCode, SiteCode, Site_cities{loc});
    fname3 = sprintf('%s/Public_Data/Time-resolved_PM2.5_20231218/Data_sharing/%s_%s_DailyEstimatedPM25.csv', ShareDir, SiteCode, Site_cities{loc});
    fname4 = sprintf('%s/Public_Data/Time-resolved_PM2.5/Data_sharing/%s_%s_DailyEstimatedPM25.csv', ShareDir, SiteCode, Site_cities{loc});
    if exist(fname1, 'file') && exist(fname3, 'file')
        databefore = readtable(fname1, 'VariableNamingRule', 'preserve');
        dataafter = readtable(fname2, 'VariableNamingRule', 'preserve');

        % make the figure
        date1 = datetime([databefore.Year, databefore.Month, databefore.Day]);
        date2 = datetime([dataafter.Year, dataafter.Month, dataafter.Day]);
        figure('Position', [10 10 1000 600])
        [N1, N2, r] = makescatplot(date1, date2, databefore.Total_Scat_Red, dataafter.Total_Scat_Red, 1, 'Red Scatter');
        [Scatter_N_before, Scatter_N_after, Scatter_r] = ...
        makescatplot(date1, date2, databefore.Total_Scat_Green, dataafter.Total_Scat_Green, 2, 'Green Scatter');
        [N1, N2, r] = makescatplot(date1, date2, databefore.Total_Scat_Blue, dataafter.Total_Scat_Blue, 3, 'Blue Scatter');

        clear databefore dataafter

        % Time-resolved PM2.5
        databefore = readtable(fname3, 'VariableNamingRule', 'preserve');
        dataafter = readtable(fname4, 'VariableNamingRule', 'preserve');

        date1 = datetime([databefore.Year_local, databefore.Month_local, databefore.Day_local]);
        date2 = datetime([dataafter.Year_local, dataafter.Month_local, dataafter.Day_local]);

        [TimeResovled_N_before, TimeResovled_N_after,TR_r1] = ...
        makescatplot(date1, date2, databefore.Value, dataafter.Value, 4, 'Time-Resolved PM_{2.5}');

        % save figure
        sfname = sprintf('%s/Public_Data/Neph_Processed/Plots/DataVersion1.1_PM25_TotalScatter_Time-resolvedPM25_%s.png', ShareDir, SiteCode);
        saveas(gcf, sfname)
        fprintf('%s saved\n', sfname)
        close

        % write the statistics to the table
        SiteCode = Site_codes(loc);
        T = readtable(TrackerFileName, 'Sheet', 'Statis_TimeResolvedPM25');
        if ~isempty(T)
            ind = find(ismember(T.SiteCode, SiteCode));
            T(ind, :) = [];
            delete_sheet(TrackerFileName, 'Statis_TimeResolvedPM25')
            writetable(T, TrackerFileName, 'Sheet', 'Statis_TimeResolvedPM25')
        end
        T = table(SiteCode,Scatter_N_before, Scatter_N_after, Scatter_r, TimeResovled_N_before, TimeResovled_N_after, TR_r1);
        writetable(T, TrackerFileName, 'Sheet', 'Statis_TimeResolvedPM25', 'WriteMode', 'Append')
    end
end
end

if S_StatisFig == 1
Colors = {'#0072BD', '#D95319', '#77AC30'};
% making bars comparing N and R
fz = 12; % font size
t = readtable(TrackerFileName, 'Sheet', 'Statis_TimeResolvedPM25');
figure('Position',[10 10 1200 400])

subplot(1,2,1)
yyaxis left
b = bar([t.Scatter_N_before, t.Scatter_N_after]);
for i = length(b)
    b(i).FaceColor = Colors{i};
end
xticks(1:length(t.SiteCode))
xticklabels(t.SiteCode)
ylabel('Number')
hold on
yyaxis right
plot(1:height(t), t.Scatter_r, '^', 'DisplayName', 'Ratio - V1 vs V2');
ylabel('Ratio')
ylim([0.4 1.4])
set(gca,'fontsize',fz)
legend({'V1 #', 'V2 #', 'Ratio - V1 vs V2'}, 'location', 'southoutside')
hold off
title('Green Scatter')

subplot(1,2,2)
yyaxis left
d = bar([t.TimeResovled_N_before, t.TimeResovled_N_after]);
for i = length(d)
    d(i).FaceColor = Colors{i};
end
xticks(1:length(t.SiteCode))
xticklabels(t.SiteCode)
ylabel('Number')
hold on
yyaxis right
plot(1:height(t), t.TR_r1, 'k^', 'DisplayName', 'Ratio - V1 vs V2');
% hold on 
% plot(1:height(t), t.TR_r2, 'k*','DisplayName','Ratio - 532 nm vs 550 nm');
ylabel('Ratio')
ylim([0.4 1.4])
set(gca, 'fontsize', fz)
legend({'V1 #', 'V2 #', 'Ratio - V1 vs V2'}, 'location', 'southoutside')
% legend({'V1 #', 'V2 #', '550 nm #', 'Ratio - V1 vs V2', 'Ratio - 532 nm vs 550 nm'}, 'location', 'southoutside')
title('Time Resolved PM2.5')

sfname = sprintf('%s/Public_Data/Neph_Processed/Plots/Statistics_all.png', ShareDir);
saveas(gcf, sfname)
fprintf('Fig saved\n')
end

function [N1, N2, ratios] = makescatplot(date1, date2, databefore, dataafter, ps, tit)
    % Titles_daily = {'Year' 'Month' 'Day' 'num_hourly_points' 'Temp_degC' 'Pressure' '%RH' 'Back_Scat_Red' 'Back_Scat_Green' ...
    %                 'Back_Scat_Blue' 'Total_Scat_Red' 'Total_Scat_Green' 'Total_Scat_Blue'};
    ind = find(databefore < 0 );
    if ~isempty(ind)
        fprintf('Negative %s found in version 1.0. Removed the following data from plotting:\n',tit)
        disp(databefore(ind))
        databefore(ind) = nan;
    end
    ind = find(dataafter < 0);
    if ~isempty(ind)
        fprintf('Negative %s found in version 2.0. Removed the following data from plotting:\n', tit)
        disp(dataafter(ind))
        dataafter(ind);
    end

    Colors = {'r', 'g', 'b', '#76ab93', 'm', 'c'};
    subplot(2, 2, ps)
    plot(date1, databefore, '.', 'color', Colors{ps}, 'MarkerSize', 15);
    hold on
    plot(date2, dataafter, '.', 'color', 'k', 'MarkerSize', 10);
    hold on
    legend({'V1.0', 'V2.0'}, 'Location', 'NorthWest', 'FontWeight', 'bold', 'FontSize', 10)

    % adding N and R
    N1 = sum(~isnan(databefore));
    N2 = sum(~isnan(dataafter));
    alldates = [min([date1; date2]):max([date1; date2])]';
    X = nan(size(alldates)); Y = nan(size(alldates));
    [~, ind] = ismember(date1, alldates);
    X(ind) = databefore; clear ind
    [~, ind] = ismember(date2, alldates);
    Y(ind) = dataafter; clear ind
    ind = find(~isnan(X) & ~isnan(Y));
    r = corrcoef(X(ind), Y(ind), 'rows', 'complete');
    r = r(1,2);
    str = sprintf('N_{V1.0} = %d  N_{V2.0} = %d  r = %4.2f',N1,N2,r);
    text(0.95, 0.95, str, 'Units', 'normalized', 'hori', 'right', 'verti', 'top');

    % cal ratio
    ratios = mean(Y./X,'omitnan');


    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    xtickangle(45)
    ylabel('Signal', 'FontWeight', 'bold', 'FontSize', 10);
    title(tit)


end

function [N1, N2, N3, ratio1, ratio2] = makescatplot2(date1, date2, date3, databefore, dataafter, dataafter2, ps, tit)
    % Titles_daily = {'Year' 'Month' 'Day' 'num_hourly_points' 'Temp_degC' 'Pressure' '%RH' 'Back_Scat_Red' 'Back_Scat_Green' ...
    %                 'Back_Scat_Blue' 'Total_Scat_Red' 'Total_Scat_Green' 'Total_Scat_Blue'};
    ind = find(databefore < 0);

    if ~isempty(ind)
        fprintf('Negative %s found in version 1.0. Removed the following data from plotting:\n', tit)
        disp(databefore(ind))
        databefore(ind) = nan;
    end

    ind = find(dataafter < 0);

    if ~isempty(ind)
        fprintf('Negative %s found in version 2.0. Removed the following data from plotting:\n', tit)
        disp(dataafter(ind))
        dataafter(ind);
    end

    ind = find(dataafter2 < 0);

    if ~isempty(ind)
        fprintf('Negative %s found in version 550. Removed the following data from plotting:\n', tit)
        disp(dataafter(ind))
        dataafter(ind);
    end

    Colors = {'r', 'g', 'b', '#76ab93', 'm', 'c'};
    subplot(2, 2, ps)
    plot(date1, databefore, '.', 'color', Colors{ps}, 'MarkerSize', 20);
    hold on
    plot(date2, dataafter, '.', 'color', 'k', 'MarkerSize', 15);
    hold on
    plot(date3, dataafter2, '.', 'color', 'r', 'MarkerSize', 10);
    hold on
    legend({'V1.0', 'V2.0','550nm'}, 'Location', 'NorthWest', 'FontWeight', 'bold', 'FontSize', 10)

    % adding N and R
    N1 = sum(~isnan(databefore));
    N2 = sum(~isnan(dataafter));
    N3 = sum(~isnan(dataafter2));
    alldates = [min([date1; date2;date3]):max([date1; date2; date3])]';
    X = nan(size(alldates)); Y = nan(size(alldates)); Y2 = nan(size(alldates));
    [~, ind] = ismember(date1, alldates);
    X(ind) = databefore; clear ind
    [~, ind] = ismember(date2, alldates);
    Y(ind) = dataafter; clear ind
    [~, ind] = ismember(date3, alldates);
    Y2(ind) = dataafter2; clear ind
    ind = find(~isnan(X) & ~isnan(Y));
    r = corrcoef(X(ind), Y(ind), 'rows', 'complete');
    r = r(1, 2); clear ind
    ind = find(~isnan(Y) & ~isnan(Y2));
    r2 = corrcoef(Y(ind), Y2(ind), 'rows', 'complete');
    r2 = r2(1, 2); clear ind
    str = sprintf('N_{V1.0} = %d  N_{V2.0} = %d  r1 = %4.2f N_{550} = %d  r2 = %4.2f', N1, N2, r, N3,r2);
    text(0.95, 0.95, str, 'Units', 'normalized', 'hori', 'right', 'verti', 'top');
    ratio1 = mean(Y ./ X, 'omitnan');
    ratio2 = mean(Y2 ./ Y, 'omitnan');

    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    xtickangle(45)
    ylabel('Signal', 'FontWeight', 'bold', 'FontSize', 10);
    title(tit)

end
