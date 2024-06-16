    

function [] = PM25filter_timeseries(PM25_data, PM25_dates, Site_cities,Site_codes, direc_output,loc)

Start_Date = datenum(PM25_dates(1,1),PM25_dates(1,2),PM25_dates(1,3));
End_Date   = datenum(PM25_dates(end,5),PM25_dates(end,6),PM25_dates(end,7));
Plotting_Dates = [Start_Date:1:End_Date];
PM_Filter_Start = datenum(PM25_dates(:,1),PM25_dates(:,2),PM25_dates(:,3));
PM_Filter_End   = datenum(PM25_dates(:,5),PM25_dates(:,6),PM25_dates(:,7));

PM_Plotting = horzcat(Plotting_Dates',NaN(size(Plotting_Dates))');

for fil = 1:size(PM25_data,1)
    strt_idx = find(Plotting_Dates == PM_Filter_Start(fil));
    end_idx  = find(Plotting_Dates == PM_Filter_End(fil));
    PM_Plotting(strt_idx:end_idx,2) = PM25_data(fil,1);
end

Site_PM_Avg = mean(PM_Plotting(:,2),1,'omitnan');
Avg_Line = repmat(Site_PM_Avg,1,size(PM_Plotting,1));
Day_num = size(PM_Plotting,1);

figure(1)
scatter(PM_Plotting(:,1),PM_Plotting(:,2),'s','filled','LineWidth',1)
if 1 == 1
    if (Day_num <= 1095)
        tick_pos = [PM_Plotting(1,1):31:PM_Plotting(end,1)];
    elseif  (Day_num> 1095)
        tick_pos = [PM_Plotting(1,1):62:PM_Plotting(end,1)];
    end
end
xticks(tick_pos)
datetick('x','mmm-yyyy','keeplimits','keepticks')
xlim([PM_Plotting(1,1) PM_Plotting(end,1)])
ax=gca;
ax.FontSize = 7;
ax.FontWeight = 'bold';
xtickangle(45)
ylabel('PM_{2.5} Concentration (\mug/m^{3})','FontWeight','bold','FontSize',7) %** check font size ok
xlabel('Date, mmm-yyyy','FontWeight','bold','FontSize',7) %** check font size ok
hold on
%plot(PM_Plotting(:,1),PM_Plotting(:,2),'-')
title(sprintf('Filter-based PM_{2.5}, %s',Site_cities{loc}),'FontSize',11) %** check font size ok
plot(PM_Plotting(:,1),Avg_Line,'--','LineWidth',1.5)
lgd = legend('Filter Data','Site Average');
lgd.FontSize = 7; %** check font size ok
grid on
set(gcf, 'outerposition', [520 260 500 288]);
% to make plot smaller in case of stetson shutdown messing with fonts/sizing,etc
% may also need to play with fonts above, see "** check font size ok" above
% set(gcf,'outerposition',[520 260 490 275]);
lgd.Location = 'best';

saveas(gcf,sprintf('%s/%s/%s_PM25_SiteAvg.png',direc_output,'Chemical_Filter_Data/Plots/Filter_PM25_plots',Site_codes{loc}))
close all

clear Avg_Line Day_num End_Date end_idx fil Plotting_Dates PM_Filter_End PM_Filter_Start
clear PM_Plotting Site_PM_Avg Start_Date strt_idx tick_pos lgd



    
    
    