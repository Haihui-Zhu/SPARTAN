function [] = Spec_barplot(Spec_data,filter_dates,filter_ids,Site_cities,Site_codes,loc,direc_output)

% --------------------

% Move water from bottom to just below OM
Spec_data = Spec_data(:, [2:end-2 1 end-1:end]) ;

Colors=[ 
%          109,207,246; % water
         237,48,41;  % red sulfate 
         240,103,166; % pink ammonium 
         245,126,32; % orange nitrate
         57,84,165; % blue sodium
         252,238,30; % yellow dust
         128,130,133; % grey TEO
         35,31,32; % black carbon
         109,207,246; % water   % Move water from bottom to just below OM
         55,98,60; % dark green OC
         80,184,72; % green residual
                        ]./255;
SpecName = {'Sulfate','Ammonium','Nitrate','Sea Salt','Dust','Trace Element Oxides','Black Carbon','Water','Organic Carbon','Residual Matter'};

avg = mean(Spec_data,1,'omitnan');

if sum(isnan(avg))<1 % no species that is all nan. 
    figure(3)
    hData = bar(Spec_data,'stacked');
    set(gca, 'xlim',[0.5,size(Spec_data,1) + .5], 'box','off','XTick',1:size(Spec_data,1));
    set(gcf,'outerposition',[190 80 1600 600]);
    set(gca,'XTickLabel',filter_ids)
    ax = gca;
    ax.FontSize = 6; %** check font size ok following stetson shutdown
    ax.FontWeight = 'bold';
    xtickangle(45)
    for sp = 1:length(SpecName)
        hData(sp).FaceColor = Colors(sp,:); hData(sp).EdgeColor = Colors(sp,:);
    end

    xlabel('Filter ID','fontsize', 8, 'fontweight', 'bold');%** check font size ok following stetson shutdown
    ylabel('Attributed Concentration (\mug m^{-3})','fontsize', 8, 'fontweight', 'bold');%** check font size ok following stetson shutdown
    hLegend = legend(SpecName,'fontsize',6,'fontweight','bold'); %** check font size ok following stetson shutdown
    hLegend.Location = 'best';
    title(sprintf('%s Chemical Speciation',Site_cities{loc}),'FontSize',10) %** check font size ok following stetson shutdown
    saveas(gcf,sprintf('%s%s/%s_PM25_Spec_all_data_bar.png',direc_output,'/Chemical_Filter_Data/Plots/Bar_spec_plots_filterlabels',Site_codes{loc}))
    
    % reset xtick labels and x-axis label and save with dates
%     set(gca,'XTick',1:4:size(merged_data_wet,1))
%     set(gca,'XTickLabel',x_dates(1:4:end,:))
    set(gca,'XTickLabel',filter_dates)
    xlabel('Filter sampling end date','fontsize', 8, 'fontweight', 'bold');% ** check font size ok following stetson shutdown
    saveas(gcf,sprintf('%s%s/%s_PM25_Spec_all_data.png',direc_output,'/Chemical_Filter_Data/Plots/Bar_spec_plots_dates',Site_codes{loc}))
    print(sprintf('%s%s/%s_PM25_Spec_all_data.eps',direc_output,'/Chemical_Filter_Data/Plots/Bar_spec_plots_dates',Site_codes{loc}),'-depsc')
    close
end


