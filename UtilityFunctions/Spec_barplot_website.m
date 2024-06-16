

function [] = Spec_barplot_website(PM25_dates,Site_codes, direc_output,loc, merged_data_wet)

% site specific instructions
% x is used by the main plotting section to place the legend in a location
% appropriate for the data from that site (e.g. as not to interfere with data, if possible)
% ----------------------------------------------------

% Move water from bottom to just below OM
merged_data_wet = merged_data_wet(:, [2:end-2 1 end-1:end]) ;

Colors=[ 237,48,41;  % red sulfate 
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


% exclude filters that don't have enough speciation data
ind = find(merged_data_wet(:,1)>0 & sum(merged_data_wet(:,2:3),2,'omitnan')>0 & merged_data_wet(:,5)>0 ); % Sulfate; Ammonium + Nitrate;  Dust

plotlength = length(ind);

if plotlength < 6 % less than 1 cartridge
    return
elseif plotlength < 50
    % as more data from these sites become available, they can be added
    % with an "elseif" statement as done for sites above
    dates(:,1) = datenum(PM25_dates(ind,5), PM25_dates(ind,6), PM25_dates(ind,7),0,0,0);
    x_dates = char(datetime(dates,'ConvertFrom','datenum', 'Format','MM-dd-yyyy'));
    figure
    hData = bar(merged_data_wet(ind,:),'stacked');
    set(gca, 'xlim',[0.5, plotlength+.5], 'box','off','XTick',1:plotlength);
    x = [0.854 0.850 0.0128 0.0403];

else
    % as more data from these sites become available, they can be added
    % with an "elseif" statement as done for sites above
    ind2 = ind(end-48:end);
    dates(1:49,1) = datenum(PM25_dates(ind2,5), PM25_dates(ind2,6), PM25_dates(ind2,7),0,0,0);
    x_dates = char(datetime(dates,'ConvertFrom','datenum', 'Format','MM-dd-yyyy'));
    figure
    hData = bar(merged_data_wet(ind2,:),'stacked');
    set(gca, 'xlim',[0.5, 49.5], 'box','off','XTick',1:49);
    x = [0.854 0.850 0.0128 0.0403];

end

% formatting
set(gca,'XTickLabel',x_dates)
ax=gca;
ax.FontSize = 6; %** check font size ok following stetson shutdown
ax.FontWeight = 'bold';
xtickangle(45)
for sp = 1:length(SpecName)
    hData(sp).FaceColor = Colors(sp,:); hData(sp).EdgeColor = Colors(sp,:);
end
Title_label = title('Example Reconstructed Fine Mass by Filter','fontsize', 9, 'fontweight', 'bold'); %** check font size ok following stetson shutdown
xlab = xlabel('Filter Sample End Date, mm-dd-yyyy','fontsize', 6, 'fontweight', 'bold'); %** check font size ok following stetson shutdown
ylab = ylabel('Attributed Concentration (\mug m^{-3})','fontsize', 6, 'fontweight', 'bold'); %** check font size ok following stetson shutdown
hLegend = legend(SpecName,'fontsize',4,'fontweight','bold'); %** check font size ok following stetson shutdown
hLegend.Position = x;
ax.YLim(1) = 0; % make sure plot does not go below attributed concentration of 0

set(gcf, 'outerposition', [520 260 500 288]);
% to make plot smaller in case of stetson shutdown messing with fonts/sizing,etc
% may also need to play with fonts above, see "** check font size ok" above
% set(gcf,'outerposition',[520 260 490 275]);
saveas(gcf,sprintf('%s/%s/%s_PM25_bar_website.png',direc_output,'Chemical_Filter_Data/Plots/Bar_spec_website',Site_codes{loc}))
print(sprintf('%s/%s/%s_PM25_bar_website.eps',direc_output,'Chemical_Filter_Data/Plots/Bar_spec_website',Site_codes{loc}),'-depsc')

close
clear hData dates x_dates hLegend xlab ylab Title_label x hLegend.Position












