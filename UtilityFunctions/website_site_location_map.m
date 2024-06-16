% Script to make a map for the website using shapefiles to define country
% borders

close all; clear all; clc

site_details = readtable('\\storage1.ris.wustl.edu\rvmartin\Active\SPARTAN-shared\\Site_Sampling\\Site_details.xlsx');
%site_details = readtable('/data1/spartan/Site_Sampling/Site_details.xlsx');
plot_status = table2array(site_details(:,13)); % status = 1 means currently sampling, 0 is a retired site, 2 is a future/planned site 
% index active sites
active_idx = find(plot_status == 1); % active sites status =1
lats_active = table2array(site_details(active_idx,5));
longs_active = table2array(site_details(active_idx,6));
% index retired sites
retired_idx = find(plot_status == 0); % retired sites status =0
lats_retired = table2array(site_details(retired_idx,5));
longs_retired = table2array(site_details(retired_idx,6));
% index future sites
future_idx = find(plot_status == 2); % future sites status =2
lats_future = table2array(site_details(future_idx,5));
longs_future = table2array(site_details(future_idx,6));

% get country boundaries based on GBD 
load('\\storage1.ris.wustl.edu\rvmartin\Active\SPARTAN-shared\Public_Data\Documents\GBD_Country_Masks.mat'); % GBD countries
%load('/misc/data1/spartan/Public_Data/Documents/GBD_Country_Masks.mat'); % GBD countries

countries_lats = []; 
countries_lons = []; 

for i = 1:length(sGBDCountries)
    test_lats = sGBDCountries(i).Lat;
    test_lons = sGBDCountries(i).Lon;
    countries_lats = [countries_lats; test_lats]; clear test_lats
    countries_lons = [countries_lons; test_lons]; clear test_lons
end

% Create figure
figure(1)
h = worldmap([-55 70],[-140 170]);
setm(gca,'Grid','off','MapProjection','miller','parallellabel','off','meridianlabel','off')
geoshow(countries_lats, countries_lons, 'DisplayType','line','color',[0.05 0.05 0.05]);
scatterm(lats_active, longs_active, 20, 'ok','filled','MarkerFaceColor','b','LineWidth',1); % currently active sites are in blue
scatterm(lats_retired, longs_retired, 20, 'ok','filled','MarkerEdgeColor','[0.5 0.5 0.5]','LineWidth',1); % retired sites are in grey
% scatterm(lats_future, longs_future, 20,'ok','filled','MarkerEdgeColor','y','LineWidth',1); % future sites are in yellow
% set(gcf,'outerposition',[520 260 490 275]); % can fiddle with this if linking directly to webiste, but for now the figure is NOT linked
print('\\storage1.ris.wustl.edu\rvmartin\Active\SPARTAN-shared\\Public_Data\\Chemical_Filter_Data\\Plots\\Website_Site_Locations.png','-dpng','-r400','-loose')
print('\\storage1.ris.wustl.edu\rvmartin\Active\SPARTAN-shared\\Public_Data\\Chemical_Filter_Data\\Plots\\Website_Site_Locations.eps','-depsc','-r400','-loose')

disp('Finished Making Map')





