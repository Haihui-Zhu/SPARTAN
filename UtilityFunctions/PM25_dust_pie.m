
function [] = PM25_dust_pie(species, Site_cities,Site_codes,loc,direcin)

SpecName = {'Al','Si','Ca','Fe','Pb','K','Mg','Na','Zn','Mn'};
% change the order:
species = species(:,[8 7 6 5 4 3 2 1 9 10]);
SpecName = SpecName([8 7 6 5 4 3 2 1 9 10]);
notes ='';

% 2022-05-15 Haihui: replace negative species (usually OM) by 0
Ind = find(species<0); 
for i = 1:length(Ind) % prepare note in the chart
    notes = sprintf('%s\n%s = %.2f and is set to 0 in the chart.',notes,SpecName{Ind(i)},species(Ind(i)));
end
species(Ind)=0; % set negative to 0;
% --------------------
Colors = [ 240,103,166; % pink Al 
             109,207,246; % light blue Si
             245,126,32; % orange Ca
             237,48,41;  % red Fe
             57,84,165; % blue Pb
             252,238,30; % yellow K
             128,130,133; % grey Mg
             35,31,32; % black Na
             80,184,72; % green Zn
             200,50,120; % unknown Mn
                        ]./255;
Colors = Colors([8 7 6 5 4 3 2 1 9 10],:); % change the order

% Colors =
% {'#A2142F','#0072BD','#D95319','k','#77AC30','#EDB120','#4DBEEE',''}; % Randall colors
           % SO4         NO3     NH4      BC    OM       Dust        SS
figure(1)
h = pie(species); % this prints fractions
% h = pie(species,SpecName); % this prints species names 
htext = findobj(h,'Type','text');
hPatch = findobj(h,'Type','patch');
for sp = 1:length(SpecName)
    hPatch(sp).FaceColor = Colors(sp,:); 
    htext(sp).String=''; 
%     htext(sp).FontSize = 20;
    hold on
end

% centering the PM value in middle
% PM25str = sprintf('%1.0f',PM25_total);
% if PM25_total < 100
%     fz = 44;
%     sc = 0.11;
% else
%     fz = 40;
%     sc = 0.12;
% end
% plot(0, 0,'o', 'MarkerSize',90,'MarkerEdgeColor','k','MarkerFaceColor','w'); % PM2.5 white circle in centre
% text(-sc*length(PM25str),0,sprintf('%1.0f',PM25_total),'fontsize',fz,'fontweight','bold')
% hold on

% 2022-05-15 Haihui: adding notes about negative number
if ~isempty(notes)
    text(-0.3, 0, notes,'units','normalized')
end
% --------------------

hold off
T = title(sprintf('%s',Site_cities{loc}),'FontSize',30,'FontWeight','bold');
set(T,'position',[-0.0039 1.05 1.001]);

saveas(gcf,sprintf('%s/Public_Data/Chemical_Filter_Data/Plots/Pie_spec_plots/%s_PM25_dust.png',direcin,Site_codes{loc}))
print(sprintf('%s/Public_Data/Chemical_Filter_Data/Plots/Pie_dust_plots/%s_PM25_dust.eps',direcin,Site_codes{loc}),'-depsc')

close all

