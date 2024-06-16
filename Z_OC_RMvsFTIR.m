%% %%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: Compare BC from HIPS and SSR

% Written by: Haihui Zhu
% Created: Oct 2023

clear ;close all ;clc
addpath('UtilityFunctions')

% Set up directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

direc_master = strcat(direc,'/Analysis_Data/Master_files');
savedir = sprintf('%s/Analysis_Data/Black_Carbon/HIPSvsSSR',direc);
savefname = sprintf('%s/Analysis_Data/Black_Carbon/HIPSvsSSR/HIPSvsSSR.xlsx',direc);


%-------------   SITE DETAILS   --------------
site_details = readtable(strcat(direc,'Site_Sampling/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));


% for loc = 25(Site_codes)
for loc = 1:length(Site_codes)

    % Find if there is a master file
    master_file = sprintf('%s/%s_master.csv',direc_master,Site_codes{loc});
    [Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, projectIDs_master,Master_hours, Master_masstype, ...
        Master_dates, Master_mass, Master_IC, Master_ICP, Master_XRF,...
        Master_carbon, Master_Method, Master_flags] = ReadMaster(master_file,Site_codes{loc});

    if isempty(Master_mass)
        continue
    else
        Filter_ID = Master_IDs;
        Cartridge_ID = CartridgeIDs_master;
        Mass_Type = Master_masstype;
        SSR = Master_carbon(:,1);
        HIPS = Master_carbon(:,2);

        T = table(Filter_ID,Cartridge_ID,Mass_Type,SSR,HIPS);

        writetable(T,savefname,'Sheet',Site_codes{loc})

        makescatter(SSR,HIPS,Mass_Type,Site_codes{loc},savedir);
        
    end
end

%% function

function fig = makescatter(SSR,HIPS,Mass_Type,Site_codes,savedir)
fontsize = 16;
scsz = 30;
colors=[0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; 0 0.4470 0.7410;0.3010 0.7450 0.9330;...
        0.8500 0.3250 0.0980; 0.4940 0.1840 0.5560; 0.6350 0.0780 0.1840];

% data preparation
validspot = find( ~isnan(SSR)& ~isnan(HIPS) & SSR>-1 & HIPS>-1);
if numel(validspot) <= 1
    fprintf('Less than 2 records, cannot make scatter plot\n')
else
 
X = SSR(validspot); Y = HIPS(validspot); C = Mass_Type(validspot);

R2 = corrcoef(X,Y,'rows','complete');
R2 = R2(1,2)^2; % R2
[m, b] = organic_regress(Y,X);
pm = '+-';
NRMSD = sqrt(mean((Y-X).^2))/mean(X);

tslope = sprintf('y = %.2fx%s%.1f',m,pm((b<0)+1),abs(b));
tRMSE = sprintf('NRMSD = %.2f',NRMSD);
tr = ['R^2 = ' sprintf('%.2f',R2)];
tN = sprintf('N = %d',numel(X));

ta1 = sprintf('%s\n%s\n%s\n%s',tr,tslope,tRMSE,tN);


% make the figure
figure 
fig=scatter(X, Y, scsz, C ,'filled');
colormap(gca,colors)
hold on

% add colorbar
c=colorbar;
caxis([-0.5 6.5]);
c.Ticks = 0:1:6;
c.TickLabels = ["Blank" "PM_{2.5}"  "PM_{10}" "TSP" "TSP-invalid" "Negative Mass" "Invalid Flow"];

% ajust axis and add axis label
rng2 = max(max(xlim),max(ylim));
rng1 = min(min(xlim),min(ylim));
set(gca,'Xlim',[rng1 rng2],'Ylim',[rng1 rng2]);
xlabel(sprintf('SSR BC (%s)','\mug'),'FontSize',fontsize)
ylabel(sprintf('HIPS BC (%s)','\mug'),'FontSize',fontsize)
pbaspect([1 1 1])

% add fitting line
xl=rng1:0.01*(rng2-rng1):rng2;
yl=m*xl + b;
plot(xl,yl,'Color','k','LineWidth',1.5)
hold on
% 1:1 line
yll=xl;
plot(xl,yll,'--k','LineWidth',1.2)


% add text
text(rng1+0.03*(rng2-rng1),rng2-0.03*(rng2-rng1),ta1, ... % statistics
    'horizontalalignment','left','verticalalignment','top','FontSize',fontsize);
text(rng2-0.03*(rng2-rng1),rng1+0.03*(rng2-rng1),Site_codes, ... % site code
    'horizontalalignment','right','verticalalignment','bottom','FontSize',fontsize);

% save figure
savefname = sprintf('%s/Scatter_%s.png',savedir,Site_codes);
saveas(gcf,savefname)
end

end


function [slope, offset] = organic_regress(Y,X) 
r = corrcoef([Y X]);
r = r(1,2);

b = regress(Y,[X ones(length(Y),1)]);
slope = b(1)/r;  
offset = b(2)+ mean(X)*b(1)*(1-1/r); 

end