%% %%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: collect reanalysis T, P for SPARTAN sites
% Haihui Zhu
% January 25, 2024

close all; clear all; clc

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set directories
% direc = '\\storage1.ris.wustl.edu\rvmartin\Active\SPARTAN-shared\';
direc_in = '/storage1/fs1/rvmartin/Active/SPARTAN-shared/';

direc_sampling = strcat(direc_in,'Site_Sampling');
direc_master = strcat(direc_in,'Analysis_Data/Master_files');
direc_met = '/storage1/fs1/rvmartin/Active/GEOS-Chem-shared/ExtData/GEOS_0.5x0.625/MERRA2/';

%-------------   SITE INFO   --------------
% Sites are listed in the file "Site_details.xlsx" in alphabetical order
% Column info:
% 1 = 4-letter site codes
% 2 = Country
% 3 = City
site_details = readtable(strcat(direc_sampling,'/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));
latitudes = table2array(site_details(:,5));
longitudes = table2array(site_details(:,6));

% the output file
sfname = sprintf('%s/MERRA-2_Surface_T_P_at_SPARTAN_2019-2023.mat',direc_sampling);

Years = 2020:2023;

for Yr = Years % do it year by year, so that it doesn't need to start over if crashed in the middle
    fprintf('processing for the year %d\n',Yr)
    DatesNum = datenum([Yr 01 01]):1:datenum([Yr 12 31]);
    D1_Dates = datevec(DatesNum);

    Ti = nan(length(DatesNum),24,length(Site_cities));
    Pi = nan(length(DatesNum),24,length(Site_cities));

    for d = 1:length(DatesNum)
        metfname = sprintf('%s/%d/%.2d/MERRA2.%d%.2d%.2d.I3.05x0625.nc4',direc_met,D1_Dates(d,1),D1_Dates(d,2),...
            D1_Dates(d,1),D1_Dates(d,2),D1_Dates(d,3));
        if exist(metfname)
            Tm = ncread(metfname,'T');  % K
            Tm = squeeze(Tm(:,:,1,:)); % only keep the surface temp
            Pm = ncread(metfname,'PS'); % Pa
            lat = ncread(metfname,'lat'); lon = ncread(metfname,'lon');
            % interp to get T & P at site
            for hr = 1:8
                Ti(d, (hr-1)*3+1 ,:) = interp2(lat',lon,Tm(:,:,hr),latitudes, longitudes);
                Pi(d, (hr-1)*3+1 ,:) = interp2(lat',lon,Pm(:,:,hr),latitudes, longitudes);
                Ti(d, (hr-1)*3+2 ,:) = Ti(d, (hr-1)*3+1 ,:); % repeat the data for another 2 hours since the reanalysis is a 3-hour mean
                Ti(d, (hr-1)*3+3 ,:) = Ti(d, (hr-1)*3+1 ,:);
                Pi(d, (hr-1)*3+2 ,:) = Pi(d, (hr-1)*3+1 ,:);
                Pi(d, (hr-1)*3+3 ,:) = Pi(d, (hr-1)*3+1 ,:);
            end
        else
            fprintf('%s not exist; skipped\n',metfname)
        end
    end
    sfname2 = sprintf('%s/MERRA-2_Surface_T_P_at_SPARTAN_%d.mat',direc_sampling,Yr);
    save(sfname2,'Pi','Ti','D1_Dates','DatesNum','Site_codes')


    if exist(sfname) == 2
        load(sfname)
        P(end+1:end+size(Pi,1),:,:) = Pi;
        T(end+1:end+size(Ti,1),:,:) = Ti;
        alldates_vec = [alldates_vec; D1_Dates ];
        alldates_num = [alldates_num; DatesNum' ];
        % sort by dates
        [alldates_num, i] = sort(alldates_num);
        alldates_vec = alldates_vec(i,:);
        P = P(i,:,:);
        T = T(i,:,:);
        save(sfname,'P','T','alldates_vec','alldates_num','Site_codes')
    else
        P = Pi;
        T = Ti;
        alldates_vec =  D1_Dates ;
        alldates_num =  DatesNum' ;
        save(sfname,'P','T','alldates_vec','alldates_num','Site_codes')
    end
    fprintf('%s saved\n',sfname)
end


