%% %%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: to read hourly-averaged nephelometer data files and PM2.5 and
% kappa values to create time-resolved PM2.5 estimates

% Written by: Crystal Weagle
% Created: 11 March 2019

% EDITS: 
% -- 2020-04-03, by Crystal Weagle: Included derivative screening to
% nephelometer data prior to calculating PM2.5. If a MAPLE site, derivative
% screening threshold is 50 Mm-1, for all other sites is 250 Mm-1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear ; clc
addpath('UtilityFunctions')

% Setup directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

% the diary function saves the processing history into an annual record
diary  (sprintf('%s/Public_Data/Data_Processing_Records/%s_TimeResolvedPM25_Record',rootdir,datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))

%% %%%%%%% USER SWITCHES %%%%%%%%%%
direc_PM25_in = sprintf('%s/Public_Data/Time-resolved_PM2.5/Input_data',rootdir); 
direc_Neph = sprintf('%s/Public_Data/Neph_Processed',rootdir); 

direc_PM25_out = sprintf('%s/Public_Data/Time-resolved_PM2.5/Data_sharing',rootdir); 
direc_plot_out = sprintf('%s/Public_Data/Time-resolved_PM2.5/Plots',rootdir); 

max_RH = 80;
deriv_screen = 1;
data_version = 'Data version 1.1';

%-------------   SITE DETAILS   --------------
% Sites are listed in the file "Site_details.xlsx" in alphabetical order based on 4-letter site code
% Column info:
% 1 = 4-letter site codes
% 2 = Country
% 3 = City
% 4 = Host Institute
% 5 = Latitude
% 6 = Longitude
% 7 = Elevation (meters)
% 8 = Sampling Mode (SPARTAN or MAIA)

%site_details = readtable('/data1/spartan/Site_Sampling/Site_details.xlsx');
site_details = readtable(sprintf('%s/Site_Sampling/Site_details.xlsx',rootdir));
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));
Site_countries = table2array(site_details(:,2));
latitudes = table2array(site_details(:,5));
longitudes = table2array(site_details(:,6));
elevations = table2array(site_details(:,7));  

%PM25_parameters = readtable('/data1/spartan/Public_Data/Chemical_Filter_Data/Sampling_Parameters_Methods.xlsx','Sheet','PM2.5 mass');
PM25_parameters = readtable(sprintf('%s/SOPs/Public SOPs/Sampling_Parameters_Methods_2.3.xlsx',rootdir),'Sheet','PM2.5 mass');

%% --------- Find and read data files from file directory ---------
for loc =  1:length(Site_codes)
      
    
    %%  ---- Read in PM2.5 data file -----

    if ~exist(sprintf('%s/PM25only_%s.csv',direc_PM25_in,Site_codes{loc}),'file') % checks if a filter masses file exists
        fprintf('WARNING: No PM2.5 file exists for %s \n', Site_codes{loc})
        continue % if no PM2.5 file exists, moves to the next site
    end
    
    fprintf('Reading PM2.5 file for %s \n', Site_codes{loc})
    PM25_file = sprintf('%s/PM25only_%s.csv',direc_PM25_in,Site_codes{loc});
    PM25_dataB = table2array(readtable(PM25_file));
    
    % PM25_data columns are:
    % 1. start date in datenum format
    % 2. end date in datenum format
    % 3. PM2.5 mass (ug/m3)
    % 4. kappa volume (kv, unitless)
    
    % need to make sure the data is ordereer by date, not by filter number as in master files
    PM25_data = sortrows(PM25_dataB); clear PM25_dataB % sorts based on values in first column, which is datenum date format
    
    %% Read in hourly Neph Bsp file
    names2 = dir(sprintf('%s/%s/%s_%s_*_hourly.csv',direc_Neph,Site_codes{loc},Site_codes{loc},Site_cities{loc}));
    
    if isempty(names2) == 1
        fprintf('WARNING: No hourly neph file exists for %s \n', Site_codes{loc})
        continue % if no hourly neph file exists, moves to the next site
    end
    
    k =0;
    for i = 1:size(names2,1)
        if names2(i).bytes > 5000
            k = k+1;
            names(k) = names2(i);
        end
    end
    files = {names.name}; 
    files = files'; clear names
    
    Neph_dataB = [];
    
    for i = 1:length(files)
        filename = files(i);
        if contains(files(i),'TSP') || contains(files(i),'PM25')
            Neph_dataC = csvread(sprintf('%s/%s/%s',direc_Neph,Site_codes{loc},char(filename)),2,0);
            Neph_dataB = [Neph_dataB; Neph_dataC];
            clear Neph_dataC
        end
    end
  
    Neph_data(:,1) = datenum(Neph_dataB(:,1), Neph_dataB(:,2), Neph_dataB(:,3), Neph_dataB(:,4),0,0);
    Neph_data(:,2:10) = Neph_dataB(:,5:end); clear Neph_dataB
    Neph_data = sortrows(Neph_data);  % sorts based on values in first column
    
    % To account for transient meteorological events a derivative-based screening protocol is applied 
    % to screen hourly scatter values where the change is greater than 50 Mm-1 hour-1 in low PM2.5 
    % environments and 250 Mm-1 hour-1 in high PM2.5 environments. Values exceeding the threshold are 
    % not used to estimate hourly PM2.5.
    if deriv_screen == 1
        deriv = diff(Neph_data(:,8:10));
        data_version = 'Data version 2.0 deriv_screen = ';
        if Site_codes{loc} == 'CAKE' | Site_codes{loc} == 'CALE' | Site_codes{loc} == 'CASH' | Site_codes{loc} == 'CADO' | Site_codes{loc} == 'CAHA'
            deriv_thres = 50;
        else
            deriv_thres =250;
        end
        for i = 1:3
            j = i + 7; % column index in Neph_data matrix
            screen_high = find(deriv(:,i) >deriv_thres) + 1;
            screen_low = find(deriv(:,i) < -deriv_thres) + 1;
            Neph_data(screen_high,j) = nan;
            Neph_data(screen_low,j) = nan;
            clear screen_high screen_low
        end
        clear i j deriv
    end
    
    Neph_data(Neph_data(:,9) < 0,9) = NaN;
    bsp_nan = find(isnan(Neph_data(:,9))); % find rows where total scatter green is Nan
    Neph_data = removerows(Neph_data,bsp_nan); clear bsp_nan
    % Neph_data columns are:
    % 1. Date + hour (datenum format)
    % 2. Temp (degC)
    % 3. Pressure
    % 4. RH (%)
    % 5. Back Scat Red (Mm-, ambient RH)
    % 6. Back Scat Green (Mm-, ambient RH)
    % 7. Back Scat Blue (Mm-, ambient RH)
    % 8. Total Scat Red (Mm-, ambient RH)
    % 9. Total Scat Green (Mm-, ambient RH)
    % 10. Total Scat Blue (Mm-, ambient RH) 
    
    %% Find overlapping sampling periods
    
    merged_file_hourly = [];
  
    % Create date columns in merged file, this is continuous from first
    % filter to end of last filter. There may be large gaps with no data
    % but those will either be filled with lower level data or removed.
    merged_file_hourly(:,1) = PM25_data(1,1):(1/24):PM25_data(end,2);
    
    if Neph_data(end,1) > PM25_data(end,2) % the nephelometer data extends beyond the filter data
        extend = PM25_data(end,2):(1/24):Neph_data(end,1); extend = extend';
        merged_file_hourly = [merged_file_hourly; extend];
        clear extend
    end
    merged_file_hourly(:,2:8) = NaN;

    for i = 1:size(PM25_data,1)
        
        %-- hourly_data/merged_file columns are:
        % 1. date + hour in datenum format
        % 2. data level > 1 = default kappa, no filter masses; 2 = default kappa with filter masses; 3 = filter-specific kappa and filter masses
        % 3. kappa
        %.4. RH (NaN if not measured or > max_RH)
        % 5. fRH, based on  kappa
        % 6. Bsp_green, ambient RH
        % 7. Bsp_green, dry (0 % RH)
        % 8. Hourly-estimated PM2.5 (35 % RH)
        
        hourly_data(:,1) = PM25_data(i,1):(1/24):PM25_data(i,2); % hours sampled by filter for indexing neph data
        hourly_data(:,2) = 3; % set data level to 3 for now as filter masses are available. Index below will set to 2 if no filter-specific kappa found
        
        if PM25_data(i,4) < 0.1
            hourly_data(:,3) = 0.1; % round kappa to 0.1
        else
            hourly_data(:,3) = PM25_data(i,4); % kappa estimated by filter measurements
        end

        [~, bsp_idx] = ismembertol(hourly_data(:,1), Neph_data(:,1),1e-10); % find neph data sampled same time as filter

        missing_bsp = find(bsp_idx == 0); % will give index if there is no scatter for a given sampled period
       
        if isempty(missing_bsp) == 0 % there is period of time with filter data but no scatter
            filler = 1:length(bsp_idx)'; % sets up a list of 1 to length of hourly_data
            bsp_idx2 = removerows(filler',missing_bsp); clear filler
            bsp_idx = removerows(bsp_idx,missing_bsp); % remove gaps in scatter data that overlaps with filter sample
            
            hourly_data(bsp_idx2,4) = Neph_data(bsp_idx,4); % ambient RH
            hourly_data(bsp_idx2,6) = Neph_data(bsp_idx,9); % Bsp_green @ ambient RH
            hourly_data(missing_bsp,[4 6]) = NaN; % setting period with no scatter to NaN
            clear bsp_idx2
        else
            hourly_data(:,4) = Neph_data(bsp_idx,4); % ambient RH
            hourly_data(:,6) = Neph_data(bsp_idx,9); % Bsp_green, ambient RH
        end
        
        hourly_data((hourly_data(:,4) > max_RH),4) = NaN; % set RH above RH_max to NaN for purposs of inferring hourly PM2.5
        noK_idx = find(isnan(hourly_data(:,3))); % find hours without a kappa from filter mesurements
        hourly_data(noK_idx,3) = 0.2; % set kappa to default 0.2 if none obtained from filter
        hourly_data(:,5) = 1+ hourly_data(:,3).*(hourly_data(:,4)./(100-hourly_data(:,4))); % calculate fRH based on filter-specific kappa and ambient RH
        hourly_data(:,7) = hourly_data(:,6)./hourly_data(:,5); % calculate dry (0% RH) Bsp_green 
        %hourly_data(:,8)=(PM25_data(i,3)/nanmean(hourly_data(:,7))).*hourly_data(:,7);
        hourly_data(:,8) = (PM25_data(i,3)/mean(hourly_data(:,7),"omitnan")).*hourly_data(:,7);   % calculate hourly PM2.5 estimates
        
        % at this point, data will be either level 2 or 3, bc filter masses are available
        % use noK_idx to set level to 2 or 3
        hourly_data(noK_idx,2) = 2; 
        
        [~, date_idx] = ismember(hourly_data(:,1), merged_file_hourly(:,1)); % find index of hourly_data dates in merged_file
        hourly_data = removerows(hourly_data,find(date_idx ==0));
        date_idx = removerows(date_idx,find(date_idx==0));
        merged_file_hourly(date_idx,2:8) = hourly_data(:,2:8);
        
        clear hourly_data noK_idx bsp_idx missing_bsp date_idx
 
    end
    

    %%  Infer PM2.5 hourly where possible with record average PM2.5/Bsp_green ratio----
    % this section allows us to estimate PM2.5 hourly from only nephelometer data prior to receiving filters
    % uses default kappa
    
    % ---- fill in merged file with Bsp data where no filter data is available ---- 
    % indexing
    nolvl_idx2 = find(isnan(merged_file_hourly(:,2))); % where there is no lvl there is no PM2.5 or Bsp filled in
    [~, bspfill_idx] = ismember(merged_file_hourly(nolvl_idx2,1), Neph_data(:,1)); % find neph data sampled where no filter data is available
    bspfill_idx = removerows(bspfill_idx,find( bspfill_idx == 0)); % remove rows where there is no Bsp data
    [~, chck_idx] = ismember(Neph_data(bspfill_idx,1),merged_file_hourly(nolvl_idx2,1)); % reindex to find merged_file index where bsp data exists and no lvl available
    nolvl_idx = nolvl_idx2(chck_idx);  clear nolvl_idx2% reindex to find merged_file index where bsp data exists and no lvl available
    % fill in Bsp, kappa, RH  in merged file
    merged_file_hourly(nolvl_idx,4) = Neph_data(bspfill_idx,4); % ambient RH
    merged_file_hourly(nolvl_idx,6) = Neph_data(bspfill_idx,9); % Bsp_green, ambient RH
    merged_file_hourly((merged_file_hourly(:,4) > max_RH),4) = NaN; % set RH above RH_max to NaN for purposs of inferring hourly PM2.5
    merged_file_hourly(nolvl_idx,3) = 0.2; % set kappa to default 0.2 if none obtained from filter
    
    merged_file_hourly(nolvl_idx,5) = 1+ merged_file_hourly(nolvl_idx,3).*(merged_file_hourly(nolvl_idx,4)./(100-merged_file_hourly(nolvl_idx,4))); % calculate fRH based on filter-specific kappa and ambient RH
    merged_file_hourly(nolvl_idx,7) = merged_file_hourly(nolvl_idx,6)./merged_file_hourly(nolvl_idx,5); % calculate dry (0% RH) Bsp_green
    % now everything is all set for the next step!
       
    % need to find rows where level is nan but bsp is available, and there
    % is data available prior to for calculation of record PM2.5/Bsp (mass scattering efficiency). 
    % nolvl_idx gives the index of rows that are missing hourly PM2.5 and have Bsp, but RH may have been to high, therefore now searching for
    % an index that finds all rows where PM2.5 hourly is missing due to lack of calculation (it could be missing due to ambient RH > RH_max, 
    % for example). REMINDER: lack of a level indicates calculation did not happen
    clear nolvl_idx
    
    test = find(~isnan(merged_file_hourly(:,7))); % finds where dry Bsp_green exists
    
    inferPM_idx =test(find(isnan(merged_file_hourly(test,8)))); clear test% indexes to find where Bsp_green_dry exists (RH < RH_max) and where there is NaN in the PM2.5 hourly column
    
    % if there is hourly PM2.5, there must be Bsp, else no PM2.5 could have been inferred, thus sum all PM2.5 and Bsp prior to
    % inferPM_idx(i) to get record mass scattering efficency (PM.25/Bsp)
    mass_scat_rec = mean(merged_file_hourly(:,8),'omitnan')/mean(merged_file_hourly(:,7),'omitnan'); % <PM2.5>/<Bsp>; <> indicates record average

    if isempty(inferPM_idx) == 0
        for i = inferPM_idx(1:end)
            before(i) = merged_file_hourly(i,8);
            merged_file_hourly(i,8) = mass_scat_rec*merged_file_hourly(i,7);
            after(i) = merged_file_hourly(i,8);
            merged_file_hourly(i,2) = 1; % set level = 1
        end
    end
    
    clear  inferPM_idx i mass_scat
    merged_file_hourly(isinf(merged_file_hourly)) = NaN;
    

    %% Create daily file 
    
        %-- merged_file_hourly columns are:
        % 1. date + hour in datenum format
        % 2. data level > 1 = default kappa, no filter masses; 2 = default kappa with filter masses; 3 = filter-specific kappa and filter masses
        % 3. kappa
        % 4. RH (NaN if not measured or > max_RH)
        % 5. fRH, based on default kappa
        % 6. Bsp_550nm, ambient RH
        % 7. Bsp_550nm, dry (0 % RH)
        % 8. Daily PM2.5 (35 % RH, average of hourly-estimated)
        
        notnan_hourly = ~isnan(merged_file_hourly(:,8));
        merged_file_hourly = merged_file_hourly(notnan_hourly,:); clear notnan_hourly
        
        dates_expand = datevec(merged_file_hourly(:,1));
        data_dates = datenum(dates_expand(:,1), dates_expand(:,2), dates_expand(:,3),0,0,0);
        unique_dates = unique(data_dates);
        
        if isempty(merged_file_hourly) == 0
            for i = 1:length(unique_dates)
                this_day = unique_dates(i); 
                merged_file_daily(i,1) = this_day; % date column (serial format)
                day_idx = find(data_dates == this_day);
                levels_day = merged_file_hourly(day_idx,2);
                merged_file_daily(i,2) = round(mean(levels_day,"omitnan")); % keep round?
                merged_file_daily(i,3) = length(day_idx); % determines number of hourly data points used to create daily values
                merged_file_daily(i,4:9) = mean(merged_file_hourly(day_idx,3:8),1, "omitnan");
                
                clear this_day levels_day
            end
        elseif isempty(merged_file_hourly) == 1
            clear bspfill_idx chck_idx data_dates dates_expand filename files i k merged_file_hourly names2 
            clear Neph_data PM25_data PM25_file unique_dates
            continue
        end
    
    
    %% Make plot 
    % only plot most recent days with data, up to and including most recent data point 
    % locate most recent data point 
    
    idx_last  = find(~isnan(merged_file_daily(:,9)), 1, 'last'); 
        
            %recent_idx = idx_last-7:idx_last;  use me if less than 30 days of recent data
            %recent_idx = idx_last-34:idx_last-4; use me if need to select
            %a particular section of data (most recent data very small and
            %far from the last section that is 30 days)
    recent_idx = idx_last-28:idx_last; %usual case

    figure
    plot(merged_file_daily(recent_idx,1), merged_file_daily(recent_idx,9),'-ob','MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b','MarkerSize',6,'Linewidth',1.1)
    set(gcf,'outerposition',[520 260 670 400]); 
    xticks([merged_file_daily(recent_idx(1),1):2:merged_file_daily(recent_idx(end),1)]) % change "2" to 1" if less than 12 points
    xlim([merged_file_daily(recent_idx(1),1) merged_file_daily(recent_idx(end),1)])
    datetick('x','mm-dd-yyyy','keeplimits','keepticks')
    ax=gca;
    ax.FontSize = 8;
    ax.FontWeight = 'bold';
    xtickangle(45)
    ylabel('PM_{2.5} \mug/m^{3}','FontWeight','bold','FontSize',9);
    xlabel('Date, mm-dd-yyyy','FontWeight','bold','FontSize',9);
    title(sprintf('Daily Estimated PM_{2.5} Concentration, %s',Site_cities{loc}),'FontSize',12)
    grid on
    
%     if Site_codes{loc} ~='VNHN' % not saving VNHN right now due to short recent data availability
        saveas(gcf,sprintf('%s/%s_EstimatedPM2.5.png',direc_plot_out,Site_codes{loc}))
%     end
    
    
    %if Site_codes{loc}~= 'ZAPR' 
      %  recent_idx = idx_last-4:idx_last; % use me if less than 30 days of recent data
   % end
    close all
    
    clear idx_last recent_idx xlab ylab
    
    %% Save site average PM.25/Bsp ratio to table 
    
    opts = detectImportOptions(sprintf('%s/MassScat_Record.csv',direc_PM25_in));
    opts.DataLine = 2;
    MassScat_Record_initial = readtable(sprintf('%s/MassScat_Record.csv',direc_PM25_in),opts); 
    
    MassScat_Record(:,1) = table2array(MassScat_Record_initial(:,1));
    MassScat_Record(:,2) = num2cell(table2array(MassScat_Record_initial(:,2)));
  
    site_idx = find(contains(MassScat_Record(:,1), Site_codes(loc)));
    MassScat_Record(site_idx,2) = num2cell(mass_scat_rec);
    
    Titles = {'SiteCode','MassScat'};
    
    MassScat_table = cell2table(MassScat_Record,'VariableNames',Titles);
    writetable(MassScat_table,sprintf('%s/MassScat_Record.csv',direc_PM25_in));
    clear site_idx MassScat_Record MassScat_table  MassScat_Record_initial opts  Titles

    
    %% Save data files
    
    % remove rows with NaN PM2.5 values
    notnan_daily =~isnan(merged_file_daily(:,9));
    merged_file_daily = merged_file_daily(notnan_daily,:); clear notnan_daily
    
    Titles_hourly = {'Site_Code' 'Country' 'City' 'Latitude' 'Longitude' 'Elevation_meters' 'Year_local' 'Month_local' 'Day_local' 'hour_local'...
        'Parameter_Code' 'Parameter_Name' 'Value' 'Units' 'Method_Code' 'Collection_Description' 'Analysis_Description' 'Conditions'};
    
    Titles_daily = {'Site_Code' 'Country' 'City' 'Latitude' 'Longitude' 'Elevation_meters' 'Year_local' 'Month_local' 'Day_local' 'Hours_sampled'...
        'Parameter_Code' 'Parameter_Name' 'Value' 'Units' 'Method_Code' 'Collection_Description' 'Analysis_Description' 'Conditions'};

    
    Hourly_public = cell(size(merged_file_hourly,1),size(Titles_hourly,2));
    Daily_public = cell(size(merged_file_daily,1),size(Titles_daily,2));
    
    % prefill info that is the same for all points
    Hourly_public(1:size(Hourly_public,1),1) = Site_codes(loc); % site code
    Hourly_public(1:size(Hourly_public,1),2) = Site_countries(loc); % country
    Hourly_public(1:size(Hourly_public,1),3) = Site_cities(loc); % city
    Hourly_public(1:size(Hourly_public,1),4) = num2cell(latitudes(loc)); % latitude
    Hourly_public(1:size(Hourly_public,1),5) = num2cell(longitudes(loc)); % longitude
    Hourly_public(1:size(Hourly_public,1),6) = num2cell(elevations(loc)); % elevation
    Hourly_public(1:size(Hourly_public,1),end) = cellstr('Ambient local');
    
    Daily_public(1:size(Daily_public,1),1) = Site_codes(loc); % site code
    Daily_public(1:size(Daily_public,1),2) = Site_countries(loc); % country
    Daily_public(1:size(Daily_public,1),3) = Site_cities(loc); % city
    Daily_public(1:size(Daily_public,1),4) = num2cell(latitudes(loc)); % latitude
    Daily_public(1:size(Daily_public,1),5) = num2cell(longitudes(loc)); % longitude
    Daily_public(1:size(Daily_public,1),6) = num2cell(elevations(loc)); % elevation
    Daily_public(1:size(Daily_public,1),end) = cellstr('Ambient local');
    
    for  i = 1:size(merged_file_hourly,1)
        
        this_day = datevec(merged_file_hourly(i,1));
        Hourly_public(i,7:10) = num2cell(this_day(1:4)); clear this_day
        Hourly_public(i,13) = num2cell(round(merged_file_hourly(i,8),1)); % PM2.5 concentration
        
        Hourly_public(i,11) = num2cell(table2array(PM25_parameters(merged_file_hourly(i,2)+5,2))); % parameter code
        Hourly_public(i,12) = table2array(PM25_parameters(merged_file_hourly(i,2)+5,1)); % parameter name
        Hourly_public(i,14) = table2array(PM25_parameters(merged_file_hourly(i,2)+5,7)); % units
        Hourly_public(i,15) = num2cell(table2array(PM25_parameters(merged_file_hourly(i,2)+5,3))); % method code
        Hourly_public(i,16) = table2array(PM25_parameters(merged_file_hourly(i,2)+5,5)); % collection description
        Hourly_public(i,17) = table2array(PM25_parameters(merged_file_hourly(i,2)+5,6)); % analysis description

    end
    
    for  i = 1:size(merged_file_daily,1)
        
        this_day = datevec(merged_file_daily(i,1));
        Daily_public(i,7:9) = num2cell(this_day(1:3)); clear this_day
        Daily_public(i,10) = num2cell(merged_file_daily(i,3));
        Daily_public(i,13) = num2cell(round(merged_file_daily(i,9),1));
        
        Daily_public(i,11) = num2cell(table2array(PM25_parameters(merged_file_daily(i,2)+5,2))); % parameter code
        Daily_public(i,12) = table2array(PM25_parameters(merged_file_daily(i,2)+5,1)); % parameter name
        Daily_public(i,14) = table2array(PM25_parameters(merged_file_daily(i,2)+5,7)); % units
        Daily_public(i,15) = num2cell(table2array(PM25_parameters(merged_file_daily(i,2)+5,3))); % method code
        Daily_public(i,16) = table2array(PM25_parameters(merged_file_daily(i,2)+5,5)); % collection description
        Daily_public(i,17) = table2array(PM25_parameters(merged_file_daily(i,2)+5,6)); % analysis description

    end
    
    
    fileID = fopen(sprintf('%s/%s_%s_HourlyEstimatedPM25.csv',direc_PM25_out,Site_codes{loc},Site_cities{loc}),'w');
    fprintf(fileID,'File Updated: %s, %s %g \n',datestr(today), data_version,deriv_thres);
    fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',Titles_hourly{1:18});
    for h = 1:size(Hourly_public,1)
        fprintf(fileID,'%s,%s,%s,%g,%g,%g,%g,%g,%g,%g,%g,%s,%g,%s,%g,%s,%s,%s\n', char(Hourly_public(h,1)),char(Hourly_public(h,2)),char(Hourly_public(h,3)),...
            cell2mat(Hourly_public(h,4:11)),char(Hourly_public(h,12)),cell2mat(Hourly_public(h,13)),char(Hourly_public(h,14)),cell2mat(Hourly_public(h,15)),...
            char(Hourly_public(h,16)),char(Hourly_public(h,17)),char(Hourly_public(h,18)));
    end
    fclose(fileID);
    
    fileID = fopen(sprintf('%s/%s_%s_DailyEstimatedPM25.csv',direc_PM25_out,Site_codes{loc},Site_cities{loc}),'w');
    fprintf(fileID,'File Updated: %s, %s %g \n',datestr(today), data_version,deriv_thres);
    fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',Titles_daily{1:18});
    for h = 1:size(Daily_public,1)
        fprintf(fileID,'%s,%s,%s,%g,%g,%g,%g,%g,%g,%g,%g,%s,%g,%s,%g,%s,%s,%s\n', char(Daily_public(h,1)),char(Daily_public(h,2)),char(Daily_public(h,3)),...
            cell2mat(Daily_public(h,4:11)),char(Daily_public(h,12)),cell2mat(Daily_public(h,13)),char(Daily_public(h,14)),cell2mat(Daily_public(h,15)),...
            char(Daily_public(h,16)),char(Daily_public(h,17)),char(Daily_public(h,18)));
    end
    fclose(fileID);

    
    clear bspfill_idx chck_idx data_dates dates_expand day_idx filename files i k merged_file_daily merged_file_hourly names2 deriv_thres
    clear Neph_data PM25_data PM25_file unique_dates mass_scat_rec Daily_public fileID h Hourly_public Titles_daily Titles_hourly data_version
    
    fprintf('Finished processing data for %s \n', Site_cities{loc})
    
    
end



