%% %%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: to read in raw nephelometer data files provided by site
% operators (unopened) and perform QA/QC prior to creating hourly and daily
% average files. Specifically:
% 1. Find new file in the Site_Sampling/$CODE_$city/Raw/ folder
% 2. Reading and format data to desired format with [neph_file_read]
% 3. QA/QC with [neph_diag], making plots refecting data conditions
% 4. Write to output files: Public_Data/Neph_Processed/$CODE/xxxx_dailyORhourly.csv
% 5. Move rejected files (didn't pass QA/AC) to Rejected, and move accepted to Archived

% Written by: Crystal Weagle
% Created: 20 April 2018

% EDITS:
% Haihui Zhu [Nov 2023]:
%   1. Added the 'Reprocess' section function to enable reprocess archive data. 
%      NOTE: Files without a CR system are skipped for reprocessing due to the 
%      uncertainties in baseline correction 
%   2. Updated neph_file_read to make it flexible for input file format. 
%      'data_title' and 'cr_title' were added to help this.
%   3. Updated Neph_diag to include all performance checks in SOP. 
%      NOTE: SOP was updated accordingly during the iterations (v3.0 to v3.1). 
%   4. Update Neph_diag to produce figures for files with sufficient data.
%   5. Added a FileTracker to monitor and summarize condition of each neph input 
%      file - model type, whether there is a CR system, are they rejected or archived,
%      flags for quality checks. 
%   6. Update data verson to 1.1.

% Haihui Zhu [July 7 2022]:
%   1. Fix a bug when defining IN_model
%   2. Fix a bug that creating repeated lines when data for one day comes
%   from different raw files
%   3. Sew daily data for the same day into one line. (They were seperated
%   due to the bug above.)
%   4. Add a check to prevent duplicated lines when reprocessing raw data
%   5. Fix a bug causes extra columns in XX_PM10_Daily.csv

% Haihui Zhu [June 9 2021]: creat corresponded output dir when adding new
% sites. Add diary function.

% Haihui Zhu [Aug. 22 2021]: Fix a bug that replace existing data with NaN
% when the existing data and the new data have overlapping dates but
% non-overlapping hours.

close all; clear;  clc
warning('off','backtrace')
addpath('UtilityFunctions')

%% %%%%%%% USER SWITCHES %%%%%%%%%%
% Setup directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

direc_output = sprintf('%s/Public_Data/Neph_Processed',direc);
direc_outputCR = sprintf('%s/Public_Data/Neph_Processed/CR_data',direc);
TrackerFileName = sprintf('%s/File_condition_tracker.xlsx', direc_output);

% if or not reprocess all the data
Reprocess = 0; % should always be 0 unless you are 100% sure you want to reprocess - it takes a long time (~5 min/site)

% The diary function saves the processing history into an annual record
diary(sprintf('%s/Public_Data/Data_Processing_Records/%s_Neph_Record', direc, datestr(today, 'yyyy-mm')))
fprintf('%s \n', datestr(today))

%%%%%%%%%%%%%%%%%%
DataVersion = '1.1';
kappa=0.2;   % calculating estimated 'dry' neph scatter using this generic value (where 0 < kappa < 1)
RHmax=100;   % maximum RH cutoff used in calculating angstrom exponent and lidar ratios
Scatmax=2500; % upper linear range of green (532nm) scatter. Signal non-linear above 2500 Mm^-1.
RHmin=0;   % lower limit of RH
%%%%%%%%%%%%%%%%%%

site_details = readtable(sprintf('%s/Site_Sampling/Site_details.xlsx',direc),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));

% model = 1 is IN101, continuous TSP, may or may not have CR system
% model = 2 is IN102, toggles between PM2.5 and PM10 size cut, may or may not have CR system
% model = 3 is prototype (rarely used, old data only usually), continuous TSP (consider skipping)

data_title = {'DateNum','Hour','RF_PMT', 'RF_REF', 'GF_PMT', 'GF_REF', 'BF_PMT', 'BF_REF', 'RB_PMT', 'RB_REF', 'GB_PMT', 'GB_REF', 'BB_PMT', 'BB_REF',... 1-14
            'Temp' ,'Amb_Pressure', 'RH', ... 15-17
            'PMT_DARK', 'Forw_Dark_Ref', 'Back_Dark_Ref', ... 18:20
            'Back_scatt_Red', 'Back_scatt_Green', 'Back_scatt_Blue', 'Total_scatt_Red', 'Total_scatt_Green', 'Total_scatt_Blue', ... 21-26
            'Fan_RPM', 'Flow_Pressure', 'Flow', 'DAC'}; % 27:30

cr_title = {'RF_CR', 'GF_CR', 'BF_CR', 'RB_CR', 'GB_CR', 'BB_CR' };

for loc =  1:length(Site_codes)

    raw_dir = sprintf ('%s/Site_Sampling/%s_%s/Nephelometer_Data/Raw',direc,Site_codes{loc},Site_cities{loc});
    archive_dir = sprintf('%s/Site_Sampling/%s_%s/Nephelometer_Data/Archived', direc, Site_codes{loc}, Site_cities{loc});
    fig_dir = sprintf('%s/Site_Sampling/%s_%s/Nephelometer_Data/Plots_Data_Checking', direc, Site_codes{loc}, Site_cities{loc});
    reject_dir = sprintf('%s/Site_Sampling/%s_%s/Nephelometer_Data/Rejected', direc, Site_codes{loc}, Site_cities{loc});

    if Reprocess == 1 
        arclist = movetoraw(archive_dir, raw_dir);
        rejlist = movetoraw(reject_dir, raw_dir);
        fprintf('Moved archived and rejected raw data to raw for reprocessing\n')

        names = dir(sprintf('%s/*', fig_dir));
        files = {names.name}; clear names
        files = files';
        if ~isempty(files)
            for ii = 1:length(files)
                if files{ii}(1) ~= '.'
                    delete(sprintf('%s/%s', fig_dir, files{ii}))
                end
            end
        end
        fprintf('Removed diagnosis figures\n')

        names = dir(sprintf('%s/*', direc_outputCR));
        files = {names.name}; clear names
        files = files';
        if ~isempty(files)
            for ii = 1:length(files)
                if files{ii}(1) ~= '.'
                    delete(sprintf('%s/%s', direc_outputCR, files{ii}))
                end
            end
        end

        % remove file records in File_condition_tracker.xlsx
        FileSum = readtable(TrackerFileName, 'Sheet', 'List');
        SiteCode = FileSum.SiteCode;
        if ~isempty(SiteCode)
            Ind = find(contains(SiteCode, Site_codes{loc}));
            FileSum(Ind, :) = [];
            delete_sheet(TrackerFileName, 'List')
            writetable(FileSum, TrackerFileName, 'Sheet', 'List')
        end
        clear FileSum SiteCode
    end

    %% -------- Search for new Neph files ---------------------------------
    names = dir(sprintf('%s/IN*',raw_dir));
    files = {names.name}; clear names
    files = files';
    IsModel3 = 0;

    if isempty(files) ==1
        names = dir(sprintf('%s/LOG*',raw_dir)); 
        files = {names.name}; clear names
        files = files';
        IsModel3 = 1; % as of 2023 Nov, LOG* files are coming from model 3. It is unlikely that we will need to reprocess them.  
    end

    if isempty(files) == 1 % if empty then skip this site
        fprintf('\nNo data files found for %s \n', Site_cities{loc})
        continue
    end

    %----------- If not empty, aggregate data -----------------------------
    fprintf('\nReading data files for %s \n', Site_cities{loc})

    for Neph=1:length(files)

        %% %%%%% Step1: Reading data files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filename = sprintf('%s/%s', raw_dir, files{Neph});
        % This section gets info for what type (e.g IN model, CR system) of file is being read.
        % could be blended to the neph_file_read2 function (Haihui)
        opts = detectImportOptions(filename);
        opts.VariableNamesLine = 1;
        opts.VariableNamingRule = 'preserve';
        rawtable = readtable(filename, opts);
        headers = rawtable.Properties.VariableNames;
        model_test = find(contains(headers, 'flow','IgnoreCase',true));
        CR_test = find(contains(headers, 'RF_CR'));
        clear opts

        FileTracker = readtable(TrackerFileName, 'Sheet', 'List', Range = "A1:S2");
        FileTracker.FileName = files{Neph};
        FileTracker.SiteCode = Site_codes{loc};

        if isempty(model_test) % regular fan, measuring TSP
            if IsModel3 == 1
                FileTracker.Model = 3;
            else
                FileTracker.Model = 1;
            end
        else
            FileTracker.Model = 2; % variable speed fan
        end

        if isempty(CR_test)
            FileTracker.CR = 0; % No CR system
            % we don't want to reprocess files without CR, but we will let it go through 
            % diag and make figures for it before 
        else
            FileTracker.CR = 1;
        end

        % skip small file and move it to rejected folder
        test_file=dir(filename);
        if test_file.bytes < 3000
            fprintf('Filesize too small. Moving to rejected. \n')
            [status,msg,msgid] = force_movefile(filename, reject_dir);
            FileTracker.TooSmall = 1;
            FileTracker.Reject = 1;
            writetable(FileTracker, TrackerFileName, 'Sheet', 'List', 'WriteMode', 'Append')
            continue
        end

        fprintf('\nReading %s from %s \n',files{Neph}, Site_cities{loc})


        % CR data collecting
        [data, CR_data, FileTracker] = neph_file_read(filename, FileTracker, data_title, cr_title);

        % remove file with less than 1 hour data 
        %   Sometimes, malfunction of Neph can procude hundreds of small files. Reject them 
        %   to prevent processing questionable data.
        if length(data) < 240
            warning('File %s rejected: less than one hour of data with fan on', files{Neph});
            force_movefile(filename, reject_dir);
            FileTracker.TooSmall = 1;
            FileTracker.Reject = 1;
            writetable(FileTracker, TrackerFileName, 'Sheet', 'List', 'WriteMode', 'Append')
            disp('Moving on to next file ')
            continue
        end

        FileTracker.StartDate = datetime(datevec(data(1, 1)));
        FileTracker.EndDate = datetime(datevec(data(end, 1)));

        % checks for any problems with PMT, reference sensors, and DAC/flow if IN102
        [data, FileTracker] = neph_diag(data, data_title, FileTracker, Scatmax, files{Neph}, fig_dir);

        % ------- Now decide what to do next --------------------------------------------------------
        % if with CR 
        %   Move rejected to the rejected folder. 
        %   No action for non-rejected
        if FileTracker.CR == 1
            if FileTracker.Reject == 1 % neph_diag determines if a file is rejected
             % for file with CR, the reject is valid
                force_movefile(filename, reject_dir);
                writetable(FileTracker, TrackerFileName, 'Sheet', 'List', 'WriteMode', 'Append')
                disp('Moving on to next file ')
                continue
            end
        else
        % if without CR
        %   When reprocessing: skip the file and move it to where it was (it should be in either 'arclist' or 'rejlist')
        %   When working on new files: Need to manually uncomment step 2 (following this) to do baseline correcton
            if Reprocess == 1 % if we are reprocessing old files, we want to move it back to where they were
                if sum(contains(arclist, files{Neph}))==1
                    force_movefile(filename, archive_dir);
                    FileTracker.Reject = 0; % need to change this to match what is being done
                    writetable(FileTracker, TrackerFileName, 'Sheet', 'List', 'WriteMode', 'Append')
                elseif sum(contains(rejlist, files{Neph}))==1 
                    force_movefile(filename, reject_dir);
                    writetable(FileTracker, TrackerFileName, 'Sheet', 'List', 'WriteMode', 'Append')
                else
                    error('file name error!')
                end
                fprintf('Reprocessed non-CR file. Skipped.\n')
                continue
            else % we don't expect to reach this point since we don't expect any new file without CR
                error ('%s is a new file without a CR system', files{Neph}) 
                % edit this part if we do have new data without CR 
                % also uncomment step 2 for baseline correction
            end
        end


        % sorts based on dates
        data = sortrows(data);
        good_dates = ~isnan(data(:, 1));
        data = data(good_dates, :);

        totscat_idx = find(contains(data_title, {'Back_scatt', 'Total_scatt', 'IgnoreCase', 'True'}));
%{
    % ----- not doing this anymore since new files are corrected by CR system ------

        %% %%%%% Step2: Correct baseline when CR_system = 0 %%%%%%%%%%%%%%%%%%%
        % needs to go after hourly-averages are created and gaps between files
        % filled with NaN    
        if FileTracker.CR == 1
            noCR_idx = [];
        else
            noCR_idx = [1:length(data)]';  
        end

        data0 = data; % a copy of data before baseline correction

        if ~isempty(noCR_idx) 

            chunk = 24000; % 100-hour chunks of time, it's important there are no large gaps in data or this does weird things (should only matter for reprocessing entire data sets from beginning to end)
            Bsp_chunks = ceil(size(noCR_idx,1)/chunk);
            Bsp_corr=msbackadj(data(noCR_idx,1),data(noCR_idx,totscat_idx),'StepSize',50);

            % removing negative points, adjusting the min in increment stages, making sure baseline values are at least 10 Mm-1 (N2 background)
            for CHK=1:Bsp_chunks
                Bsp_i=chunk*(CHK-1)+1;
                Bsp_f=min(chunk*(CHK-1)+chunk,size(noCR_idx,1));
                Min_bsp=min(Bsp_corr(Bsp_i:Bsp_f,:),[],1);
                for q=1:6
                    Bsp_corr(Bsp_i:Bsp_f,q)=Bsp_corr(Bsp_i:Bsp_f,q)-Min_bsp(:,q)+10;
                end

            end

            data(noCR_idx,totscat_idx)=Bsp_corr; %overwriting old values with baseline-corrected

            clear chunk Bsp_chunks Bsp_corr q
            % make a plot compare before and after corrections (Haihui)
            baseline_corr_plot(data0(:, totscat_idx), data(:, totscat_idx), files{Neph}, Site_codes, Site_cities, loc, fig_dir)
        end
%}
        %% %%%%% Step 2: Check for negative/infinite scatter values or negative RH %%%%%%%%%
        % if there is a "too" negative for green total scatter, set all scatter to NaN
        for ii = 1:length(totscat_idx)
            data( data(:,totscat_idx(ii)) < -2 , totscat_idx ) = NaN;
        end

        data(isinf(data)) = NaN;

        % Nan out negative RH
        % RH as low as -30 were spotted at the Melbourne site. Below are response from Chris:  
        % -30 likely means that power to the RH probe or signal from the RH probe was lost. 
        % Other aspects of the nephelometer probably still operated.
        % Those -30 values are likely better expressed as NaN.
        rh_col = find(ismember(data_title,'RH'));
        data(data(:,rh_col)<0,rh_col) = nan;


        %% %%%%% Step3: writing to output csv files %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % making site dir if doesn't exist
        SiteOutdir = sprintf('%s/%s', direc_output, Site_codes{loc});

        if ~exist(SiteOutdir, 'dir')
            mkdir(SiteOutdir)
        end

        time_now = datestr(today);

        Titles_hourly = {'Year' 'Month' 'Day' 'Hour' 'Temp_degC' 'Pressure' '%RH' 'Back_Scat_Red' 'Back_Scat_Green' ...
                         'Back_Scat_Blue' 'Total_Scat_Red' 'Total_Scat_Green' 'Total_Scat_Blue'};

        Titles_daily = {'Year' 'Month' 'Day' 'num_hourly_points' 'Temp_degC' 'Pressure' '%RH' 'Back_Scat_Red' 'Back_Scat_Green' ...
                        'Back_Scat_Blue' 'Total_Scat_Red' 'Total_Scat_Green' 'Total_Scat_Blue'};
                        
        % find the columns and rows to be printed out
        amb_cond_idx = find(ismember(data_title, {'DateNum', 'Hour', 'Temp', 'Amb_Pressure', 'RH'}));
        column_index = [amb_cond_idx totscat_idx]; % length = 11 

        % ------------------------------- TSP -------------------------------------
        if FileTracker.Model == 1 || FileTracker.Model == 3
            % ----------- regriding data to hourly and daily mean -----------------
            data = data(:, column_index);
            data_days = unique(data(:,1),'stable'); % this is a list of the dates(serial format) with data for looping
            data_hourly = NaN((24*length(data_days)), length(column_index));
            data_daily = NaN((length(data_days)), length(column_index));

            j = -23:0; j = j';
            for i = 1:length(data_days) % [daily mean]
                j = j+ 24; % j will keep increasing until 24*length(data_days)

                date_index = find(data(:,1) == data_days(i));
                data_hourly(j,1) = data_days(i);
                data_daily(i,1) = data_days(i);

                for k = 0:23 %[hourly mean]
                    hour_indexB = find(data(date_index,2) == k);
                    hour_index = date_index(hour_indexB); % column 1
                    data_hourly(j(k+1),2) = k; % column 2

    %                 ang_exp1 = log(data(hour_index,10)./data(hour_index,11))/log(457/532);
    %                 ang_exp2 = log(data(hour_index,11)./data(hour_index,12))/log(532/634);
    %                 ang_exp = nanmean([ang_exp1 ang_exp2],2);
    %                 bsp_550 =data(hour_index,11).*((550/532).^ang_exp);

                    data_hourly(j(k+1),3:end) = nanmean(data(hour_index,3:end)); % columns 3 - 12
    %                 data_hourly(j(k+1),end-2) = nanmean(bsp_550); % column 13
    %                 data_hourly(j(k+1),end-1) = nanmean(ang_exp);% column 14

                    clear hour_index hour_indexB ang_exp1 ang_exp2 ang_exp bsp_550
                end

                num_points = 24 - sum(isnan(data_hourly(j,11))); % determines number of hourly data points per day based on total green scatter

                data_daily(i,2) = num_points;
                data_daily(i,3:end) = nanmean(data_hourly(j,3:end));

                clear date_index  num_points
            end

            clear data_days data data column_index

            % ----------------------------- writing output -------------------------------
            % ------ Hourly TSP --------
            fname = sprintf('%s/%s/%s_%s_NephTSP_hourly.csv', direc_output, Site_codes{loc}, Site_codes{loc}, Site_cities{loc});
            datatype = 'hourly';
            write_scatter(fname, data_hourly, Titles_hourly, datatype, RHmax, Scatmax, time_now, DataVersion)

            % -------- Daily TSP ---------
            fname = sprintf('%s/%s/%s_%s_NephTSP_daily.csv', direc_output, Site_codes{loc}, Site_codes{loc}, Site_cities{loc});
            datatype = 'daily';
            write_scatter(fname, data_daily, Titles_daily, datatype, RHmax, Scatmax, time_now, DataVersion)

            % -------- CR ---------
            % always save CR in a new file named with start and end dates  
            if FileTracker.CR == 1 && ~isempty(CR_data)
                write_cr (CR_data, direc_outputCR, Site_codes, loc, time_now, DataVersion)
            end

            fprintf('Saved hourly and daily average files for %s \n \n', Site_cities{loc})

        end
        
        % ---------------------------------- PM2.5 and PM10 ----------------------------------
        if FileTracker.Model == 2
            % ----------- regriding data to hourly and daily mean -----------------
            data_PM10 = data(data(:,29) <= 2.0,column_index);
            data_PM25 = data(data(:,29) >= 4.5,column_index);

            data_days = unique(data(:,1),'stable'); % this is a list of the dates(serial format)  with data for looping
            data_PM25_hourly = NaN((24*length(data_days)), length(column_index));
            data_PM10_hourly = NaN((24*length(data_days)), length(column_index));
            data_PM25_daily = NaN((length(data_days)), length(column_index));
            data_PM10_daily = NaN((length(data_days)), length(column_index));

            j = -23:0; j = j';
            for i = 1:length(data_days) % loop through days
                j = j+ 24;

                date_index_PM25 = find(data_PM25(:,1) == data_days(i));
                data_PM25_hourly(j,1) = data_days(i);
                data_PM25_daily(i,1) = data_days(i);

                date_index_PM10 = find(data_PM10(:,1) == data_days(i));
                data_PM10_hourly(j,1) = data_days(i);
                data_PM10_daily(i,1) = data_days(i);

                for k = 0:23 % loop through hours
                    %%%% PM25 %%%%
                    hour_index_PM25B = find(data_PM25(date_index_PM25,2) == k);
                    hour_index_PM25 = date_index_PM25(hour_index_PM25B);
                    data_PM25_hourly(j(k+1),2) = k; % column 2

        %             fRH_PM25=1+kappa*data_PM25(hour_index_PM25,5)./(100-data_PM25(hour_index_PM25,5));
        %             ang_exp1_PM25 = log(data_PM25(hour_index_PM25,10)./data_PM25(hour_index_PM25,11))/log(457/532);
        %             ang_exp2_PM25 = log(data_PM25(hour_index_PM25,11)./data_PM25(hour_index_PM25,12))/log(532/634);
        %             ang_exp_PM25 = nanmean([ang_exp1_PM25 ang_exp2_PM25],2);
        %             bsp_550_PM25 =data_PM25(hour_index_PM25,11).*((550/532).^ang_exp_PM25);

                    data_PM25_hourly(j(k+1),3:end) = mean(data_PM25(hour_index_PM25,3:end),'omitnan'); % columns 3 - 12
        %             data_PM25_hourly(j(k+1),end-2) = nanmean(bsp_550_PM25); % column 13
        %             data_PM25_hourly(j(k+1),end-1) = nanmean(ang_exp_PM25);% column 14
        %             data_PM25_hourly(j(k+1),end) = nanmean(fRH_PM25); % column 13

                    %%%% PM10 %%%
                    hour_index_PM10B = find(data_PM10(date_index_PM10,2) == k);
                    hour_index_PM10 = date_index_PM10(hour_index_PM10B);
                    data_PM10_hourly(j(k+1),2) = k;

        %             fRH_PM10=1+kappa*data_PM10(hour_index_PM10,5)./(100-data_PM10(hour_index_PM10,5));
        %             ang_exp1_PM10 = log(data_PM10(hour_index_PM10,10)./data_PM10(hour_index_PM10,11))/log(457/532);
        %             ang_exp2_PM10 = log(data_PM10(hour_index_PM10,11)./data_PM10(hour_index_PM10,12))/log(532/634);
        %             ang_exp_PM10 = nanmean([ang_exp1_PM10 ang_exp2_PM10],2);
        %             bsp_550_PM10 =data_PM10(hour_index_PM10,11).*((550/532).^ang_exp_PM10);

                    data_PM10_hourly(j(k+1),3:end) = mean(data_PM10(hour_index_PM10,3:end),'omitnan'); % columns 3-12
        %             data_PM10_hourly(j(k+1),end-2) = nanmean(bsp_550_PM10); % column 13
        %             data_PM10_hourly(j(k+1),end-1) = nanmean(ang_exp_PM10);% column 14
        %             data_PM10_hourly(j(k+1),end) = nanmean(fRH_PM10); % column 13

                    clear hour_index_PM25B hour_index_PM25 hour_index_PM10 hour_index_PM10B fRH_PM25 fRH_PM10 ang_exp1_PM25 ang_exp2_PM25 ang_exp_PM25 fRH_PM10 ang_exp1_PM10 ang_exp2_PM10 ang_exp_PM10 bsp_550_PM25 bsp_550_PM10
                end

                num_PM25points = 24 - sum(isnan(data_PM25_hourly(j,11))); % determines number of hourly data points per day based on total green scatter
                num_PM10points = 24 - sum(isnan(data_PM10_hourly(j,11))); % determines number of hourly data points per day based on total green scatter

                data_PM25_daily(i,2) = num_PM25points;
                data_PM25_daily(i,3:end) = mean(data_PM25_hourly(j,3:end),'omitnan');

                data_PM10_daily(i,2) = num_PM10points;
                data_PM10_daily(i,3:end) = mean(data_PM10_hourly(j,3:end),'omitnan');

                clear date_index_PM25 date_index_PM10  num_PM25points num_PM10points

            end

            clear data_days data_PM25 data_PM10 data column_index

            % remove 'empty' rows: when processing PM10 and PM2.5 together,
            % there will be sampling dates for one dataset that contain no
            % data. Need to remove those rows to prevent rows with full 'NaN'
            data_PM25_hourly(isnan(data_PM25_hourly(:,11)), :) = [];
            data_PM10_hourly(isnan(data_PM10_hourly(:,11)), :) = [];
            data_PM25_daily(isnan(data_PM25_daily(:, 11)), :) = [];
            data_PM10_daily(isnan(data_PM10_daily(:, 11)), :) = [];

            % ----------------------------- writing output -------------------------------
            % Hourly PM25
            fname = sprintf('%s/%s/%s_%s_NephPM25_hourly.csv', direc_output, Site_codes{loc}, Site_codes{loc}, Site_cities{loc});
            datatype = 'hourly';
            write_scatter(fname, data_PM25_hourly, Titles_hourly, datatype, RHmax, Scatmax, time_now, DataVersion)
            
            % Daily PM25
            fname = sprintf('%s/%s/%s_%s_NephPM25_daily.csv', direc_output, Site_codes{loc}, Site_codes{loc}, Site_cities{loc});
            datatype = 'daily';
            write_scatter(fname, data_PM25_daily, Titles_daily, datatype, RHmax, Scatmax, time_now, DataVersion)

            % Hourly PM10
            fname = sprintf('%s/%s/%s_%s_NephPM10_hourly.csv',direc_output,Site_codes{loc},Site_codes{loc},Site_cities{loc});
            datatype = 'hourly';
            write_scatter(fname, data_PM10_hourly, Titles_hourly, datatype, RHmax, Scatmax, time_now, DataVersion)

            % Daily PM10
            fname = sprintf('%s/%s/%s_%s_NephPM10_daily.csv', direc_output, Site_codes{loc}, Site_codes{loc}, Site_cities{loc});
            datatype = 'daily';
            write_scatter(fname, data_PM10_daily, Titles_daily, datatype, RHmax, Scatmax, time_now, DataVersion)
            %%%%%  CR data  %%%%%%
            if FileTracker.CR == 1 && ~isempty(CR_data)
                write_cr (CR_data, direc_outputCR, Site_codes, loc, time_now, DataVersion)
            end
        end

        
        %% %%%%% Step4: move raw file to archive %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if FileTracker.Reject == 0
            fprintf('Moving %s to archived\n', files{Neph})
            [status, msg, msgid] = force_movefile(sprintf('%s/%s', raw_dir, files{Neph}), archive_dir);
        end
        writetable(FileTracker, TrackerFileName, 'Sheet', 'List', 'WriteMode', 'Append')
                
        clear FileTracker CR_data  CR_test data_PM10_daily data_PM10_hourly data_PM25_daily data_PM25_hourly
        clear file_size fileID filename  headers  j k   model_test noCR_idx test_file time_now  
    end 
    clear Neph  Titles_daily Title_hourly files IsModel3
        
    %% ------ Update Nephelometer File Summary -------------------------------------------
    FileTracker = readtable(TrackerFileName, 'Sheet', 'List');
    FileTracker(ismember(FileTracker.FileName, 'TEMPLATE'),:) = [];
    SiteCode = unique(FileTracker.SiteCode);
    AllSiteCodes = FileTracker.SiteCode;
    headers = FileTracker.Properties.VariableNames;
    headers = headers(5:end);
    FileTrackerMat = table2array(FileTracker(:, 5:end));

    for loc = 1:length(SiteCode)
        ind = find(contains(AllSiteCodes, SiteCode{loc}));
        RecNum = length(ind); % number of records for this site
        CR_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'CR'))) ./ RecNum;
        Reject_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'Reject'))) ./ RecNum;
        TooSmall_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'TooSmall'))) ./ RecNum;
        FanOff25Perc_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'FanOff25Perc'))) ./ RecNum;
        PoorFlow_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'PoorFlow'))) ./ RecNum;
        UnstableRef_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'UnstableRef')) > 0) ./ RecNum;
        LowRef_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'LowRef')) > 0) ./ RecNum;
        DarkTooHigh_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'DarkTooHigh'))) ./ RecNum;
        UnstableDark_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'UnstableDark'))) ./ RecNum;
        LowPMT_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'LowPMT'))) ./ RecNum;
        Insensitive_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'Insensitive')) > 0) ./ RecNum;
        TempCorrDark_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'TempCorrDark')) > 0) ./ RecNum;
        PressCorrDark_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'PressCorrDark')) > 0) ./ RecNum;
        RHCorrDark_rate(loc, 1) = 100 * sum(FileTrackerMat(ind, ismember(headers, 'RHCorrDark')) > 0) ./ RecNum;
    end

    SiteCode{end + 1} = 'Total';
    loc = loc + 1;
    RecNum = length(FileTracker.CR); % number of records for this site
    CR_rate(loc, 1) = 100 * sum(FileTracker.CR) ./ RecNum;
    Reject_rate(loc, 1) = 100 * sum(FileTracker.Reject) ./ RecNum;
    TooSmall_rate(loc, 1) = 100 * sum(FileTracker.TooSmall) ./ RecNum;
    FanOff25Perc_rate(loc, 1) = 100 * sum(FileTracker.FanOff25Perc) ./ RecNum;
    PoorFlow_rate(loc, 1) = 100 * sum(FileTracker.PoorFlow) ./ RecNum;
    UnstableRef_rate(loc, 1) = 100 * sum(FileTracker.UnstableRef > 0) ./ RecNum;
    LowRef_rate(loc, 1) = 100 * sum(FileTracker.LowRef > 0) ./ RecNum;
    DarkTooHigh_rate(loc, 1) = 100 * sum(FileTracker.DarkTooHigh) ./ RecNum;
    UnstableDark_rate(loc, 1) = 100 * sum(FileTracker.UnstableDark) ./ RecNum;
    LowPMT_rate(loc, 1) = 100 * sum(FileTracker.LowPMT) ./ RecNum;
    Insensitive_rate(loc, 1) = 100 * sum(FileTracker.Insensitive > 0) ./ RecNum;
    TempCorrDark_rate(loc, 1) = 100 * sum(FileTracker.TempCorrDark > 0) ./ RecNum;
    PressCorrDark_rate(loc, 1) = 100 * sum(FileTracker.PressCorrDark > 0) ./ RecNum;
    RHCorrDark_rate(loc, 1) = 100 * sum(FileTracker.RHCorrDark > 0) ./ RecNum;

    Summary = table(SiteCode, CR_rate, Reject_rate, TooSmall_rate, FanOff25Perc_rate, PoorFlow_rate, UnstableRef_rate, LowRef_rate, DarkTooHigh_rate, ...
        UnstableDark_rate, LowPMT_rate, Insensitive_rate, TempCorrDark_rate, PressCorrDark_rate, RHCorrDark_rate);
    writetable(Summary, TrackerFileName, 'Sheet', 'Summary')
    clear Summary *_rate RecNum loc FileTrackerMat FileTracker headers


end


diary off

%% FUNCTION 
function files = movetoraw(archive_dir, rawdir)
    names = dir(sprintf('%s/IN*', archive_dir));
    files = {names.name}; clear names
    files = files';

    if ~isempty(files)

        for Neph = 1:length(files)
            filename = sprintf('%s/%s', archive_dir, files{Neph});
            [status, msg, msgid] = force_movefile(filename, rawdir);
        end

    end

    % LOG* are not reprocessed since they are model 3 (TSP) data
    % names = dir(sprintf('%s/LOG*', archive_dir));
    % files = {names.name}; clear names
    % files = files';

    % if ~isempty(files)

    %     for Neph = 1:length(files)
    %         filename = sprintf('%s/%s', archive_dir, files{Neph});
    %         [status, msg, msgid] = force_movefile(filename, rawdir);
    %     end

    % end

end

function write_cr (CR_data, direc_outputCR, Site_codes, loc, time_now, DataVersion)

    Dates_CR = datevec(CR_data(:, 1));
    start_date = Dates_CR(1, 1) * 100 + Dates_CR(1, 2);
    end_date = Dates_CR(end, 1) * 100 + Dates_CR(end, 2);

    Merged_CR = [Dates_CR(:, 1:3) CR_data(:, 2:end)];

    Titles_CR = {'Year' 'Month' 'Day' 'RF_CR' 'GF_CR' 'BF_CR' 'RB_CR' 'GB_CR' 'BB_CR'};

    fileID = fopen(sprintf('%s/%s_%06d_to_%06d_NephCRdata.csv', direc_outputCR, Site_codes{loc}, start_date, end_date), 'w');
    fprintf(fileID, 'Created on %s. Data Version = %s\n', time_now, DataVersion);
    fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s\n', Titles_CR{1, 1:9});
    fprintf(fileID, '%04d, %02d, %02d, %6.2f,%6.2f,%6.2f, %6.2f,%6.2f,%6.2f\n', Merged_CR.');
    fclose(fileID);

end

function baseline_corr_plot(data0, data, file, Site_codes, Site_cities, loc, fig_dir)
    Colors = {'r', 'g', 'b', 'm', '#76ab93', 'c'};
    figure('Position',[10 10 1100 600])

    subplot(2,1,1)
    for ii = 1:3
        plot(1:length(data0), data0(:,ii)',':','color',Colors{ii});
        hold on
        plot(1:length(data), data(:, ii)', '-', 'color', Colors{ii});
        hold on
    end

    % xlabel('Data point # ', 'FontWeight', 'bold', 'FontSize', 10);

    % set(gcf, 'outerposition', [520 260 1000 550]);
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    xtickangle(45)
    ylabel('Backward Signal', 'FontWeight', 'bold', 'FontSize', 10);

    legend({'Before','After'}, 'Location', 'best', 'FontWeight', 'bold', 'FontSize', 10)

    title('Baseline correction - Backward', 'FontSize', 12)

    subplot(2, 1, 2)
    for ii = 4:6
        plot(1:length(data0), data0(:, ii)', ':', 'color', Colors{ii});
        hold on
        plot(1:length(data), data(:, ii)', '-', 'color', Colors{ii});
        hold on
    end

    xlabel('Data point # ', 'FontWeight', 'bold', 'FontSize', 10);

    % set(gcf, 'outerposition', [520 260 1000 550]);
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    xtickangle(45)
    ylabel('Forward Signal', 'FontWeight', 'bold', 'FontSize', 10);

    legend({'Before', 'After'}, 'Location', 'best', 'FontWeight', 'bold', 'FontSize', 10)

    title('Baseline correction - Forward', 'FontSize', 12)

    % save figure
    filename = file(1:end-4);
    sfname = sprintf('%s/%s_baseline_corr.png', fig_dir, filename);
    saveas(gcf,sfname)
end

function write_scatter(fname, data, Titles, datatype, RHmax, Scatmax, time_now, DataVersion)

    Dates = datevec(data(:, 1));

    if exist(fname, 'file') % checks if a filter masses file exists
        % read existing neph data
        Neph_exist = csvread(fname, 2, 0);
        Neph_exist(isnan(Neph_exist(:, 8)), :) = []; % remove NaN rows
        
        % format the new data
        data_new = [Dates(:, 1:3) data(:, 2:end)];

        % SEWING OVERLAPING DATES [Haihui 2021-08-22]
        date_num_new = data(:, 1);
        date_num_existing = datenum(Neph_exist(:, 1), Neph_exist(:, 2), Neph_exist(:, 3), 0, 0, 0);
        overlap_dates = find(ismember(date_num_existing, date_num_new));
        if ~isempty(overlap_dates)
            switch datatype
                case 'hourly'
                    if length(overlap_dates) <= 25 % one hour is splited into two Neph files
                        for rp = 1:length(overlap_dates) 
                            Ind = find(ismember(date_num_new, date_num_existing(overlap_dates(rp))) & ismember(data_new(:, 4), Neph_exist(overlap_dates(rp), 4)));
                            if Neph_exist(overlap_dates(rp), 5) ~= data_new(Ind, 5) % the split hour
                                data_new(Ind, 5:end) = 0.5 * (Neph_exist(overlap_dates(rp), 5:end) + data_new(Ind, 5:end));
                                % assuming the average is representative, since it is difficult figure out how to weight them
                            % else 
                                % other hours: the new and old data are the same, no action needed.
                            end
                        end
                    % else 
                        % reprocessing old files, simply remove the value in the existing one.
                    end
                    Neph_exist(overlap_dates, :) = []; % remove overlaped dates in the existing data since they are also availale in the data_new 

                case 'daily' 
                    if length(overlap_dates) == 1 % one of the hours splited to two raw files.
                        Ind = find(ismember(date_num_new, date_num_existing(overlap_dates)));
                        tot_hour = data_new(Ind, 4) + Neph_exist(overlap_dates, 4); % total hour = sum of both
                        data_new(Ind, 5:end) = ((data_new(Ind, 5:end) * data_new(Ind, 4) + ...
                            Neph_exist(overlap_dates, 5:end) * Neph_exist(overlap_dates, 4))) ./ ...
                            tot_hour; % weighted mean
                        if tot_hour == 25
                            data_new(Ind, 4) = 24; % Ignore the repeated hour
                        elseif tot_hour <= 24
                            data_new(Ind, 4) = tot_hour; 
                        end
                    % elseif length(overlap_dates) > 1 % reprocessing old files, simply remove the value in the existing one.
                    end
                    Neph_exist(overlap_dates, :) = [];
            end
        end

        Merged_hourly = [Neph_exist; data_new];
        clear data_new Neph_exist
    else
        Merged_hourly = [Dates(:, 1:3) data(:, 2:end)];
    end
    
    % sort data chronologically
    dates_sorting = datenum(Merged_hourly(:, 1), Merged_hourly(:, 2), Merged_hourly(:, 3), Merged_hourly(:, 4), 0, 0);
    [~, index] = sortrows(dates_sorting);
    Merged_hourly = Merged_hourly(index, :);

    fileID = fopen(fname, 'w');
    fprintf(fileID, 'Updated %s. Data version %s, RHmax=%d,Scatmax=%d \n', time_now, DataVersion,RHmax, Scatmax);
    fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', Titles{1, 1:13});
    fprintf(fileID, '%04d, %02d, %02d, %02d, %6.1f, %6.2f, %6.0f, %6.1f, %6.1f, %6.1f,%6.1f,%6.1f, %6.1f\n', Merged_hourly.');
    fclose(fileID);
    fprintf('Data saved to: %s\n', fname)
end

