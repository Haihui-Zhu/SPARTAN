%% %%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: read and perform QA/QC on BC data as measured by the SSR in
% triplicate and perform QA/QC on the data

% Written by: Crystal Weagle
% Created: 21 November 2018

% EDITS:
% 2019-01-13, Crystal Weagle: added reading of cartridge number (i.e. XXXX-###) and Lot
% number (i.e. L#####) from filter masses spreadsheets. Lot number is
% important because if a filter has a lot number it is a stretched teflon
% and needs to be flagged since the value of q below may change (waiting on comparison with HIPS)

% 2019-09-27, Crystal Weagle: added loop to correct new stretched teflon BC
% based on the measured reflectance of blank (usually between 75 and 85).

% 2020-10-07, Emmie Le Roy: fixed bug where BC values were not being
% filtered properly; also corrected the new stretched teflon correction loop
% for edge case where cartridges don't have a blank (because the mass_type was set
% to 6 instead of 0 due to flow issues). In this case the correction is
% based on the average surface reflectance of 80, rather than the
% corresponding field blank

% 2021-1-21, Chris Oxford: Repaired the errors created by xlsread. xlsread
% was not placing the numbers in the proper row,column. Changed to
% readcell and made variable conversions where required. Also added many variables
% (deleted one) to clear commands prior to the continue commands throughout
% the document.

% 2021-11-06, Haihui Zhu: adding ReadMaster and WriteToMaster functions

% 2023-01-26, Haihui Zhu:
% 1. Change the structure so that SSR BC is read cartridge by cartridge
% (instead of chunk by chunk).
% 2. Add a 'clean up master' section to prepare master files for
% reprocessing all variables mentioned above.

close all; clear ; clc
addpath('UtilityFunctions')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% USER SWITCHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Clean_Up_Master = 0;

% Set up directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

direc_master = strcat(direc,'/Analysis_Data/Master_files');
direc_SSR = strcat(direc,'/Analysis_Data/Black_Carbon/SSR_by_site');
direc_sampling = strcat(direc,'/Site_Sampling');

% The diary function saves the processing history into an annual record
diary(sprintf('%s/Public_Data/Data_Processing_Records/%s_BC_SSR_Record',direc,datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SSR_threshold = [20 90]; % averages above 90 or below 20 are outside linear range to interpret into BC mass
SSR_std = 1.5; % the standard deviation between triplicate measurements must be 1.5 or less
SA = 3.142;  % 25 mm filter surface area in cm2
ABS_coeff = 0.06; % default absorption cross section in cm2/ug (Snider et al. 2016, ACP)
R_100 = 100; % initial reflectance (blank filter reflectance)
R_avg_stretch = 85;
q = 1/1.5; % accounts for thickness of mesh teflon filter - * needs to be updated for stretched teflon, waiting for HIPS

%-------------   SITE DETAILS   --------------
site_details = readtable(strcat(direc_sampling,'/Site_details.xlsx'));
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));
R_toolow_table = zeros(size(Site_codes));
%% --------- Find and read data file in SSR and Master directories ---------

% for loc = 25(Site_codes)
for loc = 1:length(Site_codes)

    % Check if SSR BC file exists
    if exist(sprintf('%s/%s_SSR.xlsx',direc_SSR,Site_codes{loc})) == 0
        fprintf('WARNING: No SSR file exists for %s \n', Site_codes{loc})
        % If no SSR BC measurement file exists, moves to the next site
        continue
    end

    % Find if there is a master file
    master_file = sprintf('%s/%s_master.csv',direc_master,Site_codes{loc});
    [Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, projectIDs_master,Master_hours, Master_masstype, ...
        Master_dates, Master_mass, Master_IC, Master_ICP, Master_XRF,...
        Master_carbon, Master_Method, Master_flags] = ReadMaster(master_file,Site_codes{loc});

    if isempty(Master_mass)
        continue
    else

        %%%%%%% Clean Up Master %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Clean up the dates, hours, volumes, masstype in master file to prepare
        % it for reprocessing
        if Clean_Up_Master ==1
            Master_carbon(:,1) = NaN; % set SSR BC to NaN.
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [Haihui, 2023, 1, 26 ]

        % finds Master data lines that do not have a BC measurement based on NaN
        index = find(isnan(Master_carbon(:,1)) == 1);

        if isempty(index) == 1
            fprintf('All samples in %s master file have BC values \n', Site_codes{loc})
            disp('Moving on to next site')
            clear CartridgeIDs_master  LotIDs_master Master_*  projectIDs_master Master_mass
            continue
        end

        % Read data from SSR file
        filename = sprintf('%s/%s_SSR.xlsx',direc_SSR,Site_codes{loc});
        fprintf('Reading %s \n', filename);
        %[num,txt,raw] = xlsread(filename);
        txtandnum=readcell(filename);
        txt=string(txtandnum);
        num=txtandnum;
        % Create analysis IDs for consistency with chemical analysis labeling
        [I_id, J_id] = findIndexContaining('Sample ID', txt);
        labels = txt(I_id+1:end,J_id);

        % Format filter ID in SSR sheet
        analysis_IDs = cell(length(labels),1);
        for i = 1:length(labels)
            label = char(deblank(labels(i,:))); % the deblank function removes any spaces in the string
            if size(label,2)==13 % old labeling system (e.g. 15001-KRSE-1T)
                tlabel(1,1:4) = char(Site_codes{loc});
                tlabel(1,5:6) = '-0';
                tlabel(1,7:9) = label(3:5);
            elseif size(label,2) == 11 % new labeling system (e.g. KRSE-0001-1)
                tlabel(1,1:4) = char(Site_codes{loc});
                tlabel(1,5) = '-';
                tlabel(1,6:9) = label(6:9);
            else
                disp('WARNING: filterID label not size 13 or 11, check for mis-labeling');
            end
            analysis_IDs{i} = tlabel;
            clear lab label
        end

        % Read all the SSR BC average and STD
        [I,J] = ind2sub(size(txt),find(contains(txt,'Mean Reading'))); % J is the col
        foof=num(:,J);% this section added so that xlsread command, which is causing errors due to assumptions. Chris Oxford 1/25/2021
        testchar='NaN'; % for testing for NaN
        k=0; % k is the number of deleted rows in used in the temporary variable foof.
        for i=1:length(foof(:,1)) % foof should contain a header that needs to be removed.
            if (ischar(foof{i,1}) || ismissing(foof{i,1})) && ~(string(foof{i,1})==string(testchar)) % does this row need to be deleted?
                foof{i,:}=[]; % delete row;
                k=k+1; % add 1 to the count
            end
            if (string(foof{i,1})==string(testchar)) % does this row contain the string 'NaN'?
                foof{i,1}=str2num(foof{i,1}); % change the string to a number.
            end
        end
        % start_row = find(~isnan(num(:,J)));
        % R_avgB = num(start_row(1):end, J); % select cells with mean readings
        R_avgB=cell2mat(foof);
        clear I J


        [I,J] = ind2sub(size(txt),find(contains(txt,'Std.'))); %  J is the col. This is the original command
        foof=num(:,J);
        for i=1:k % for every row deleted above. This line could throw an error as it assumes that only the first rows need to be deleted. Do not put strings in the mean reading column.
            foof{i,:}=[]; % remove same rows that were removed in mean reading.
        end
        for i=1:length(foof(:,1))
            if (string(foof{i,1})==string(testchar)) % is the cell equal to the string 'NaN'?
                foof{i,1}=str2num(foof{i,1}); % convert the string to a number
            end
        end
        R_stdB=cell2mat(foof); % this assumes that every average has a standard deviation
        % end of this section added so that xlsread command, which is causing errors.
        %R_stdB = num(start_row(1):end, J); % select cells with standard deviation of triplicate readings
        clear I J start_row

        % sanity check
        if length(analysis_IDs) ~= length(R_avgB)
            error('analysis_IDs does not match R_avgB')
        end

        % Now looping through all cartridges that need to add SSR
        cartridge_IDs = unique(CartridgeIDs_master(index));

        for ii = 1:length(cartridge_IDs)

            % info needed from master
            Master_ind = ismember(CartridgeIDs_master,cartridge_IDs{ii});
            tLot_IDs  = cellstr(LotIDs_master(Master_ind));
            mass_type =       Master_masstype(Master_ind);
            tFilters  =            Master_IDs(Master_ind);
            % info needed from SSR
            SSR_ind      = find(ismember(analysis_IDs,tFilters));
            if ~isempty(SSR_ind)
                tFilters_SSR = analysis_IDs(SSR_ind); % tFilters_SSR should be the same as tFilters, but in case it is not, compare them when necessary
                R_avg        = R_avgB(SSR_ind); % standard deviation of triplicate R for this cartridge
                R_std        = R_stdB(SSR_ind); % standard deviation of triplicate R for this cartridge

                % ====== Calculate BC (ug) and perform QA on data =============
                % QA includes checking that the standard deviation of triplicate SSR
                % measurements is 1.5 or less, and that the average of triplicate
                % measurements are between 20 and 90 for linear conversion from
                % reflectance to BC using ABS_coeff = 0.075 cm2/ug

                BC_standard = ((q*SA)/ABS_coeff)*log(R_100./R_avg);

                % ------ blank correction for stretched teflon filters ------
                if strcmp(tLot_IDs(1),"Unknown") == 0 % if Lot_ID is known, it is stretched teflon then we need to apply a correction to account to lower reflectance on unsampled filters
                    % finds the blanks for the cartridge
                    Ind_blanks_master = find(mass_type == 0);
                    Ind_blanks_SSR = find(ismember(tFilters_SSR,tFilters(Ind_blanks_master))==1);
                    if ~isempty(Ind_blanks_master) % if cartridge blanks exist, correct R based on mean R of all blanks in that cartridge
                        % If R of the blank is greater than all other R in
                        % catridges, correct based on the blank
                        if R_avg(Ind_blanks_SSR) >= max(R_avg)
                            BC_blankCorr = ((q*SA)/ABS_coeff)*log(mean(R_avg(Ind_blanks_SSR))./R_avg);
                            % Otherwise, corrected based on R_max or R_avg
                        else
                            BC_blankCorr = ((q*SA)/ABS_coeff)*log(R_avg_stretch./R_avg);
                        end

                        % If cartridge blank doesn't exist, use mean R of Teflon filters
                    elseif isempty(Ind_blanks_master)==1
                        BC_blankCorr = ((q*SA)/ABS_coeff)*log(R_avg_stretch./R_avg);
                    end
                else % Lot_ID is a nan, therefore a mesh teflon,  use normal BC calculation with R = 100
                    BC_blankCorr = ((q*SA)/ABS_coeff)*log(R_100./R_avg);
                end
                clear Ind_blanks_*

                % %%%%% Line plot to look at changes %%%%%
                x = 1:length(BC_blankCorr);

                figure(1)
                plot(x,BC_standard,'-ob')
                hold on
                plot(x,BC_blankCorr,'-*k')
                xlabel('Filter number')
                ylabel('Estimated black carbon,   \mug/m^3')
                title(sprintf('%s',Site_cities{loc}))
                legend('Standard ','With blank correction')

                SaveDir = sfmkdir (sprintf('%s/Plots_BCtest/%s',direc_SSR,Site_codes{loc}));
                saveas(gcf,sprintf('%s/%s-%s.png',SaveDir,Site_cities{loc},cartridge_IDs{ii}))

                clear BC_standard
                close

                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % check average reflectance measurement range
                R_toohigh = find(R_avg > SSR_threshold(2));
                R_toolow = find(R_avg < SSR_threshold(1));
                Rstd_toohigh = find(R_std > SSR_std);

                if isempty(R_toohigh)==0
                    for h = 1:length(R_toohigh)
                        if mass_type( ismember(tFilters, tFilters_SSR{R_toohigh(h)}) ) ==1
                            fprintf('WARNING: Reflectivity too high for filter %s despite normal flows\n',tFilters_SSR{R_toohigh(h)})
                        end
                    end
                    clear h
                end
                if isempty(R_toolow)==0
                    for h = 1:length(R_toolow)
                        if mass_type( ismember(tFilters, tFilters_SSR{R_toolow(h)}) ) ==1
                            fprintf('WARNING: Reflectivity too low for filter %s despite normal flows \n',tFilters_SSR{R_toolow(h)})
                        end
                    end
                    clear h
                end
                if isempty(Rstd_toohigh)==0
                    for h = 1:length(Rstd_toohigh)
                        fprintf('WARNING: Standard deviation too high for filter %s \n',tFilters_SSR{Rstd_toohigh(h)})
                    end
                    clear h
                end
                
                % no longer mask out R_toolow. This measurement will be
                % used to estimate HIPS_BC
%                 BC_blankCorr(R_toolow) = -899; % set to invalid becasuse below 20 means excessive BC but not reliable measurement
                BC_blankCorr(R_toohigh) = 0; % set to zero because reflectance is above 90, thus negligable BC via this method
                BC_blankCorr(Rstd_toohigh) = -899;
                R_toolow_table(loc)= R_toolow_table(loc)+length(R_toolow);

                for jj = 1:length(BC_blankCorr)
                    tfilter  = tFilters_SSR{jj};
                    master_ind = ismember(Master_IDs,tfilter);
                    Master_carbon(master_ind,1) = BC_blankCorr(jj);
                    clear tfilter master_ind
                end
                clear R_toohigh R_toolow Rstd_toohigh R_avg R_std tFilters_SSR tFilters BC_blankCorr BC_standard Master_ind
            end
            
        end

        Master_carbon(Master_masstype==3,1)= -899; % sets all nuclepore filters to -899; SSR can not measure reflectance on nuclepore filters
        Master_carbon(Master_masstype==4,1)= -899; % sets all saturated nuclepore and associated teflon to -899
        
        %% Write to Master data file
        
        WriteToMaster (  Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, projectIDs_master,...
            Master_hours, Master_masstype, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
            Master_carbon, Master_Method, Master_flags,...
            direc_master, Site_codes{loc})

        clear filename  i I_id analysis_IDs tLot_IDs mass_type num  txt  Master_mass ...
            J_id labels   Master_* foof ID_index_master_A x CartridgeIDs_master...
            master_file  cartridge_IDs LotIDs_master projectIDs_master

    end
end
t = table(Site_codes,R_toolow_table);
writetable(t,sprintf('%s/SSR_too_low.xlsx',direc_SSR))

fprintf('END of SSR BC processing for %s \n', datestr(today))

diary off

%% Functions

function [row, col] = findIndexContaining(searchstring, file)
% find row and column index of string in file

index = find(contains(file,searchstring));
[row, col] = ind2sub(size(file),index); clear index
end











