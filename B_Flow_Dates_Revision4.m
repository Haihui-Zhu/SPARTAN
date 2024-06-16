%% %%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: read and perform QA/QC on excel sheets that contain sampling
% times and flow data for SPARTAN cartridges

% Sheets are named as: [SITE_Code]_dates_flows.xlsx

% Written by: Crystal Weagle
% Created: 24 October 2018

% EDITS:
% --2023-01-26, Haihui Zhu
% 1. Change the structure so that dates, volume, sampled hour are read
% filter by filter (instead of chunk by chunk). Also mass type and method
% are determined filter by filter
% 2. Add a 'clean up master' section to prepare master files for
% reprocessing all variables mentioned above.

% --2021-09-06, Haihui Zhu, applied ReadMaster and WriteToMaster functions,
% added a section to pick up flag in XXX_flow_dates.xlxs

% --2021-06-29, Haihui Zhu, adding lines (around 244-258) to pick up
% nucleopore filter date

% --2020-11-06, Emmie Le Roy, changed directory to storage1, added
% constants section at the top, added warning messages for flow data
% failures, fixed bug where mass_type of blanks (0) was still being overwritten
% as 6 (invalid flow), added loop to ignore ALL partly filled data in
% dates_flows sheet

% -- 2020-05-19, Crystal Weagle, added short loop to check for NaN in
% maa_type column. If NaN's are found, loops through and check flow rates
% to assign a mass type of 1 (PM2.5), 2 (PM10), or 0 (blank)

% -- 2019-09-24, Crystal Weagle, fixed bug that was overwriting mass_type of
% blanks from 0 (blank) to 6 (invalid flow), approx. lines 172-174

% --2019-05-21, Crystal Weagle, added ability to determine method index and
% reference to Site_details.xlsx spreadsheet

close all; clear all; clc
addpath('UtilityFunctions')

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Clean_Up_Master = 0; % set to 1 if want to reprocess all dates and volumes

% Setup directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

direc_sampling = strcat(direc,'Site_Sampling');
direc_master = strcat(direc,'Analysis_Data/Master_files');

% The diary function saves the processing history into an annual record
diary(sprintf('%sPublic_Data/Data_Processing_Records/%s_Flow_Processing_Record',direc,datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))

%-------------   SITE INFO   --------------
% Sites are listed in the file "Site_details.xlsx" in alphabetical order
% Column info:
% 1 = 4-letter site codes
% 2 = Country
% 3 = City
site_details = readtable(strcat(direc_sampling,'/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));

%% Read dates, hours, and volumes. Determine masstype and methods
for loc = 1:length(Site_codes)

    filename = sprintf('%s/%s_%s/Cartridge_Data/%s_dates_flows.xlsx',...
        direc_sampling,Site_codes{loc},Site_cities{loc},Site_codes{loc});
    % Check if a site flow sheet exists
    if exist(filename,'file')==0
        fprintf('WARNING: No flow file exists for %s \n', Site_codes{loc})
        % If no filter flow file exists, moves to the next site
        continue
    end

    fprintf('Reading file for %s \n', Site_codes{loc})

    % ---- make sure volume column is read as double ----
    opts = detectImportOptions(filename);
            FlowTitle = opts.VariableNames;
            col_double = find(contains(FlowTitle,{'volume','start_','stop_','Flow_','hour','SSe_ID'},'IgnoreCase',true)); % columns that need to be read as double 
            for ii = 1:length(col_double)
            opts.VariableTypes{col_double(ii)} = 'double'; end 
            clear ii

    flowraw = readtable(filename,opts);
    % ---- make sure filter ID is in standard format (XXXX_NNNN) ----
    flow_IDs = flowraw.Analysis_ID;
    for i = 1:size(flowraw,1)
        tID = flow_IDs{i};
        if length(tID)<9
            tID2(1,1:5) = tID(1:5);
            tID2(1,6) = '0';
            tID2(1,7:9) = tID(6:8);
        elseif length(tID)>9
            error('filter ID %d too long',tID)
        else 
            tID2 = tID;
        end
        flow_IDs{i}=tID2;
        clear tID tID2
    end
    flowraw.Analysis_ID = flow_IDs; clear flow_IDs
   

    % ---- read master files ----
    master_file = sprintf('%s/%s_master.csv',direc_master,Site_codes{loc});

    [Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, ProjectIDs_master,Master_hours, Master_masstype, ...
        Master_dates, Master_mass, Master_IC, Master_ICP, Master_XRF,...
        Master_carbon, Master_Method, Master_flags] = ReadMaster(master_file,Site_codes{loc});

    if isempty(Master_mass) % skip if master file not exist
        continue
    end

    %%%%%%% Clean Up Master %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Clean up the dates, hours, volumes, masstype in master file to prepare
    % it for reprocessing
    if Clean_Up_Master ==1
        Master_dates     = NaN.*Master_dates;
        Master_mass(:,2) = NaN; % set only volume to NaN, keep mass.
        Master_hours     = NaN.* Master_hours;
        Master_masstype  = NaN.*Master_masstype;
        Master_Method    = NaN.*Master_Method;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [Haihui, 2023, 1, 26 ]
  
    % Finds Master data lines that do not have a start year, therefore missing flow and date data
    index = find(isnan(Master_dates(:,1)));
    if isempty(index) == 1
        fprintf('All samples in %s master file have dates \n', Site_codes{loc})
        disp('Moving on to next site')
        continue
    end

    for ii = 1:length(index)

        tMaster_ID  = Master_IDs{index(ii)};
        tflow_ind   = find(ismember(flowraw.Analysis_ID, tMaster_ID)==1); % find the corresponding filter position in flow data

        if ~isempty(tflow_ind)

            tflow_dates         = table2array(flowraw(tflow_ind,contains(FlowTitle,{'start_','stop_'})));
            tflow_hours_sampled = table2array(flowraw(tflow_ind,contains(FlowTitle,{'hours_sampled','Hours_sampled'})));
            tflow_volume        = table2array(flowraw(tflow_ind,contains(FlowTitle,{'volume'})));
            tflow_rates         = table2array(flowraw(tflow_ind,contains(FlowTitle,{'Flow_'})));
            tflow_SSe_ID        = table2array(flowraw(tflow_ind,contains(FlowTitle,{'SSe_ID'})));
            tflow_filter_ID     = table2array(flowraw(tflow_ind,contains(FlowTitle,{'Filter_ID'})));
            
            tMaster_Mass      = Master_mass(index(ii),1); % col 1 is mass, 2 and 3 are volume and SSR_BC
            tMaster_LotID     = LotIDs_master{index(ii)};
            tMaster_ProjectID = ProjectIDs_master{index(ii)};

            % Check flow rates
            mass_type  = Master_masstype(index(ii));
            % 0 = blank
            % 1 = PM2.5
            % 2 = PM10
            % 3 = PMcoarse (nuclepore)
            % 4 = unknown/void, nuclepore filter saturated
            % 5 = negative mass
            % 6 = invalid flow rates

            %%%%% Pick up nucleopore data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Before setting mass type according to flow rate, check if it is from
            % a nucleopore filter (although unlikely).
            FilterType = tflow_filter_ID{1}(end);
            if FilterType == 'N'
                mass_type = 3; % PMc data from a nucleopore filter
                % if FilterType is a number, it is a SS5 sampler. Set mass_type
                % based on flow rate as shown below
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [Haihui, 2021, 6, 29]

            %%%%% Set expected mass type based on external start flow rate %%%%%%%%
            col_ext_start = 2;
            targets = [5 1.5 0]; % possible target flow rates
            types   = [1  2  0];
            if isnan(mass_type)
                % Grab external start flow rates
                exFlow = tflow_rates(col_ext_start);
                % Find the closest target flow rate
                [~,idx] = min(abs(targets - exFlow)); % idx = idx of min
                % Set mass type based on closest target
                mass_type = types(idx); % sets mass type

                clear idx exFlow
            end


            %%%%% Set Method Index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            method_index = NaN;
            % If SS4 instrument, PM2.5 and PM10 target flow is 4 LPM
            if tflow_SSe_ID/5000 < 1 % SS4
                % Exception for SGSU
                if Site_codes{loc} == 'SGSU' % will give <1 if SSe_ID is SS4 and OLD SGSU sampler
                    target = 5; % SGSU SS4 has target flow rate of 5 LPM
                    method_index = 2;
                    start_flow_ex = tflow_rates(col_ext_start); % external start flow rate
                    r_SF     = start_flow_ex.*[0.9 1.1]; % 10% range of external start flow rate
                    r_target = [target*0.9, target*1.1];  % range of acceptable start flow rates (withing 10% of target flow rate)
                    % Regular SS4 Instruments
                else
                    target = 4;
                    method_index = 1;
                    start_flow_ex = tflow_rates(col_ext_start); % external start flow rate
                    r_SF     = start_flow_ex.*[0.9 1.1]; % 10% range of external start flow rate
                    r_target = target.*[0.9 1.1];  % range of acceptable start flow rates (withing 10% of target flow rate)
                end

                % If SS5, PM10 target flow is 1.5 LPM
            elseif mass_type == 2
                start_flow_ex = tflow_rates(col_ext_start); % external start flow rate
                r_SF = start_flow_ex.*[0.8 1.2]; % range of acceptable flow rates
                r_target = [1.2 1.8]; % 10 % of target flow rate
                if strcmp(tMaster_LotID,'Unknown') == 1 % no LotID, therefore mesh teflon
                    method_index = 1;
                elseif strcmp(tMaster_ProjectID,'S') == 1 % there is a LotID, therefore stretched teflon, now must determine based on project_ID whether SPARTAN or MAIA sampling protocol
                    method_index = 2;
                elseif strcmp(tMaster_ProjectID,'M') == 1
                    method_index = 3;
                end
                % If SS5 and not PM10, target flow rate is 5 LPM
            elseif tflow_SSe_ID/5000 >= 1
                target = 5;
                start_flow_ex = tflow_rates(col_ext_start); % external start flow rate
                r_SF     = start_flow_ex.*[0.9 1.1]; % 10% range of external start flow rate
                r_target = target.*[0.9 1.1];   % range of acceptable start flow rates (withing 10% of target flow rate)
                if strcmp(tMaster_LotID,'Unknown') == 1 % no LotID, therefore mesh teflon
                    method_index = 3;
                elseif strcmp(tMaster_ProjectID,'S') == 1 % there is a LotID, therefore stretched teflon, now must determine based on project_ID whether SPARTAN or MAIA sampling protocol
                    method_index = 4;
                elseif strcmp(tMaster_ProjectID,'M') == 1
                    method_index = 5;
                end
                % If there is no listed SSe_ID, likely did not sample
            elseif isnan(tflow_SSe_ID)
                r_target = [0 0];
                r_SF = [0 0];
                method_index = NaN;
                fprintf('No SSe ID found for flows of %s \n',tMaster_ID)
            end

            %%  ------------------ Flow Rate Validation ------------------
            % after indexing above, write mass type as 6 if flow is invalid.
            % this section will override mass_type = 5 (negative mass) since
            % flow problems are causal

            % Skip blanks
            if tflow_dates(1) == 0 && tflow_volume == 0
                mass_type = 0;

            % Skip neg masses and nuclepore filters
            elseif mass_type == 5 || mass_type == 3
                
            else
                % CHECK: External Start Flow within 10% of Target Flow Rate
                if tflow_rates(2) > r_target(2) || tflow_rates(2) < r_target(1)
                    fprintf('WARNING: external start flow rate not within 10%% of target flow rate for %s \n', tMaster_ID)
                    mass_type = 6; % set mass type to 6 for invalid flow rates
%                     tflow_volume = NaN; % because invalid flow rates no volume is recognized, keeps dates
                end

                % CHECK: External End Flow within 10% of External Start Flow Rate
                if tflow_rates(4) > r_SF(2) || tflow_rates(4) < r_SF(1)
                    fprintf('WARNING: external end flow rate not within 10%% of external start flow rate, indicates clogging for %s \n', tMaster_ID)
                    mass_type = 6;  % set mass type to 6 for invalid flow rates
%                     tflow_volume = NaN; % because invalid flow rates no volume is recognized, keeps dates
                end

                % CHECK: Invalid internal flow reported by instrument
                if tflow_rates(1) ~= -899 % internal flow rate for CASH never work, skip CASH to do the check
                    if tflow_rates(1) > r_SF(2) || tflow_rates(1) < r_SF(1) % internal start flow rate must be within 10 % of external start flow rate
                        mass_type = 6;  % set mass type to 6 for invalid flow rates
%                         tflow_volume = NaN; % because invalid flow rates no volume is recognized, keeps dates
                        fprintf('WARNING: internal start flow rate not within 10%% of external start flow rate for %s \n', tMaster_ID)
                    end
                end

            end
            % -------------------------------------------------

            Master_hours(index(ii))    = tflow_hours_sampled;
            Master_dates(index(ii),:)  = tflow_dates; % writes sampling dates to master file
            Master_masstype(index(ii)) = mass_type; % writes new mass type bc new flag = 6 (invalid flow)
            Master_mass(index(ii),1)   = tMaster_Mass; % writes updated Mass to master file
            Master_mass(index(ii),2)   = tflow_volume; % writes volumes to master file
            Master_Method(index(ii))   = method_index;
                
            clear invflow start_flow_ex r_SF r_target  target method_index
          
        end

    end

    %% Check & Add Flags
    Flow_flag  = table2array(flowraw(:,contains(FlowTitle,'flags')));
    Flow_AnaID = table2array(flowraw(:,contains(FlowTitle,'Analysis_ID')));
    if isa(Flow_flag,'cell') == 1 % if Flow_flag is a mat, there isn't any flag. Otherwise, Flow_flag will be a cell
        for ii = 1:length(Flow_flag) % Even Flow_flag{ii} is empty, the AddFlag fxn will reformat flags (e.g. removing extra blank)
            jj = find(ismember(Master_IDs,Flow_AnaID{ii}));
            if ~isempty(jj)
                Master_flags{jj} = AddFlag(Master_flags{jj},Flow_flag{ii});
            end
        end
    end
    % add flag for negative sampled mass and invalid flow rate
    ind = find(Master_masstype==5); % negative mass
    for ii = 1:length(ind)
        Master_flags{ind(ii)} = AddFlag(Master_flags{ind(ii)},'Negative Mass');
    end
    ind = find(Master_masstype==6); % invalid flow
    for ii = 1:length(ind)
        Master_flags{ind(ii)} = AddFlag(Master_flags{ind(ii)},'Invalid Flow');
    end

    %% Write to Master Data File
    WriteToMaster(Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, ProjectIDs_master,Master_hours, Master_masstype, ...
        Master_dates, Master_mass, Master_IC, Master_ICP, Master_XRF,...
        Master_carbon, Master_Method, Master_flags,...
        direc_master, Site_codes{loc})

    clear Master_* filename tflow_* mass_type flags_flow  master_file method_index  ...
          CartridgeIDs_master LotIDs_master target Master_Barcodes ProjectIDs_master Flow_flag Flow_AnaID

end
disp('Finished reading flow dates');
diary off

