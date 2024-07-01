%% %%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: Read  filter masses from MTL weighing system output and write to
% site-specific file for further processing by other scripts. 

% FILE REQUIREMENTS: Filter mass data from MTL should be exported as a Test 
% Report ("pre-weights" or "pre- & post-weights") for each individual cart-
% ridge and the comma-delimited text file should be saved as a '.xlsx' file. 
% The full filterID (i.e.'CAHA-0105-3') should be used in the MTL system.

% To rerun a raw MTL file: 
% Simply move the raw file to direc_new
% If you do not want to see warnings about replacing existing weights to 
% new values, set replace_warning to 0. 


% Written by: Haihui Zhu, Crystal Weagle, Emme Le Roy, and Chris Oxford
% Created: 12 May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear ; clc
addpath('UtilityFunctions')
warning('off')
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%% USER SWITCHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);
replace_warning = 1; % if you expect to have many values replaced and don't want to see the warnings, turn this off (= 0). 

direc_new=strcat(direc,'Analysis_Data/Filter_Masses/NEW_MTL');
direc_output=strcat(direc,'Analysis_Data/Filter_Masses/Masses_by_site/MTL_weighing_WashU'); 
direc_archive=strcat(direc, 'Analysis_Data/Archived_Filter_Data/MTL_masses');

% Diary
diary(sprintf('%s/Public_Data/Data_Processing_Records/%s_MTLmasses_Record',direc,datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))

% Title of the output MTL masses table
titles = {'FilterID','AnalysisID','Filter_Barcode','CartridgeID','LotID', ...
    'Preweigh_Date','Preweight_avg_mg','Preweight_std_ug','Preweight_NumReps',...
    'Postweigh_Date','Postweight_avg_mg','Postweight_std_ug','Postweight_NumReps', 'Net_Weight_ug'};

%% %%%%%%%%%%%%%%%%%%%%%%%% FIND & READ DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section will find the sample IDs in the MTL report, isolate the
% Site IDs for indexing to site MTL mass files

% Read Acceptance Testing Summary File (to fill in CartridgeIDs and LotIDs)
Sum_lots = readtable((strcat(direc,'Analysis_Data/Acceptance_Testing/Acceptance_Testing_Summary_WashU.xlsx'))); 

%  Read all MTL filter masses test reports in the 'MTL_NEW' directory
files = getFiles(direc_new);

for data_file = 1:length(files)
    %% ------------------------READ MTL TEST REPORT-------------------------
     
    % Do not read hidden files
    char_file = char(files(data_file));
    if char_file(1)=="."  || char_file(1)== "~" 
        continue
    end
    
    filename = sprintf('%s/%s',direc_new,files{data_file,1});
    fprintf('Reading %s \n', files{data_file,1});
    opts = detectImportOptions(filename);

    % Read Pre-weight and Post-weight Dates
    opts.VariableTypes(1:3) = {'char'}; % make sure barcode is read as char, not numbers
    opts.VariableTypes(4:13) = {'double'}; 
    opts.DataRange = 'A1'; % reading additional info at the top of the table
    DateInfo = readtable(filename,opts); % DateInfo + raw_file complete the table
    pre_date = table2array(DateInfo(contains(table2array(DateInfo(:,1)),'Pretest Date'), 2)); 
    post_date = table2array(DateInfo(contains(table2array(DateInfo(:,1)),'Posttest Date'), 2));
    
    if isempty(char(post_date)) == 1 % no post-weigh data
        postweights_avail = 0; % will be used later
        post_date = {'Unknown'};
    else % There are post-weights
        postweights_avail = 1;
    end

    % Read filter mass
    DataVarNameInd = find(contains(table2array(DateInfo(:,1)),'Position'));
    opts.VariableNamesRange = ['A',num2str(DataVarNameInd)];
    opts.DataRange = ['A',num2str(DataVarNameInd+1)]; 
    raw_file = readtable(filename,opts); 
    raw_file = raw_file(1:20,1:14);
    raw_file.PretestStandardDeviation = round(raw_file.PretestStandardDeviation,2);
    raw_file.PosttestStandardDeviation = round(raw_file.PosttestStandardDeviation,2);
    raw_titles = raw_file.Properties.VariableNames; 

    clear DateInfo opts 

    % check if standard deviation too high:
    if sum(find(raw_file.PretestStandardDeviation > 2.5)) > 0
        error('Standard devaition of preweight mass too high! Measure the mass again.')
    elseif sum(find(raw_file.PosttestStandardDeviation > 2.5)) > 0
        error('Standard devaition of postweight mass too high! Measure the mass again.')
    end

    % Get row # of the filters from MTL Test Report 
    samples_rowstart = find(contains( raw_file.Position,'Phase #1'));
    sample_rows = samples_rowstart(1):samples_rowstart(2)-2; % Looks for the second time it says Phase 1 and subtracts 2. This gets you to Phase #8, or filter 8 of the cartridge
    
    % IDs 
    filter_IDs = raw_file.FilterID(sample_rows); 
    filter_IDs_char = char(filter_IDs);
    site_IDs = unique(filter_IDs_char(:,1:4),'rows');
    analysis_IDs = cellstr(filter_IDs_char(:,1:9)); clear filter_IDs_char
    
    % Find barcode
    if sum(contains(raw_titles,'MTLFilterID')) == 1 % the MTL Filter ID (barcode) column exist
        filter_barcodes = raw_file.MTLFilterID(sample_rows);  

    elseif sum(contains(raw_titles,'MTLFilterID')) == 0 
        error("Missing Filter Barcodes! Please fill in manually in the MTL Mass File or re-export data!")
    
    else 
        error('More than 1 barcode columns! Check the MTL Mass File!')
    end


    %% --------------- READ & WRITE TO MTL MASS SHEET ---------------------
    % Search for MTL mass sheet for each site with data in the MTL report
    % One file for each site. 
    for site = 1:size(site_IDs,1)
       
        % ---- Collecting pre (and post weight) for this site ----

        % Index of samples for site_ID(i) in report files
        samples_site_idx = find(contains(filter_IDs,site_IDs(site,:)));
        site_filterIDs = filter_IDs(samples_site_idx);

        site_analysisIDs = analysis_IDs(samples_site_idx);
        site_barcodes = filter_barcodes(samples_site_idx);

        [site_lotIDs, site_cartridgeIDs] = find_lot_cart(site_filterIDs,Sum_lots);

        % Get mass information for this site from MTL report file
        preweigh_date(1:size(samples_site_idx,1),1) = pre_date;
        preweight_avgs = raw_file.CorrectedPretest(sample_rows(samples_site_idx)); % units = mg
        preweight_stds = raw_file.PretestStandardDeviation(sample_rows(samples_site_idx)); % units = ug
        preweight_reps = raw_file.PretestRepetitions(sample_rows(samples_site_idx));  % unitless

        % possible that first file found has both preweights and postweights
        if postweights_avail == 1
            postweight_date(1:size(samples_site_idx,1),1) = post_date;
            postweight_avgs = raw_file.CorrectedPosttest(sample_rows(samples_site_idx)); % units = mg
            postweight_stds = raw_file.PosttestStandardDeviation(sample_rows(samples_site_idx)); % unit = ug
            postweight_reps = raw_file.PosttestRepetitions(sample_rows(samples_site_idx)); % unitless
            netweight = round((postweight_avgs - preweight_avgs)*1000,1);  % multiply by 1000 to get in units of ug
            % no postweights, therefore set columns to NaN
        elseif postweights_avail ==0
            postweight_date(1:size(preweight_avgs,1),1) = post_date;
            postweight_avgs(1:size(preweight_avgs,1),1) = NaN;
            postweight_stds(1:size(preweight_avgs,1),1) = NaN;
            postweight_reps(1:size(preweight_avgs,1),1) = NaN;
            netweight(1:size(preweight_avgs,1),1) = NaN;
        end


        % ---- Reading MTL MASS SHEET  ----
        % Three Scenarios based on completeness of MTL mass sheet:
        % 1. No MTL mass file exist - create one. OutputType = 1. 
        % 2. MTL mass file exists, but doesn't include any filters in the
        % report - adding the pre (and possibly post) weight data to file. 
        % OutputType = 2.
        % 3. MTL mass file exists, and includes all filters in the new report -
        % overwrite the existing value with the new values. OutputType = 3.

        mass_file = sprintf('%s/%s_MTL_masses.csv',direc_output,site_IDs(site,:));
        
        OutputType = 1; 

        if exist(mass_file,'file') % Scenario 3: MTL mass file exist. Read it.

            OutputType = 2; 

            opts = detectImportOptions(mass_file);
            opts.ExtraColumnsRule = 'ignore'; 
            opts.VariableTypes(1:3) = {'char'};
            mass_file_table = readtable(mass_file,opts);
                
            % find filter IDs NOT in MTL mass file
            newIDs_idx = find(ismember(site_filterIDs,mass_file_table.FilterID) == 0); 
            % note: if there isempty(newIDs_idx) is not true, we assume 
            % there are at least 8 filters (a cartridge) not in
            % mass_file_table. [Senario 2]

            if isempty(newIDs_idx) % Scenario 3: all filters included. 
                % Need to compare with the existing value. Give a warning if
                % existing value is not 'nan' or not the same as the new value. 
                OutputType = 3;

                for ii = 1:length(site_filterIDs)
                    % find index of filters in mass_file_table
                    idx_massfile = find(ismember(mass_file_table.FilterID, site_filterIDs(ii)) == 1);

                    mass_file_table = check_and_fill(mass_file_table,'AnalysisID', idx_massfile, site_analysisIDs(ii),site_filterIDs(ii), replace_warning);
                    mass_file_table = check_and_fill(mass_file_table,'Filter_Barcode', idx_massfile, site_barcodes(ii),site_filterIDs(ii), replace_warning);
                    mass_file_table = check_and_fill(mass_file_table,'CartridgeID', idx_massfile, site_cartridgeIDs(ii),site_filterIDs(ii), replace_warning);
                    mass_file_table = check_and_fill(mass_file_table,'LotID', idx_massfile, site_lotIDs(ii),site_filterIDs(ii), replace_warning);

                    mass_file_table  = check_and_fill(mass_file_table,'Preweigh_Date', idx_massfile, preweigh_date(ii),site_filterIDs(ii), replace_warning);
                    mass_file_table  = check_and_fill(mass_file_table,'Preweight_avg_mg', idx_massfile, preweight_avgs(ii),site_filterIDs(ii), replace_warning);
                    mass_file_table  = check_and_fill(mass_file_table,'Preweight_std_ug', idx_massfile, preweight_stds(ii),site_filterIDs(ii), replace_warning);
                    mass_file_table  = check_and_fill(mass_file_table,'Preweight_NumReps', idx_massfile, preweight_reps(ii),site_filterIDs(ii), replace_warning);
                    
                    if postweights_avail == 1
                        mass_file_table = check_and_fill(mass_file_table,'Postweigh_Date', idx_massfile, postweight_date(ii),site_filterIDs(ii), replace_warning);
                        mass_file_table = check_and_fill(mass_file_table,'Postweight_avg_mg', idx_massfile, postweight_avgs(ii),site_filterIDs(ii), replace_warning);
                        mass_file_table = check_and_fill(mass_file_table,'Postweight_std_ug', idx_massfile, postweight_stds(ii),site_filterIDs(ii), replace_warning);
                        mass_file_table = check_and_fill(mass_file_table,'Postweight_NumReps', idx_massfile, postweight_reps(ii),site_filterIDs(ii), replace_warning);

                        mass_file_table = check_and_fill(mass_file_table,'Net_Weight_ug', idx_massfile, netweight(ii),site_filterIDs(ii), replace_warning);
                    end
                end

            end
                       
        end


        % ---- Writing to MTL MASS SHEET  ----

        % combine all variable to a table in the same format as the mass file:
        newtable = table(site_filterIDs,site_analysisIDs,site_barcodes,site_cartridgeIDs,site_lotIDs,...
            preweigh_date,preweight_avgs,preweight_stds,preweight_reps,...
            postweight_date,postweight_avgs,postweight_stds,postweight_reps,netweight,'VariableNames',titles);

        if OutputType == 1
            fprintf('No MTL mass file exists for %s. Creating new MTL mass file.\n', site_IDs(site,:))
            mass_file_table = newtable;
        elseif OutputType == 2
            % add the new rows to the original table:
            mass_file_table(height(mass_file_table)+1 : height(mass_file_table)+height(newtable),:) = newtable;
            
        % elseif OutputType == 3 % New info in mass_file_table already. Nothing need to be done.             
        end


        %  WRITE TO MTL MASS SHEET
        outfile = sprintf('%s/%s_MTL_masses.csv',direc_output,site_IDs(site,:));
        writetable(mass_file_table,outfile)
        fprintf('Finished writing to %s MTL mass file \n', site_IDs(site,:))

        % Copy file to archived MTL masses folder for future reference
        file_destination = strcat(direc_archive,sprintf('/%s/',site_IDs(site,:)));
        status = CopyFile(filename,file_destination);

        if isempty(find(status(:)==0,1))==1
            delete(filename)
            disp('File has been moved from NEW_MTL folder to Archives')
        else
            disp('Cannot delete file from NEW_MTL data folder as it has not been moved to all archived site folders')
        end

        % Clear variables
        clear site_filterIDs_new site_barcodes samples_site_idx...
            site_cartridgeIDs site_lotIDs preweigh_date ...
            preweight_avgs preweight_stds preweight_reps ...
            postweight_date postweight_avgs postweight_stds ...
            postweight_reps netweight site_filterIDs ...
            site_analysisIDs site_barcodes site_cartridgeIDs site_lotIDs ...
            site_preweigh_dates status file_destination...
            site_postweigh_dates netweight mass_file idx_massfile ...
            mass_file_table compare_IDs newIDs_idx rowidx
    end
     
    % clear variables specific to each data file so not to interfere with next file processing
    clear mass_file samples_site_idx char_file  filename raw_file ...
          raw_file samples_rowstart sample_rows IDs_columns filter_IDs site_IDs ...
          analysis_IDs filter_barcodes pre_date_row pre_date_col post_date_row ...
          post_date_col preweigh_date postweight_date postweights_avail ...
          ans
            
end

function [site_lotIDs, site_cartridgeIDs] = find_lot_cart(site_filterIDs,Sum_lots)

for ii = 1:length(site_filterIDs)
    for col = 1:size(Sum_lots,2) % looping through column
        tdata = table2array(Sum_lots(:,col));
        if isa(tdata,'cell')
            ind = find( contains(tdata,site_filterIDs{ii}) );
            if ~isempty(ind)
                % LotIDs are 2 columns to left of the filterID column
                site_lotIDs(ii,1) = table2array(Sum_lots(ind, col-2));
                % CartridgeIDs are 1 column to left of the filterID column
                site_cartridgeIDs(ii,1) = table2array(Sum_lots(ind, col-1));
            end
        end
    end

    if isempty(site_lotIDs(ii,1)) || isempty(site_cartridgeIDs(ii,1))
        error('Lot ID or cartridge ID not found for %s!',site_filterIDs{ii})
    end
end

end

function mass_file_table = check_and_fill(mass_file_table, var, idx_massfile, new_value, fileterID, replace_warning)
titles = mass_file_table.Properties.VariableNames;
tcolumn = table2array( mass_file_table(:, ismember(titles,var)) );

if isa(new_value,'double')
   if ~isnan(new_value) && tcolumn(idx_massfile) ~= new_value
       tcolumn(idx_massfile) = new_value;
       if replace_warning == 1
            warning('%s %s old value %.2d replaced with %.2d\n',fileterID, var, tcolumn(idx_massfile), new_value)
       end
   end
    mass_file_table(:, ismember(titles,var)) = table(tcolumn);

else % not a double format => filter info 
    tcolumn{idx_massfile} = new_value{1};
    mass_file_table(:, ismember(titles,var)) = tcolumn;
end


end
