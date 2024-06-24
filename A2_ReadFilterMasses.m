%% %%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: read  filter masses from site-specific file to create or add to 
% the master data file. Perform QC on filter masses.
% Before a filter ID is added to a Master file using this script, no other
% data can be added for that filter ID

% Written by: Crystal Weagle
% Created: 14 May 2020

% -----EDITS (yyyy-mm-dd, author name): 
% Revision 2: 2020-May-14, Crystal Weagle 
% Revamped Revision1 script to include checking for MTL mass
% summary spreadsheets for each site. There is no longer reading in of all 
% filter masses in the spreadsheets if a Master file already exists, rather 
% the script will check what filters are already in the master file and only
% read in data for new post-weighed filters. Thus, I removed the section that
% checks for nuclepore filters. Other than setting mass type to 5 if a filter
% mass is negative, this script no longer assigns mass type. The master file
% is still created in a similar manner. 

% Revision 3: 
% 2020-Oct-27, Emmie Le Roy:
% Changed directories for move to Storage1, Fixed bug where filter IDs with
% no post-weights were being added to the master file when "pre-weight
% only" files were being processed. This prevented future post-weights from
% being written to the master file. Now, only MTL data with net-weights are
% written to the master file. Also fixed bug where Project IDs were blank,
% now Project IDs pull from the Site_Details sheet and select the first 
% character of the Sampling_Mode (i.e. 'SPARTAN' = 'S'; 'MAIA' = 'M').

% 2020-Dec-05, Emmie Le Roy:
% (maybe not)Re-added mass_type assignment for nuclepore filters to re-run scripts
% from scratch.

% Revision 4:
% 2021-Aug-25 Haihui Zhu
% 1. Applied the ReadMaster function and the WriteToMaster function 
% 2. Replaced readcell by readtable(simplified the script)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ; close all;  clc
addpath('UtilityFunctions')

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%% USER SWITCHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

direc_new = strcat(direc,'Analysis_Data/Filter_Masses/Masses_by_site'); % filter mass data files
direc_master = strcat(direc,'Analysis_Data/Master_files'); % master data files

%Read the table
site_details = readtable(sprintf('%s/Site_Sampling/Site_details.xlsx',direc),'PreserveVariableNames',true);%,opts,'ReadVariableNames', true);
Site_codes = site_details.Site_Code;
Site_cities = site_details.City;
Sampling_Mode = site_details.Sampling_Mode;

% The diary function saves the processing history into an annual record
diary(sprintf('%s/Public_Data/Data_Processing_Records/%s_FilterMasses_Record',direc,datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
label_length_new = 11; % i.e. 'AEAZ-0010-1'
label_length_old = 13; % i.e. '10001-AEAZ-1T'

neg_threshold = -0.3; 
                            
%% %%%%%%%%%%%%%%%%%%%%%%%% FIND & READ DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section will read both the master files and mass files and determine
% which filters need to be added to the master file
% for loc = 27:27 %(Site_codes)
for loc = 1:length(Site_codes)
 
    % Check if the master file for the site already exists
      master_file = sprintf('%s/%s_master.csv',direc_master,Site_codes{loc});
      
     [Titles,Master_IDs, Master_Barcodes, Master_CartridgeIDs, Master_LotIDs, Master_ProjectIDs,Master_hours, Master_masstype, ...
      Master_dates, Master_mass, Master_IC, Master_ICP, Master_XRF,...
      Master_carbon, Master_Method, Master_flags] = ReadMaster(master_file,Site_codes{loc});
     
    if ~isempty(Master_CartridgeIDs) 
        master_exist = 1; % will read exist master file
        master_cartridges = unique(Master_CartridgeIDs); 
    else
        master_exist = 0; % will create a new file
    end
    clear master_file
 
    %% Check if manual filter weighing mass sheet exists for the site 
    % and whether it contains post-weights
    filename = sprintf('%s/Manual_weighing_Dal/%s_masses.xlsx',direc_new,Site_codes{loc});
    if exist(filename,'file') % file found
       fprintf('Reading manual weighing file for %s \n', Site_codes{loc})
        
       % Read manual data   
       table_dal = readtable(filename);
       titles_dal = table_dal.Properties.VariableNames;
       % Find dal columns
       dal_postweight_col = find(contains(titles_dal,'Post_Sampling_Weight')==1);
       dal_netweight_col  = find(contains(titles_dal,'Actual_Weight_of_Sample_PM')==1);
       dal_PMerror_col    = find(contains(titles_dal,'Error_in_PM')==1);
       dal_cartridgeID_col = find(contains(titles_dal,'Cartridge_ID')==1);
       dal_filterID_col = find(contains(titles_dal,'Filter_ID')==1);
       dal_barcode_col  = find(contains(titles_dal,'Barcode')==1);
       dal_lotID_col    = find(contains(titles_dal,'Lot_ID')==1);

        % If postweights exist in Dal manual filter mass file
        if ~isempty(dal_postweight_col) 
            fprintf('Post-weights found in %s manual weighing file \n',Site_codes{loc})
            
            % Find indices where net-weights exist (screen NaNs)
            netweight_dal = table2array(table_dal(:,dal_netweight_col));
            masses_idx = find(~isnan(netweight_dal));
            
            % If there are non-NaN post-weights in file and master file exists
            if isempty(masses_idx)==0 && master_exist==1                
                % Compare cartridges and find index of new cartridges to be added to master file
                cartridges_dal = table2array(table_dal(:,dal_cartridgeID_col)); 
                newIDs_idx2 = find(ismember(cartridges_dal,master_cartridges)==0); % 0 = not in the master file; find index of new masses in "Masses" spreadsheet. The same numbers will be new rows in the Master data file
                newIDs_idx = newIDs_idx2(~isnan(netweight_dal(newIDs_idx2))); % new net weight 
              
                if isempty(newIDs_idx) == 0
                    % Formating filterID and mass_type
                    Labels_dal = table2array(table_dal(newIDs_idx,dal_filterID_col));
                    Labels_dal_char = char(deblank(Labels_dal));

                    AnalysisIDs_dal = cell(length(Labels_dal),1);
                    Filter_type = char();
                    MassType_dal = NaN(size(AnalysisIDs_dal));

                    for h = 1:length(Labels_dal)
                        if length(Labels_dal{h})==label_length_new
                            AnalysisIDs_dal(h,:) = cellstr(Labels_dal_char(h,1:9));
                            Filter_type(h,:) = 'T';
                        elseif length(Labels_dal{h})==label_length_old
                            AnalysisIDs_dal(h,:) = cellstr(strcat(Labels_dal_char(h,7:10),'-',Labels_dal_char(h,3:5)));
                            Filter_type(h,:) = Labels_dal_char(h,13);
                        else
                            fprintf('WARNING: label size not 11 or 13 for %s \n',Labels_dal{h} );
                        end
                    end

                    MassType_dal(Filter_type=='N') = 3; % sets nuclepore filters to mass type 3
                    Barcodes_dal = table2array(table_dal(newIDs_idx,dal_barcode_col)); 
                    CartridgeIDs_dal = table2array(table_dal(newIDs_idx,dal_cartridgeID_col));
                    LotIDs_dal = table2array(table_dal(newIDs_idx,dal_lotID_col));
                    Masses_dal = table2array(table_dal(newIDs_idx,dal_netweight_col));

                    Barcodes_dal = Reformat(Barcodes_dal);
                    LotIDs_dal = Reformat(LotIDs_dal);
                else
                    fprintf('No post-weights in %s manual weighing file to be added to existing master file \n', Site_codes{loc})
                end
            
            % If there are non-NaN post-weights in file but there is no existing master file    
            elseif isempty(masses_idx) == 0 && master_exist == 0
                % we will not compare IDs in the two files, all post-weighs will need to be added
                % Formating filterID and mass_type
                    Labels_dal = table2array(table_dal(masses_idx,dal_filterID_col));
                    Labels_dal_char = char(deblank(Labels_dal));

                    AnalysisIDs_dal = cell(length(Labels_dal),1);
                    Filter_type = char();
                    MassType_dal = NaN(size(AnalysisIDs_dal));

                for h = 1:length(Labels_dal)
                    if length(Labels_dal{h})==label_length_new
                        AnalysisIDs_dal(h,:) = cellstr(Labels_dal_char(h,1:9));
                        Filter_type(h,:) = 'T';
                    elseif length(Labels_dal{h})==label_length_old
                        AnalysisIDs_dal(h,:) = cellstr(strcat(Labels_dal_char(h,7:10),'-',Labels_dal_char(h,3:5)));
                        Filter_type(h,:) = Labels_dal_char(h,13);
                    else
                        fprintf('WARNING: label size not 11 or 13 for %s \n',Labels_dal{h} );
                    end    
                end    
                MassType_dal(Filter_type=='N') = 3; % sets nuclepore filters to mass type 3
                Barcodes_dal = table2array(table_dal(masses_idx,dal_barcode_col));
                CartridgeIDs_dal = table2array(table_dal(masses_idx,dal_cartridgeID_col));
                LotIDs_dal = table2array(table_dal(masses_idx,dal_lotID_col));
                Masses_dal = table2array(table_dal(masses_idx,dal_netweight_col));

                Barcodes_dal = Reformat(Barcodes_dal);
                LotIDs_dal = Reformat(LotIDs_dal);
            end

        else
            fprintf('No post-weights in %s manual weighing file \n', Site_codes{loc})        
        end
    else
        fprintf('No manual weighing file exists for %s \n', Site_codes{loc})
    end
    
    clear filename masses_idx cartridges_dal table_dal master_cartridges Labels_dal Labels_dal_char Filter_type dal_*
    
    %% Check if MTL automatic weighing mass sheet exists for the site 
    % and whether it contains post-weights
    filename = sprintf('%s/MTL_weighing_WashU/%s_MTL_masses.csv',direc_new,Site_codes{loc});
    if exist(filename,'file') % file found
        fprintf('Reading MTL weighing file for %s \n', Site_codes{loc})

        MTL_data_initial = readtable(filename,'Format','%s %s %s %s %s %s %n %n %n %s %n %n %n %n'); 
        
        % Find indices where net-weights exist (screen NaNs)
        netweights_mtl = MTL_data_initial.Net_Weight_ug;
        masses_idx = find(~isnan(netweights_mtl)); % will be empty if all netweights are NaN

        % If there are non-NaN post-weights in the file and the master file exists
        if isempty(masses_idx)==0 && master_exist == 1
            
            % Compare cartridges and find index of new cartridges to be added to master file
            analysisIDs_mtl = MTL_data_initial.AnalysisID;
            newIDs_idx2 = find(ismember(analysisIDs_mtl,Master_IDs)==0); % 0 = not in the master file; Find index of new masses in "Masses" spreadsheet
            newIDs_idx = newIDs_idx2(~isnan(netweights_mtl(newIDs_idx2))); % screen for NaN netweights (i.e. pre-weight only data)
            clear analysisIDs_mtl 
            
            % If there are cartridges in the MTL sheet that need to be added to master file
            if isempty(newIDs_idx)==0
                % Select the new MTL data to be added
                AnalysisIDs_mtl = MTL_data_initial.AnalysisID(newIDs_idx);
                Barcodes_mtl = MTL_data_initial.Filter_Barcode(newIDs_idx);
                CartridgeIDs_mtl = MTL_data_initial.CartridgeID(newIDs_idx);
                LotIDs_mtl = MTL_data_initial.LotID(newIDs_idx);
                Masses_mtl = MTL_data_initial.Net_Weight_ug(newIDs_idx);
                MassType_mtl = NaN(size(AnalysisIDs_mtl));

                Barcodes_mtl = Reformat(Barcodes_mtl);
                LotIDs_mtl = Reformat(LotIDs_mtl);
            else
                fprintf('No post-weights in %s MTL weighing file to be added to existing master file \n', Site_codes{loc})
            end

        % If there are non-NaN post-weights in file but there is no existing master file    
        elseif isempty(masses_idx) == 0 && master_exist == 0
            % we will not compare IDs in the two files, all post-weighs will need to be added
            AnalysisIDs_mtl = MTL_data_initial.AnalysisID(masses_idx);
            Barcodes_mtl    = MTL_data_initial.Filter_Barcode(masses_idx);
            CartridgeIDs_mtl= MTL_data_initial.CartridgeID(masses_idx);
            LotIDs_mtl   = MTL_data_initial.LotID(masses_idx);
            Masses_mtl   = MTL_data_initial.Net_Weight_ug(masses_idx);
            MassType_mtl = NaN(size(AnalysisIDs_mtl));

            Barcodes_mtl = Reformat(Barcodes_mtl);
            LotIDs_mtl = Reformat(LotIDs_mtl);

        else % there are no non-NaN post-weights in the file
            fprintf('No post-weights in %s MTL weighing file \n', Site_codes{loc})
        end

    else
        fprintf('WARNING: No MTL weighing file exists for %s \n', Site_codes{loc})
    end
    
    clear masses_idx netweights_mtl newIDs_idx MTL_data_initial
    
%% %%%%%%%%%%%%%%%% Concatenate MTL and Dal weighing data %%%%%%%%%%%%%%%%%

    if exist('Masses_mtl','var') == 1 && exist('Masses_dal','var') == 1 % MTL and Dal masses are available
        % if both exist, dal will have "older" samples, therefore put them first in filter order

        new_AnalysisIDs = vertcat(AnalysisIDs_dal, AnalysisIDs_mtl);
        new_Barcodes = vertcat(Barcodes_dal, Barcodes_mtl);
        new_CartridgeIDs = vertcat(CartridgeIDs_dal, CartridgeIDs_mtl);
        new_LotIDs = vertcat(LotIDs_dal, LotIDs_mtl);
        char_SamplingMode = char(Sampling_Mode(loc));
        new_ProjectIDs = repelem({char_SamplingMode(1)},length(new_AnalysisIDs))'; clear char_SamplingMode 
        new_Masses2 = vertcat(Masses_dal, Masses_mtl); % concatenate collected mass information (ug)
        new_MassType = vertcat(MassType_dal, MassType_mtl);
        clear *_dal *_mtl

    elseif exist('Masses_mtl','var')==1 && exist('Masses_dal','var')==0 % Just MTL masses available

        new_AnalysisIDs = AnalysisIDs_mtl;
        new_Barcodes = Barcodes_mtl;
        new_CartridgeIDs = CartridgeIDs_mtl;
        new_LotIDs = LotIDs_mtl;
        char_SamplingMode = char(Sampling_Mode(loc));
        new_ProjectIDs = repelem({char_SamplingMode(1)},length(new_AnalysisIDs))'; clear char_SamplingMode 
        new_Masses2 = Masses_mtl; % units = ug
        new_MassType = MassType_mtl; % set mass_type to NaN, will be filled in below if negative mass, otherwise will be filled in by B_Flow_Dates script
        clear *_mtl

    elseif exist('Masses_mtl','var')==0 && exist('Masses_dal','var')==1 %Just Dal masses available

        new_AnalysisIDs = AnalysisIDs_dal;
        new_Barcodes = Barcodes_dal;
        new_CartridgeIDs = CartridgeIDs_dal;
        new_LotIDs = LotIDs_dal;
        char_SamplingMode = char(Sampling_Mode(loc));
        new_ProjectIDs = repelem({char_SamplingMode(1)},length(new_AnalysisIDs))'; clear char_SamplingMode
        new_Masses2 = Masses_dal; % units = ug
        new_MassType = MassType_dal;  
        clear *_dal 
    end
    
    if exist('new_Masses2','var')
        % Search for negative masses and set to NaN if below threshold (and
        % set mass type to 5 if neg mass)
        [new_Masses, new_MassType] = screenNegMasses(new_Masses2, neg_threshold, new_MassType); clear new_Masses2
    else % no new post-weights to add, move on to next site
        fprintf('No post-weights available for %s \nMoving on to next site \n\n', Site_codes{loc})
        clear Master_* 
        continue % no post-weights available from MTL or Dal weighing
    end
 
    
%% %%%%%%%%%%%%%%%%%% Write to master data csv files %%%%%%%%%%%%%%%%%%%%%%
        % need to concatenate the columns before writing to data file
        rows = size(Master_IDs,1)+1:size(Master_IDs,1)+ size(new_AnalysisIDs,1);
        Master_IDs(rows,:) = new_AnalysisIDs;
        Master_Barcodes(rows,:) = new_Barcodes;
        Master_CartridgeIDs(rows,:) = new_CartridgeIDs;
        Master_LotIDs(rows,:) = new_LotIDs;
        Master_ProjectIDs(rows,:) = new_ProjectIDs;
        Master_masstype(rows,:) = new_MassType;
        Master_mass(rows,1) = new_Masses;

        for t = rows(1):rows(end)
            Master_flags{t,1}= '';
        end

        % Fill in all numerical master data with NaN
        Master_mass(rows,2) = NaN;
        Master_hours(rows,:) = NaN;
        Master_dates(rows,:) = NaN;
        Master_IC(rows,:) = NaN;
        Master_ICP(rows,:) = NaN;
        Master_XRF(rows,:) = NaN;
        Master_carbon(rows,:) = NaN;
        Master_Method(rows,:) = NaN;

        clear new_*

        % ---- WRITE FILE -----
        WriteToMaster(Titles,Master_IDs, Master_Barcodes, Master_CartridgeIDs, Master_LotIDs, Master_ProjectIDs,...
            Master_hours, Master_masstype, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
            Master_carbon, Master_Method, Master_flags,direc_master,Site_codes{loc})

        clear fileID Master_* master_exist master_file h
end

fprintf('END of filter mass processing for %s \n', datestr(today))

diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [masses, mass_type] = screenNegMasses(masses, threshold, mass_type)
    % Screens negative masses and turns to 0 if above threhold and turns
    % to NaN if below threshold, threshold should be a negative value
    % also updates the mass_type to 5 if negative
    neg_mass_ind = find(masses<threshold);
    if isempty(neg_mass_ind) == 0 % there are negative masses below accepance 
        masses(neg_mass_ind) = NaN; % set neg values below threshold to NaN
        mass_type(neg_mass_ind) = 5; % set mass type to 5 for neg masses
        fprintf('WARNING: There are %1.0f out of %1.0f new filter masses that are negative \n', length(neg_mass_ind), length(masses))
    end
    clear neg_mass_ind
end  

function DataOut = Reformat(DataIn)

    if isa(DataIn,'double') == 1
        DataOut = num2cell(DataIn);
    else
        DataOut = DataIn;
    end
    for i = 1:length(DataOut)
        if ismember(DataOut(i),'NaN')
            DataOut{i}= 'Unknown';
        end
    end
end