%% %%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: read and sort ICP-MS data obtained from the water studies lab
% with no pre-processing; place ICP-MS data and relevant flags in the
% Master data file for each SPARTAN site that has data contained in the
% ICP-MS file; move ICP-MS files from "New" to relevant  folder

% Written by: Crystal Weagle
% Created: 2 December 2018

% EDITS:
% Haihui Zhu 2021-06-01

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear ; clc
fprintf('%s \n', datestr(today))
addpath('UtilityFunctions')

%% %%%%%%% USER SWITCHES %%%%%%%%%%
% Setup directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

site_details = readtable(sprintf('%s/Site_Sampling/Site_details.xlsx',direc));
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));

Elog_filename = sprintf('%s/Analysis_Data/E-Logs/ICP-MS E-Log.xlsx',direc); % contains extraction volumes for individual filters
[elog_status, Sheets] = xlsfinfo(Elog_filename); 

diary(sprintf('%s/Public_Data/Data_Processing_Records/%s_ICPMS_Record',direc,datestr(today,'yyyy-mm')))


%% --------- Find and read data for writing to master file ---------
% files in the raw file directory should be striaght from the IC system, unmodified
% this will read all files in the raw directory

names = dir(sprintf('%s/*.xlsx',direc_ICPMS)); names2 = struct2cell(names)';
files = {names.name}; files_bytesB = string(names2(:,4)); files_bytes = double(files_bytesB);
clear names names2 files_bytesB
files = files';

% column indices for metal data
metals = {'7Li' '24Mg' '27Al' '31P' '47Ti' '51V' '52Cr' '55Mn' '56Fe' '59Co' '60Ni' '65Cu'...
          '66Zn' '75As' '82Se' '107Ag' '111Cd' '121Sb' '137Ba' '140Ce' '208Pb'};

for ICP_file = 1:length(files)
    if files_bytes(ICP_file) < 10000 % do not read small files
        continue
    end
    
    filename = sprintf('%s/%s',direc_ICPMS,files{ICP_file,1});
    try
    [num,txt,raw] = xlsread(filename,'Sheet1');
    catch
        [num,txt,raw] = xlsread(filename,'raw');
    end
    fprintf('Reading %s \n', files{ICP_file,1});
    
    test = 0; % zero means the file does not have any test data contained, there is a loop below that changes this to 1 if test data is detected and will save the file in a designated folder
    
    % find indices for data in file (this is a time consuming loop - improvements to be made?)
    label_indexB = []; sites = [];
    for i = 1:length(Site_codes)
        index = find(contains(txt(:,3),Site_codes{i}));
        if isempty(index) == 0
            label_indexB = [label_indexB; index];
            sites = [sites; Site_codes(i)];
        end
        clear index
    end
    clear i
    
    LB = ["LB","BL"];  % label for (lab) blank [Revision1 removes 'BL' and 'LB' seperately. This change make the script more concise]

    % indexing labels and mean values from data file
    labelsB = txt(label_indexB,3); % index of labels in txt matrix
    LB_indB = find(contains(txt(label_indexB,3),LB)); % index for lab blank labels in label_indexB
    if isempty(LB_indB) == 0 
        labelsC = char(removerows(labelsB,LB_indB));
        LB_ind = label_indexB(LB_indB); % index of LB labels in txt matrix
        labels = labelsC(:,1:10);
        label_index = removerows(label_indexB,LB_indB);
    else
        labelsC = char(labelsB);
        labels = labelsC(:,1:10);
        label_index = label_indexB;
    end
    
    clear labelsC label_indexB labelsB
    
    
    % find row indices for mean trace element data
    x_idx = find(contains(txt(:,1), 'x'));
    
    for i = 1:length(label_index)
        belowrange = find(x_idx < label_index(i) - 1);
        aboverange = find(x_idx > label_index(i) + 4);
        bad_idx = vertcat(belowrange,aboverange);
        index_means(i,1) = removerows(x_idx,bad_idx);
        clear belowrange aboverange bad_idx
    end
    
    metal_row = find(contains(txt(:,4), '7Li'));
    missing_idx = [];
    
    index_metals = zeros(1,length(metals));
    for i = 1:length(metals)
        if isempty(find(contains(txt(metal_row,:),metals{i}))) == 0
            index_metals(i) = find(contains(txt(metal_row,:),metals{i}));
        else % if the metal in the metals list is not contained in the file
            missing_idx = [missing_idx; i];
        end
    end
    clear i
    
    if isempty(missing_idx) == 1
        ICP_data = cell2mat(raw(index_means,index_metals));
%         ICP_data = num(index_means, index_metals);
    else
        index_metals2 = removerows(index_metals', missing_idx);
        ICP_data = cell2mat(raw(index_means,index_metals2));
%         ICP_data = num(index_means, index_metals2); % actual data from the file. Index_metals2 is teh index to teh column with data in the file, not the metal in the metals list
    end
    
    clear index_metals2 
    
    %% Generate proper sample labels
    
    data_labels(length(label_index),:) = blanks(8); % create blank matrix for writing proper data labels 
    test_index = [];
    sites = char(sites);
    for i = 1:size(labels,1)
        temp_label = labels(i,(find(~isspace(labels(i,:)))));
        
        for j = 1:size(sites,1) % finds the location of the site index in the label
            idx = find(ismember(temp_label,sites(j,:)));
            if isempty(idx) == 0 && length(idx) == 4
                site = temp_label(idx);
                break
            end
        end

        if contains(temp_label,'X') == 1 || contains(temp_label,'x') == 1 || contains(temp_label,'-A')
                test_index = [test_index; i];
                test = 1; % set test == 1 to make sure file is saved in test directory for future reference
                continue
        end
        
        data_labels(i,1:4) = site;
        
        if length(temp_label) == 8 % likely XXXX-###, all the numbers and the dash are there!
            data_labels(i,5:8) = temp_label(5:8);
            
        elseif length(temp_label) == 7 % likely XXXX### OR XXXX-##
            if contains(temp_label,'-') == 1 % XXXX-## OR XXXX-0#
                data_labels(i,5) = temp_label(5);
                if temp_label(6) == 0 % XXXX-0#
                    data_labels(i,6:7) = '00';
                    data_labels(i,8) = temp_label(7);
                else % XXXX-##
                    data_labels(i,6) = '0';
                    data_labels(i,7:8) = temp_label(6:7);
                end
            else % XXXX###
                data_labels(i,5) = '-';
                data_labels(i,6:8) = temp_label(5:7);
            end
            
        elseif  length(temp_label) == 6 % likely  XXXX## OR XXXX-#
            if contains(temp_label,'-') == 1 % XXXX-#
                data_labels(i,5) = temp_label(5);
                data_labels(i,6:7) = '00';
                data_labels(i,8) = temp_label(6);
            else % XXXX##
                if temp_label(5) == 0 % XXXX0#
                    data_labels(i,5:7) = '-00';
                    data_labels(i,8) = temp_label(6);
                else %XXXX##
                    data_labels(i,5:6) = '-0';
                    data_labels(i,7:8) = temp_label(5:6);
                end
            end
            
        elseif length(temp_label) == 5 % likely XXXX#
            data_labels(i,5:7) = '-00';
            data_labels(i,8) = temp_label(5);
            
        elseif length(temp_label) == 9 % this would be an unusual label
            if idx(1) == 6 % 1####SITE
                data_labels(i,1:4) = temp_label(6:9); % site ID is last 4 in ICP-MS data label
                data_labels(i,5) = '-';
                data_labels(i,6:8) = temp_label(3:5);
            else
                test_index = [test_index; i];
                test = 1; % set test == 1 to make sure file is saved in test directory for future reference
            end
            
        elseif length(temp_label) == 10 % this would be an unusual label
            if idx(1) == 7 % 1####-SITE
                data_labels(i,1:4) = temp_label(7:10); % site ID is last 4 in ICP-MS data label
                data_labels(i,5) = '-';
                data_labels(i,6:8) = temp_label(3:5);
            else
                test_index = [test_index; i];
                test = 1; % set test == 1 to make sure file is saved in test directory for future reference
            end
        end
        clear temp_label idx
        
    end
    clear i
    
    if isempty(test_index) ==0
        data_labelsB = data_labels; clear data_labels
        data_labels = char(removerows(cellstr(data_labelsB), test_index));
        clear data_labelsB
        
        ICP_dataB = ICP_data; clear ICP_data
        ICP_data = removerows(ICP_dataB, test_index);
        clear ICP_dataB
    end

    if isempty(data_labels) == 0
    site_IDs = unique(data_labels(:,1:4),'rows'); % IDs of sites with data in data file
    data_labels_cell = cellstr(data_labels);
    elseif test == 1 && isempty(data_labels) ==1
        file_destination = ('/data1/weagle/SPARTAN/Filter_data_raw/Metal_data/ICP-MS/Archived/Test_data/');
        status= copyfile(sprintf('%s',filename), sprintf('%s',file_destination));
        disp('Test data detected, file has been saved to test data folder, no actual samples found');
        disp('Moving on to next file');
        delete(filename); disp('File deleted from New data folder');
        clear status file_destination
        continue
    else
        disp('No samples found in this file')
        disp('Moving on to next file')
        continue
    end
    %% Find and write to master data files

    for i = 1:size(site_IDs,1)
        master_file = sprintf('%s/%s_master.csv',direc_master,site_IDs(i,:));
        if exist(master_file,'file') == 0
            fprintf('WARNING: No master file file exists for %s \n', site_IDs(i,:))
            fprintf('Looking for master file for next site in IC data file \n')
            continue
        elseif exist(master_file,'file') == 2
            fprintf('Master file found for %s \n', site_IDs(i,:))
            samples_ICPfile_ind = find(contains(data_labels_cell,site_IDs(i,:))); % index of samples for site_ID(i) in ICP-MS file
            
            Master_data_initial = readtable(master_file); % master data table with all values prior to manipulating in this script
            % columns used in this script 
            Master_IDs = cell(table2array(Master_data_initial(:,1)));
            Master_data = table2array(Master_data_initial(:,7:end-1)); % creates a master data matrix to add new mass data to
            % no need untill writing master file
            Master_Barcode = table2array(Master_data_initial(:,2));
            CartridgeIDs_master = table2array(Master_data_initial(:,3));
            LotIDs_master = table2array(Master_data_initial(:,4));
            projectID = table2array(Master_data_initial(:,5));
            hours_sampled = table2array(Master_data_initial(:,6));
            Master_flags = char(table2array(Master_data_initial(:,end)));
            
            
            masterID_sample_ind = [];
            for k = 1:length(samples_ICPfile_ind)
                if isempty(find(contains(Master_IDs, data_labels_cell(samples_ICPfile_ind(k))))) == 0
                    masterID_sample_ind(k) = find(ismember(Master_IDs, data_labels_cell(samples_ICPfile_ind(k)))==1); % finds the row index in the master file for samples in ICP file
                    samples_ICP_Master(k) = samples_ICPfile_ind(k);
                else
                    continue
                end
            end
            
            if isempty(masterID_sample_ind) == 0
                
                
                Elog_SheetName = char(Sheets(contains(Sheets,site_IDs(i,:)))); % find out the sheet for current site ( site_IDs(i) )
                Elograw = readtable(Elog_filename,'Sheet', Elog_SheetName); 
                ElogAnalysisID = table2array(Elograw(:,3)); % Analysis ID (filter ID) used to match Master ID
                ElogVolume = table2array(Elograw(:,6)); % Volume data to calculate mass;
                
                extraction_volumes_idxB = find(contains(ElogAnalysisID, Master_IDs(masterID_sample_ind))); % analysis id 
                extraction_volumes_idx = removerows(extraction_volumes_idxB,find(contains(ElogAnalysisID(extraction_volumes_idxB), 'X'))); clear extraction_volumes_idxB
                extraction_volumes = ElogVolume(extraction_volumes_idx);
                
                if length(extraction_volumes)>length( samples_ICP_Master)
                    extraction_volumes = extraction_volumes(1:length(samples_ICP_Master));
                elseif isempty(extraction_volumes)
                    fprintf('WARNING: No extraction volumes found for %s data\n', site_IDs(i,:))
                    fprintf('Moving on to next site, no data written to %s master file \n', site_IDs(i,:))
                    continue
                end
                    
                
                Master_ICP_totalcolumns = 26:46;
                
                if exist('missing_idx','var') == 1
                    ICP_columns = removerows(Master_ICP_totalcolumns',missing_idx);
                else 
                    ICP_columns = Master_ICP_totalcolumns;
                end
                
                Master_data(masterID_sample_ind, ICP_columns) = ICP_data(samples_ICP_Master,:).*extraction_volumes;
                
                clear  index 
                
     Titles = {'FilterID','Filter_Barcode','CartridgeID','LotID','projectID',...
              'hours_sampled','Mass_type','start_year','start_month','start_day',...
              'start_hour','stop_year','stop_month','stop_day','stop_hour',...
              'mass_ug','Volume_m3','SSR_BC_ug',...
              'IC_F_ug','IC_Cl_ug','IC_NO2_ug','IC_Br_ug','IC_NO3_ug','IC_PO4_ug',...
              'IC_SO4_ug','IC_Li_ug','IC_Na_ug','IC_NH4_ug','IC_K_ug','IC_Mg_ug',...
              'IC_Ca_ug',...
              'Li_ICP_ng','Mg_ICP_ng','Al_ICP_ng','P_ICP_ng','Ti_ICP_ng',...
              'V_ICP_ng','Cr_ICP_ng','Mn_ICP_ng','Fe_ICP_ng','Co_ICP_ng',...
              'Ni_ICP_ng','Cu_ICP_ng','Zn_ICP_ng','As_ICP_ng','Se_ICP_ng',...
              'Ag_ICP_ng','Cd_ICP_ng','Sb_ICP_ng','Ba_ICP_ng','Ce_ICP_ng',...
              'Pb_ICP_ng',...
              'Al_XRF_ng','Si_XRF_ng','S_XRF_ng','Cl_XRF_ng','K_XRF_ng',...
              'Ca_XRF_ng','Ti_XRF_ng','V_XRF_ng','Fe_XRF_ng','Zn_XRF_ng',...
              'Ce_XRF_ng','Pb_XRF_ng','As_XRF_ng','Co_XRF_ng','Cr_XRF_ng',...
              'Cu_XRF_ng','Mg_XRF_ng','Mn_XRF_ng','Ni_XRF_ng','Sb_XRF_ng',...
              'Rb_XRF_ng','Sr_XRF_ng','Cd_XRF_ng','Se_XRF_ng','Sn_XRF_ng',...
              'BC_HIPS_ug','EC_FTIR_ug','OC_FTIR_ug','method_index','Flags'};


                        fileID = fopen(sprintf('%s/%s_master.csv',direc_master,site_IDs(i,:)),'w');
                        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,      %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,     %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,  %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,   %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,   %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,   %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,   %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,   %s,%s\n',Titles{1,1:length(Titles)});

                if isa(LotIDs_master,'cell') == 1
                    for h=1:size(Master_data,1)
                        fprintf(fileID,'%s,%s,%s,%s,%s,%f,%6.1f,%f,%f,%f,   %f,%f,%f,%f,%f,%6.2f,%f,%f,%f,%f,  %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,  %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,   %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,   %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,   %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,   %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,   %f,%s\n',...
                    char(Master_IDs(h,:)), char(Master_Barcode(h,:)),char(CartridgeIDs_master(h,:)), char(LotIDs_master(h,:)), char(projectID(h)), hours_sampled(h), Master_data(h,:),Master_flags(h,:));
                    end
                    
                elseif isa(LotIDs_master,'double') == 1
                    for h=1:size(Master_data,1)
                        fprintf(fileID,'%s,%s,%s,%s,%s,%f,%6.1f,%f,%f,%f,   %f,%f,%f,%f,%f,%6.2f,%f,%f,%f,%f,  %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,  %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,   %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,   %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,   %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,   %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,   %f,%s\n',...
                    char(Master_IDs(h,:)), char(CartridgeIDs_master(h,:)), LotIDs_master(h,:), char(projectID(h)), hours_sampled(h), Master_data(h,:),Master_flags(h,:));
                    end
                end
                
                fclose(fileID);
                fprintf('Finished writing ICP-MS data to %s master data file \n', site_IDs(i,:))
                
            else
                fprintf('Samples in ICP-MS data not found in %s Master file \n', site_IDs(i,:));
                disp('Moving on to next site')
            end
            
        end
        
            clear Master_data Master_data_initial txt_elog num_elog raw_elog extraction_volumes
            clear Master_IDs masterID_sample_ind samples_ICPfile_ind samples_ICP_Master
            clear Master_flags Master_flags_initial master_file hours_sampled CartridgeIDs_master 
            clear LotIDs_master hours_sampled ICP_columns Master_ICP_totalcolumns
        
    end
    
    for i = 1:size(site_IDs,1)
        file_destination = sprintf('%s/Archived_Filter_Data/ICP-MS/%s/',direc,site_IDs(i,:));
        status(i)= copyfile(sprintf('%s',filename), sprintf('%s',file_destination));
        clear file_destination
    end
    
    if test == 1
        file_destination = (sprintf('%s/Archived_Filter_Data/ICP-MS/Test_data/',direc));
        status(i)= copyfile(sprintf('%s',filename), sprintf('%s',file_destination));
        clear file_destination
        disp('Test data detected, file has been saved to test data folder')
    end
    
    if isempty(find(status == 0)) == 1
        delete(filename)
        disp('File has been deleted from new ICP-MS data folder')
    elseif test == 1
        delete(filename)
    else
        disp('Cannot delete file from new ICP-MS data folder as it has not been moved to all site folders')
    end
    
    clear CartridgeIDs_master data_index data_labels data_labels_cell fileID filename h i ICP_data ICP_file index_means index_metals
    clear k labels test LotIDs_master test_index num num_sz raw site_IDs site status Titles txt txt_sz 
    clear raw_elog num_elog txt_elog metal_row first_sample extraction_volumes_idx label_index sites missing_idx
    
end

 fprintf('END of ICP-MS data processing for %s \n', datestr(today))
 diary off

disp('Program finished')







