%% %%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: read and perform QA/QC on XRF data straight from the XRF
% system with no pre-processing; put XRF data and relevant flags in the
% Master data file for each SPARTAN site that has data contained in the XRF
% file; move XRF files from "Raw" to relevant site folder

% NOTE: Not filter flag should be added from XRF-NEW file. Filter flag
% should be always added by Script_B from xxxx_dates_flows.xlsx

% Written by: Crystal Weagle
% Created: 16 January 2020

% EDITS:

% 25-Aug-2021 Haihui Zhu:
%   1. Implement ReadMaster function and WriteToMaster function so that 
% future change in the master file can be automatically applied without 
% making change in this script.
%   2. Replace xlsread to Readtable
%   3. Use force_movefile function instead of 'moveFile' to move processed files.

% 11-Dec-2020 Emmie Le Roy: changed directories to storage1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear ; clc
addpath('UtilityFunctions')

%% Directories and User Switches

% Setup directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

direc_new = strcat(direc,'/Analysis_Data/XRF/NEW');
direc_archive = strcat(direc,'/Analysis_Data/Archived_Filter_Data/XRF_data/');
direc_master = strcat(direc,'/Analysis_Data/Master_files');
direc_sampling = strcat(direc,'/Site_Sampling');

% The diary function saves the processing history into an annual record
diary(sprintf('%s/Public_Data/Data_Processing_Records/%s_XRF_Record',direc,datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))


%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% column indices for metal data (element of interest is in brackets)
metals = {'C (Na)' 'C (Al)' 'C (Si)' 'C (S)' 'C (Cl)' 'C (K)' 'C (Ca)' 'C (Ti)' 'C (V)' 'C (Fe)' 'C (Zn)' 'C (Ce)' 'C (Pb)' 'C (As)' ...
    'C (Co)' 'C (Cr)' 'C (Cu)' 'C (Mg)' 'C (Mn)' 'C (Ni)' 'C (Sb)' 'C (Rb)' 'C (Sr)' 'C (Cd)' 'C (Se)' 'C (Sn)'};

% From Feb19 2020 calculation
% metal_accLimit = [0.008 0.026 0.003 0.049 0.101 0.041 0.031 0.001 0.008 0.003 0.003 0.006 0.002 ...
%   0.001 0.002 0.007 0.077 0.002 0.002 0.058 0.006 0.008 0.015 0.001 0.082]; % acceptance limit for metals in XRF data files, units ug/cm2
% BW commented this out 07-09-21 as we are no longer setting values to
% zero if below acceptance limit, also accetance limits/MDL values and
% way to calculate this changing

r = (21.2/10)/2; % diameter is 21.2 mm, need to convert to cm and then to radius
filter_area = pi*r^2;

%-------------   SITE DETAILS   --------------
site_details = readtable(strcat(direc_sampling,'/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));

%% --------- Find and read data files in "NEW" file directory ---------
% files in the raw file directory should be straight from the XRF system, unmodified
% this section will find the sample IDs in the XRF file, isolate the
% Site IDs for indexing to master files, and select out the data in the XRF file

files = getFiles(direc_new);

for XRF_file = 1:length(files)
    
    % Do not read hidden files
    char_file = char(files(XRF_file));
    if char_file(1)=="."
        continue
    end
    clear char_file
    
    filename = sprintf('%s/%s',direc_new,files{XRF_file,1});
    fprintf('Reading %s \n', files{XRF_file,1});
    
%     [num,txt,raw] = xlsread(filename);
    % update to readtable [Haihui]
    initialdata = readtable(filename,'PreserveVariableNames',true);
    XRFtitle = initialdata.Properties.VariableNames;
    
    rows = find(~isnan(initialdata.Pos));
    xrf_IDs = initialdata.Ident(rows);  % filter ID
    xrf_IDs_char = char(xrf_IDs);
    site_ID = unique(xrf_IDs_char(:,1:4),'rows'); clear xrf_IDs_char
    
    if size(site_ID,1) > 1
        warning('File skipped: There are more than 1 site included in %s. Please split.',files{XRF_file,1})
        continue
    end
    
    % ------- Find index for metals in xrf file -----
    index_metals = zeros(1,length(metals));
    for i = 1:length(metals)
        index_metals(i) = find(contains(XRFtitle,metals{i}));
    end
    clear i
    
    xrf_dataB = table2array(initialdata(rows, index_metals)); % isolate XRF data, units of ug/cm2
    
    % BW commented out July 12 2021 as we are no longer setting values below
    % the acceptance limit to zero for XRF. A section will be added where
    % the MDL and analytical detection limits will be recorded.
    % for i = 1:length(index_metals)
    %    below_mdl_idx = find(xrf_dataB(:,i) < metal_accLimit(i));
    %    if isempty(below_mdl_idx) == 0
    %        xrf_dataB(below_mdl_idx,i) = 0;
    %    end
    %    clear below_mdl_idx
    % end
    
    xrf_data = (xrf_dataB.*filter_area)*1000; % multiply by filter surface area to get in units of ug, multiply by 1000 to get in final units of ng
    
    clear xrf_dataB below_mdl_idx index_metals num txt raw
    
    %% The ReadMaster function checks for existance of a master file and read it  
    master_file = sprintf('%s/%s_master.csv',direc_master,site_ID);
    [Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, projectIDs_master,...
        Master_hours, Mass_type, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
        Master_carbon, Master_Method, Master_flags] = ReadMaster ( master_file, site_ID );
    
    % writing new XRF data to Master_XRF
    master_idx = find(contains(Master_IDs, xrf_IDs));
    Master_XRF(master_idx, :) = xrf_data;

    % Write master file 
    WriteToMaster(Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, projectIDs_master,...
        Master_hours, Mass_type, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
        Master_carbon, Master_Method, Master_flags,direc_master,site_ID)
    
    % move raw data to archive
    file_destination = strcat(direc_archive,sprintf('%s/',site_ID));
%     status = moveFile(filename,file_destination);
    [status,msg,msgid] = force_movefile(filename, file_destination);
    clear file_destination
    

    %% adding those filter to MDL_Unc_Ref.mat
    MDLfname = sprintf('%s/Analysis_Data/XRF/MDL_Unc_Ref.mat',direc);
    load(MDLfname,'FilterGroup','Ref_Dates','Ref_Values')
    today = datetime('today');
    uniRefDates = unique(Ref_Dates); uniRefDates(end+1) = 999999;

    SampleDate = floor(convertTo(table2array(initialdata(:,contains(XRFtitle,'Time'))),'datenum'));
    FilterID = table2array(initialdata(:,contains(XRFtitle,'Ident')));

    for i = 1:numel(uniRefDates)-1
        if i == 1 % Everything before the second Reference MDL refers to the first Reference
            Ind = find(SampleDate<uniRefDates(i+1));
        else
            Ind = find(SampleDate<uniRefDates(i+1) & SampleDate>=uniRefDates(i));
        end
        % exclude rows don't contain filterID
        Ind(contains(FilterID(Ind,1),'?')) = [];

        % make sure the filters do not exist in FilterGroup
        tfilters = FilterGroup{2,i};
        overlap_ind = find(ismember(FilterID(Ind,1),tfilters)==1);
        Ind(overlap_ind) = [];
        rows = size(FilterGroup{2,i},1)+1:size(FilterGroup{2,i},1)+numel(Ind);
        FilterGroup{2,i}(rows,1) = FilterID(Ind,1);
        clear tfilters
    end

    date_mark = today;
    save(MDLfname,'FilterGroup','Ref_Dates','Ref_Values','date_mark')


    clear xrf_IDs XRF_file xrf_data Master_* initialdata 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('END of XRF data processing for %s \n', datestr(today))
diary off
disp('Program finished')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function files = getFiles(directory)
% retrieve all files in specified directory
names = dir(sprintf('%s/*.xlsx',directory));
files = {names.name}; files = files';
end


function [row, col] = findIndexContaining(string, file)
% find row and column index of string in file
index = find(contains(file,string));
[row, col] = ind2sub(size(file),index); clear index
end
