%% %%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: to read master data files and examine mass and ions 
% Data validation steps include: 
% 1. Apply a default blank correction for metals determined by ICP-MS
% 2. Search for any flag or signs of invalid sampling (pinholes, bad flow
% rate, etc) and remove the filter for further analysis 
% 3. Calculate BC using SSR when HIPS measurement is not available. 
% 4. Uses flags from IC and the 'Bad_sodium' list to remove invalid IC measurements
% 5. Check sum of measured species vs. mass collected and mark filter as invalid if sum of species is greater than collected mass
% 6. Reconstruct fine mass using measured constituents
% 5. Determine filter-specific kappa values based on measured constituents
% 6. Generate data products for research and posting on SPARTAN website.
%    Data products generated include:
%     - all measured components (dry, 0% RH), with any necessary corrections applied for PM2.5, PM10, and when applicable, PMcoarse
%     - reconstructed fine mass with residual matter and PBW at 35 %RH
% 
% Written by: Crystal Weagle
% Created: 22 January 2019
%
% EDITS (sorted from most to least recent):
% -- 2023-03-19, by Haihui Zhu:
%     1. Adding the calculation for Unc in findMDLUnc function.
%     2. Adding section to calculate HIPS BC using the HISP-SSR curve
%     3. Update calculation of PBW that matches GEOS-Chem and halves SNA mass growth
%     4. Adding lines to remove bad Na according to Bad_sodium_data_in_2017.xlsx
%     4. Rearrange sections for readability and conciseness
%    *5. Update public data to Version 2.3
% 
% -- 2022-08-18, by Haihui Zhu:
%     1. Apply a 'copyflag' function so that certian filter flag will be 
%     added to the public files. A list of flags can be found the [flags] 
%     sheet in sampling_parameter_methods.xlsx 
%     2. Replace SSR in the public files (PM25 speciation and RCFM) with 
%     HIPS whenever HIPS is available.
%
% -- 2022-06-22, by Haihui Zhu:
%    *1. Update public data to Version 2.2
%     2. Apply Xuan Liu's equation to calculate fine soil [Liu, X. 2022 JGR]
%     3. Update TEO calculation when using XRF metals.
%     4. Use updated sampling_parameter_methods.xlsx (for soil and TEO). 
%     path: SOPs/Public SOPs/
%     5. Replace SS_wet with 2 if they are higher than 2 ug/m3. (Previously
%     replaced by NaN that prevent OM calculation)  
%
% -- 2022-05-03, by Haihui Zhu:
%     1. Add lines to exclude filters labeled as 'Need to be excluded' 
% 
% -- 2022-02-15, by Haihui Zhu: 
%     1. Simplify the script for readability.
%     2. Fix the error when Na_orig does not exist.
% 
% -- 2021-10-24, by Haihui Zhu: 
%    *1. Update public data to Version 2.1
%     2. Applied ReadMaster function for easy reading master files. This 
%        funciton also provide the title of the master file, which help to
%        identify data columns. Using number istead of title as column 
%        identifier can lead to error.
%     3. Add Flag (VIH, FS) and filter ID to the public files
%     4. Add MDL and Unc for XRF to PMxx_Speciation  
% 
% -- 2021-07-19, by Haihui Zhu: update master_data/XX_Public_date column identifiers
% 
% -- 2020-04-08, by Crystal Weagle: Setting threshold for Na concentration to
% meet 2 ug/m3 RCFM sea salt threshold. Any value above 2 ug/m3 is set to
% site-average. Change is applied to all SPARTAN data, historical and
% recent --
% 
% -- 2020-02-25, by Crystal Weagle: Setting threshold for sea salt
% concentration to 2 ug/m3 in RCFM --
% 
% -- 2020-02-25, by Crystal Weagle: Added section to create site average
% contribution of major components --
% 
% -- 2020-02-24, by Crystal Weagle: Updated indexing to reference correct
% columns in accordance with new master_file.csv layout (includes, Project
% ID, XRF columns, and UV-Vis and FTIR columns) --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc
addpath('UtilityFunctions')

%% %%%%%%% USER SWITCHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataVersion = 'Data version 2.4'; % will be in the first line of the public files

% Setup directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

direc_master = strcat(direc,'/Analysis_Data/Master_files');
direc_output =  strcat(direc,'/Public_Data');

% optional figure switch:
additional_pie = 0;

% ---- the diary function saves the processing history into monthly records ----
diary(sprintf('%s/Public_Data/Data_Processing_Records/%s_MasterFile_Processing_Record.txt',direc, datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))

% ---- Load Site Details and Sampling Methods ----
Sampling_Parameters_Methods = strcat(direc,'/SOPs/Public SOPs/Sampling_Parameters_Methods_2.3.xlsx');

PM25_mass_para = readtable(Sampling_Parameters_Methods,'Sheet','PM2.5 mass','PreserveVariableNames',true);
PM25_IC_para   = readtable(Sampling_Parameters_Methods,'Sheet','PM2.5 water-soluble ions','PreserveVariableNames',true);
PM25_metals_para = readtable(Sampling_Parameters_Methods,'Sheet','PM2.5 trace elements','PreserveVariableNames',true);
PM25_carbon_para = readtable(Sampling_Parameters_Methods,'Sheet','PM2.5 BC_OC_EC','PreserveVariableNames',true);

PM10_mass_para = readtable(Sampling_Parameters_Methods,'Sheet','PM10 mass','PreserveVariableNames',true);
PM10_IC_para   = readtable(Sampling_Parameters_Methods,'Sheet','PM10 water-soluble ions','PreserveVariableNames',true);
PM10_metals_para = readtable(Sampling_Parameters_Methods,'Sheet','PM10 trace elements','PreserveVariableNames',true);
PM10_carbon_para = readtable(Sampling_Parameters_Methods,'Sheet','PM10 BC_OC_EC','PreserveVariableNames',true);

PMc_parameters  = readtable(Sampling_Parameters_Methods,'Sheet','PMc','PreserveVariableNames',true);
RCFM_parameters = readtable(Sampling_Parameters_Methods,'Sheet','PM2.5 RCFM','PreserveVariableNames',true);

Flag_parameters = readtable(Sampling_Parameters_Methods,'Sheet','Flags','PreserveVariableNames',true);

site_details = readtable(strcat(direc,'/Site_Sampling/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));
Site_countries = table2array(site_details(:,2));
latitudes = table2array(site_details(:,5));
longitudes = table2array(site_details(:,6));
elevations = table2array(site_details(:,7));

% ---- addtional output initialization ----
HIPS_Mesured_BC = nan(size(Site_codes));
HIPS_Est_BC = nan(size(Site_codes));
HIPS_total_BC = nan(size(Site_codes));
SSR_total_BC = nan(size(Site_codes));

% calculate reportable data #:
DataNum_tab_rows = {'City','PM2.5','Water','SO4','NH4','NO3','SS','Dust','TEO','BC','OC','RM'};
tSiteDataNum = nan(length(Site_codes),length(DataNum_tab_rows)-1);

% ---- Load the MDL & Unc matrix for XRF ----  
MDLfname = sprintf('%s/Analysis_Data/XRF/MDL_Unc_Ref.mat',direc);
load(MDLfname,'Ref_Dates','Ref_Values','FilterGroup')


% ---- Load the list of bad sodium cartridges ----
BadNa = sprintf('%s/Analysis_Data/Ion_Chromatography/Bad_sodium_data_in_2017.xlsx',direc);
opts = detectImportOptions(BadNa);
opts.VariableNames = {'CartID'};
BadNa = table2array(readtable(BadNa,opts));

% Mass type variable meaning:
% 0 = cartridge (traveling) blank
% 1 = PM2.5
% 2 = PM10
% 3 = PMcoarse
% 4 = saturated nuclepore, filter measurements invalid
% 5 = negative mass, filter measurements invalid
% 6 = invalid flow, filter measurements invalid

% fileID = fopen('Zero_values.txt', 'w');
% if fileID == -1
%     error('File could not be opened for appending');
% end
% fprintf(fileID,'Label    S_year S_month S_day S_hour E_year E_month E_day E_hour Parameter Value  Description\n');
% fclose(fileID);

%% Read master files and process public files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for loc = 35:35 %(Site_codes)
 for loc =  37:numel(Site_codes)%[1:33 35]
    
    % Read Master Files
    master_file = sprintf('%s/%s_master.csv',direc_master,Site_codes{loc});
    [Titles,Master_IDs,   Master_Barcodes, Master_CartridgeIDs, Master_LotIDs, Master_projectIDs,  Master_hours,  Master_masstype, ...
            Master_dates, Master_mass,     Master_IC,           Master_ICP,    Master_XRF,...
            Master_carbon, Master_Method,  Master_flags] = ReadMaster(master_file,Site_codes{loc} );
    
    Vol_Col = 2; % 2nd col in mass matrix is volume
    IC_title = Titles(contains(Titles,'IC_')); 
    ICP_title = Titles(contains(Titles,'ICP_')); 
    XRF_title = Titles(contains(Titles,'XRF_')); 

    if isempty(Master_mass)
        continue
    end
    
    % Read IC MDL Files
    ic_mdl_file = sprintf('%s/Analysis_Data/Ion_Chromatography/MDL/%s_IC_MDL.xlsx',direc,Site_codes{loc});
    if exist(ic_mdl_file,'file')
        ic_mdl = readtable(ic_mdl_file);
        ic_mdl_filterid = ic_mdl.Filter_ID;
        ic_mdl_titles = ic_mdl.Properties.VariableNames;
        ic_mdl_titles(6:end) = {'Fluoride' 'Chloride' 'Nitrite'  'Bromide' 'Nitrate' 'Phosphate' 'Sulfate'...
            'Lithium' 'Sodium' 'Ammonium' 'Potassium' 'Magnesium' 'Calcium'};
    end
    %% Data validation step 1: Corrections for ICP-MS measurement
    %  Prcedure:
    %  1. sampled filters - blank = corrected trace elements;
    
    Cart_ALL = unique(Master_CartridgeIDs); % IDs of cartridges in master file
    
    for i = 1:length(Cart_ALL)
        Cart_IND = find(strcmp(Master_CartridgeIDs, Cart_ALL(i)));
        Cart_BL = find(Master_masstype(Cart_IND) == 0); % finds the blanks for the cartridge
        if ~isempty(Cart_BL)
            if length(Cart_IND)  == 16 && isempty(Cart_BL) == 0 % for cartridges with 2 blanks, one teflon and one nuclepore, in positions 8 and 16
                Cart_BL = [8 16];
                Master_ICP(Cart_IND(Cart_BL(1)-7):Cart_IND(Cart_BL(1) - 1),:) = Master_ICP(Cart_IND(Cart_BL(1)-7):Cart_IND(Cart_BL(1) - 1),:) - Master_ICP(Cart_IND(Cart_BL(1)),:); % corrects PM2.5
                Master_ICP(Cart_IND(Cart_BL(2)-7):Cart_IND(Cart_BL(2) - 1),:) = Master_ICP(Cart_IND(Cart_BL(2)-7):Cart_IND(Cart_BL(2) - 1),:) - Master_ICP(Cart_IND(Cart_BL(2)),:); % corrects PMcoarse
                
            elseif ismember(Site_codes{loc} ,{'SGSU'}) && Cart_BL == 8 % will only be one blank (teflon) and in position 8
                Master_ICP(Cart_IND(Cart_BL-7):Cart_IND(Cart_BL - 1),:) = Master_ICP(Cart_IND(Cart_BL-7):Cart_IND(Cart_BL - 1),:) - Master_ICP(Cart_IND(Cart_BL),:); % only PM2.5 at this site
                
            elseif length(Cart_IND) == 8 && isempty(Cart_BL) == 0 % SS5 station with blanks in filter slot 7
                Master_ICP(Cart_IND(Cart_BL(end)-6):Cart_IND(Cart_BL(end) - 1),:) = Master_ICP(Cart_IND(Cart_BL(end)-6):Cart_IND(Cart_BL(end) - 1),:) - Master_ICP(Cart_IND(Cart_BL(end)),:); % corrects PM2.5 filters
                Master_ICP(Cart_IND(Cart_BL(end)+1),:) = Master_ICP(Cart_IND(Cart_BL(end)+1),:) - Master_ICP(Cart_IND(Cart_BL(end)),:); % corrects PM10 filter
                
            end
        end
        clear Cart_BL Cart_IND
    end
    
    clear Cart_ALL

%     for i = 1:size(Master_ICP,2) % ICP-MS trace metal correction
%         Master_ICP(Master_ICP(:,i) <0,i) = 0; % set any negative metals to zero (due to correction above)
%     end

    %% Data validation step 2: Exclude bad filters or filter with invalid sampling hour, method, or volume
    Ind = find(contains(Master_flags,{'Need to be excluded'})==1);
    if ~isempty(Ind)
        fprintf('Following filters excluded for public files as indicated by filter flags:\n')
        disp(Master_IDs(Ind))

          [ Master_IDs, Master_Barcodes, Master_CartridgeIDs, Master_LotIDs, Master_projectIDs,...
            Master_hours, Master_masstype, Master_dates, Master_mass, Master_IC, ...
            Master_ICP, Master_XRF, Master_carbon, Master_Method, Master_flags] = ...
            RemoveRows (Ind, ...
            Master_IDs, Master_Barcodes, Master_CartridgeIDs, Master_LotIDs, Master_projectIDs,...
            Master_hours, Master_masstype, Master_dates, Master_mass, Master_IC, ...
            Master_ICP, Master_XRF, Master_carbon, Master_Method, Master_flags);
    end

    % find valid cols: vol exist and sampled hr ~=0
    vol_idx = find(Master_mass(:,Vol_Col) ~= 0 & ~isnan(Master_mass(:,Vol_Col)) & Master_hours ~=0 & ~isnan(Master_hours));
     
    [Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,   Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method, Master_flags]...
     = SortData (vol_idx,...
    Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,   Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method, Master_flags);

    clear vol_idx

    % remove method is NaN (SSe_ID is NaN, refer to script B)
    non_nan_idx = find(~isnan(Master_Method(:,1)));
 
    [Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,   Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method, Master_flags]...
     = SortData (non_nan_idx,...
    Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,   Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method, Master_flags);

    clear non_nan_idx

    %% Data validation step 3: Calculate BC when HIPS-BC not available
    %  Set SSR Bc = -899 to NaN, set negative BC to 0.
    Master_carbon(Master_carbon(:,1) == -899,1) = NaN;
    Master_carbon(Master_carbon(:,1) < 0, 1) = 0; % set any negative BC to zero due to using blank R in calculation [rare]
    
    Ind = find( isnan(Master_carbon(:,2)) &  Master_masstype<3 & Master_carbon(:,1)>0); % HIPS not available & Teflon filters & SSR-BC available
    HIPS_Mesured_BC(loc,1) = sum(~isnan(Master_carbon(:,2)));
    HIPS_Est_BC (loc,1)    = length(Ind);
    % apply the curve:
    a1=2163.54158502065;
    b1=0.00198619072508564;
    c1=-0.00467334844964624;
    a2=262.938831079902;
    b2=0.0123767218647736;
    c2=3.10315712581732;
    Master_carbon(Ind,2) =  a1*sin(b1* Master_carbon(Ind,1) + c1) + a2*sin(b2* Master_carbon(Ind,1) + c2); % SSR-BC unit = ug; HIPS-BC unit = ug;
    % add flag to those filters
    for ii = 1:length(Ind)
        Master_flags{Ind(ii)} = AddFlag (Master_flags{Ind(ii)}, 'HIPS-BC estimated using SSR-BC');
    end

    HIPS_total_BC(loc,1) = sum( ~isnan(Master_carbon(:,2)) );
    SSR_total_BC(loc,1)  = sum( Master_masstype<3 & Master_carbon(:,1)>0 ); 

    %% Data validation step 4: Mask out invalid IC data: bad correlation coefficient or bad Na conc in blanks
    
    % ----- Mark data as NaN when flags indicate invalid data -----
    ions = {'F','Cl','NO2','Br','NO3','PO4','SO4','Li','Na','NH4','K','Mg','Ca'};
    
    for k = 1:length(ions)
        for j = 1:size(Master_IC,1) 
            if contains(Master_flags(j,:),sprintf('S-%s',ions{k})) == 1 % the correlation coefficient for species i was below the required limit.
                Master_IC(j,k) = NaN; %set value to NaN because violated QA protocol
            end
            if contains(Master_flags(j,:),sprintf('Sx-%s',ions{k})) == 1 % The correlation coefficient for species i was missing from IC data file
                Master_IC(j,k) = NaN; %set value to NaN because violated QA protocol
            end
%            if contains(Master_flags(j,:),sprintf('qc1-%s',ions{i})) == 1% The average quality control standard check for species i exceeded +/- 10 % of specified concentration.
%                Master_IC(j,k) = NaN; %set value to NaN because violated QA protocol
%            end
        end
    end
    
    % find cartridges with bad Na (Na too high in blanks)
    Ind = find(contains(BadNa,Site_codes{loc}));
    if ~isempty(Ind)
        NaInd = find(contains(IC_title,'IC_Na'));
        for ii = 1:length(Ind)
            tCartridge = BadNa{Ind(ii)};
            MasterInd = find(ismember(Master_CartridgeIDs,tCartridge));
            Master_IC(MasterInd,NaInd) = NaN; % instead of NaN, using 0 so that no SS will be reported but ASO4 won't be affected.
        end
    end
    clear j k

    
    %% Data validation step 5: Sort data according to dates  
    
    % sort dates appropriately 
    temp_dates = datenum(Master_dates(:,1),Master_dates(:,2),Master_dates(:,3),0,0,0); % create datenum for sorting data based on start dates
    [dates_sort, idx] = sort(temp_dates);

    % Divide data by mass type
    [Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,  Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method,Master_flags]...
     = SortData (idx,...
    Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,  Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method,Master_flags);

    clear temp_dates idx dates_sort

    %% Data validation step 6: Convert mass to concentration (ug/m3); change unit of metals from ng to ug 
    if isempty(Master_mass) == 1
        disp('No Mass data detected')
        continue
    else
        % converting mass to concentration
        Master_dates = round(Master_dates); % need to round to nearest hour
        Master_mass(:,1) = Master_mass(:,1) ./Master_mass(:,Vol_Col);
        Master_IC = Master_IC./Master_mass(:,Vol_Col);
        Master_ICP = Master_ICP./Master_mass(:,Vol_Col);
        Master_XRF = Master_XRF./Master_mass(:,Vol_Col);
        Master_carbon = Master_carbon./Master_mass(:,Vol_Col);
    end
    % convert units of metals from ng to ug
    Master_XRF = Master_XRF./1000; Master_ICP = Master_ICP./1000;


    %% Data validation step 7: Compare mass to sum of measured species 
    % remove filter if sum of species greater than collected mass by more than 10%
    % Mass =  Metals + IC (without F, Cl, K, Mg, Ca) + BC 
    IC_col = find(contains(IC_title,{'IC_NO2_ug','IC_Br_ug','IC_NO3_ug',...
             'IC_PO4_ug','IC_SO4_ug','IC_Li_ug','IC_Na_ug','IC_NH4_ug'}));
     
    k = 0;
    
    for i = 1:length(Master_IDs)
        massi = Master_mass(i,1); 

        Master_sumspec = SumSpec(i,Master_XRF,Master_ICP,Master_IC,Master_carbon,IC_col);
        
        if Master_sumspec > massi
            fprintf('Sum of species (%.1f %s) for filter %s greater than collected mass (%.1f %s) \n', Master_sumspec,'ug/m3',Master_IDs{i},massi,'ug/m3')
        end
        
        if (Master_sumspec - massi) > 0.1*massi % test if difference between sum of species and collected mass is > 10% of collected mass
            k = k +1;
            fprintf('Filter mass for %s is invalid: exceeding 10%% allowance.\n', Master_IDs{i})
            filter_index(k) = i;
        end
    end
    
    if exist('filter_index','var') ~=0 % remove rows containing invalid mass data 
        [Master_IDs, Master_Barcodes,Master_CartridgeIDs, Master_LotIDs, Master_projectIDs,...
         Master_hours,  Master_masstype,Master_dates,        Master_mass,   Master_IC,  ...
         Master_ICP,    Master_XRF,     Master_carbon,       Master_Method, Master_flags] = ...
         RemoveRows (filter_index, ...
         Master_IDs, Master_Barcodes,Master_CartridgeIDs, Master_LotIDs, Master_projectIDs,...
         Master_hours,  Master_masstype,Master_dates,        Master_mass,   Master_IC,  ...
         Master_ICP,    Master_XRF,     Master_carbon,       Master_Method, Master_flags);
    end

    clear Master_sumspec massi i k filter_index 
    
    %% Formating: Public data, change the unit of metals back to ng 

    SNAcol  =  [find(contains(IC_title,'IC_SO4_ug')==1);
                find(contains(IC_title,'IC_NO3_ug')==1);
                find(contains(IC_title,'IC_NH4_ug')==1)];

    otherICcol=[find(contains(IC_title,'IC_Na_ug')==1);
                find(contains(IC_title,'IC_PO4_ug')==1);
                find(contains(IC_title,'IC_NO2')==1);
                find(contains(IC_title,'IC_Br')==1);
                find(contains(IC_title,'IC_K_ug')==1);
                find(contains(IC_title,'IC_Mg')==1);
                find(contains(IC_title,'IC_Ca')==1)];

    % Rearrange data to write to public (exclude IC_F, IC_Cl, IC_Li)
    Public_data = [ Master_mass(:,1), Master_IC(:,SNAcol), Master_carbon(:,2), Master_IC(:,otherICcol),...
                    Master_ICP, Master_XRF, Master_carbon(:,4), Master_Method];
    Public_data(:,[1:4 6:10 12]) = round(Public_data(:,[1:4 6:10 12]),2); % mass and water-soluble ions minus Mg round to 2 digits after decimal
    Public_data(:,5) = round(Public_data(:,5),2); % BC round to 1 digit after decimal
    Public_data(:,11) = round(Public_data(:,11)*1000,1); % convert Mg to nanograms/m3 and round to one digit after the decimal
    Public_data(:,13:end-2) = round(Public_data(:,13:end-2)*1000,2); % convert trace elements back to nanograms/m3 and round to one digit after the decimal (at ~ line 219 it was converted from ng to ug)
    
    % Titles for the public data:
    Title_public = {'PM25','IC_SO4_ug','IC_NO3_ug','IC_NH4_ug','BC_ug','IC_Na_ug','IC_PO4_ug','IC_NO2_ug','IC_Br_ug','IC_K_ug','IC_Mg_ng','IC_Ca_ug'};
    % column 13 to later will change if there is change in columns in master file
    Title_public(end+1:end+length(ICP_title)) = ICP_title;
    Title_public(end+1:end+length(XRF_title)) = XRF_title;
    Title_public(end+1:end+2) = {'OC_FTIR_ug','method_index'}; % 'BC_SSR_ug','EC_FTIR_ug' are removed from public data
    
    
    %% Formating: Sort into PM2.5, PMcoarse, PM10, and PMc 
    
    % index for mass types
    PM25_index = find(Master_masstype(:,1) == 1);
    blank_index = find(Master_masstype(:,1) == 0);
    PM10_index = find(Master_masstype(:,1) == 2);
    PMc_index = find(Master_masstype(:,1) == 3);
    
    if isempty(PM25_index) % No PM2.5 data then no need to go further 
        fprintf('No valid PM25 data for %s\n',Site_codes{loc})
        continue
    end

    % PM2.5
    [PM25_labels, PM25_Barcodes, PM25_cartridgeIDs, PM25_LotIDs, PM25_projectIDs,...
     PM25_hours,  PM25_masstype, PM25_dates,        PM25_Mass,   PM25_IC,   ...
     PM25_ICP,    PM25_XRF,      PM25_Carbon,       PM25_Method, PM25_flags]...
     = SortData (PM25_index,...
    Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,  Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method,Master_flags);
    
    PM25_data_public = Public_data(PM25_index,:);

    % PM10
    if ~isempty(PM10_index) 
        [PM10_labels, PM10_Barcodes, PM10_cartridgeIDs, PM10_LotIDs, PM10_projectIDs,...
         PM10_hours,  PM10_masstype, PM10_dates,        PM10_Mass,   PM10_IC,   ...
         PM10_ICP,    PM10_XRF,      PM10_Carbon,       PM10_Method, PM10_flags]...
        = SortData (PM10_index,...
        Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
        Master_hours, Master_masstype,   Master_dates,        Master_mass,  Master_IC,    ...
        Master_ICP,   Master_XRF,        Master_carbon,       Master_Method,Master_flags);
        
        PM10_data_public = Public_data(PM10_index,:);
    end
    
    % PMc
    if ~isempty(PMc_index) 
        [PMc_labels, PMc_Barcodes, PMc_cartridgeIDs, PMc_LotIDs, PMc_projectIDs,...
         PMc_hours,  PMc_masstype, PMc_dates,        PMc_Mass,   PMc_IC,   ...
         PMc_ICP,    PMc_XRF,      PMc_Carbon,       PMc_Method, PMc_flags]...
        = SortData (PMc_index,...
        Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
        Master_hours, Master_masstype,   Master_dates,        Master_mass,  Master_IC,    ...
        Master_ICP,   Master_XRF,        Master_carbon,       Master_Method,Master_flags);
        
        PMc_data_public = [ PMc_Mass(:,1), PMc_IC(:,SNAcol), PMc_IC(:,otherICcol), PMc_ICP ];
    end
    
    clear Public_data
    
    %% Write to public file 1: XXXX_speciation.csv
    % components contained in PM25_data_public matrix for searching the parameters sheet
    spec_components = {'Filter PM2.5 mass'...
        'Sulfate' 'Nitrate' 'Ammonium' 'Equivalent BC' 'Sodium' ... % SNA + HIPS_BC + Na 2-6
        'Phosphate' 'Nitrite' 'Bromide' 'Potassium','Magnesium' 'Calcium'... % IC 7-12
        'Lithium' 'Magnesium' 'Aluminum' 'Phosphorus' 'Titanium' ... % ICP 13-17
        'Vanadium' 'Chromium' 'Manganese' 'Iron' 'Cobalt' ... % ICP 18-22
        'Nickel' 'Copper' 'Zinc' 'Arsenic' 'Selenium' 'Gold' ... % ICP 23-28
        'Cadmium' 'Antimony' 'Barium' 'Cerium' 'Lead' ... % ICP 29-33
        'Sodium' 'Aluminum' 'Silicon' 'Sulfur' 'Chlorine' 'Potassium' ... % XRF 34-39
        'Calcium' 'Titanium' 'Vanadium' 'Iron' 'Zinc' 'Cerium' 'Lead' ... % XRF 40-46
        'Arsenic' 'Cobalt' 'Chromium' 'Copper' 'Magnesium' 'Manganese'... % XRF 47-52
        'Nickel' 'Antimony' 'Rubidium' 'Strontium' 'Cadmium' 'Selenium' 'Tin' ... % XRF 53-59
        'OC'}; % FTIR 60
    % PMc_data_public contains 4 more cols [2 carbons 1 method]

    % Titles for the xxxx_specieation.csv
    Titles = {'Site_Code' ... % col 1:3
        'Latitude' 'Longitude' 'Elevation_meters' ... % col 4:6
        'Filter_ID' 'Start_Year_local' 'Start_Month_local' 'Start_Day_local' 'Start_hour_local' 'End_Year_local' 'End_Month_local' 'End_Day_local' 'End_hour_local' 'Hours_sampled'... % col 7:16
        'Parameter_Code' 'Parameter_Name' 'Value' 'Units'... % col 17:20
        'Analytical_MDL' 'MDL' 'UNC' 'Method_Code' ... % col 21:24
        'Collection_Description' 'Analysis_Description' 'Conditions' 'Flag'};% col 25:28
    

        %% --------------- PM25 -------------- 
    Cart_ALL = unique(PM25_cartridgeIDs);

    for i = 1:length(Cart_ALL) % loop through each cartridge

        cart_idx = find(ismember(PM25_cartridgeIDs, Cart_ALL(i)) == 1);

        % will change for each component, loops over chemical components (k) by column
        for k = 1:length(spec_components)
            if i ==1 && k==1
                idx_spec = 1:length(cart_idx); % idx_spec is the row index for new data file
            else
                idx_spec = size(RFM_value,1) +1 : size(RFM_value,1) +length(cart_idx);
            end

            % speices specific information
            if k == 1 % PM2.5 mass
                parameters = PM25_mass_para; % first 5 are PM2.5 mass, therefore no need for further detailed indexing
                index = (1:5)';
                comp_idx = index(PM25_data_public(cart_idx,end));
                clear index
            elseif k == 5 % black carbon: HIPS or SSR
                parameters = PM25_carbon_para; 
                comp_idx = PM25_data_public(cart_idx,end); % range from 1 to 5, will point to SSR-based equivalent BC. This is updated to HIPS in the following lines
                for h = 1:length(comp_idx)
                    if PM25_data_public(cart_idx(h),end) == 5 && contains(PM25_flags{cart_idx(h)},'HIPS-BC estimated using SSR-BC') % HIPS for MAIA filters
                        comp_idx(h)= 9;
                    elseif PM25_data_public(cart_idx(h),end) == 5 % HIPS for MAIA filters
                         comp_idx(h)= 8;
                    elseif contains(PM25_flags{cart_idx(h)},'HIPS-BC estimated using SSR-BC') % HIPS for SPARTAN filters
                        comp_idx(h)= 7;
                    else
                        comp_idx(h)= 6;
                    end
                end
            elseif k == 60 % FTIR OC
                parameters = PM25_carbon_para; 
                comp_idx = zeros(size(cart_idx)); % will be assigned values of either 12 or 13 as a pointer to parameter info
                for h = 1:length(comp_idx)
                    if PM25_data_public(cart_idx(h),end) == 5 % MAIA filter
                        comp_idx(h)= 13;
                    else % SPARTAN filter
                        comp_idx(h)= 12;
                    end
                end

            elseif k >=13 && k <= 33 %trace elements from ICPMS are in columns 13 to 33;
                parameters = PM25_metals_para;
                parameter_names = table2array(parameters(:,1));
                parameter_analysis = table2array(parameters(:,6));
                indexA = find(contains(parameter_names, spec_components(k)));
                indexB = find(contains(parameter_analysis(indexA),'ICP-MS'));
                index = indexA(indexB);
                comp_idx = index(PM25_data_public(cart_idx,end));
                clear index indexA indexB parameter_analysis parameter_names
            elseif k > 34 % exclude Na_xrf
                parameters = PM25_metals_para;
                parameter_names = table2array(parameters(:,1));
                parameter_analysis = table2array(parameters(:,6));
                indexA = find(contains(parameter_names, spec_components(k)));
                indexB = find(contains(parameter_analysis(indexA),'XRF'));
                index = indexA(indexB);

                % find if the cartridge belongs to SPARTAN or MAIA
                tProject = unique(PM25_projectIDs(cart_idx));
                if length(tProject)>1
                    error('one cartridge assigned to %d projects.',length(tProject))
                else
                    if char(tProject) == 'S'
                        comp_idx = index(1);
                    elseif char(tProject) == 'M'
                        comp_idx = index(2);
                    end
                end

                
                % just a check:
                if size(PM25_Mass,1)~= size(PM25_data_public,1)
                    error('PM25_data_public length does not match PM25_Mass!')
                end
                % Read XRF MDL and Unc from XRF_MDL_UNC table
                [tAnaMDL,tMDL,tUNC ] = findMDLUnc (PM25_labels(cart_idx,:),spec_components(k),PM25_data_public(cart_idx,k),PM25_Mass(cart_idx,2),Ref_Dates,Ref_Values,FilterGroup);
                
                clear index indexA indexB parameter_analysis parameter_names

            elseif k == 34
                continue

            else % IC components
                parameters = PM25_IC_para;
                parameter_names = table2array(parameters(:,1));
                index = find(contains(parameter_names, spec_components(k)));
                comp_idx = index(PM25_data_public(cart_idx,end));
                clear index

                % Read IC MDL
                if exist('ic_mdl')==1
                    tfiilters = PM25_labels(cart_idx,:);
                    indi = find(ismember(ic_mdl_filterid,tfiilters));
                    indj = find(contains(ic_mdl_titles, spec_components(k)));
                    if ~isempty(indi)
                        tMDL = table2array(ic_mdl(indi,indj));
                        tAnaMDL = nan.*tMDL;
                        tUNC = nan.*tMDL;
                    end
                    clear indi indj
                end

            end

            % Common information for all species
            RFM_labels(idx_spec,1)  = PM25_labels(cart_idx,:); % Filter ID
            RFM_dates(idx_spec,1:8) = PM25_dates(cart_idx,1:8); % dates
            RFM_hour(idx_spec,1) =    PM25_hours(cart_idx); % hours sampled
            RFM_para_code(idx_spec,1) = table2array(parameters(comp_idx,2)); % parameter code
            RFM_para_name(idx_spec,1) = table2array(parameters(comp_idx,1)); % parameter name
            RFM_value(idx_spec,1) = PM25_data_public(cart_idx,k); % value (e.g PM2.5 mass, sulfate, trace element)
            RFM_unit(idx_spec,1) = table2array(parameters(comp_idx,7)); % units

            if exist('tMDL','var')
                RFM_method(idx_spec,1) = tAnaMDL; % Analytical MDL
                RFM_method(idx_spec,2) = tMDL; % MDL
                RFM_method(idx_spec,3) = tUNC; % UNC
            else
                RFM_method(idx_spec,1) = NaN; % Analytical MDL
                RFM_method(idx_spec,2) = NaN; % MDL
                RFM_method(idx_spec,3) = NaN; % UNC
            end

            RFM_method(idx_spec,4) = table2array(parameters(comp_idx,3)); % method code

            RFM_descrip(idx_spec,1) = table2array(parameters(comp_idx,5)); % collection description
            RFM_descrip(idx_spec,2) = table2array(parameters(comp_idx,6)); % analysis description

            % Adding Filter Flags (Ignore IC flags)
            RFM_flag(idx_spec,1) = copyflag(PM25_flags(cart_idx),Flag_parameters);


%             % print out RFM_labels, spec name and RFM_dates if there is a 0 in value
%             find_zeros_print(RFM_labels(idx_spec,1),RFM_dates(idx_spec,1:8),RFM_para_name(idx_spec,1),RFM_value(idx_spec,1), RFM_descrip(idx_spec,2) )

            clear  idx_spec parameters comp_idx tAnaMDL tMDL tUNC flag_vih
        end

        clear idx cart_idx k para_idx idx_spec cart_idxB
    end

    % do not allow nan for reported value AND do not allow nan for dates, this checks the start year
    nan_idx = find(  isnan(RFM_value )  | isnan(RFM_dates(:,1))   );
    temp =  zeros(size(RFM_unit));

   [RFM_labels, RFM_para_code, RFM_para_name, RFM_hour,     RFM_dates,  ...
    RFM_value,  RFM_unit,      RFM_method,    RFM_descrip,  RFM_flag, ~,~,~,~,~] ...
    =  RemoveRows (nan_idx, ...
    RFM_labels, RFM_para_code, RFM_para_name, RFM_hour,     RFM_dates,  ...
    RFM_value,  RFM_unit,      RFM_method,    RFM_descrip,  RFM_flag, temp,temp,temp,temp,temp); clear temp

    % ---- write table for sharing on website ----
    fileID = fopen(sprintf('%s/%s/%s_PM25_speciation.csv',direc_output,'Chemical_Filter_Data/PM25',Site_codes{loc}),'w');

    fprintf(fileID,'File Updated: %s \n%s \nSite:" %s, %s" \n',datestr(today),DataVersion,    Site_cities{loc},Site_countries{loc});
    fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',Titles{1:end});
    for h = 1:size(RFM_value,1)
        fprintf(fileID,'%s,   %g,%g,%g, %s, %g,%g,%g,%g,%g,%g,%g,%g,%g,%g, %s, %g, %s, %g,%g,%g,%g,  %s,%s,%s,%s\n', ...
            Site_codes{loc},          ... % col 1:3
            latitudes(loc),          longitudes(loc),        elevations(loc),...   % col 4:6
            RFM_labels{h,1},         RFM_dates(h,1:8),       RFM_hour(h,1), ...    % col 7:16
            RFM_para_code(h,1),      RFM_para_name{h,1},     RFM_value(h,1),      RFM_unit{h,1},... % col 17:20
            RFM_method(h,1:4),       RFM_descrip{h,1},       RFM_descrip{h,2}, ... % col 21:26
            'Ambient local',         RFM_flag{h,1} ); % col 25:28
    end
    fclose(fileID);

    clear nan_idx RFM_* Cart_ALL

        %% --------------- PM10 -------------- 

    if exist('PM10_data_public','var')

        spec_components{1} = 'PM10 mass';

        Cart_ALL = unique(PM10_cartridgeIDs);

        for i = 1:length(Cart_ALL) % loop through each cartridge

            cart_idx = find(ismember(PM10_cartridgeIDs, Cart_ALL(i)) == 1); % index for values for cartridge data in PM10_data_public

            % will change for each component, loops over chemical components (k) by column
            for k = 1:length(spec_components)
                if i ==1 && k==1
                    idx_spec = 1:length(cart_idx); % idx_spec is the row index for new data file
                else
                    idx_spec = size(RFM_value,1)+1 : size(RFM_value,1)+length(cart_idx);
                end

                if k == 1 % PM10 mass
                    parameters = PM10_mass_para;
                    index = (1:3)';
                    comp_idx = index(PM10_data_public(cart_idx,end));
                    clear index
                elseif k == 5 % SSR black carbon
                    parameters = PM10_carbon_para; % first 5 are equivalent BC (SSR-based), therefore no need for further detailed indexing
                    comp_idx = PM10_data_public(cart_idx,end); % range from 1 to 5, will point to SSR-based equivalent BC. This is updated to HIPS in the following lines
                    for h = 1:length(comp_idx)
                        if PM10_data_public(cart_idx(h),end) == 3 && contains(PM25_flags{cart_idx(h)},'HIPS-BC estimated using SSR-BC') % HIPS for MAIA filters
                            comp_idx(h)= 7;
                        elseif PM10_data_public(cart_idx(h),end) == 3 % HIPS for MAIA filters
                             comp_idx(h)= 6;
                        elseif contains(PM25_flags{cart_idx(h)},'HIPS-BC estimated using SSR-BC') % HIPS for SPARTAN filters
                            comp_idx(h)= 5;
                        else
                            comp_idx(h)= 4;
                        end
                    end

                elseif k == 60 % FTIR OC
                    parameters = PM10_carbon_para;
                    comp_idx = zeros(size(cart_idx)); % will be assigned values of either 12 or 13 as a pointer to parameter info
                    for h = 1:length(comp_idx)
                        if PM25_data_public(cart_idx(h),end) == 3 % MAIA filter
                            comp_idx(h)= 10;
                        else % SPARTAN filter
                            comp_idx(h)= 11;
                        end
                    end

                elseif k >=13 && k <= 33 %trace elements from ICPMS are in columns 13 to 33
                    parameters = PM10_metals_para;
                    parameter_names = table2array(parameters(:,1));
                    parameter_analysis = table2array(parameters(:,6));
                    indexA = find(contains(parameter_names, spec_components(k)));
                    indexB = find(contains(parameter_analysis(indexA),'ICP-MS'));
                    index = indexA(indexB);
                    comp_idx = index(PM10_data_public(cart_idx,end));
                    clear index indexA indexB parameter_analysis parameter_names
                elseif k > 34 % XRF
                    parameters = PM10_metals_para;
                    parameter_names = table2array(parameters(:,1));
                    parameter_analysis = table2array(parameters(:,6));
                    indexA = find(contains(parameter_names, spec_components(k)));
                    indexB = find(contains(parameter_analysis(indexA),'XRF'));
                    index = indexA(indexB);
                    if char(PM10_projectIDs(i)) == 'S'
                        comp_idx = index(1);
                    elseif char(PM10_projectIDs(i)) == 'M'
                        comp_idx = index(2);
                    end
                    % Read XRF MDL and Unc from XRF_MDL_UNC table
                        % just a check:
                        if size(PM10_Mass,1)~= size(PM10_data_public,1)
                            error('PM25_data_public length does not match PM25_Mass!')
                        end
                    [tAnaMDL,tMDL,tUNC ] = findMDLUnc (PM10_labels(cart_idx,:),spec_components(k),PM10_data_public(cart_idx,k),PM10_Mass(cart_idx,2), Ref_Dates,Ref_Values,FilterGroup);

                    clear index indexA indexB parameter_analysis parameter_names

                elseif k == 34
                    continue
                else % IC components
                    parameters = PM10_IC_para;
                    parameter_names = table2array(parameters(:,1));
                    index = find(contains(parameter_names, spec_components(k)));
                    comp_idx = index(PM10_data_public(cart_idx,end));
                    clear index

                    % Read IC MDL
                    if exist('ic_mdl')==1
                        tfiilters = PM10_labels(cart_idx,:);
                        indi = find(ismember(ic_mdl_filterid,tfiilters));
                        indj = find(contains(ic_mdl_titles, spec_components(k)));
                        if ~isempty(indi)
                            tMDL = table2array(ic_mdl(indi,indj));
                            tAnaMDL = nan.*tMDL;
                            tUNC = nan.*tMDL;
                        end
                        clear indi indj
                    end

                end

                RFM_labels(idx_spec,1)  = PM10_labels(cart_idx,:); % Filter ID
                RFM_dates(idx_spec,1:8) = PM10_dates(cart_idx,1:8); % dates
                RFM_hour(idx_spec,1) =    PM10_hours(cart_idx); % hours sampled
                RFM_para_code(idx_spec,1) = table2array(parameters(comp_idx,2)); % parameter code
                RFM_para_name(idx_spec,1) = table2array(parameters(comp_idx,1)); % parameter name
                RFM_value(idx_spec,1) = PM10_data_public(cart_idx,k); % value (e.g PM2.5 mass, sulfate, trace element)
                RFM_unit(idx_spec,1) = table2array(parameters(comp_idx,7)); % units

                if exist('tAnaMDL','var')
                    RFM_method(idx_spec,1) = tAnaMDL; % Analytical MDL
                    RFM_method(idx_spec,2) = tMDL; % MDL
                    RFM_method(idx_spec,3) = tUNC; % UNC
                else
                    RFM_method(idx_spec,1) = NaN; % Analytical MDL
                    RFM_method(idx_spec,2) = NaN; % MDL
                    RFM_method(idx_spec,3) = NaN; % UNC
                end
                RFM_method(idx_spec,4) = table2array(parameters(comp_idx,3)); % method code

                RFM_descrip(idx_spec,1) = table2array(parameters(comp_idx,5)); % collection description
                RFM_descrip(idx_spec,2) = table2array(parameters(comp_idx,6)); % analysis description

                % Adding Filter Flags (Ignore IC flags)
                RFM_flag(idx_spec,1) = copyflag(PM25_flags(cart_idx),Flag_parameters);

                clear  idx_spec parameters comp_idx tAnaMDL tMDL tUNC flag_vih
            end

            clear idx cart_idx k para_idx idx_spec cart_idxB
        end

        % do not allow nan for reported value AND do not allow nan for dates, this checks the start year
        nan_idx = find(  isnan(RFM_value )  | isnan(RFM_dates(:,1))   );
        temp =  zeros(size(RFM_unit));

        [RFM_labels, RFM_para_code, RFM_para_name, RFM_hour,     RFM_dates,  ...
            RFM_value,  RFM_unit,      RFM_method,    RFM_descrip,  RFM_flag, ~,~,~,~,~] ...
            =  RemoveRows (nan_idx, ...
            RFM_labels, RFM_para_code, RFM_para_name, RFM_hour,     RFM_dates,  ...
            RFM_value,  RFM_unit,      RFM_method,    RFM_descrip,  RFM_flag, temp,temp,temp,temp,temp); clear temp

        fileID = fopen(sprintf('%s/%s/%s_PM10_speciation.csv',direc_output,'Chemical_Filter_Data/PM10',Site_codes{loc}),'w');

        fprintf(fileID,'File Updated: %s \n%s \nSite:" %s, %s" \n',datestr(today),DataVersion,    Site_cities{loc},Site_countries{loc});
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',Titles{1:end});
        for h = 1:size(RFM_value,1)
            fprintf(fileID,'%s,   %g,%g,%g, %s, %g,%g,%g,%g,%g,%g,%g,%g,%g,%g, %s, %g, %s, %g,%g,%g,%g,  %s,%s,%s,%s\n', ...
                Site_codes{loc},          ... % col 1:3
                latitudes(loc),          longitudes(loc),        elevations(loc),...   % col 4:6
                RFM_labels{h,1},         RFM_dates(h,1:8),       RFM_hour(h,1), ...    % col 7:16
                RFM_para_code(h,1),      RFM_para_name{h,1},     RFM_value(h,1),      RFM_unit{h,1},... % col 17:20
                RFM_method(h,1:4),       RFM_descrip{h,1},       RFM_descrip{h,2}, ... % col 21:26
                'Ambient local',         RFM_flag{h,1} ); % col 25:28
        end
        fclose(fileID);

        clear nan_idx RFM_* 
    end

        %% --------------- PMc --------------- 
    if exist('PMc_data_public','var')  && ~isempty(PMc_data_public) == 0

        spec_components = {'PMc mass'...
            'Sulfate' 'Nitrate' 'Ammonium' 'Sodium' 'Phosphate' 'Nitrite' 'Bromide' 'Potassium'...
            'Magnesium' 'Calcium' 'Lithium' 'Magnesium' 'Aluminum' 'Phosphorus' 'Titanium' ...
            'Vanadium' 'Chromium' 'Manganese' 'Iron' 'Cobalt' 'Nickel' 'Copper'...
            'Zinc' 'Arsenic' 'Selenium' 'Gold' 'Cadmium' 'Antimony' 'Barium' 'Cerium' 'Lead'};

        if length(spec_components)~=size(PMc_data_public)
            error('PMc public data size wrong')
        end

        Titles = {'Site_Code' 'Country' 'City' ...
            'Latitude' 'Longitude' 'Elevation_meters' ...
            'Filter_ID' 'Start_Year_local' 'Start_Month_local' 'Start_Day_local' 'Start_hour_local' 'End_Year_local' 'End_Month_local' 'End_Day_local' 'End_hour_local' 'Hours_sampled' ...
            'Parameter_Code' 'Parameter_Name' 'Value' 'Units' ...
            'Method_Code' ...
            'Collection_Description' 'Analysis_Description' 'Conditions'};

        Cart_ALL = unique(PMc_cartridgeIDs);

        for i = 1:length(Cart_ALL) % loop through each cartridge

            cart_idx = find(ismember(PMc_cartridgeIDs, Cart_ALL(i)) == 1); % index for values for cartridge data in PMc_data

            % will change for each component
            for k = 1:length(spec_components)
                if i ==1 && k==1
                    idx_spec = 1:length(cart_idx);
                else
                    idx_spec = size(RFM_value,1) +1:size(RFM_value,1) +length(cart_idx);
                end

                RFM_labels(idx_spec,1) = PMc_labels(cart_idx,:); % Filter ID
                RFM_dates(idx_spec,1:8) = PMc_dates(cart_idx,1:8); % dates
                RFM_hour(idx_spec,1) = PMc_hours(cart_idx); % hours sampled

                RFM_para_code(idx_spec,1) = table2array(PMc_parameters(k,2)); % parameter code
                RFM_para_name(idx_spec,1) = table2array(PMc_parameters(k,1)); % parameter name
                RFM_value(idx_spec,1) = PMc_data_public(cart_idx,k); % value (e.g PM2.5 mass, sulfate, trace element)
                RFM_unit(idx_spec,1) = table2array(PMc_parameters(k,7)); % units
                RFM_method(idx_spec,1) = table2array(PMc_parameters(k,3)); % method code
                RFM_descrip(idx_spec,1) = table2array(PMc_parameters(k,5)); % collection description
                RFM_descrip(idx_spec,2) = table2array(PMc_parameters(k,6)); % analysis description
                clear  idx_spec
            end
            clear idx cart_idx k para_idx idx_spec cart_idxB
        end

        %---- write table for sharing on website----
        nan_idx = find(  isnan(RFM_value )  | isnan(RFM_dates(:,1))   );
        temp =  zeros(size(RFM_unit));

        [RFM_labels, RFM_para_code, RFM_para_name, RFM_hour,     RFM_dates,  ...
            RFM_value,  RFM_unit,      RFM_method,    RFM_descrip,  ~, ~,~,~,~,~] ...
            =  RemoveRows (nan_idx, ...
            RFM_labels, RFM_para_code, RFM_para_name, RFM_hour,     RFM_dates,  ...
            RFM_value,  RFM_unit,      RFM_method,    RFM_descrip,  temp, temp,temp,temp,temp,temp); clear temp

        fileID = fopen(sprintf('%s/%s/%s_PMc_speciation.csv',direc_output,'Chemical_Filter_Data/PMc',Site_codes{loc}),'w');
        
        fprintf(fileID,'File Updated: %s \n%s \nSite:" %s, %s" \n',datestr(today),DataVersion,    Site_cities{loc},Site_countries{loc});
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',Titles{1:end});
        for h = 1:length(RFM_value)
            fprintf(fileID,'%s, %g,%g,%g, %s, %g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%s,%g,%s,%g,%s,%s,%s\n', ...
                Site_codes{loc},         ... % col 1:3
                latitudes(loc),          longitudes(loc),        elevations(loc),...   % col 4:6
                RFM_labels{h,1},         RFM_dates(h,1:8),       RFM_hour(h,1), ...    % col 7:16
                RFM_para_code(h,1),      RFM_para_name{h,1},     RFM_value(h,1),      RFM_unit{h,1},... % col 17:20
                RFM_method(h,1),         RFM_descrip{h,1},       RFM_descrip{h,2}, 'Ambient local' ); % col 21:24
        end
        fclose(fileID);
        clear nan_idx Titles RFM_*
    end

    %% Formating: Reconstruct PM2.5  
    
     missing_PM25 = find( isnan(PM25_data_public(:,1)) | PM25_data_public(:,1)<=0 );
    
     temp = NaN.*PM25_hours;
    [PM25_labels, PM25_Barcodes,PM25_cartridgeIDs, PM25_LotIDs, PM25_projectIDs,...
     PM25_hours,  PM25_dates,   PM25_flags,  PM25_data_public,  ...
     ~,  ~,    ~,     ~,   ~,    ~ ] = ...
     RemoveRows (missing_PM25, ...
     PM25_labels, PM25_Barcodes,PM25_cartridgeIDs, PM25_LotIDs, PM25_projectIDs,...
     PM25_hours,  PM25_dates,   PM25_flags,  PM25_data_public, ...
     temp, temp, temp, temp, temp, temp );

           % Species  = [ASO4 NH4NO3 RM   OC NaCl Soil BC  TEO  Na2SO4]
              kappa   = [0.61 0.61 0.10 0.10 1.5  0.01 0.01 0.01 0.68];   % GC: need to double check 
              Density = [1.70 1.70 1.30 0.10 2.2  2.2  1.8  2.5  2.67];   % GC 
              RH = 35;
              MGrowth = 1+kappa./Density.*RH/(100-RH); % mass growth factor
              MGrowth (1:2) = 1.1 ; % halve growth factor of SNA to be consistent with GC simulation. Ref:
              % http://wiki.seas.harvard.edu/geos-chem/index.php/Particulate_matter_in_GEOS-Chem

    % ==== Na2SO4 (sea-salt associated sulfate) ====
    NaSO4_dry=0.18.*PM25_data_public(:,6); % 0.18*[Na+]
    NaSO4_wet=NaSO4_dry*MGrowth(8);
    
     % === ANO3 ====
    ANO3_dry = nan(size(PM25_data_public,1),1);
    Ind = PM25_data_public(:,4)>0.29*PM25_data_public(:,3); % if NH4 is sufficient
    ANO3_dry(Ind) = 1.29*PM25_data_public(Ind,3); % [ANO3(dry)]=1.29*[NO3]
    % NH4 insufficient or one of them is NaN
    ANO3_dry(~Ind) = sum(PM25_data_public(~Ind,[3 4]),2,'omitnan');  % [ANO3(dry)]= [NH4] + [NO3]
    % mask out when both NH4 and NO3 are NaN
    ANO3_dry(sum(isnan(PM25_data_public(:,[3 4])),2)==2) = NaN;
    ANO3_wet=ANO3_dry*MGrowth(2);
    
    % ==== ASO4 ====
    SO4inNaSO4 = 0.12*PM25_data_public(:,6); % "0.12" factor from Henning et al 2003 JGR 108(D1) pp4030, doi:10.1029/2002JD002439
    SO4inNaSO4(isnan(SO4inNaSO4)) = 0;
    % SO4_NSS = [SO4] - 0.12*[Na+]
    SO4_NSS  = sum( [PM25_data_public(:,2) -SO4inNaSO4],2); clear SO4inNaSO4
    SO4_NSS(SO4_NSS<0)=0; % SO4_NSS should never be negative
    %[ASO4(dry)] = [SO42-] + [NH4+]-(18/62*[NO3-])
    NH4_remain = sum([PM25_data_public(:,4) -0.29*PM25_data_public(:,3)],2,'omitnan'); % remaining NH4 after bonded with NIT; if no NO3 (NaN value means did not pass QA/QC), will give NaN
    NH4_remain(NH4_remain<0) = 0; % remaining NH4 should never be negative. In a NIT exceeding event, HNO3 exist.
    ASO4_dry = SO4_NSS+NH4_remain; % low NIT condition, SO4 exists as (NH4)2SO4; high NIT condition, it can be NH4HSO4 or even H2SO4.
    ASO4_wet = ASO4_dry*MGrowth(1); 
    clear NH4_remain
    
    % ==== Sea Salt ====
    % Sea salt(dry) = (2.54*[Na+] - 0.1[Al]) or (1.65[Cl] - 0.1[Al])
    for i = 1:size(PM25_data_public,1)
        if ~isnan(PM25_data_public(i,contains(Title_public,'Al_XRF')))  % check for presence of XRF data first, column 35 is aluminum
            SSalt_dry(i,1) = 2.54*PM25_data_public(i,6)-0.1*(PM25_data_public(i,contains(Title_public,'Al_XRF'))/1000);
        elseif ~isnan(PM25_data_public(i,contains(Title_public,'Al_ICP'))) % check for presence of Al from ICP-MS
            SSalt_dry(i,1) = 2.54*PM25_data_public(i,6)-0.1*(PM25_data_public(i,contains(Title_public,'Al_ICP'))/1000);
        else % no Al from ICP or XRF, therefore do not correct with Al
            SSalt_dry(i,1) = 2.54*PM25_data_public(i,6);
        end
    end
    clear i
    SSalt_dry(SSalt_dry<-0.3)=NaN;  % nan out if too negative

    SSalt_wet = SSalt_dry*MGrowth(4);
     
    % ==== Dust & TEO ====
    % if XRF data available, define soil (aka mineral dust) by combination of {Al, Si, Ca, Ti, Fe} following Xuan Liu's eqn:
    % [SOIL] = [1.89Al×(1+MAL)+2.14Si+1.40Ca+1.36Fe+1.67Ti]×CF

    % if no XRF, only ICP-MS data available, soil is defined as:
    % [SOIL] = 10*([Al] + [Fe] + [Mg]), if no Fe use: [SOIL] = 10*([Mg] + 30*[Ti] + [Al])
    
    % soil is assumed to be hydrophobic, therefore no water uptake
    xrf_soil_factors = [1.89 2.14 1.40 1.36 1.67];
    Soil = nan(size(PM25_data_public,1),1);

    % Regional varying MAL and CF values:
    MAL_CF={'AEAZ'	'AUMN'	'ARCB'	'BDDU'	'BIBU'	'CADO'	'CAHA'	'CAKE'	'CALE'	'CASH'	'CHTS'	'CLST'	'CODC'	'ETAD'	'IDBD'	'ILHA'	'ILNZ'	'INDH'	'INKA'	'KRSE'	'KRUL'	'MXMC'	'NGIL'	'PHMO'	'PRFJ'	'SGSU'	'TWKA'	'TWTA'	'USBA'	'USBO'	'USMC'	'USNO'	'USPA'	'VNHN'	'ZAJB'	'ZAPR';
            0.72 	0.24 	0.62 	0.62 	0.62 	0.62 	0.62 	0.62 	0.62 	0.62 	0.59 	0.62 	0.27 	0.62 	0.62 	0.72 	0.72 	0.62 	0.62 	0.59 	0.59 	0.27 	0.27 	0.62 	0.27 	0.62 	0.62 	0.62 	0.27 	0.27 	0.27 	0.27	0.66 	0.62 	0.62 	0.62 ;
            1.14 	1.05 	1.02 	1.02 	1.02 	1.02 	1.02 	1.02 	1.02 	1.02 	1.11 	1.02 	1.05 	1.02 	1.02 	1.14 	1.14 	1.02 	1.02 	1.11 	1.11 	1.05 	1.05 	1.02 	1.05 	1.02 	1.02 	1.02 	1.05 	1.05 	1.05 	1.05 	1.14 	1.02 	1.02 	1.02 };
    
    Ind = find(ismember(MAL_CF(1,:),Site_codes{loc}));
    if ~isempty(Ind)
        MAL = cell2mat(MAL_CF(2,Ind)); 
        CF = cell2mat(MAL_CF(3,Ind));
    else
        error('MAL and CF not found for %s',Site_codes{loc})
    end

    % TEO ratios if ICP-MS = (1.47[V] + 1.27[Ni] + 1.25[Cu] + 1.24[Zn] + 1.32[As] + 1.2[Se] + 1.07[Ag] + 1.14[Cd] + 1.2[Sb] + 1.12[Ba] + 1.23[Ce] + 1.08[Pb]);
    TEO_ratio_ICP = [1.47 1.27 1.25 1.24 1.32 1.2 1.07 1.14 1.2 1.12 1.23 1.08];
    TEO_title_ICP =[find(contains(Title_public,'V_ICP')==1),find(contains(Title_public,'Ni_ICP')==1),find(contains(Title_public,'Cu_ICP')==1),...
                    find(contains(Title_public,'Zn_ICP')==1),find(contains(Title_public,'As_ICP')==1),find(contains(Title_public,'Se_ICP')==1),...
                    find(contains(Title_public,'Ag_ICP')==1),find(contains(Title_public,'Cd_ICP')==1),find(contains(Title_public,'Sb_ICP')==1),...
                    find(contains(Title_public,'Ba_ICP')==1),find(contains(Title_public,'Ce_ICP')==1),find(contains(Title_public,'Pb_ICP')==1)];
    % TEO ratios if XRF = (1.79[V] + 1.69[Cr] + 1.63[Mn] + 1.34[Co] + 1.27[Ni] + 1.25[Cu] + 1.24[Zn] + 1.43[As] + 1.41[Se] + 1.09[Rb] + 1.18[Sr] + 1.14[Cd] + 1.20[Sn] + 1.26[Sb] + 1.20[Ce] + 1.12[Pb]);
    TEO_ratio_XRF = [1.79 1.69 1.63 1.34 1.27 1.25 1.24 1.43 1.41 1.09 1.18 1.14 1.20 1.26 1.20 1.12];
    TEO_title_XRF =[find(contains(Title_public,'V_XRF')==1),find(contains(Title_public,'Cr_XRF')==1),find(contains(Title_public,'Mn_XRF')==1),...
                    find(contains(Title_public,'Co_XRF')==1),find(contains(Title_public,'Ni_XRF')==1),find(contains(Title_public,'Cu_XRF')==1),...
                    find(contains(Title_public,'Zn_XRF')==1),find(contains(Title_public,'As_XRF')==1),find(contains(Title_public,'Se_XRF')==1),...
                    find(contains(Title_public,'Rb_XRF')==1),find(contains(Title_public,'Sr_XRF')==1),find(contains(Title_public,'Cd_XRF')==1),...
                    find(contains(Title_public,'Sn_XRF')==1),find(contains(Title_public,'Sb_XRF')==1),find(contains(Title_public,'Ce_XRF')==1),...
                    find(contains(Title_public,'Pb_XRF')==1)];
    TEO = nan(size(PM25_data_public,1),1);
    
    
    % Start to calculate soil and TEO
    no_metals = 0;
    method_metals = ones(size(PM25_data_public,1),1); % default is ICP-MS (ICP-MS = 1, XRF = 2)
    
    for i = 1:size(PM25_data_public,1)

        if ~isnan(PM25_data_public(i,34)) % there is XRF data
            tSoil = xrf_soil_factors.*PM25_data_public(i,[35 36 40 41 43]);
            tSoil(1) = tSoil(1).*(1+MAL); 
            tSoil = tSoil.*CF;
            Soil(i,1) = (sum(tSoil,2,'omitnan'))/1000; clear tSoil % convert to ug/m3 from ng/m3
            
            TEO(i,1) = (sum(TEO_ratio_XRF.*PM25_data_public(i,TEO_title_XRF),2,'omitnan'))/1000;
            
            method_metals(i) = 2;

        elseif ~isnan(PM25_data_public(i,15)) % no XRF, but there is ICP-MS data
            if PM25_data_public(i,21) > 0 % Fe is available in ICP-MS data
                Soil(i,1) =(10*(PM25_data_public(i,14) + PM25_data_public(i,21) + PM25_data_public(i,15)))/1000; % convert to ug/m3 from ng/m3
            elseif PM25_data_public(i,21) == 0 || isnan(PM25_data_public(i,21)) % no Fe in ICP-MS data
                Soil(i,1) = (10*(PM25_data_public(i,14) +  30*PM25_data_public(i,17) + PM25_data_public(i,15)))/1000 ;% convert to ug/m3 from ng/m3
            end
            
            TEO(i,1)=(sum(TEO_ratio_ICP.*PM25_data_public(i,TEO_title_ICP),2,'omitnan'))/1000; % divide by 1000 to convert to ug/m3 from ng/m3
        
        else % no ICP-MS or XRF data
            no_metals = no_metals+1; % counts number of filters that do not have ICP-MS of XRF data available
            Soil(i,1) = nan;
            TEO(i,1)= nan;
        end
    end
    
    if no_metals > 0
        fprintf('There are %d filters out of %d that do not have metal data \n',no_metals, size(PM25_data_public,1))
    end
    clear no_metals
    % ==== BC ====
%     BC = PM25_data_public(:,5); % non-hygroscopic
    BC = PM25_data_public(:,5)*0.06./0.1; % assume sigma 0.06 is too low. scale to sigma = 0.10

    % ==== OC ====
    % Need to convert to OM?
    OC_dry = PM25_data_public(:,end-1); % FTIR
    OC_wet = OC_dry.*MGrowth(3);
   
    % ==== RESIDUAL MASS (Organic Matter) ====
    merged_data_wet = [PM25_data_public(:,1) ASO4_wet ANO3_wet SSalt_wet Soil BC OC_wet TEO NaSO4_wet];
    colmask = any(isnan(merged_data_wet(:,[1 2 3 5])),2); % if any of [ASO4_wet ANO3_wet Soil] is nan, do not calculate RM
    
    RM_wet=merged_data_wet(:,1)-nansum(merged_data_wet(:,2:end),2);
    RM_wet(colmask == 1) = NaN; % if missing input inorganic, do not calculate RM
    RM_dry=RM_wet./MGrowth(3); % Dry residual matter
    
    Too_neg_OM=find(RM_wet./merged_data_wet(:,1)<-0.1); % set to NaN if negative by more than 10% of total wet (35% RH) PM2.5 mass
    RM_dry(Too_neg_OM,1)=NaN;
    RM_wet(Too_neg_OM,1)=NaN;
    clear too_neg_OM mask colmask index_OM merged_data_wet
    
    % ==== PBW ====
    % PBW = particle-bound water (mass conc., ug/m3) [currently not reported in RCFM files]
    PBW = sum( [ ASO4_wet -ASO4_dry ANO3_wet -ANO3_dry RM_wet -RM_dry OC_wet -OC_dry SSalt_wet -SSalt_dry NaSO4_wet -NaSO4_dry ] ,2,'omitnan');
    
    % ==== Kappa mix ====
    % dry volume
    Vol_mix = [ASO4_dry  ANO3_dry  RM_dry  OC_dry  SSalt_dry  Soil  BC  TEO  NaSO4_dry]./Density;
    
    kappa(4) = kappa(4)/2; % reducing kappa for NaCl by half to retain PBW at 0% RH
    
    for i=1:size(Vol_mix,1)
        kappa_mix_spec(i,:)=Vol_mix(i,:).*kappa./sum(Vol_mix(i,:),2,'omitnan');
    end
    kappa_mix=sum(kappa_mix_spec,2,'omitnan');
    kappa_mix(kappa_mix==0)=NaN;
    kappa_mix(kappa_mix>0.9)=NaN;
    kappa_mix(isnan(kappa_mix_spec(:,2))) = NaN; % if missing SO4 we are missing that and something else, therefore set to NaN
    
    % ------------ Merge RCFM --------------------------------------------
    merged_data_wet = [PM25_data_public(:,1) ASO4_wet ANO3_wet SSalt_wet Soil BC TEO OC_wet RM_wet kappa_mix NaSO4_wet PM25_data_public(:,end) method_metals]; 
    % columns are:
    % 1 = PM2.5
    % ... see variable names
    % 11 = method index
    % 12 = method_index for xrf vs. icp
    
    %% Write to PM2.5_only file with kappa
    % ----- File to use with nephelometer data to estimate time-resolved PM2.5 -----
    % this only outputs total PM2.5 mass concentration, sampling dates, and kappa
    Titles={'Start_Date' 'End_Date' 'PM2.5_mass_ug/m3' 'kappa_mix'};
    
    fileID = fopen(sprintf('%s/Public_Data/Time-resolved_PM2.5/Input_data/PM25only_%s.csv',direc,Site_codes{loc}),'w');
    fprintf(fileID,'%s,%s,%s,%s\n',Titles{1,1:length(Titles)});
    
    for i = 1:numel(PM25_dates(:,1))
        start_date = datenum(PM25_dates(i,1),PM25_dates(i,2), PM25_dates(i,3), PM25_dates(i,4),0,0);
        end_date = datenum(PM25_dates(i,5),PM25_dates(i,6), PM25_dates(i,7), PM25_dates(i,8),0,0);
        fprintf(fileID,'%.10f,%.10f,%6.2f,%6.2f \n',start_date, end_date, PM25_data_public(i,1),kappa_mix(i));
    end
    
    fclose(fileID);
    
    %% Write to public file 3: RCFM files 
    merged_data_wet(:,7) = merged_data_wet(:,7)*1000; % converting TEO from ug/m3 to ng/m3
    merged_data_wet(:,1:10) = round(merged_data_wet(:,1:10),2);
    
    rcfm_components = {'Filter PM2.5 mass' 'Ammoniated Sulfate' 'Ammonium Nitrate' 'Sea Salt' 'Fine Soil' 'Equivalent BC' 'Trace Element Oxides'...
        'Organic Carbon' 'Residual Matter' 'kappa'};
    
    Titles = {'Site_Code' 'Latitude' 'Longitude' 'Elevation_meters' 'Filter_ID'...
        'Start_Year_local' 'Start_Month_local' 'Start_Day_local' 'Start_hour_local' 'End_Year_local' 'End_Month_local' 'End_Day_local' 'End_hour_local'...
        'Hours_sampled' 'Parameter_Code' 'Parameter_Name' 'Value' 'Units' 'Method_Code' 'Collection_Description'...
        'Analysis_Description' 'Conditions' 'Flag'};
    
    Cart_ALL = unique(PM25_cartridgeIDs);
    for i = 1:length(Cart_ALL) % loop through each cartridge
        
        cart_idx = find(ismember(PM25_cartridgeIDs, Cart_ALL(i)) == 1); % index for values for cartridge data in merged_data_wet
        
        % will change for each component
        for k = 1:size(merged_data_wet,2)-3 % -3 bc last 3 columns are Na2SO4, PM2.5/BC method index, and metals method index, respectively
            if i ==1 && k==1
                idx_spec = 1:length(cart_idx);
            else
                idx_spec = size(RFM_value,1) +1 : size(RFM_value,1)+length(cart_idx);
            end
            
            if k == 1 % PM2.5 mass
                parameters = PM25_mass_para; % first 5 are PM2.5 mass, therefore no need for further detailed indexing
                index = (1:5)';
                comp_idx = index(PM25_data_public(cart_idx,end));
                clear index
            elseif k == 6 % black carbon: HIPS or SSR
                parameters = PM25_carbon_para; 
                comp_idx = PM25_data_public(cart_idx,end); % range from 1 to 5, will point to SSR-based equivalent BC. This is updated to HIPS in the following lines
                for h = 1:length(comp_idx)
                    if PM25_data_public(cart_idx(h),end) == 5 && contains(PM25_flags{cart_idx(h)},'HIPS-BC estimated using SSR-BC') % HIPS for MAIA filters
                        comp_idx(h)= 9;
                    elseif PM25_data_public(cart_idx(h),end) == 5 % HIPS for MAIA filters
                         comp_idx(h)= 8;
                    elseif contains(PM25_flags{cart_idx(h)},'HIPS-BC estimated using SSR-BC') % HIPS for SPARTAN filters
                        comp_idx(h)= 7;
                    else
                        comp_idx(h)= 6;
                    end
                end

            elseif k == 8 % FTIR OC
                parameters = PM25_carbon_para; 
                comp_idx = zeros(size(cart_idx)); % will be assigned values of either 12 or 13 as a pointer to parameter info
                for h = 1:length(comp_idx)
                    if PM25_data_public(cart_idx(h),end) == 5 % MAIA filter
                        comp_idx(h)= 13;
                    else % SPARTAN filter
                        comp_idx(h)= 12;
                    end
                end
            elseif k == 2 || k == 3 || k == 9 || k ==10 % 2 = ASO4, 3 = NH4NO3, 9 = RM, 10 = kappa
                parameters = RCFM_parameters;
                parameter_names = table2array(parameters(:,1));
                index = find(contains(parameter_names, rcfm_components(k)));
                comp_idx = index;
                clear index parameter_names
            else % remaining RCFM parameters that need distinguishing between ICP-MS and XRF
                parameters = RCFM_parameters;
                parameter_names = table2array(parameters(:,1));
                index = find(contains(parameter_names, rcfm_components(k)));
                comp_idx = index(merged_data_wet(cart_idx,end));
                clear index parameter_names
            end
            
            RFM_labels(idx_spec,1) = PM25_labels(cart_idx); % hours sampled
            RFM_dates(idx_spec,1:8)= PM25_dates(cart_idx,1:8); % dates
            RFM_hour(idx_spec,1) = PM25_hours(cart_idx); % hours sampled
            RFM_para_code(idx_spec,1) = table2array(parameters(comp_idx,2)); % parameter code
            RFM_para_name(idx_spec,1) = table2array(parameters(comp_idx,1)); % parameter name
            RFM_value(idx_spec,1) = merged_data_wet(cart_idx,k); % value (e.g PM2.5 mass, sulfate, trace element)
            RFM_unit(idx_spec,1)  = table2array(parameters(comp_idx,7)); % units
            RFM_method(idx_spec,1) = table2array(parameters(comp_idx,3)); % method code
            RFM_descrip(idx_spec,1)  = table2array(parameters(comp_idx,5)); % collection description
            RFM_descrip(idx_spec,2) = table2array(parameters(comp_idx,6)); % analysis description

%             % print out RFM_labels, spec name and RFM_dates if there is a 0 in value
%             find_zeros_print(RFM_labels(idx_spec,1),RFM_dates(idx_spec,1:8),RFM_para_name(idx_spec,1),RFM_value(idx_spec,1), RFM_descrip(idx_spec,2) )

            % Adding Filter Flags (Ignore IC flags) 
            RFM_flag(idx_spec,1) = copyflag(PM25_flags(cart_idx),Flag_parameters);

            clear  idx_spec parameters comp_idx flag_vih flag_fs
        end
        clear idx cart_idx k idx_spec cart_idxB
    end
    
    nan_idx = find(  isnan(RFM_value )  | isnan(RFM_dates(:,1))   ); 
    temp =  zeros(size(RFM_unit));

    [RFM_labels, RFM_para_code, RFM_para_name, RFM_hour,     RFM_dates,  ...
     RFM_value,  RFM_unit,      RFM_method,    RFM_descrip,  RFM_flag, ~,~,~,~,~] ...
     =  RemoveRows (nan_idx, ...
     RFM_labels, RFM_para_code, RFM_para_name, RFM_hour,     RFM_dates,  ...
     RFM_value,  RFM_unit,      RFM_method,    RFM_descrip,  RFM_flag, temp,temp,temp,temp,temp); clear temp
 
    % write table for sharing on website
    fileID = fopen(sprintf('%s/%s/%s_PM25_RCFM.csv',direc_output,'RCFM',Site_codes{loc}),'w');
    fprintf(fileID,'File Updated: %s \n%s \nSite:" %s, %s" \n',datestr(today),DataVersion,    Site_cities{loc},Site_countries{loc});
    fprintf(fileID,'%s,%s,%s, %s,%s,%s,%s,%s, %s,%s,%s,%s,%s, %s,%s,%s,%s,%s, %s,%s,%s,%s,%s\n',Titles{1:end});
    for h = 1:size(RFM_value,1)
        fprintf(fileID,'%s,%g,%g,%g,%s,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%s,%g,%s,%g,%s,%s,%s,%s\n', ...
                Site_codes{loc},        ... % col 1:3
                latitudes(loc),          longitudes(loc),        elevations(loc),...   % col 4:6
                RFM_labels{h,1},         RFM_dates(h,1:8),       RFM_hour(h,1), ...    % col 7:16 
                RFM_para_code(h,1),      RFM_para_name{h,1},     RFM_value(h,1),      RFM_unit{h,1},... % col 17:20
                RFM_method(h,1),         RFM_descrip{h,1},       RFM_descrip{h,2}, ... % col 21:23
                'Ambient local',         RFM_flag{h,1} ); % col 24:25
    end
    fclose(fileID);

    %% Figure 1: Time series of cartridge average and site average PM2.5 for website
    PM25filter_timeseries(PM25_Mass(:,1), PM25_dates, Site_cities, Site_codes, direc_output,loc)

    species_plot = [PBW PM25_data_public(:,[2 4 3 6])  Soil TEO BC OC_dry RM_dry]; % order: {'SO4','NH4','NO3','SS','Dust','TEO','BC','OC','RM'};
    
    % calculate reportable data for later use
    tSiteDataNum(loc,1) = sum(~isnan(PM25_data_public(:,1)));
    tSiteDataNum(loc,2:end) = sum(~isnan(species_plot));
       
    %% Figure 2: Bar chart - all data with dates or filter IDs as labels  
    dates = datenum(PM25_dates(:,5), PM25_dates(:,6), PM25_dates(:,7),0,0,0);
    x_dates = char(datetime(dates,'ConvertFrom','datenum', 'Format','MM-dd-yyyy'));
    Spec_barplot(species_plot,x_dates,PM25_labels,Site_cities,Site_codes,loc,direc_output);

    
    %% Figure 3: Bar chart - most recent and complete data for website   
    Spec_barplot_website(PM25_dates, Site_codes, direc_output,loc, species_plot)
     

    %% Figure 4: RCFM Pie chart
    PM25_total = mean(PM25_data_public(:,1),'omitnan');
    % if there is a missing species (all NaNs in one species), no pie will be made.
    Pie_made = PM25_RCFMavg_pie(PM25_total, species_plot, Site_cities,loc); % function for making pie chart
    if Pie_made == 1
        % save figure
        fname = sprintf('%s/Public_Data/Chemical_Filter_Data/Plots/Pie_spec_plots/%s_PM25_RCFMavg',direc,Site_codes{loc});
        saveas(gcf,sprintf('%s.png',fname))
        print(sprintf('%s.eps',fname),'-depsc')
        close all

        % save data to a xlsx
        PM25 = PM25_data_public(:,1);
        SO4 = PM25_data_public(:,2);
        NH4 = PM25_data_public(:,4);
        NO3 = PM25_data_public(:,3);
        Na = PM25_data_public(:,6);
        Cl =  PM25_IC(:,contains(IC_title,'IC_Cl'));
        savedata_dry(PM25_labels, PM25, BC, TEO, Soil, Na, NO3, NH4, SO4, PBW, OC_dry, RM_dry, Cl,fname)
        clear PM25 SO4 NH4 NO3 Na Cl
    end

    %% Figure 5: RCFM Pie chart for 2019 and forward
    year = 2020;
    ind = find(PM25_dates(:,5) >= year);

    PM25_total = mean(PM25_data_public(ind,1),'omitnan');

    Cl =  mean(PM25_IC(ind,contains(IC_title,'IC_Cl')),'omitnan');
    K =  mean(PM25_IC(ind,contains(IC_title,'IC_K')),'omitnan');
    N = length(ind); % number of filters included

    Pie_made = PM25_RCFMavg_pie_with_Cl(PM25_total, species_plot(ind,:), Cl,K, Site_cities,loc); % function for making pie chart
    if Pie_made == 1
        % save figure
        savedir = sfmkdir(sprintf('%s/Public_Data/Chemical_Filter_Data/Plots/Pie_spec_plots/MTL_filters',direc));
        fname = sprintf('%s/%s_PM25_RCFMavg',savedir,Site_codes{loc});
        saveas(gcf,sprintf('%s.png',fname))
        print(sprintf('%s.eps',fname),'-depsc')
        close all

        % save data to a xlsx
        PM25 = PM25_data_public(ind,1);
        SO4 = PM25_data_public(ind,2);
        NH4 = PM25_data_public(ind,4);
        NO3 = PM25_data_public(ind,3);
        Na = PM25_data_public(ind,6);
        Cl =  PM25_IC(ind,contains(IC_title,'IC_Cl'));
        savedata_dry(PM25_labels(ind,:), PM25, BC(ind,:), TEO(ind,:), Soil(ind,:), Na, NO3, NH4, SO4, PBW(ind,:), OC_dry(ind,:), RM_dry(ind,:),Cl,fname)
        clear PM25 SO4 NH4 NO3 Na Cl
    end

       


    %% Additional Pie Charts: Cartridge specific pie charts
    if additional_pie == 1

    infname =  sprintf('%s/Public_Data/Chemical_Filter_Data/Plots/Pie_spec_plots/ByCartridge/Cartridges_need_pie_charts.xlsx',direc);
    cartlist = readtable(infname,'Sheet','PRFJ');
    cartlist = table2array(cartlist);
    ind = find(ismember(PM25_cartridgeIDs,cartlist));

    if ~isempty(ind)
        PM25_total = mean(PM25_data_public(ind,1),'omitnan');
        Cl =  mean(PM25_IC(ind,contains(IC_title,'IC_Cl')),'omitnan');
        K =  mean(PM25_IC(ind,contains(IC_title,'IC_K')),'omitnan');
        N = length(ind); % number of filters included

        Pie_made = PM25_RCFMavg_pie(PM25_total, species_plot(ind,:) , Site_cities,loc); % function for making pie chart
        if Pie_made == 1
            % add note of filter number
            notes = sprintf('Num of Filters = %d',N);
            text(-0.3, -0.05, notes,'units','normalized')

            fname = sprintf('%s/Public_Data/Chemical_Filter_Data/Plots/Pie_spec_plots/ByCartridge/%s_PM25_RCFMavg',direc,Site_codes{loc});
            saveas(gcf,sprintf('%s.png',fname))
            print(sprintf('%s.eps',fname),'-depsc')
            close all
        end

        % making another pie chart that includes Cl
        Pie_made = PM25_RCFMavg_pie_with_Cl(PM25_total, species_plot(ind,:) , Cl,K, Site_cities,loc); % function for making pie chart
        if Pie_made == 1
            % add note of filter number
            notes = sprintf('Num of Filters = %d',N);
            text(-0.1, -0.05, notes,'units','normalized')
            saveas(gcf,sprintf('%s_with_Cl.png',fname))
            print(sprintf('%s_with_Cl.eps',fname),'-depsc')
            close all


            % save data to a xlsx
            PM25 = PM25_data_public(ind,1);
            SO4 = PM25_data_public(ind,2);
            NH4 = PM25_data_public(ind,4);
            NO3 = PM25_data_public(ind,3);
            Na = PM25_data_public(ind,6);
            Cl =  PM25_IC(ind,contains(IC_title,'IC_Cl'));
            savedata_dry(PM25_labels(ind,:), PM25, BC(ind,:), TEO(ind,:), Soil(ind,:), Na, NO3, NH4, SO4, PBW(ind,:), OC_dry(ind,:), RM_dry(ind,:),Cl,fname)
            savedata_wet(PM25_labels(ind,:), PM25, BC(ind,:), TEO(ind,:), Soil(ind,:), ANO3_wet(ind,:), ASO4_wet(ind,:), NaSO4_wet(ind,:), OC_wet(ind,:), RM_wet(ind,:),SSalt_wet(ind,:),fname)
            clear  PM25 SO4 NO3 NH4 Na Cl

        end
    end

    end

    %%
    close all
    clear *_dry *_wet BC *_index  colmask nan_idx  Titles  RFM_* Master_* SO4_NSS Soil SSalt_water start_date end_date TEO  Vol_mix 
    clear fileID i IN_tot_wet index_OM kappa_mix kappa_mix_spec PBW Too_neg_OM missing_PM25 PM*_index  PM25_dates Cart_ALL
    clear merged_data_wet species  PM*_data_public Master_CartridgeIDs h method_metals ic_mdl ic_mdl_*
    clear PM*_labels PM*_Barcodes PM*_cartridgeIDs PM*_LotIDs PM*_projectIDs PM*_masstype PM*_dates PM*_Mass PM*_Carbon PM*_Method PM*_flags PM*_ICP PM*_XRF PM*_IC
    
    fprintf('Finished processing data for %s \n\n', Site_cities{loc})
    
 end

%% Print out data summary tables
% print out the BC availability table
BC_table = table(Site_codes, HIPS_Mesured_BC, HIPS_Est_BC, HIPS_total_BC, SSR_total_BC);
writetable(BC_table,sprintf('%s/Analysis_Data/Black_Carbon/BC_availability.xlsx',direc));

%  print out reportable data num 
T = table(Site_cities, tSiteDataNum(:,1),tSiteDataNum(:,2),tSiteDataNum(:,3),tSiteDataNum(:,4),tSiteDataNum(:,5),...
                       tSiteDataNum(:,6),tSiteDataNum(:,7),tSiteDataNum(:,8),tSiteDataNum(:,9),tSiteDataNum(:,10),'VariableNames',DataNum_tab_rows);
writetable(T,sprintf('%s/Analysis_Data/Reportable_Data_Count_New.xlsx',direc));


fprintf('END of Master File processing for %s \n', datestr(today))

diary off

%% FUNCTIONS
function savedata_wet(PM25_labels, PM25, BC, TEO, Soil, ANO3_wet, ASO4_wet, NaSO4_wet, OC_wet, RM_wet,SSalt_wet,fname)
        T = table(PM25_labels, PM25, BC, TEO, Soil, ANO3_wet, ASO4_wet, NaSO4_wet, OC_wet, RM_wet, SSalt_wet);
        tablename = sprintf('%s.xlsx',fname);
        writetable(T,tablename,'Sheet','wet')
end

function savedata_dry(PM25_labels, PM25, BC, TEO, Soil, Na, NO3, NH4, SO4, PBW, OC_dry, RM_dry,Cl,fname) 
        T = table(PM25_labels, PM25, BC, TEO, Soil, Na, NO3, NH4, SO4, PBW,Cl,OC_dry, RM_dry);
        tablename = sprintf('%s.xlsx',fname);
        delete(tablename) % delete existing file 
        writetable(T,tablename,'Sheet','dry')
end

function [tAnaMDL,tMDL,tUNC ]= findMDLUnc (tPM25_labels,ThisElement,Mass_conc, Vol, Ref_Dates,Ref_Values,FilterGroup)

    Titles = Ref_Values.Properties.VariableNames;
    Ref_Elements = table2array(Ref_Values(:,contains(Titles,'Element')));
%     Vol          = table2array(Ref_Values(1,contains(Titles,'Volume')));
    Unc_Pro_Vol  = 0.01*table2array(Ref_Values(1,contains(Titles,'Unc_proportion_volume')));

    for ii = 1:size(tPM25_labels,1)
        ThisFilter = tPM25_labels(ii,:);
        DateInd = zeros(size(FilterGroup,2)-1,1);

        for dd = 1:length(DateInd)
            DateInd(dd) = sum(contains(FilterGroup{2,dd},ThisFilter));
        end
        if sum(DateInd)>0
            DateInd(DateInd>1)=1; % raw data are duplicated in Archive
            MDL_Ind = find( contains(Ref_Elements,ThisElement) & Ref_Dates == FilterGroup{1,DateInd==1} );
            tAnaMDL(ii,1) = table2array(Ref_Values(MDL_Ind,contains(Titles,'Analytical_MDL'))) ./Vol(ii);
            tMDL(ii,1)    = table2array(Ref_Values(MDL_Ind,ismember(Titles,'MDL (ng)')))       ./Vol(ii);
            UncA          = table2array(Ref_Values(MDL_Ind,contains(Titles,'Unc_add')));
            UncP          = 0.01*table2array(Ref_Values(MDL_Ind,contains(Titles,'Unc_proportion_mass')));
            tUNC(ii,1) = sqrt( (UncA./Vol(ii))^2 + (UncP^2 + Unc_Pro_Vol^2)*Mass_conc(ii)^2 );

        else
            tAnaMDL(ii,1) = NaN;
            tMDL(ii,1) = NaN;
            tUNC(ii,1) = NaN;
        end
    end
    tAnaMDL = round(tAnaMDL,2);
    tMDL = round(tMDL,2);
    tUNC = round(tUNC,2);
end

function [PM25_labels,PM25_Barcodes,PM25_cartridgeIDs,PM25_LotIDs,PM25_projectIDs,...
        PM25_hours,PM25_masstype,PM25_dates,PM25_mass,PM25_IC,PM25_ICP,PM25_XRF,...
        PM25_carbon,PM25_Method,PM25_flags] = ...
        SortData (PM25_index, Master_IDs, Master_Barcodes,CartridgeIDs_master, LotIDs, ...
        projectIDs_master,Master_hours, Master_masstype, Master_dates,  ...
        Master_mass,  Master_IC,  Master_ICP,  Master_XRF, Master_carbon, ...
        Master_Method,Master_flags)

        PM25_labels   = Master_IDs(PM25_index,:);
        PM25_Barcodes = Master_Barcodes(PM25_index,:);
        PM25_cartridgeIDs = CartridgeIDs_master(PM25_index);
        PM25_LotIDs     = LotIDs(PM25_index,:);
        PM25_projectIDs = projectIDs_master(PM25_index);
        PM25_hours = Master_hours(PM25_index,:);
        
        PM25_masstype  = Master_masstype(PM25_index,:);
        PM25_dates  = Master_dates(PM25_index,:);
        PM25_mass  = Master_mass(PM25_index,:);
        PM25_IC  =   Master_IC(PM25_index,:);
        PM25_ICP  =  Master_ICP(PM25_index,:); 
        PM25_XRF  =  Master_XRF(PM25_index,:); 
        PM25_carbon  = Master_carbon(PM25_index,:);
        PM25_Method  = Master_Method(PM25_index,:);
        
        PM25_flags = Master_flags(PM25_index,:);
end

function [PM25_labels,PM25_Barcodes,PM25_cartridgeIDs,PM25_LotIDs,PM25_projectIDs,...
        PM25_hours,PM25_masstype,PM25_dates,PM25_mass,PM25_IC,PM25_ICP,PM25_XRF,...
        PM25_carbon,PM25_Method,PM25_flags] = ...
        RemoveRows (filter_index, ...
        PM25_labels,PM25_Barcodes,PM25_cartridgeIDs,PM25_LotIDs,PM25_projectIDs,...
        PM25_hours,PM25_masstype,PM25_dates,PM25_mass,PM25_IC,PM25_ICP,PM25_XRF,...
        PM25_carbon,PM25_Method,PM25_flags)
        
        PM25_hours = removerows(PM25_hours, filter_index);
        PM25_labels = removerows(PM25_labels, filter_index);
        PM25_cartridgeIDs = removerows(PM25_cartridgeIDs, filter_index);
        PM25_projectIDs = removerows(PM25_projectIDs, filter_index);
        PM25_Barcodes = removerows(PM25_Barcodes, filter_index);
        PM25_LotIDs   = removerows(PM25_LotIDs, filter_index);
        
        PM25_masstype  = removerows(PM25_masstype, filter_index);
        PM25_dates  = removerows(PM25_dates, filter_index);
        PM25_mass  = removerows(PM25_mass, filter_index);
        PM25_IC  =   removerows(PM25_IC, filter_index);
        PM25_ICP  =  removerows(PM25_ICP, filter_index);
        PM25_XRF  =  removerows(PM25_XRF, filter_index);
        PM25_carbon  = removerows(PM25_carbon, filter_index);
        PM25_Method  = removerows(PM25_Method, filter_index);
        
        PM25_flags = removerows(PM25_flags, filter_index);
        
end

function sumspec = SumSpec(i,XRF,ICP,IC,carbon,IC_col)

    carbon_sum = sum(carbon(i,[1 4]),2,'omitnan'); % Sum of HIPS BC and FTIR OC (if available)

    if ~isnan(XRF(i,2))  % check for presence of XRF data (Al here) first and sum only XRF metals
        % Metals + IC (without F, Cl, K, Mg, Ca) + BC
        spec = [sum(XRF(i,[2:3,5:end]),2,'omitnan')   sum(IC(i,IC_col),2,'omitnan')   carbon_sum];
        % currently don't use Na, S from XRF; ==> use IC_Na IC_SO4 already

    elseif ~isnan(ICP(i,3))  % check for presence of Al from ICP-MS and sum only ICP-MS metals
        spec = [sum(ICP(i,:),2,'omitnan')   sum(IC(i,IC_col),2,'omitnan')   carbon_sum];
     
    else % no Al from ICP or XRF, therefore sum xrf as default
        spec = [sum(XRF(i,[2:3,5:end]),2,'omitnan')   sum(IC(i,IC_col),2,'omitnan')   carbon_sum];
    
    end
    sumspec = sum(spec, 2,'omitnan');
end

function out = copyflag(cart_flags,Flag_parameters)
    
    publicflags = Flag_parameters.Flag;
    out = cart_flags;
    % remove flags
    for i = 1:length(cart_flags)
        out{i}='';
    end
    % add filter flags
    for i = 1:length(publicflags)
        flag_ind = find(contains(cart_flags,publicflags{i})==1);

        for j = 1:length(flag_ind)
            out{flag_ind(j)} = AddFlag( out{flag_ind(j)} , publicflags{i} );
        end

    end

end
