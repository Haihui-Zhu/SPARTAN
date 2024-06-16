function [Titles,Master_IDs, Master_Barcodes, Master_CartridgeIDs, Master_LotIDs, Master_ProjectIDs,Master_hours, Master_masstype, ...
          Master_dates, Master_mass, Master_IC, Master_ICP, Master_XRF,...
          Master_carbon, Master_Method, Master_flags] = ReadMaster(master_file,site_ID)

% This function read the master file when provided the filename (with 
% absolute path) and the Site ID. Site IS is just for error message. 

% Outputs:
% Titles = 83 x 1 cell 
% Master_IDs = N x 1 cell   
% Master_Barcodes = N x 1 cell
% Master_CartridgeIDs = N x 1 cell
% Master_LotIDs = N x 1 cell
% Master_ProjectIDs = N x 1 cell
% Master_hours = N x 1 mat; hour sampled
% Master_masstype = N x 1 mat 
% Master_dates = N x 8 mat; Start year month day and End year month day
% Master_mass = N x 2 mat; 'mass_ug','Volume_m3'
% Master_IC = N x 13 mat
% Master_ICP = N x 21 mat
% Master_XRF = N x 26 mat
% Master_carbon = N x 4 mat; 'BC_SSR_ug', 'BC_HIPS_ug','EC_FTIR_ug','OC_FTIR_ug'
% Master_Method = N x 1 mat
% Master_flags = N x 1 cell

% Written by: Haihui Zhu
% Created: 2021-08-25   

Titles = {'FilterID','Filter_Barcode','CartridgeID','LotID','projectID','hours_sampled','Mass_type',... % 7 cols
    'start_year','start_month','start_day','start_hour','stop_year','stop_month','stop_day','stop_hour',... % 8 cols
    'mass_ug','Volume_m3',... % 2 cols
    'IC_F_ug','IC_Cl_ug','IC_NO2_ug','IC_Br_ug','IC_NO3_ug','IC_PO4_ug',...
    'IC_SO4_ug','IC_Li_ug','IC_Na_ug','IC_NH4_ug','IC_K_ug','IC_Mg_ug',...
    'IC_Ca_ug',... % IC = 13 cols
    'Li_ICP_ng','Mg_ICP_ng','Al_ICP_ng','P_ICP_ng','Ti_ICP_ng',...
    'V_ICP_ng','Cr_ICP_ng','Mn_ICP_ng','Fe_ICP_ng','Co_ICP_ng',...
    'Ni_ICP_ng','Cu_ICP_ng','Zn_ICP_ng','As_ICP_ng','Se_ICP_ng',...
    'Ag_ICP_ng','Cd_ICP_ng','Sb_ICP_ng','Ba_ICP_ng','Ce_ICP_ng',...
    'Pb_ICP_ng',... % ICP = 21 cols
    'Na_XRF_ng','Al_XRF_ng','Si_XRF_ng','S_XRF_ng','Cl_XRF_ng','K_XRF_ng',...
    'Ca_XRF_ng','Ti_XRF_ng','V_XRF_ng','Fe_XRF_ng','Zn_XRF_ng',...
    'Ce_XRF_ng','Pb_XRF_ng','As_XRF_ng','Co_XRF_ng','Cr_XRF_ng',...
    'Cu_XRF_ng','Mg_XRF_ng','Mn_XRF_ng','Ni_XRF_ng','Sb_XRF_ng',...
    'Rb_XRF_ng','Sr_XRF_ng','Cd_XRF_ng','Se_XRF_ng','Sn_XRF_ng',... % 26 cols
    'BC_SSR_ug','BC_HIPS_ug','EC_FTIR_ug','OC_FTIR_ug','method_index','Flags'}; % 6 cols

if exist(master_file,'file')
    fprintf('Master file found for %s \n', site_ID)

    opts = detectImportOptions(master_file);
    opts.ExtraColumnsRule = 'ignore'; 
    opts.VariableTypes{2} = 'char'; % make sure barcode read as char
    Master_data_initial = readtable(master_file,opts); clear opts

    Orig_title = Master_data_initial.Properties.VariableNames; % reading old title but it won't be output. 
    % the output title will be the 'Title' variable. If 'Orig_title' isn't
    % same as 'Title', this script will add column(s).



    if ~isempty(Master_data_initial)
    Master_IDs = cell(table2array(Master_data_initial(:,1)));
    Master_Barcodes = table2array(Master_data_initial(:,2)); % cell won't work when it is full of NaN
    Master_CartridgeIDs = table2array(Master_data_initial(:,3));  
    Master_LotIDs = table2array(Master_data_initial(:,4));
    Master_ProjectIDs = table2array(Master_data_initial(:,5)); 

    for i = 1:length(Master_IDs)
        tID = Master_IDs{i};
        if length(tID)<9
            tID2(1,1:5) = tID(1:5);
            tID2(1,6) = '0';
            tID2(1,7:9) = tID(6:8);
        elseif length(tID)>9
            error('filter ID %d too long',tID)
        else 
            tID2 = tID;
        end
        Master_IDs{i}=tID2;
    end
    
    if isa(Master_LotIDs,'double') == 1
        Ind = find(isnan(Master_LotIDs));
        Ind2 = find(~isnan(Master_LotIDs));
        Master_LotIDs = num2cell(Master_LotIDs);
        for i = 1:numel(Ind)
            Master_LotIDs{Ind(i)} = 'NaN';
        end
        for i = 1:numel(Ind2)
            Master_LotIDs{Ind2(i)} = num2str(Master_LotIDs{Ind2(i)});
        end
    end
%     if isa(Master_Barcodes,'double') == 1
%         Ind = find(isnan(Master_Barcodes));
%         Ind2 = find(~isnan(Master_Barcodes));
%         Master_Barcodes = num2cell(Master_Barcodes);
%         for i = 1:numel(Ind)
%             Master_Barcodes{Ind(i)} = 'NaN';
%         end
%         for i = 1:numel(Ind2)
%             Master_Barcodes{Ind2(i)} = num2str(Master_Barcodes{Ind2(i)});
%         end
%         
%     end
    
    Master_hours = table2array(Master_data_initial(:,6));
    Master_masstype = table2array(Master_data_initial(:,7));
    Master_dates = table2array(Master_data_initial(:,8:15));
    Master_mass = table2array(Master_data_initial(:,16:17)); % mass Volume_m3
    Master_IC  = table2array(Master_data_initial(:,contains(Orig_title,'IC_')));
    Master_ICP = table2array(Master_data_initial(:,contains(Orig_title,'ICP')));
    Master_XRF = table2array(Master_data_initial(:,contains(Orig_title,'XRF'))); 
    Master_carbon = table2array(Master_data_initial(:,end-5:end-2));
    Master_Method = table2array(Master_data_initial(:,contains(Orig_title,'method_index')));
    Master_flags = table2array(Master_data_initial(:,end));
    
    if isa(Master_flags,'double') == 1
        Ind = find(isnan(Master_flags));
        Ind2 = find(~isnan(Master_flags));
        Master_flags = num2cell(Master_flags);
        for i = 1:numel(Ind)
            Master_flags{Ind(i)} = 'NaN';
        end
        for i = 1:numel(Ind2)
            Master_flags{Ind2(i)} = num2str(Master_flags{Ind2(i)});
        end
        
    end
    
    % adding new Na_XRF if wasn't in Orig_title
    if  numel(Orig_title) == 82 && size(Master_XRF,2) == 25
        Master_XRF(:,2:end+1) = Master_XRF;
        Master_XRF(:,1) = NaN; % adding col for Na_XRF
    end
    
    else

    fprintf('WARNING: Master file for %s has no content\n\n', site_ID)
    Master_IDs = cell(0, 1);
    Master_Barcodes = cell(0,1);
    Master_CartridgeIDs = cell(0,1);
    Master_LotIDs = cell(0,1);
    Master_ProjectIDs = cell(0,1);
    Master_hours = NaN.*zeros(0,1);
    Master_masstype = NaN.*zeros(0,1);
    Master_dates = NaN.*zeros(0,8);
    Master_mass = NaN.*zeros(0,2);
    Master_IC  = NaN.*zeros(0,13);
    Master_ICP = NaN.*zeros(0,21);
    Master_XRF = NaN.*zeros(0,26);
    Master_carbon = NaN.*zeros(0,4);
    Master_Method = NaN.*zeros(0,1);
    Master_flags = cell(0, 1);
        
    end

else
    fprintf('WARNING: No master file exists for %s \n\n', site_ID)
    Master_IDs = cell(0, 1);
    Master_Barcodes = cell(0,1);
    Master_CartridgeIDs = cell(0,1);
    Master_LotIDs = cell(0,1);
    Master_ProjectIDs = cell(0,1);
    Master_hours = NaN.*zeros(0,1);
    Master_masstype = NaN.*zeros(0,1);
    Master_dates = NaN.*zeros(0,8);
    Master_mass = NaN.*zeros(0,2);
    Master_IC  = NaN.*zeros(0,13);
    Master_ICP = NaN.*zeros(0,21);
    Master_XRF = NaN.*zeros(0,26);
    Master_carbon = NaN.*zeros(0,4);
    Master_Method = NaN.*zeros(0,1);
    Master_flags = cell(0, 1);

end
