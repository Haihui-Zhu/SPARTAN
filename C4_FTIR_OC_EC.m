%% %%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: read OC data from FTIR reports provided by Ann Dillner from UCD

% Written by: Haihui Zhu
% Created: May 22 2024

close all; clear ; clc
addpath('UtilityFunctions')
warning('backtrace','off')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% USER SWITCHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

Redo_All_Archieve = 0; % set to 1 if want to re-process all archived data

direc_master = strcat(direc,'/Analysis_Data/Master_files');
direc_FTIR = strcat(direc,'/Analysis_Data/FTIR/raw_data_and_reports/'); 
direc_sampling = strcat(direc,'/Site_Sampling'); 
direc_archive = sfmkdir(strcat(direc,'/Analysis_Data/Archived_Filter_Data/FTIR/'));

% The diary function saves the processing history into an annual record
diary(sprintf('%s/Public_Data/Data_Processing_Records/%s_FTIR_Record.txt',direc,datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))

%-------------   SITE DETAILS   --------------
site_details = readtable(strcat(direc_sampling,'/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));

%% --------- Pull out Archived if needed ---------
if Redo_All_Archieve == 1
    
    % read the list of files that do not need to reprocess:
    fname = strcat(direc_archive,'No_need_to_reprocess.txt');
    fileID = fopen(fname, 'r');
    if fileID == -1
        error('File could not be opened for appending');
    end
    SkipList = {};
    while ~feof(fileID)
        line = fgetl(fileID);
        SkipList{end+1} = line;
    end
    fclose(fileID);

    % moving archive data to the new dir
    for loc = 1:length(Site_codes)

        files = getFiles(direc_archive);

        for fid = 1:length(files)
            if contains(SkipList,files{fid,1})
                fprintf('%s skipped\n',files{fid,1})
                continue
            else
                tfile = strcat(direc_archive,'/',files{fid,1});
                [SUCCESS,MESSAGE,MESSAGEID] = copyfile(tfile, direc_FTIR);
                if SUCCESS ~= 1
                    disp(MESSAGE)
                end
            end
        end
       
    end
    fprintf('\nDone copying archived FTIR file to NEW.\n\n')
end
%% --------- Find and read data file in the FTIR raw directory ---------

allfname = getFiles(direc_FTIR);

for ff = 1:length(allfname)
    if allfname{ff}(1) =='.' ||  allfname{ff}(1) =='~' 
        continue
    end

    % read the raw data
    fname = strcat(direc_FTIR,allfname{ff});
    opts = detectImportOptions(fname);
    opts.VariableNamesRange = 'A2';
    rawdata = readtable(fname,opts);
    title = rawdata.Properties.VariableNames;
    sample_name = table2array(rawdata(:,contains(title,'sample','IgnoreCase',true)));

    % select only Tensor II_SN151 data
    ind = find(contains(sample_name,'Tensor II_SN151'));
    sample_name = sample_name(ind);
    rawdata = rawdata(ind,:);

    % read other needed info
    site_code = table2array(rawdata(:,contains(title,'site','IgnoreCase',true)));
    volume = table2array(rawdata(:,contains(title,'volume','IgnoreCase',true))); % m3
    if sum(contains(title,'FTIR_OC'))>0
        OC = volume .* rawdata.FTIR_OC; % ug
        EC = volume .* rawdata.FTIR_EC; % ug
    else 
        OC = volume .* rawdata.OC; % ug
        EC = volume .* rawdata.EC; % ug
    end
    note = table2array(rawdata(:,contains(title,'note','IgnoreCase',true)));

    % extract filter ID from the sample_names
    filter_id = cell(size(sample_name));
    for ss = 1:length(sample_name)
        t_sample_name = sample_name{ss};
        filter_id{ss} = extract_filter_id(t_sample_name);
    end

    % find the sites included in this report
    unique_site = unique(site_code);

    % loop through the sites and read the corresponding master files
    for loc = 1:length(unique_site)
        
        master_file = sprintf('%s/%s_master.csv',direc_master,unique_site{loc});
        [Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, projectIDs_master,...
            Master_hours, Mass_type, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
            Master_carbon, Master_Method, Master_flags] = ReadMaster ( master_file, unique_site{loc});
        
        ind = find(contains(filter_id,unique_site{loc}));
        filter_id_site = filter_id(ind);

        if ~isempty(Master_IDs)
            % add FTIR data to Master_carbon
            % Master_carbon = N x 4 mat; 'BC_SSR_ug', 'BC_HIPS_ug','EC_FTIR_ug','OC_FTIR_ug'
            for ii = 1:length(ind)
                master_ind = find(contains(Master_IDs,filter_id{ind(ii)}));
                if ~isempty(master_ind) % found filter in the master file
                    % regardless of if current master file contains any
                    % FTIR data, write the new data in.

                    % Give a warning if overwriting with a new value:
                    if ~isnan(Master_carbon(master_ind,4)) && abs(Master_carbon(master_ind,4) - OC(ind(ii))) > 0.1
                        warning('Overwriting %s FTIR EC and OC in master: %.2f and %.2f.\n                                             New values: %.2f and %.2f',...
                            Master_IDs{master_ind},...
                            Master_carbon(master_ind,3), Master_carbon(master_ind,4),...
                            EC(ind(ii)), OC(ind(ii)))
                    end
                    Master_carbon(master_ind,3) = EC(ind(ii));
                    Master_carbon(master_ind,4) = OC(ind(ii));
                    
                else
                    fprintf('%s not found in master file\n',Master_IDs{master_ind})
                end

            end
            % Write master file
            WriteToMaster(Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, projectIDs_master,...
                Master_hours, Mass_type, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
                Master_carbon, Master_Method, Master_flags,direc_master,unique_site{loc})

        else
            fprintf('FTIR data not written to %s_master.csv\n',unique_site{loc})
        end

    end


    % copy/move file to the corresponding archive
    [status,msg,msgid] = force_movefile(fname, direc_archive );
    if status == 1
        fprintf('Successfully moved file to Archived: %s\n',fname)
    end

end
disp('Done!')


%% Functions
function filter_id = extract_filter_id(sample_name)
        % example sample name: Tensor II_SN151_S_13034_AEAZ_2T_PM25_05_05_2019.0.csv
        loc1 = find(sample_name == '_'); % filter name after third '_'
        if contains(sample_name , {'PM25','PM2_5','PM2.5'},'IgnoreCase',true)
            loc2 = find(sample_name == 'P'); % filter name before 'PM25'
        elseif contains(sample_name , 'FB') % cannot find 'PM25' could be a field blank 'FB'
            loc2 = find(sample_name == 'F');
            loc2 = loc2(end); % the last F is in 'FB'
        end
        if exist('loc2','var')
            FilterID = sample_name(loc1(3)+1 : loc2(end)-2);
        else
            FilterID = sample_name(loc1(3)+1 : loc1(5)-1);
        end

        label_length_new = 11; % i.e. 'AEAZ-0010-1'
        label_length_old = 13; % i.e. '10001-AEAZ-1T'

        if length(FilterID)==label_length_new
            filter_id= strcat(FilterID(1:4),'-',FilterID(6:9));

        elseif length(FilterID)==label_length_old
            filter_id = strcat(FilterID(7:10),'-0',FilterID(3:5));

        elseif length(FilterID)== 9 % standard name
            filter_id = FilterID;

        % filter id with minor typos:
        elseif length(FilterID)==29
            filter_id = FilterID(3:11);

        elseif length(FilterID) == label_length_new-1 % % i.e. 'AEAZ-010-1'
            filter_id= strcat(FilterID(1:4),'-0',FilterID(6:8));


        elseif ~isempty(FilterID)
            fprintf('WARNING: %s filter name cannot be standardized.\n',FilterID );
        end

end