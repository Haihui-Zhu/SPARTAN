%% %%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: read and perform QC on BC data as measured by the HIPS 

% Written by: Haihui Zhu
% Created: Oct 16 2022

close all; clear ; clc
addpath('UtilityFunctions')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% USER SWITCHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

Redo_All_Archieve = 1; % set to 1 if want to re-process all archived data

direc_master = strcat(direc,'/Analysis_Data/Master_files');
direc_HIPS = strcat(direc,'/Analysis_Data/HIPS/Reports'); 
direc_sampling = strcat(direc,'/Site_Sampling'); 
direc_archive = strcat(direc,'/Analysis_Data/Archived_Filter_Data/HIPS/');

% The diary function saves the processing history into an annual record
diary(sprintf('%s/Public_Data/Data_Processing_Records/%s_HIPS_BC_Record.txt',direc,datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))

%-------------   SITE DETAILS   --------------
site_details = readtable(strcat(direc_sampling,'/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));

%------------- define absorption coefficient --------------
sigma = 0.06;

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

        ardir = strcat(direc_archive,Site_codes{loc});
        files = getFiles(ardir);

        for fid = 1:length(files)
            if contains(SkipList,files{fid,1})
                fprintf('%s skipped\n',files{fid,1})
                continue
            else
                tfile = strcat(ardir,'/',files{fid,1});
                [SUCCESS,MESSAGE,MESSAGEID] = copyfile(tfile, direc_HIPS);
                if SUCCESS ~= 1
                    disp(MESSAGE)
                end
            end
        end
       
    end
    fprintf('\nDone copying archived HIPS file to NEW.\n\n')
end

%% --------- Find and read data file in the HIPS raw directory ---------

allfname = getFiles(direc_HIPS);

for ff = 1:length(allfname)
    if allfname{ff}(1) =='.' ||  allfname{ff}(1) =='~' 
        continue
    end
    
    if ismember(allfname{ff}(end-3:end), '.csv')
        tfile = sprintf('%s/%s',direc_HIPS,allfname{ff});
        rawdata = readtable(tfile,'PreserveVariableNames',true);
    else
        tfile = sprintf('%s/%s',direc_HIPS,allfname{ff});
        rawdata = readtable(tfile,'sheet','HIPS Results','PreserveVariableNames',true);
    end
        tabeltitles  = rawdata.Properties.VariableNames;
        fprintf('Reading %s\n',direc_HIPS)
    % useful columns (no need to read volume; conc. will be calculated in Script F)
    col = find(contains(tabeltitles,'ID','IgnoreCase',true)==1);
    if length(col)>1
        col = contains(tabeltitles,{'filter ID','filterID'},'IgnoreCase',true);
    end
    FilterID = table2array(rawdata(:,col));
    
    col = contains(tabeltitles,{'Filter Type','FilterType'},'IgnoreCase',true);
    FilterType = table2array(rawdata(:,col));
    
    col = contains(tabeltitles,'tau','IgnoreCase',true);
    tau = table2array(rawdata(:,col));
    
    col = contains(tabeltitles,{'Deposit Area','DepositArea'},'IgnoreCase',true);
    depArea = table2array(rawdata(:,col));

    col = find(contains(tabeltitles,'comment','IgnoreCase',true)==1);
    if isempty(col)
       comments = [];
    elseif length(col)==1 
        comments = table2array(rawdata(:,col));
    else
        col2 = contains(tabeltitles(col),'filter','IgnoreCase',true);
        comments = table2array(rawdata(:,col(col2)));
    end
    if isa(comments,"double") % sometimes comments are made under the 'lab comment' column
        col2 = contains(tabeltitles(col),'lab','IgnoreCase',true);
        comments = table2array(rawdata(:,col(col2)));
    end
    if isa(comments,"double") % if still 'double', there are no text comments
        comments=[];
    end

    % remove lab blank / field blank
    if ~isempty(FilterType)
        ind1 = find( contains(FilterType,{'LB','FB'},'IgnoreCase',true)== 1); 
        ind2 = find( contains(FilterID,{'LB','FB'},'IgnoreCase',true)  == 1);
        ind = [ind1;ind2];
    else
        ind = find(contains(FilterID,{'LB','FB'},'IgnoreCase',true)==1);
    end 

    for ii = 1:length(ind)
        FilterID{ind(ii)} = '';
        tau(ind(ii)) = NaN;
        depArea(ind(ii)) = NaN;
        if ~isempty(comments)
           comments{ind(ii)} = '';
        end
    end

    % format filter ID
    label_length_new = 11; % i.e. 'AEAZ-0010-1'
    label_length_old = 13; % i.e. '10001-AEAZ-1T'
    format_filter = cell(size(FilterID));
    for h = 1:length(FilterID)
        format_filter{h}='';

        if length(FilterID{h})==label_length_new
            format_filter{h}= strcat(FilterID{h}(1:4),'-',FilterID{h}(6:9));
        
        elseif length(FilterID{h})==label_length_old
            format_filter{h} = strcat(FilterID{h}(7:10),'-',FilterID{h}(3:5));

        elseif length(FilterID{h})==29
            format_filter{h} = FilterID{h}(3:11);

        elseif ~isempty(FilterID{h})
            fprintf('WARNING: label size not 11 or 13 for %s \n',FilterID{h} );
        end
    end


    % calculate BC mass (ug)
    bcmass = depArea.*tau./sigma;
%     bcmass(bcmass<0) = 0; % negative tau


%% ---------- Write HIPS BC to the master files -------------------------
format_filter_char = char(format_filter);
site_IDs = unique(format_filter_char(:,1:4),'rows'); clear format_filter_char

for loc = 1:size(site_IDs,1)
    if ~isempty(deblank(site_IDs(loc,:)))

        master_file = sprintf('%s/%s_master.csv',direc_master,site_IDs(loc,:));
        [Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, projectIDs_master,...
            Master_hours, Mass_type, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
            Master_carbon, Master_Method, Master_flags] = ReadMaster ( master_file, site_IDs(loc,:) );

        if ~isempty(Master_IDs)
            Master_HIPS = Master_carbon(:,2);

            siteidx = find( contains(format_filter,site_IDs(loc,:)) == 1 );
            tbcmass = bcmass(siteidx);

            % writing new HIPS data & flags to Master_carbon

            for ii = 1:length(siteidx)
                filteridx = find( contains(Master_IDs,format_filter{siteidx(ii)}) == 1);
                if isempty(filteridx)
                    % sometimes filters named XXXX_0XXX in the master file but
                    % XXXX_XXX in the raw data file, or the other way. Double check here:
                    if length(format_filter{siteidx(ii)})==8
                        newname = strcat(format_filter{siteidx(ii)}(1:5),'0',format_filter{siteidx(ii)}(6:8));
                    else
                        newname = strcat(format_filter{siteidx(ii)}(1:5),format_filter{siteidx(ii)}(7:9));
                    end
                    filteridx = find( contains(Master_IDs,newname) == 1);
                end

                if isempty(filteridx)
                    warning('%s not found in master file',format_filter{siteidx(ii)})
                else
                    tMasterHIPS = Master_HIPS(filteridx);

                    if isnan(tMasterHIPS) || isempty(tMasterHIPS)  % no HIPS for this filter in master file
                        Master_HIPS(filteridx) = tbcmass(ii);
                    elseif  abs(tMasterHIPS-tbcmass(ii)) < 0.3 % there is already HIPS data but is slightly different from the new HIPS data 
                        Master_HIPS(filteridx) = tbcmass(ii);
                    elseif ~isnan(tbcmass(ii))
                        Master_HIPS(filteridx) = tbcmass(ii);
                        warning ('HIPS data exists for %s (%.2f), new data = %.2f.',format_filter{siteidx(ii)},tMasterHIPS,tbcmass(ii))
                    end
                    
                    % Add flags
                    if ~isempty(comments)
                        tflag = '';
                        tcomment = comments{siteidx(ii)};

                        % VIH
                        if contains(tcomment,{'weight','inhomogenous','concentrated'},'IgnoreCase',true)
                            Master_flags{filteridx} = AddFlag(Master_flags{filteridx},'VIH');
                            tflag = AddFlag(tflag,'VIH');
                        end

                        % FS
                        if contains(tcomment,{'stretch','wrinkle'},'IgnoreCase',true)
                            Master_flags{filteridx} = AddFlag(Master_flags{filteridx},'FS');
                            tflag = AddFlag(tflag,'FS');
                        end

                        % PH
                        if contains(tcomment,{'Pin','hole'},'IgnoreCase',true)
                            Master_flags{filteridx} = AddFlag(Master_flags{filteridx},'PH');
                            tflag = AddFlag(tflag,'PH');
                        end

                        % VCT
                        if contains(tcomment,{'contaminat','streak','particle','flecks'},'IgnoreCase',true)
                            Master_flags{filteridx} = AddFlag(Master_flags{filteridx},'VCT');
                            tflag = AddFlag(tflag,'VCT');
                        end

                        % SGP
                        if contains(tcomment,{'screen pattern','grid'},'IgnoreCase',true)
                            Master_flags{filteridx} = AddFlag(Master_flags{filteridx},'SGP');
                            tflag = AddFlag(tflag,'SGP');
                        end

                        % MTL
                        if contains(tcomment,{'metallic','reflect'},'IgnoreCase',true)
                            Master_flags{filteridx} = AddFlag(Master_flags{filteridx},'MTL');
                            tflag = AddFlag(tflag,'MTL');
                        end

                        % SUS
                        if contains(tcomment,{'suspect','compromise','QUESTIONABLE'},'IgnoreCase',true)
                            Master_flags{filteridx} = AddFlag(Master_flags{filteridx},'SUS');
                            tflag = AddFlag(tflag,'SUS');
                        end

                        if  isempty(tflag) && ~isempty(tcomment) && bcmass(siteidx(ii))>0
                            fprintf('%s comment not added to flags: %s\n',format_filter{siteidx(ii)},tcomment)
                        end
                    end

                end
            end

            % Write master file
            Master_carbon(:,2) = Master_HIPS;
            WriteToMaster(Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs_master, projectIDs_master,...
                Master_hours, Mass_type, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
                Master_carbon, Master_Method, Master_flags,direc_master,site_IDs(loc,:))

        else
            fprintf('HIPS data not written to %s_master.csv\n',site_IDs(loc,:))
        end


        % copy/move file to the corresponding archive
        file_destination = strcat(direc_archive,sprintf('%s/',site_IDs(loc,:)));
        [status(loc),~,~]= copyfile(tfile, file_destination );
        clear file_destination


    end
end
if sum(status == 0) == 0
    delete(tfile)
else
    Warning ('File moving not successful. Not deleted: %s\n',tfile)
end
clear format_filter_char format_filter comments FilterID tau depArea FilterType

end
fprintf('END of HIPS BC processing for %s \n', datestr(today))

diary off


