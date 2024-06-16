% This script goes through all the Neph new and archived files to re
% generate CR data files with the period before CR. This will help
% diagnose. 

% Written by: Haihui Zhu
% Created: 17 Jan 2022

close all; clear; 
clc
warning('off','backtrace')


%% %%%%%%% USER SWITCHES %%%%%%%%%%
% rootdir = '/Volumes/rvmartin/Active/haihuizhu/6.SPARTAN'; % for Haihui
rootdir = '/Volumes/rvmartin/Active/SPARTAN-shared'; % for Haihui
% rootdir = '\\storage1.ris.wustl.edu\rvmartin\Active\SPARTAN-shared';
direc_outputCR = sprintf('%s/Public_Data/Neph_Processed/CR_data_redo',rootdir);

% The diary function saves the processing history into an annual record
diary(sprintf('%s/Public_Data/Data_Processing_Records/%s_Neph_CR_redo',rootdir,datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))
time_now = datestr(today);

site_details = readtable(sprintf('%s/Site_Sampling/Site_details.xlsx',rootdir),'PreserveVariableNames',true); 
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));

%% %%%%%%% Reading data files %%%%%%%%%%%%%%%%%%

for loc =  [1 5 15 34]%:length(Site_codes)
    
    % find files from both Archived and Rejected folders 
    for DIR = 1:3
        if DIR == 1
            Folder = 'Rejected';
        elseif DIR == 2
            Folder = 'Archived';
        else
            Folder = 'Raw';
        end
    
    
    direc = sprintf ('%s/Site_Sampling/%s_%s/Nephelometer_Data/%s',rootdir,Site_codes{loc},Site_cities{loc},Folder);
        
    names = dir(sprintf('%s/IN*',direc));
    files = {names.name}; clear names
    files = files';
    
    if isempty(files) ==1
        names = dir(sprintf('%s/LOG*',direc));
        files = {names.name}; clear names
        files = files';
        fprintf('\nNo data files found for %s %s-folder\n\n', Site_cities{loc},Folder)
        continue
    end
    
    fprintf('\nReading and aggregating data files for %s %s-folder\n', Site_cities{loc},Folder)
    
    data = []; CR_data = []; RemovedData = [];
    
    reject_index = [];
    
    for Neph=1:length(files)
        reject = 0; % set to 1 in neph_diag if file is rejected 
        
        fid=fopen(sprintf('%s/%s',direc,files{Neph}));
        filename = sprintf('%s/%s',direc,files{Neph});
        
        test_file=dir(filename);
        file_size=test_file.bytes;
        
        if file_size < 5000
            fprintf('Filesize too small.\n')
            continue
        end
        
        fprintf('reading %s from %s \n',files{Neph}, Site_cities{loc})
        
        
        % This section gets info for what type (e.g IN model, CR system) of file is being read.
        
        headers2 =textscan(fid, '%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,\r\n]', 1);
        headers = string(headers2); clear headers2
        model_test = find(contains(headers,["Flow","FLOW","flow"]));
        IN_model = [];
        CR_test = find(contains(headers,'RF_CR'));
        
        if isempty(model_test) == 1
            if isempty(IN_model) == 1
                IN_model = 1; % regular fan
            end
        else
            IN_model = 2; % variable speed fan
        end
        
        if isempty(CR_test) == 1
            CR_system = 0; % No CR system
        else
            CR_system = 1; %CR_system == 1 indicates output in last 6 columns
        end
        
        [RemovedData0, CR_datab, neph_num] = neph_file_CR_redo(fid, CR_system, IN_model);
        if ~isempty(RemovedData0)
            RemovedDatab = [RemovedData0(:,1) RemovedData0(:,end-15:end)]; % Back scatt Red to CR
        else
            continue
        end 
        clear RemovedData0 

        if CR_system == 1
            CR_data = [CR_data; CR_datab];
            RemovedData = [RemovedData;RemovedDatab]; % only collect when CR == 1
        end
    end
    
    % write to new CR file
    if ~isempty(CR_data)
        Dates_CR = datevec(CR_data(:,1));
        start_date = Dates_CR(1,1)*100 + Dates_CR(1,2);
        end_date = Dates_CR(end,1)*100 + Dates_CR(end,2);
        
        Merged_CR = [Dates_CR(:,1:3) CR_data(:,2:end)];
        
        Titles_CR = {'Year' 'Month' 'Day' 'RF_CR' 'GF_CR' 'BF_CR' 'RB_CR' 'GB_CR' 'BB_CR'};
        
        fileID = fopen(sprintf('%s/%s_%06d_to_%06d_NephCRdata_%s.csv',direc_outputCR,Site_codes{loc},start_date,end_date,Folder),'w');
        fprintf(fileID,'Created on %s \n',time_now);
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n',Titles_CR{1,1:9});
        fprintf(fileID,'%04d, %02d, %02d, %6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f\n', Merged_CR.');
        fclose(fileID);
        
        clear start_date end_date Merged_CR Dates_CR 
        
    end

    if ~isempty(RemovedData)
        Dates = datevec(RemovedData(:,1));
        start_date = Dates(1,1)*100 + Dates(1,2);
        end_date = Dates(end,1)*100 + Dates(end,2);
        
        Merged = [Dates(:,1:3) RemovedData(:,2:end)];
        
        Titles = {'Year' 'Month' 'Day' 'Back scatt Red' 'Back scatt Green' 'Back scatt Blue' 'Tot scatt Red' 'Tot scatt Green' 'Tot scatt Blue' 'Fan RPM' 'Pres' 'Flow' 'DAC' 'RF_CR' 'GF_CR' 'BF_CR' 'RB_CR' 'GB_CR' 'BB_CR'};
        
        fileID = fopen(sprintf('%s/%s_%06d_to_%06d_RM_Nephdata_%s.csv',direc_outputCR,Site_codes{loc},start_date,end_date,Folder),'w');
        fprintf(fileID,'Created on %s \n',time_now);
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',Titles{1,1:end});
        fprintf(fileID,'%04d, %02d, %02d, %6.2f,%6.2f,%6.2f, %6.2f,%6.2f,%6.2f, %6.2f,%6.2f,%6.2f, %6.2f, %6.2f,%6.2f,%6.2f, %6.2f,%6.2f,%6.2f\n', Merged.');
        fclose(fileID);
        
        clear start_date end_date Merged Dates 
        
    end
    
    fprintf('Saved CR-additional data for %s %s\n \n', Site_cities{loc},Folder)
    pause(5)
    end
end


%% Functions

% this function reads in nephelometer files
function [RemovedData, CR_datab, neph_num] = neph_file_CR_redo(fid, CR_system, IN_model)

%% IN101 models that AUTOMATICALLY account for CR correction
%%%%%%%%%%%%%%%%%%%

CR_datab = [];

if IN_model == 1 || IN_model == 3
    
    if CR_system == 1
        
        % this file format has a date format with dashes
        A = textscan(fid,    '%s  %n-%n-%n  %n:%n:%n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1);
        if isempty(A{1,8}) % did not read first value after date/time column
            clear A
            % this file format has a date format with slashes
            A = textscan(fid,'%s  %n/%n/%n  %n:%n:%n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1);
        end
        % To delete incomplete last raw due to (maybe) an SD card memory issue [Haihui 2021-6-6]
        for i = 1:size(A,2)
            if length(A{i}) ~= length(A{end})
                A{i}(end) = []; end
        end
        %%%%%
        neph_num = A{1,1};
        dataB = cell2mat(A(2:end)); clear A
        index = find(dataB(:,3) < 2000);
        dataB(index,3) = dataB(index,3) + 2000; clear index
        Dates = datenum(dataB(:,3), dataB(:,1), dataB(:,2), 0,0,0); % creating a serial number for the date only (YYYY, MM, DD)
        
        data_fileB = [Dates dataB(:,4) dataB(:,7:end)];
        
        % work with CR system
        fan_offB = find(data_fileB(:,27) == 0);
        
        ind = find(isnan(dataB(:,32)) == 0); % looks for non-NaN values in RF_CR column
        CR = [Dates(ind) dataB(ind, 32:37)];
        CR_datab = [CR_datab; CR];
        clear Dates dataB ind CR
        
        if isempty(fan_offB) == 1
            RemovedData = data_fileB;
        elseif length(fan_offB) < 200
            % this section checks that there is enough space before and after
            % CR system is on to remove +/- 10 more data points
            if fan_offB(1) < 10
                start = 1;
            else
                start = fan_offB(1) - 9;
            end
            if fan_offB(end) + 10 > length(data_fileB(:,26))
                stop = length(data_fileB(:,26));
            else
                stop = fan_offB(end) + 10;
            end
            
            fan_off = start:stop; fan_off = fan_off'; clear fan_offB
            
            RemovedData = data_fileB(fan_off,:);
%             [data_file, ~] = removerows(data_fileB,fan_off); clear fan_off
        else
            fan_off = [];
            ind_B = find(diff(fan_offB) > 5);
            ind = [1; ind_B+1];
            for i = 1:length(ind)
                if fan_offB(ind(i)) < 10
                    start = 1;
                else
                    start = fan_offB(ind(i)) - 9;
                end
                if i == length(ind)
                    if  fan_offB(end) + 10 > length(data_fileB(:,26))
                        stop = length(data_fileB(:,26));
                    else
                        stop = fan_offB(end) + 10;
                    end
                else
                    stop = fan_offB(ind(i+1)-1) + 10;
                end
                fan_offC = start:stop; fan_offC = fan_offC';
                fan_off = [fan_off; fan_offC]; clear fan_offC
            end
%             [data_file, ~] = removerows(data_fileB,fan_off);

            RemovedData = data_fileB(fan_off,:);
            percent_off = (size(fan_off,1)/size(data_fileB,1))*100; clear fan_off
            if percent_off > 25
                warning('Fan off for > 25 % of file')
            end
            
            clear fan_off i ind ind_B
        end
        
        
        %% %%%%  No CR system in place %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    elseif CR_system == 0
        A = [];
        while feof(fid) ==0
            % this file format has a date format with dashes
            B = textscan(fid,'%s  %n-%n-%n  %n:%n:%n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1,'CollectOutput',1);
            if size(B,2) < 31 % did not read all of file for some reason
                % this file format has a date format with slashes
                B = textscan(fid,'%s  %n/%n/%n  %n:%n:%n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1,'CollectOutput',1);
            end
            % To delete incomplete last raw due to (maybe) SD card memory issue [Haihui 2021-6-6]
            for i = 1:size(B,2)
                if length(B{i}) ~= length(B{end})
                B{i}(end) = []; end
            end
            B2 = cell2mat(B(2:end));
            A = [A;B2]; 
            neph_num = B{1,1}; clear B B2
        end
        dataB = A; clear A
        index = find(dataB(:,3) < 2000);
        dataB(index,3) = dataB(index,3) + 2000; clear index
        Dates = datenum(dataB(:,3), dataB(:,1), dataB(:,2), 0,0,0); % creating a serial number for the date only (YYYY, MM, DD)
        data_fileB = [Dates dataB(:,4) dataB(:,7:end)];
        
        data_file = data_fileB;
        
        fan_off1 = find(data_file(:,27) == 0 | isinf(data_file(:,27)));
        fan_off2 = find(data_file(:,27) > 50000); % when the fan is dead/dying the rpm often gets displayed as a really high number rather than zero
        fan_off = [fan_off1; fan_off2];
%         [data_file, ~] = removerows(data_fileB,fan_off);
        RemovedData = data_fileB(fan_off,:);
        percent_off = (size(fan_off,1)/size(data_fileB,1))*100; clear fan_off
        if percent_off > 25
            warning( 'Fan off for > 25 % of file')
        end
        
    end
    
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% IN102 models with variable speed fan and AUTOMATICALLY accounts for CR correction
    
    if IN_model == 2
        
        if CR_system == 1
            
            % this file format has a date format with dashes
            A = textscan(fid,    '%s  %n/%n/%n  %n:%n:%*n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1);
            if isempty(A{1,8}) % did not read first value after date/time column
                clear A
                % this file format has a date format with slashes
                A = textscan(fid,'%s  %n-%n-%n  %n:%n:%*n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1);
            end
            if isempty(A{1,8}) % did not read first value after date/time column, no seconds in file time
                clear A
                % this file format has a date format with slashes
                A = textscan(fid,'%s  %n-%n-%n  %n:%n      %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1);
            end
            if isempty(A{1,8}) % has one 'DATE_TIME' column and 6 addtional dates and time columns 
                clear A
                A = textscan(fid,'%s %s %n %n %n %n %n %n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %s %n  %n','delimiter',',','headerlines',1);
            end
            if length(A) > 40 % raw data with additional dates and time columns
                A([2 8 24 25 26 46:end]) = []; % clear the additional columns 
            end
            % To delete incomplete last raw due to (maybe) SD card memory issue [Haihui 2021-6-6]
            for i = 1:size(A,2)
                if length(A{i}) ~= length(A{end})
                A{i}(end) = []; end
            end
            %%%%%%
            neph_num = A{1,1};
            dataB = cell2mat(A(2:end)); clear A 
            if dataB(1,3) >12 & dataB(1, 1:3) < 2000 % therefore third column must be year column
                index = find(dataB(:,3) < 2000);
                dataB(index,3) = dataB(index,3) + 2000; clear index
                Dates = datenum(dataB(:,3), dataB(:,1), dataB(:,2), 0,0,0); % creating a serial number for the date only (YYYY, MM, DD)
            elseif dataB(1,1) >12 % therefore first column must be year column
                index = find(dataB(:,1) < 2000);
                dataB(index,1) = dataB(index,1) + 2000; clear index
%                 % ----- site-specific correction----- % this is for when the date/time in the neph file is not correct (off by some set amount of time)
%                 Dates_datenum = datenum(dataB(:,1), dataB(:,2), dataB(:,3), dataB(:,4), dataB(:,5),0);
%                 Dates_correct = Dates_datenum+ 16.6750; 
%                 Dates_datevec = datevec(Dates_correct); clear Dates_correct Dates_datenum
%                 dataB(:,1:5) = Dates_datevec(:,1:5);
%                 % ----------------
%                 dataB(:,5) =  dataB(:,5) - 80; % this is a weird exception, minutes is off by 80...
                  % --------
                Dates = datenum(dataB(:,1), dataB(:,2), dataB(:,3), 0,0,0); % creating a serial number for the date only (YYYY, MM, DD)
            end
            data_fileB = [Dates dataB(:,4) dataB(:,6:end)]; % [dates, hours, Scatter data]
            
% if time in neph file is not local, use this short section to correct as
% well as add in fraction of time differeence to lines 160 and 164 above 
% e.g. Dates = datenum(dataB(:,1), dataB(:,2), dataB(:,3), dataB(:,4),0,0) + 11/24; is  +11 hour time difference 
% test = datevec(Dates); 
% Dates2 = datenum(test(:,1), test(:,2), test(:,3),0,0,0);
%             data_fileB = [Dates2 test(:,4) dataB(:,6:end-6)]; clear test Dates
%             Dates = Dates2; clear Dates2
            
            % work with CR system
            ind = find(isnan(dataB(:,35)) == 0);
            CR = [Dates(ind) dataB(ind, 34:39)];
            CR_datab = [CR_datab; CR];
            clear Dates dataB ind CR
            % remove Fan off data
            fan_onB = find(data_fileB(:,27) > 0); 
            fan_offB = find(data_fileB(:,27) == 0); 
            threshold = 5; % 120 = CR; 100 could be a dying fan; larger than 5 the flowing could be 
                
            if isempty(fan_offB) == 1 % no fan_off data; nothing needs to be removed
                RemovedData = data_fileB;
                
            elseif isempty(fan_onB) == 1 
                RemovedData =[];
            else % good chance we have a CR or a dead or dying fan
                fan_off = []; % will collect all the fan off data
                
                ind_B = find(diff(fan_onB) >= threshold);  
                % ind_B find fan-off longer than 75 sec, fan_onB(ind_B+1) 
                % locating the first data after the off period 

                if fan_onB(1) >= threshold % for file starting with a CR or dying fan
                    fan_restart(1) = fan_onB(1);
                    fan_restart(2:numel(ind_B)+1) = fan_onB(ind_B+1); 
                    % fan_restart = all the data right after a CR or dying fan
                else
                    fan_restart = fan_onB(ind_B+1); 
                end
                
                for i = 1:length(fan_restart)
                    if fan_restart(i) <= 130 
                        start = 1; 
                        stop  = fan_restart(i)+10;
                    else
                        start = fan_restart(i) - 130;  % - 130 is to ensure all the fan off are included 
                        stop  = fan_restart(i) + 9;    % the first few (~10) data discarded 
                    end 
                    if stop > length(data_fileB(:,26))
                        stop = length(data_fileB(:,26));
                    end
                    fan_offC = [start:stop]'; 
                    fan_off = [fan_off; fan_offC]; clear fan_offC
                end
                fan_off = [fan_off; fan_offB]; % Sometimes there are random 
                % fan_off data that not included in ind_B. Should be also removed.
                fan_off = unique(fan_off); % remove repeating rows 
                RemovedData = data_fileB(fan_off,:);
%                 [data_file, ~] = removerows(data_fileB,fan_off);
                percent_off = (size(fan_off,1)/size(data_fileB,1))*100;
                if percent_off > 25
                    warning( 'Fan off for > 25 % of file')
                end
                
                clear fan_off i fan_onB percent_off
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        elseif CR_system == 0
            
            try
                % this file format has a date format with dashes
                A = textscan(fid,'%s  %n-%n-%n  %n:%n:%n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1);
                neph_num = A{1,1};
                dataB = cell2mat(A(2:end)); clear A
                index = find(dataB(:,3) < 2000);
                dataB(index,3) = dataB(index,3) + 2000; clear index
                Dates = datenum(dataB(:,3), dataB(:,1), dataB(:,2), 0,0,0); % creating a serial number for the date only (YYYY, MM, DD)
                data_fileB = [Dates dataB(:,4) dataB(:,7:end-1)];
            catch
                % this file format has a date format with slashes
                A = textscan(fid,'%s  %n/%n/%n  %n:%n:%n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1);
                neph_num = A{1,1};
                dataB = cell2mat(A(2:end)); clear A
                index = find(dataB(:,3) < 2000);
                dataB(index,3) = dataB(index,3) + 2000; clear index
                Dates = datenum(dataB(:,3), dataB(:,1), dataB(:,2), 0,0,0); % creating a serial number for the date only (YYYY, MM, DD)
                data_fileB = [Dates dataB(:,4) dataB(:,7:end-1)];
            end
            
            data_file = data_fileB;
            
            fan_off1 = find(data_file(:,27) == 0 | isinf(data_file(:,27)));
            fan_off2 = find(data_file(:,27) > 50000); % when the fan is dead/dying the rpm often gets displayed as a really high number rather than zero
            fan_off = [fan_off1; fan_off2];

            RemovedData = data_fileB(fan_off,:);
%             [data_file, ~] = removerows(data_fileB,fan_off);
            percent_off = (size(fan_off,1)/size(data_fileB,1))*100; clear fan_off
            if percent_off > 25
                warning( 'Fan off for > 25 % of file')
            end
            
        end
        
    end
    
end




