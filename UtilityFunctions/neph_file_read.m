% this function reads in nephelometer files
% ---- Output ----
% data_file:
% N x 30 in size
% Contains datenum(1), hour(2), detailed scatter data (3:14), met data(15:17), dark reference(18:20), back and total scatter (21:26) and Fan RPM (27), Pressure (28), Flow (29), DAC (30).
% Datenum is a number, including year month date info
% The fan off periods are excluded. 
%
% CR_datab:
% Clean Air Reference scatter data for red, green, and blue light, forward and backward. Therefore a total of 6 columns

function [data_file, CR_data, FileTracker] = neph_file_read(filename, FileTracker, data_title, cr_title)
CR_system = FileTracker.CR;
IN_model = FileTracker.Model;
% data_title = {'DateNum','Hour','RF_PMT', 'RF_Ref', 'GF_PMT', 'GF_REF', 'BF_PMT', 'BF_REF', 'RB_PMT', 'RB_REF', 'GB_PMT', 'GB_REF', 'BB_PMT', 'BB_REF',... 1-14
%     'Temp' ,'Amb_Pressure', 'RH', ... 15-17
%     'PMT_DARK', 'Forw_Dark_Ref', 'Back_Dark_Ref', ... 18:20
%     'Back_scatt_Red', 'Back_scatt_Green', 'Back_scatt_Blue', 'Total_scatt_Red', 'Total_scatt_Green', 'Total_scatt_Blue', ... 21-26
%     'Fan_RPM', 'Flow_Pressure', 'Flow', 'DAC'}; % 27:30

% cr_title = {'RF_CR', 'GF_CR', 'BF_CR', 'RB_CR', 'GB_CR', 'BB_CR' };

% NOTE:
% Titles lised above is just a reference. The actual data_title is passed from script G. Should update it here if changes are made in G.
% PMT = photomultiplier tube. Measurement of scatter with aerosol 
% There are two pressures - one ambient, one for the flow, not all file have both available. The ReadWithTitile takes care of this.

%% start reading

if CR_system == 0

    % ReadWithTitilen function reads Neph data and standardizes the format
    [data_fileB, CR_data] = ReadWithTitile(filename, data_title, cr_title, IN_model, CR_system);

    % remove rows when the fan is off
    fan_off1 = find(data_fileB(:,27) == 0 | isinf(data_fileB(:,27)));
    fan_off2 = find(data_fileB(:,27) > 50000); % when the fan is dead/dying the rpm often gets displayed as a really high number rather than zero
    fan_off = [fan_off1; fan_off2];
    [data_file, ~] = removerows(data_fileB,fan_off);

    percent_off = (size(fan_off,1)/size(data_fileB,1))*100; clear fan_off
    if percent_off > 25
        warning( 'Fan off for > 25 % of file')
        FileTracker.FanOff25Perc = 1;
    end

elseif CR_system == 1

    %% IN101 models that AUTOMATICALLY account for CR correction
    %%%%%%%%%%%%%%%%%%%
    if IN_model == 1 || IN_model == 3


        % ReadWithTitilen function reads Neph data and standardizes the format
        [data_fileB, CR_data] = ReadWithTitile(filename, data_title, cr_title, IN_model, CR_system);


        % remove rows when the fan is off
        fan_offB = find(data_fileB(:,27) == 0);

        fan_off = [];
        if isempty(fan_offB) == 1
            data_file = data_fileB;
        elseif length(fan_offB) < 200 % one fan off contains about 120 rows of data. < 200 means there is only one fan off
            % this section checks that there is enough space before and after CR system is on to remove +/- 10 more data points
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

            [data_file, ~] = removerows(data_fileB,fan_off);

            percent_off = (size(fan_off,1)/size(data_fileB,1))*100;
            if percent_off > 25
                warning('Fan off for > 25 % of file')
            end

            clear fan_off

        else % there are more than 1 fan off
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
            [data_file, ~] = removerows(data_fileB,fan_off);
            clear i ind ind_B


            percent_off = (size(fan_off,1)/size(data_fileB,1))*100;
            if percent_off > 25
                warning('Fan off for > 25 % of file')
            end

            clear fan_off
        end

    end

    %% IN102 models with variable speed fan and AUTOMATICALLY accounts for CR correction
    if IN_model == 2

        % ReadWithTitilen function reads Neph data and standardizes the format
        [data_fileB, CR_data] = ReadWithTitile(filename, data_title, cr_title, IN_model, CR_system);
    
        % remove rows when the fan is off
        fan_onB = find(data_fileB(:,27) > 0);
        fan_offB = find(data_fileB(:,27) == 0);

        threshold = 5; % 120 = CR; 100 could be a dying fan; larger than 5 the flowing could be

        if isempty(fan_offB) == 1 % no fan_off data; nothing needs to be removed
            data_file = data_fileB;

        elseif isempty(fan_onB) == 1
            data_file =[];
            warning('Fan off for > 25 % of file')
            FileTracker.FanOff25Perc = 1;
        else
            % good chance we have a CR or a dead or dying fan
            fan_off = []; % will collect all the fan off data to be removed
            fan_on_gap = diff(fan_onB);
            ind_B = find(fan_on_gap >= threshold); clear fan_on_gap

            if fan_onB(1) >= threshold % for file starting with a CR or dying fan
                fan_restart(1) = fan_onB(1);
                fan_restart(2:numel(ind_B)+1) = fan_onB(ind_B+1); % all the data right after a CR or dying fan
            else
                fan_restart = fan_onB(ind_B+1);
            end

            for i = 1:length(fan_restart)
                if fan_restart(i) <= 130
                    start = 1;
                    stop  = fan_restart(i)+10;
                else
                    start = fan_restart(i) - 130;
                    stop  = fan_restart(i) + 9;
                end
                if stop > size(data_fileB,2)
                    stop = size(data_fileB,2);
                end
                fan_offC = start:stop; fan_offC = fan_offC';
                fan_off = [fan_off; fan_offC]; clear fan_offC
            end
            fan_off = [fan_off; fan_offB]; % adding random fan_off data;
            fan_off = unique(fan_off); % remove repeating rows
            [data_file, ~] = removerows(data_fileB,fan_off);

            percent_off = (size(fan_off,1)/size(data_fileB,1))*100;
            if percent_off > 25
                warning( 'Fan off for > 25 % of file')
                FileTracker.FanOff25Perc = 1;
            end
            clear fan_off i fan_onB percent_off
        end

    end

end
clear data_fileB  PS fan_off start stop fan_offB
end



%% functions call in this function
function [alldata, crdata] = ReadWithTitile(filename, data_title, cr_title, IN_model, CR_system)

% reading the csv file for the first time
opts = detectImportOptions(filename);
opts.VariableNamesLine = 1; % sometimes readtable fails to find read the variable names, likely because of specieal sympols (# etc)
opts.VariableNamingRule = 'preserve'; 
A = readtable(filename, opts);
headers = A.Properties.VariableNames;

% make sure that most data types are read as double 
vartypes = opts.VariableTypes;
indt = 1:length(headers);
inds = find(contains(headers, 'SELECT', 'IgnoreCase', true));
ind = find(contains(headers, {'date_time', 'date/time', 'UTC_Time','date'}, 'IgnoreCase', true));
indt([1 ind inds]) = []; % also column 1 should be model ID.
vartypes(indt) = {'double'};

% now reading the csv file again
opts.VariableNames = headers;
opts.VariableTypes = vartypes;
A = readtable(filename, opts);

%% catch date + time info
if length(ind) ==1
    tdatevec = datevec(table2array(A(:,ind)));
    if tdatevec(:, 1) < 2000
        tdatevec(:, 1) = tdatevec(:, 1) + 2000;
    end
    tdatenum = datenum(tdatevec(:,1),tdatevec(:,2),tdatevec(:,3),...
        zeros(length(tdatevec),1),zeros(length(tdatevec),1),zeros(length(tdatevec),1));
    Hour = tdatevec(:,4);

elseif isempty(ind) % not Date_time column. There should be other columns containting year,month, date, and hour info
    datetitles = {'year','month','date','hour'};
    for ii = 1:length(datetitles)
        ind = find(contains(headers,datetitles{ii},'IgnoreCase',true));
        if length(ind) == 1
            tdatevec(:,ii) = table2array(A(:,ind));
        else
            error('%s data not found! Model type = %d',datetitles{ii},IN_model)
        end
    end
    tdatenum = datenum(tdatevec(:,1),tdatevec(:,2),tdatevec(:,3),...
        zeros(length(tdatevec),1),zeros(length(tdatevec),1),zeros(length(tdatevec),1));
    Hour = tdatevec(:,4);
else
    ind = find(contains(headers, {'date_time', 'date/time', 'UTC_Time'}, 'IgnoreCase', true));
    if length(ind) == 1
        tdatevec = datevec(table2array(A(:, ind)));
        if tdatevec(:, 1) < 2000
            tdatevec(:, 1) = tdatevec(:, 1) + 2000;
        end
        tdatenum = datenum(tdatevec(:, 1), tdatevec(:, 2), tdatevec(:, 3), ...
            zeros(length(tdatevec), 1), zeros(length(tdatevec), 1), zeros(length(tdatevec), 1));
        Hour = tdatevec(:, 4);
    else
        disp(headers(ind))
        error('More than 1 column of date time infomation from the input. Titles are shown above.')
    end
end

alldata(:,1) = tdatenum;
alldata(:,2) = Hour;

clear tdatevec Hour

%% catch other info specified in the data_title
if IN_model == 1
    jj = 3:length(data_title)-3; % no flow for model 1
else
    jj = 3:length(data_title); % start from 3 since column 1 is 'datenum'; 2 is Hour; both obtained already. 
end
for ii = jj 

    ind = find(matches(headers,data_title{ii},'IgnoreCase',true));

    if length(ind)==1
        alldata(:,ii) = table2array(A(:,ind));

    elseif isempty(ind) % If there aren't exact match, find the column contains the target title

        ind = find(contains(headers,data_title{ii},'IgnoreCase',true));

        if length(ind) == 1
            alldata(:,ii) = table2array(A(:,ind));

        elseif length(ind) > 1
            % list solutions here
            disp(headers(ind))
            error('More than 1 column found for %s. Titles are shown above.',data_title{ii})

        else 

            solution_exist = 0; % will turn to 1 if there is a relevant solution
            
            % ----------- list solutions here ---------------------------------------------------
            % Headers use space instead of '_'
            ttitle = data_title{ii};
            ind = find(ttitle == '_');
            ttitle(ind) = ' '; % remove '_' which sometimes is a space in the raw data
            ind = find(contains(headers, ttitle, 'IgnoreCase', true)); % expect to at least one
            if length(ind) == 1
                alldata(:, ii) = table2array(A(:, ind));
                solution_exist = 1;
            end

            switch data_title{ii}
                % Amb_pressure and flow_pressure not found
                case 'Amb_Pressure' 
                    solution_exist = 1;
                    % if reach here, temp was found and there was only one
                    % temp column. Ambient pressure usually follow right
                    % after ambient temperature
                    ind1 = find(contains(headers,'pressure','IgnoreCase',true)); % expect to at least one
                    ind2 = find(contains(headers,'temp','IgnoreCase',true)); % expect to at least one
                    ind = find(ismember(ind1, ind2 + 1)); % expect pressure comes right after temp
                    if length(ind)==1
                        alldata(:,ii) = table2array(A(:,ind1(ind)));
                    elseif length(ind) == 2 % there are 'temp' and 'temp_raw', delete the ind leading to temp_raw
                        ind3 = find(contains(headers,'temp_raw','IgnoreCase',true)); % expect to have only one match
                        if length(ind3) == 1
                            ind2(ind2 == ind3) =[];
                            ind = find(ismember(ind1, ind2 + 1)); % expect pressure comes right after temp
                            alldata(:, ii) = table2array(A(:, ind1(ind)));
                        else % there isn't temp_raw, not sure what is going no. need manual check
                            error('Ambient pressure column cannot be identified.')
                        end
                    else
                        error('Ambient pressure column cannot be identified.')
                    end
                case 'Flow_Pressure'
                    solution_exist = 1;
                    % if reach here, Fan_RPM was found and there was only one
                    % Fan_RPM column. flow pressure usually follow right
                    % after Fan_RPM.
                    ind1 = find(contains(headers,'pressure','IgnoreCase',true)); % expect to at least one
                    ind2 = find(contains(headers, {'FAN_RPM','FAN RPM'}, 'IgnoreCase', true)); % expect to at least one
                    ind = find(ismember(ind1, ind2 + 1)); % expect pressure comes right after Fan RPM
                    if length(ind)==1
                        alldata(:,ii) = table2array(A(:,ind1(ind)));
                    else % sometimes there isn't flow pressure, print a note
                        fprintf('Flow pressure column not found.')
                    end
            end 

            % for cases not listed above, give an error message 
            if solution_exist == 0
                error('No column found for %s',data_title{ii})
            end
            clear solution_exist
        end

        % unlikely to reach here - this means to columns in the input file have identical title
    else
        disp(headers(ind))
        error('More than 1 column match for %s. Titles are shown above.',data_title{ii})
    end
end

%% output CR data with date (don't need hour info)
if CR_system == 1
    for ii = 1:length(cr_title)
        ind = find(contains(headers,cr_title{ii},'IgnoreCase',true));

        if length(ind)==1
            crdata(:,ii) = table2array(A(:,ind));

        elseif isempty(ind)
            error('No column found for %s',cr_title{ii})
        else
            disp(ind)
            error('More than 1 column match for %s. Titles are shown above.',cr_title{ii})
        end
    end
    % keep only rows of cr date with values
    ind = find(~isnan(crdata(:,1)));
    crdata = [tdatenum(ind,:),crdata(ind,:)];
else
    crdata = [];
end

%% for files provide the 'SELECT' column, use it to select valid data (similar to fan off removal but more robust!)
if ~isempty(inds)
    select = table2array(A(:,inds));
    valid_data = find(contains(select, 'DATA', 'IgnoreCase', true));
    alldata = alldata(valid_data,:); % only deliver the valid data 
end

end




