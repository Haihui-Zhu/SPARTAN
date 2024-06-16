%% %%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: Examine the uncertainty in air volume in the filter sampling process 
% Haihui Zhu
% January 23, 2024

close all; clear all; clc
addpath('/storage1/fs1/rvmartin/Active/haihuizhu/1.code')
warning('off')

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_Compile_data = 1; % if 0 , do not clear up the work space
S_TimeSeries = 1; 
% Set directories
% direc = '\\storage1.ris.wustl.edu\rvmartin\Active\SPARTAN-shared\';
direc_in = '/storage1/fs1/rvmartin/Active/SPARTAN-shared/';
direc_out = '/storage1/fs1/rvmartin/Active/haihuizhu/6.SPARTAN/';

direc_sampling = strcat(direc_in,'Site_Sampling');
direc_master = strcat(direc_in,'Analysis_Data/Master_files');
direc_met = '/storage1/fs1/rvmartin/Active/GEOS-Chem-shared/ExtData/GEOS_0.5x0.625/MERRA2/';

% The diary function saves the processing history into an annual record
diary(sprintf('%sPublic_Data/Data_Processing_Records/%s_volume_calculation.txt',direc_in,datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))

%-------------   SITE INFO   --------------
% Sites are listed in the file "Site_details.xlsx" in alphabetical order
% Column info:
site_details = readtable(strcat(direc_sampling,'/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));
latitudes = table2array(site_details(:,5));
longitudes = table2array(site_details(:,6));

target = find(ismember(Site_codes,'USPA'));


% met data 
metfname = sprintf('%s/MERRA-2_Surface_T_P_at_SPARTAN_2019-2023.mat',direc_sampling);
load(metfname) % 'P','T','alldates_vec','alldates_num','Site_codes'
            
T1 = 22.5+273; % K
% P1 = 101325; % Pa % use USBO ambient pressure
usbo = find(ismember(Site_codes,'USBO'));


for loc = target %1:length(Site_codes)
if S_Compile_data == 1
    % Read the flowrate_dates file
    % this file contains data that we want to examine
    flowrate_fname = sprintf('%s/%s_%s/Cartridge_Data/%s_dates_flows.xlsx',direc_sampling,Site_codes{loc},Site_cities{loc},Site_codes{loc});
    flowrate_raw = readtable(flowrate_fname);
    filterid = flowrate_raw.Filter_ID; % with _N filter location
    analysisid = flowrate_raw.Analysis_ID; % no filter location - match master file
    startdate = [flowrate_raw.start_year flowrate_raw.start_month flowrate_raw.start_day flowrate_raw.start_hour];
    stopdate = [flowrate_raw.stop_year flowrate_raw.stop_month flowrate_raw.stop_day flowrate_raw.stop_hour];
    intstart = flowrate_raw.Flow_internal_start;
    intstop = flowrate_raw.Flow_internal_end;
    samphr = flowrate_raw.hours_sampled;

    volume1 = flowrate_raw.volume_m3; % current volume
    clear flowrate_raw

    % volume based on other method
    volume2 = nan.*volume1; % prepare for external to internal 
    volume3 = nan.*volume1; % prepare for integration
    volume4 = nan.*volume1; % prepare for calibration uncertainty

    % Read master file to connect filter ID to cartridge ID
    masterfname = sprintf('%s/%s_master.csv',direc_master,Site_codes{loc});
    masterraw = readtable(masterfname);
    m_cartridgeid = masterraw.CartridgeID;
    m_filterid = masterraw.FilterID;

    cartridgeid = cell(size(filterid));
    for ii = 1:length(m_filterid)
        ind = find(ismember(analysisid,m_filterid(ii)));
        if ~isempty(ind)
            cartridgeid(ind) = m_cartridgeid(ii);
        end
    end
    for ii = 1:length(cartridgeid)
        if isa(cartridgeid{ii},'double') 
            cartridgeid{ii} = '';
        end
    end
    unique_cartr = unique(m_cartridgeid);
    clear m_filterid m_cartridgeid masterraw masterfname

    % now read the instantanuous flow rate
%     cartridgeid = cartridgeid(~cellfun('isempty',cartridgeid))  ;
    

    for cid = 1:length(unique_cartr)

        tcartr = unique_cartr{cid};
        tcartr = tcartr(6:8); % read the last three characters 

        for flid = [1:6 8]
            % find the location of this filter 
            ind = find(contains(cartridgeid,unique_cartr{cid}) & contains(filterid,sprintf('-%d',flid)));

            if ~isnan(startdate(ind,1))
                % read the raw flow rate file
                rawfdir = dir(sprintf('%s/%s_%s/Cartridge_Data/%s/FC0*-%d.csv',direc_sampling,Site_codes{loc},Site_cities{loc},...
                            unique_cartr{cid},flid)); %%% need to find an easier way
                if ~isempty(rawfdir)
                    rawfname = sprintf('%s/%s',rawfdir.folder, rawfdir.name);
                    rawdata = import_raw_flowrate(rawfname);
                    samplingtime = [rawdata.YEAR rawdata.MONTH rawdata.DATE rawdata.HOUR rawdata.MINUTE rawdata.SECOND];
                    if flid == 8 % PM10
                        flowrate_raw = rawdata.FLOW_B_LPM; % unit = Litter/minute;
                    else
                        flowrate_raw = rawdata.FLOW_A_LPM; % unit = Litter/minute;
                    end


                    % cut non-sampling period
                    ind1 = find(samplingtime(:,1) == startdate(ind,1) & samplingtime(:,2) == startdate(ind,2) & ...
                        samplingtime(:,3) == startdate(ind,3) & samplingtime(:,4) == startdate(ind,4) );
                    ind2 = find(samplingtime(:,1) == stopdate(ind,1) & samplingtime(:,2) == stopdate(ind,2) & ...
                        samplingtime(:,3) == stopdate(ind,3) & samplingtime(:,4) == stopdate(ind,4)-1 ); % -1 because in flow_dates, the end hour is the hour after it actually ends (8:59:59 labeled as 9:00:00)
                    if isempty(ind2)
                        ind2 = find(samplingtime(:,1) == stopdate(ind,1) & samplingtime(:,2) == stopdate(ind,2) & ...
                            samplingtime(:,3) == stopdate(ind,3) & samplingtime(:,4) == stopdate(ind,4) ); % -1 because in flow_dates, the end hour is the hour after it actually ends (8:59:59 labeled as 9:00:00)
                        if isempty(ind2) % still not found, skip
                            ind2 = find(samplingtime(:,1) == stopdate(ind,1) & samplingtime(:,2) == stopdate(ind,2) & ...
                                samplingtime(:,3) == stopdate(ind,3) & samplingtime(:,4) == ceil(stopdate(ind,4)) ); % using 'ceil' beause BDDU label end hour as 5.5 when start time is 6. 'ceil' might sample more time than it actually run but the flow rate will be nearly 0. If use 'floor', likely introduce too much underestimation.
                            if isempty(ind2)
                                fprintf('%s flow sampling stop time not found\n',filterid{ind})
                                continue
                            end
                        end
                    end
                    ind3 = ind1(1):ind2(end); % it works for USPA but not necessarily others
                    samplingtime = samplingtime(ind3,:);
                    flowrate_raw = flowrate_raw(ind3,:);
                    ind4 = isnan(samplingtime(:,1)); % sometimes there are pauses in a file.  (9 day sampling method?)
                    samplingtime(ind4,:) = [];
                    flowrate_raw(ind4,:) = [];
                    ind4 = isnan(flowrate_raw<0.5); % not sampling
                    samplingtime(ind4,:) = [];
                    flowrate_raw(ind4,:) = [];

                    % now read the ambient T P
                    T2 = nan(size(flowrate_raw));
                    P2 = nan(size(flowrate_raw));
                    P1 = nan(size(flowrate_raw)); % USBO pressure
                    % The hard way:
                    %{
                oldfname = '';
                for d = 1:length(samplingtime)
                    metfname = sprintf('%s/%d/%.2d/MERRA2.%d%.2d%.2d.I3.05x0625.nc4',direc_met,samplingtime(d,1),samplingtime(d,2),...
                        samplingtime(d,1),samplingtime(d,2),samplingtime(d,3));

                    if matches(metfname,oldfname) % merra-2 T P exist
                        % find applicable rows
                        T2(d) = Ti(samplingtime(d,4)+1);
                        P2(d) = Pi(samplingtime(d,4)+1);

                    else
                        % reading new file for T P
                        fprintf('reading %s\n',metfname)
                        Tm = ncread(metfname,'T');  % K
                        Tm = squeeze(Tm(:,:,1,:)); % only keep the surface temp
                        Pm = ncread(metfname,'PS'); % Pa
                        lat = ncread(metfname,'lat'); lon = ncread(metfname,'lon');
                        % interp to get T & P at site
                        Ti = nan(24,1); Pi = nan(24,1);
                        for hr = 1:8
                            Ti( (hr-1)*3+1:(hr-1)*3+3 )= interp2(lat',lon,Tm(:,:,hr),latitudes(loc), longitudes(loc));
                            Pi( (hr-1)*3+1:(hr-1)*3+3 ) = interp2(lat',lon,Pm(:,:,hr),latitudes(loc), longitudes(loc));
                        end
                        % find applicable rows
                        T2(d) = Ti(samplingtime(d,4)+1);
                        P2(d) = Pi(samplingtime(d,4)+1);
                    end
                    oldfname = metfname;
                end
                    %}

                    % The easy way:
                    for d = 1:length(samplingtime)
                        dateind = find(alldates_vec(:,1)==samplingtime(d,1) & alldates_vec(:,2)==samplingtime(d,2) & alldates_vec(:,3)==samplingtime(d,3));
                        if isempty(dateind) % for years before 2019, use 2020 met data (leap year)
                            dateind = find(alldates_vec(:,1)==2020 & alldates_vec(:,2)==samplingtime(d,2) & alldates_vec(:,3)==samplingtime(d,3));
                        end
                        T2(d) = T(dateind,samplingtime(d,4)+1,loc);
                        P2(d) = P(dateind,samplingtime(d,4)+1,loc);
                        P1(d) = P(dateind,samplingtime(d,4)+1,usbo);
                    end
                    %}

                    % convert sampling time to time bin for integration
                    timebin = diff(datenum(samplingtime))*24*60; % unit = minute
                    timebin(timebin>1) = nan;

                    % now calculate ambient flow rate
                    flowrate_new = flowrate_raw./T1./P2.*T2.*P1; % for V4

                    volume2(ind) = (intstart(ind) + intstop(ind))/2.* samphr(ind)*60 /1000; % m3
                    volume3(ind) = sum((flowrate_raw(2:end) + flowrate_raw(1:end-1))/2.*timebin,'omitnan') /1000; % m3
                    volume4(ind) = sum((flowrate_new(2:end) + flowrate_new(1:end-1))/2.*timebin,'omitnan') /1000;  % need to convert unit
                end
            else % if startdate(ind,:) are nans, plot out the instantaneous flow rate for reference
               % read the raw flow rate file
                rawfdir = dir(sprintf('%s/%s_%s/Cartridge_Data/%s/FC0*-%d.csv',direc_sampling,Site_codes{loc},Site_cities{loc},...
                            unique_cartr{cid},flid)); %%% need to find an easier way
                if ~isempty(rawfdir)
                    rawfname = sprintf('%s/%s',rawfdir.folder, rawfdir.name);
                    rawdata = import_raw_flowrate(rawfname);
                    samplingtime = [rawdata.YEAR rawdata.MONTH rawdata.DATE rawdata.HOUR rawdata.MINUTE rawdata.SECOND];
                    flowrate_raw = rawdata.FLOW_A_LPM; % unit = Litter/minute;

                    figure('Position',[10 10 1200 900])
                    X = datetime(samplingtime(:,1:3));
                    fontsize = 12;
                    plot(X,flowrate_raw)

                    sfname = sprintf('%s/Temp/Instantaneous_flow_%s_%s_%d.png',direc_out,Site_codes{loc},tcartr,flid);
                    saveas(gcf,sfname)
                    fprintf('%s saved\n',sfname)
                end
            end
        end
    end
end
%%
if S_TimeSeries == 1
    % PM10 are reported in another column. skip for now
    volume1b = volume1;% once only!!! And comment out 'clear'!!
    volume2b = volume2;
    volume3b = volume3;
    volume4b = volume4;
    ind = find(volume1b<4 | volume3b>10 | volume4b<4 | isnan(volume4b));
    volume1=volume1b; volume1(ind) = nan; 
    volume2=volume2b; volume2(ind) = nan;
    volume3=volume3b; volume3(ind) = nan;
    volume4=volume4b; volume4(ind) = nan;

    volume1(volume1==0) = NaN;
    volume2(volume2==0) = NaN;


    figure('Position',[10 10 1200 900])
    X = datetime(startdate(:,1:3));
    subplot(2,2,1)
    subtimeseries(X,volume1,volume2,'volume1','volume2')
    subplot(2,2,2)
    subtimeseries(X,volume2,volume3,'volume2','volume3')
    subplot(2,2,3)
    subtimeseries(X,volume3,volume4,'volume3','volume4')
    subplot(2,2,4)
    subtimeseries(X,volume1,volume4,'volume1','volume4')

    sfname = sprintf('%s/Temp/Volume_Calculation_%s.png',direc_out,Site_codes{loc});
    saveas(gcf,sfname)
    fprintf('%s saved\n',sfname)
end
end

%% FUNCTION
function subtimeseries(X,volume1,volume2,name1, name2)
    
    fontsize = 16;
    plot(X,volume1,X,volume2)
    ylabel(sprintf('volume (m^3)'))
    [ta1,m,b,~] = getStatics(volume1,volume2);

    m1 = mean(volume1,'omitnan');
    m2 = mean(volume2,'omitnan');
    ratio = m2/m1;
    
    ta1 = sprintf('%s \nmean1 = %4.2f \nmean2 = %4.2f\nratio = %4.3f',ta1, m1, m2,ratio);

    text(0.05,0.05,ta1, 'Units','normalized', ...
        'horizontalalignment','left','verticalalignment','bottom','color','k','fontsize',fontsize);
    legend({name1,name2},'location','northeast','fontsize',fontsize)
    set(gca,'fontsize',fontsize)

end


function outdata = import_raw_flowrate(fname)

    %% Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 25);

    % Specify range and delimiter
    opts.DataLines = [9, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["SN", "DATE_TIME", "YEAR", "MONTH", "DATE", "HOUR", "MINUTE", "SECOND", "FLOW_A_RAW", "FLOW_A_LPM", "FLOW_B_RAW", "FLOW_B_LPM", "TEMP_RAW", "TEMP_C", "VACUUM_A_RAW", "VACUUM_A_KPA", "VACUUM_B_RAW", "VACUUM_B_KPA", "ADC_1_RAW", "ADC_1_UNIT", "ADC_2_RAW", "ADC_2_UNIT", "ADC_3_RAW", "ADC_3_UNIT", "STATUS"];
    opts.VariableTypes = ["double", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, "STATUS", "EmptyFieldRule", "auto");
    opts = setvaropts(opts, "DATE_TIME", "InputFormat", "yyyy-MM-dd'T'HH:mm:ss");
    opts = setvaropts(opts, "SN", "TrimNonNumeric", true);
    opts = setvaropts(opts, "SN", "ThousandsSeparator", ",");

    % Import the data
    outdata = readtable(fname, opts);

    year = outdata.YEAR; 

if isnan(year(1)) % wrong format
    clear opts

    opts = delimitedTextImportOptions("NumVariables", 20);

    % Specify range and delimiter
    opts.DataLines = [8, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["SN", "DATES", "TIME", "FLOW_A_RAW", "FLOW_A_LPM", "FLOW_B_RAW", "FLOW_B_LPM", "TEMP_RAW", "TEMP_C", "VACUUM_A_RAW", "VACUUM_A_KPA", "VACUUM_B_RAW", "VACUUM_B_KPA", "ADC_1_RAW", "ADC_1_UNIT", "ADC_2_RAW", "ADC_2_UNIT", "ADC_3_RAW", "ADC_3_UNIT", "STATUS"];
    opts.VariableTypes = ["char", "datetime", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, "STATUS", "EmptyFieldRule", "auto");
    opts = setvaropts(opts, "DATES", "InputFormat", "yyyy-MM-dd");
    opts = setvaropts(opts, "TIME", "InputFormat", "HH:mm:ss");
%     opts = setvaropts(opts, "SN", "TrimNonNumeric", true);
%     opts = setvaropts(opts, "SN", "ThousandsSeparator", ",");

    % Import the data
    outdata = readtable(fname, opts);

    % adding YEAR, MONTH, DATE, HOUR, MINUTE, SECOND
    a = datevec(outdata.DATES);
    outdata.YEAR = a(:,1);
    outdata.MONTH = a(:,2);
    outdata.DATE = a(:,3);
    a = datevec(outdata.TIME);
    outdata.HOUR = a(:,4);
    outdata.MINUTE = a(:,5);
    outdata.SECOND = a(:,6);

end

end