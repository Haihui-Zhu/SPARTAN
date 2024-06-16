function [Ref_Dates,Ref_Values,FilterGroup] = Load_MDL_Unc_XRF(direc)
% Purpose: group filters according to the XRF running date so that to
% find the corresponded MDL and UNC

MDLfname = sprintf('%s/Analysis_Data/XRF/MDL_Unc_Ref.mat',direc);
load(MDLfname,'Ref_Dates','Ref_Values','FilterGroup')
% today = datetime('today');

% comment out the following since a new section is added to script E to
% update MDL_Unc_Ref.mat
%{
if today-date_mark > 168 % if the MDL_Unc_Ref is out dated for 30 days, update it. 
   
    % make a new MDL ref file
    fprintf('Updating XRF MDL\n')
    XRF_MDL_UNC  = readtable(strcat(direc,'/Analysis_Data/XRF/SPARTAN_XRF_MDL_Unc.xlsx'),'PreserveVariableNames',true);
    Ref_Dates = floor(convertTo(XRF_MDL_UNC.Date,'datenum'));
    Ref_Values = XRF_MDL_UNC(:,2:end);
%     Ref_AnMDL = table2array(XRF_MDL_UNC(:,3)); % analytical MDL
%     Ref_MDL = table2array(XRF_MDL_UNC(:,4));% MDL
%     Ref_UncA = table2array(XRF_MDL_UNC(:,5));% Unc_Additive
%     Ref_UncP = table2array(XRF_MDL_UNC(:,6));% Unc_proportion_mass
%     Ref_Vol = table2array(XRF_MDL_UNC(1,7));
%     Unc_Prop_Vol = table2array(XRF_MDL_UNC(1,8));
    
    direc_sampling = strcat(direc,'/Site_Sampling');
    site_details = readtable(strcat(direc_sampling,'/Site_details.xlsx'),'PreserveVariableNames',true);
    Site_codes = table2array(site_details(:,1));

    uniRefDates = unique(Ref_Dates); uniRefDates(end+1) = 999999;% the last date set to 999999
    for i = 1:numel(uniRefDates)
        FilterGroup{1,i} = uniRefDates(i);
        FilterGroup{2,i} = {};
    end

    % read from Archive
    fprintf('Processing Archived XRF Files\n');
    direc_archive = strcat(direc,'/Analysis_Data/Archived_Filter_Data/XRF_data/');
    for loc = 1:length(Site_codes)
        file_destination = strcat(direc_archive,sprintf('%s/',Site_codes{loc}));
        arc_fname = getFiles(file_destination);
        for fID = 1:numel(arc_fname)
            % Do not read hidden files
            char_file = char(arc_fname(fID));
            if char_file(1)=="."
                continue
            elseif char_file(1)=="~"
                continue
            end
            clear char_file

            filename = sprintf('%s/%s',file_destination,arc_fname{fID,1});
            initialdata = readtable(filename,'PreserveVariableNames',true);
            XRFtitle = initialdata.Properties.VariableNames;
            SampleDate = floor(convertTo(table2array(initialdata(:,contains(XRFtitle,'Time'))),'datenum'));
            FilterID = table2array(initialdata(:,contains(XRFtitle,'Ident')));

            for i = 1:numel(uniRefDates)-1
                if i == 1 % Everything before the second Reference MDL goes to the first Reference
                    Ind = find(SampleDate<uniRefDates(i+1));
                else
                    Ind = find(SampleDate<uniRefDates(i+1) & SampleDate>=uniRefDates(i));
                end
                % exclude rows don't contain filterID
                Ind(contains(FilterID(Ind,1),'?')) = [];

                rows = size(FilterGroup{2,i},1)+1:size(FilterGroup{2,i},1)+numel(Ind);
                FilterGroup{2,i}(rows,1) = FilterID(Ind,1);
            end
        end
    end

    % read from New
    direc_new = strcat(direc,'/Analysis_Data/XRF/NEW');
    files = getFiles(direc_new);
    fprintf('Processing New XRF Files\n');

    for fID = 1:numel(files)
        % Do not read hidden files
        char_file = char(files(fID));
        if char_file(1)=="."
            continue
        elseif char_file(1)=="~"
            continue
        end
        clear char_file

        filename = sprintf('%s/%s',direc_new,files{fID,1});
        initialdata = readtable(filename,'PreserveVariableNames',true);
        XRFtitle = initialdata.Properties.VariableNames;
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

            rows = size(FilterGroup{2,i},1)+1:size(FilterGroup{2,i},1)+numel(Ind);
            FilterGroup{2,i}(rows,1) = FilterID(Ind,1);
        end
    end

%     % remove previous files
%     previousfiles = dir(sprintf('%s/Analysis_Data/XRF/MDLRe*.mat',direc));
%     for pf = 1:length(previousfiles)
%         thisfile = sprintf('%s/Analysis_Data/XRF/%s',direc,previousfiles{pf}.name);
%         delete(thisfile)
%     end
    date_mark = today;
    save(MDLfname,'FilterGroup','Ref_Dates','Ref_Values','date_mark')

    
end

%}


% function files = getFiles(directory)
% % retrieve all files in specified directory
% names = dir(sprintf('%s/*.xlsx',directory));
% files = {names.name}; files = files';
% end

end
