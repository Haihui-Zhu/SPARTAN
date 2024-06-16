function [Ref_Dates,Ref_Element,Ref_AnMDL,Ref_MDL,Ref_Unc,MDLRef] = SummarizeMDLRef(direc)
% Purpose: group filters according to the XRF running date so that to
% find the corresponded MDL and UNC

MDLfname = sprintf('%s/Analysis_Data/XRF/MDLRef_%s.mat',direc,datestr(today,'yyyy-mm'));

if exist(MDLfname,'file')
    load(MDLfname)
else
    XRF_MDL_UNC  = readtable(strcat(direc,'/Analysis_Data/XRF/SPARTAN_XRF_MDL_Unc.xlsx'),'PreserveVariableNames',true);
    Ref_Dates = floor(convertTo(XRF_MDL_UNC.Date,'datenum'));
    Ref_Element = XRF_MDL_UNC.Element;
    Ref_AnMDL = table2array(XRF_MDL_UNC(:,3));
    Ref_MDL = table2array(XRF_MDL_UNC(:,4));
    Ref_Unc = table2array(XRF_MDL_UNC(:,5));

    direc_sampling = strcat(direc,'/Site_Sampling');
    site_details = readtable(strcat(direc_sampling,'/Site_details.xlsx'),'PreserveVariableNames',true);
    Site_codes = table2array(site_details(:,1));

    uniRefDates = unique(Ref_Dates); uniRefDates(end+1) = 999999;% the last date set to 999999
    for i = 1:numel(uniRefDates)
        MDLRef{1,i} = uniRefDates(i);
        MDLRef{2,i} = {};
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
                if i == 1 % Everything before the second Reference MDL refers to the first Reference
                    Ind = find(SampleDate<uniRefDates(i+1));
                else
                    Ind = find(SampleDate<uniRefDates(i+1) & SampleDate>=uniRefDates(i));
                end
                % exclude rows don't contain filterID
                Ind(contains(FilterID(Ind,1),'?')) = [];

                rows = size(MDLRef{2,i},1)+1:size(MDLRef{2,i},1)+numel(Ind);
                MDLRef{2,i}(rows,1) = FilterID(Ind,1);
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

            rows = size(MDLRef{2,i},1)+1:size(MDLRef{2,i},1)+numel(Ind);
            MDLRef{2,i}(rows,1) = FilterID(Ind,1);
        end
    end
    previousfiles = dir(sprintf('%s/Analysis_Data/XRF/MDLRe*.mat',direc));
    previousfiles = {previousfiles.name}; 
    previousfiles = previousfiles';
    for pf = 1:length(previousfiles)
        thisfile = sprintf('%s/Analysis_Data/XRF/%s',direc,previousfiles{pf});
        delete(thisfile)
    end
    
    save(MDLfname,'MDLRef','Ref_Dates','Ref_Element','Ref_AnMDL','Ref_MDL','Ref_Unc')
end


    function files = getFiles(directory)
    % retrieve all files in specified directory
    names = dir(sprintf('%s/*.xlsx',directory));
    files = {names.name}; files = files';
    end

end