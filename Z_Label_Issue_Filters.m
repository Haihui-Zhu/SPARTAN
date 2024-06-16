
% Label issue filter in Master files

% RootDir= '/Volumes/rvmartin/Active/haihuizhu/6.SPARTAN/';
RootDir= '/Volumes/rvmartin/Active/SPARTAN-shared/';

Infile = strcat(RootDir,'/Analysis_Data/XRF/IssueSampleList.xlsx');
direc_master = strcat(RootDir,'/Analysis_Data/Master_files');
direc_sampling = strcat(RootDir,'/Site_Sampling');


% the diary function saves the processing history into an annual record
diary(sprintf('%s/Public_Data/Data_Processing_Records/%s_Label_issue_filters.txt',RootDir, datestr(today,'yyyy-mm')))
fprintf('%s \n', datestr(today))


InTable = readtable(Infile);
SampleID = InTable.Sample;
Issue = InTable.Issue;
Notes = InTable.Notes;

for i = 1:length(SampleID)
    
    SiteCode = SampleID{i}(1:4);
    MasterFile = sprintf('%s/%s_master.csv',direc_master,SiteCode);
    
    [Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs, projectIDs_master,Master_hours, Master_masstype, ...
     Master_dates, Master_mass, Master_IC, Master_ICP, Master_XRF,...
     Master_carbon, Master_Method, Master_flags] = ReadMaster( MasterFile, SiteCode );
    
    if isempty(Master_mass)
        continue
    end
    
    Ind = ismember(Master_IDs,SampleID{i});
    if ismember(Notes{i},'Need to be excluded')
        NewFlag = sprintf('%s:%s',Notes{i},Issue{i});
    else 
        NewFlag = Notes{i};
    end

    Master_flags{Ind} = AddFlag(Master_flags{Ind}, NewFlag);
    fprintf('Adding %s to %s flags\n',NewFlag,SampleID{i})

    % finish this sample
    WriteToMaster( Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs, projectIDs_master,...
        Master_hours, Master_masstype, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
        Master_carbon, Master_Method, Master_flags,...
        direc_master, SiteCode)

end

diary off
