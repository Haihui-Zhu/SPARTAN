
% this script will make a copy of current master files to the backup folder
clear; close all; clc
addpath('UtilityFunctions')

% Setup directories 
debug_mode = 0;
SharedDir = find_root_dir(debug_mode);

% destination
DESTINATION = sprintf('%s/Analysis_Data/Master_files/Backups/Backup_%s',SharedDir,datestr(today,'yyyy-mm-dd'));
mkdir(DESTINATION)
% Source
direc_master = sprintf('%s/Analysis_Data/Master_files',SharedDir);

% additional info
site_details = readtable(strcat(SharedDir,'/Site_Sampling/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));

% start backing up
for loc = 1:length(Site_codes)

    MasterFile = sprintf('%s/%s_master.csv',direc_master,Site_codes{loc});
    [SUCCESS,MESSAGE,MESSAGEID] = copyfile(MasterFile, DESTINATION);
    if SUCCESS ~= 1
        disp(MESSAGE)
    end

end
fprintf('Backup completed. \n')