function [] = WriteToMaster( Titles,Master_IDs, Master_Barcodes, Master_CartridgeIDs, Master_LotIDs, Master_ProjectIDs,...
                             Master_hours, Master_masstype, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
                             Master_carbon, Master_Method, Master_flags,...
                             direc_master,site_ID)
% This function will write all the input data into the master file 
% Inputs:
% direc_master should be the dir for the master file
% Site ID shoule be the 4 letter Site Code.
% Other inputs (master data) should be same as listed in the ReadMaster
% function

% Written by: Haihui Zhu
% Created: 2021-08-25   


% If the script reaches here, the test write is successful. Can start to write to master files:

    fileID = fopen(sprintf('%s/%s_master.csv',direc_master,site_ID),'w');
    if fileID == -1
        error('Cannot open master file: %s',site_ID)
    end
    % Print Titles
    printTitles(fileID, Titles)

    % Print Data
    for h=1:size(Master_IDs,1) % 7 cols                 % 8 cols               % 2 cols                     % 13 cols                                                   % 21 cols                                                                              % 26 cols                              % 6 cols
        fprintf(fileID,'%s,%s,%s,%s,%s,%f,%6.1f,   %f,%f,%f,%f,%f,%f,%f,%f,   %6.2f,%f,  %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,   %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,    %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,%f,%f,%f, %f,%f,%f,%f,%f,%s\n',...
            char(Master_IDs(h,:)), char(Master_Barcodes(h,:)), char(Master_CartridgeIDs(h,:)), char(Master_LotIDs(h,:)), char(Master_ProjectIDs(h,:)),...
            Master_hours(h),       Master_masstype(h,:),       Master_dates(h,:),        Master_mass(h,:),         Master_IC(h,:), ...
            Master_ICP(h,:),       Master_XRF(h,:),            Master_carbon(h,:),       Master_Method(h,:), char(Master_flags(h,:)));
    end

    fclose(fileID);
    fprintf('Finished writing to %s master file \n\n', site_ID)

end
