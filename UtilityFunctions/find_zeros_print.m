function find_zeros_print(RFM_labels,RFM_dates,RFM_para_name,RFM_value,RFM_descrip)

ind = find(RFM_value == 0 & RFM_dates(:,1)>2019); % Historical data should be omitted

if ~isempty(ind)
    % Append new lines to the file
    fileID = fopen('Zero_values.txt', 'a');
    if fileID == -1
        error('File could not be opened for appending');
    end
    for h = 1:length(ind)
        fprintf(fileID, '\n%s %4d %4d %4d %4d %4d %4d %4d %4d %s %4.2d %s',...
            RFM_labels{ind(h),:}, RFM_dates(ind(h),:), RFM_para_name{ind(h),:}, RFM_value(ind(h),:), RFM_descrip{ind(h),:});
    end
    fclose(fileID);
end

end