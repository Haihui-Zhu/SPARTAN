function printTitles(fileID, titles)
    % write specified titles in specified fileID
    rep_string = repmat('%s,',1,length(titles));
    char_string = rep_string(1:end-1);
    print_string = strcat(char_string, '\n');
    fprintf(fileID,print_string,titles{1,1:length(titles)});
end  