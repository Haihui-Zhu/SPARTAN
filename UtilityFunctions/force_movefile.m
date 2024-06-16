function [status,msg,msgid] = force_movefile(filename, file_destination)
    [status,msg,msgid] = movefile(filename, file_destination,'f');
    if status == 1
        % fprintf('Successfully moved file: %s\n',filename)
    else
        fprintf('Moving not successful: %s \n msg:%s \n',filename,msg);
    end
end