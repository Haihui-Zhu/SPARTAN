function status = CopyFile(filename,destination)
    % creat destination folder if not exist
    if ~exist(destination,'dir')
        mkdir(destination)
    end
    % move specified filename to specified destination
    [status,~,~]= copyfile(filename, destination );
end