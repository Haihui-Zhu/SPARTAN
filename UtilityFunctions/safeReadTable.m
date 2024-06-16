function table = safeReadTable(file)
    % reads table while ignoring extra columns created by extra delimiters
    % at end of lines
    opts = detectImportOptions(file);
    opts.ExtraColumnsRule = 'ignore'; 
    table= readtable(file, opts); clear opts
end