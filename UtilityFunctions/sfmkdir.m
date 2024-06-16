function  dirname = sfmkdir(dirname)

if ~exist(dirname,'dir')
    mkdir(dirname)
end

end