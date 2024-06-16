function direc = find_root_dir(debug_mode)
    currentdir = pwd;

    % For Haihui 
    if contains(currentdir,'Volumes') && debug_mode == 0
        direc='/Volumes/rvmartin/Active/SPARTAN-shared/';
        
    elseif contains(currentdir,'Volumes') && debug_mode == 1
        direc='/Volumes/rvmartin/Active/haihuizhu/6.SPARTAN/'; 
    
    elseif contains(currentdir,'storage1/fs1') && debug_mode == 0
        direc='/storage1/fs1/rvmartin/Active/SPARTAN-shared'; 

    elseif contains(currentdir,'storage1/fs1') && debug_mode == 1
        direc='/storage1/fs1/rvmartin/Active/haihuizhu/6.SPARTAN/'; 
        
        
    % for Chris or anyone using Windows OS
    elseif contains(currentdir,{'\\storage1.ris.wustl.edu','Z:'}) && debug_mode == 0 
        direc='\\storage1.ris.wustl.edu\rvmartin\Active\SPARTAN-shared\';

    else
        error('Cannot identify the root directory')
    end

end
