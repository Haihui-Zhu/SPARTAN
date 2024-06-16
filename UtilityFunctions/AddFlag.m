function flagOut = AddFlag(orig_flags,new_flag)

orig_flags = strtrim(orig_flags); % remove leading and tailing blank 

% ----- preventing repeating flags ----------------
if ~isempty(new_flag)
    while contains(orig_flags, new_flag)
        orig_flags = erase( orig_flags, new_flag); % prevent repeated part
    end
end

% ------ a special check for 'qc2-FSx-Li' -------------------------------
% when there is 'qc2-FSx-Li' and need to add 'FS' previous section will
% erroenously remove 'FS' in 'qc2-FSx-Li'.
if contains(orig_flags,'qc2-x-Li') 
    replace(orig_flags,'qc2-x-Li','qc2-FSx-Li')
end

% ------ adding new flag -------------------------------------------------
if isempty(new_flag) % just in case new flag is empty
    flagOut = orig_flags;
else
    flagOut = [orig_flags,'; ',new_flag];
end

% ------ Formatting flag -------------------------------------------------
while contains(flagOut, "; ;")
    flagOut = replace(flagOut,"; ;",";"); % remove extra ';' in the middle of a flag
end
while contains(flagOut, "  ")
    flagOut = replace(flagOut,"  "," "); % remove extra blank in the middle of a flag
end

if ~isempty(flagOut)
    while flagOut(1) == ';' || flagOut(1) == ' ' % prevent starting with ';' or a blank
        flagOut = flagOut(2:end);
    end
end

flagOut = strtrim(flagOut); % Remove the leading and trailing whitespace

end