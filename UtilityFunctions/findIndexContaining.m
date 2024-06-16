function [row, col] = findIndexContaining(string, file)
     % find row and column index of string in file
     index = find(contains(file,string));
     [row, col] = ind2sub(size(file),index); clear index
end