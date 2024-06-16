% this script delete a sheet in excel function
function delete_sheet(filename, sheetname)
    newFilename = filename(1:end - 5);
    newFilename(end + 1:end + 6) = '2.xlsx';

    sheets = sheetnames(filename);

    if sum(strcmp(sheets, sheetname)) == 0
        error('%s not found in file!', sheetname)
    else

        for i = 1:length(sheets)

            if ~strcmp(sheets{i}, sheetname)
                % Read data from the sheet
                data = readtable(filename, 'Sheet', sheets{i});

                % Write data to the new file
                % Note: 'WriteMode', 'append' is used to add sheets to the file
                writetable(data, newFilename, 'Sheet', sheets{i});
            end
        end

        % now delete the old file and rename the new one
        delete(filename)
        for i = 1:length(sheets)
            if ~strcmp(sheets{i}, sheetname)
                data = readtable(newFilename, 'Sheet', sheets{i});
                writetable(data, filename, 'Sheet', sheets{i});
            end
        end
        delete(newFilename)

    end

end
