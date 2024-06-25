
function find_ic_zeros(data_samples, Master_IDs, ic_zero_fname, data_type)

if data_type == 1
    species = {'Fluoride','Chloride','Nitrite','Bromide','Nitrate','Phosphate','Sulfate'};
elseif data_type == 2
    species = {'Lithium','Sodium','Ammonium','Potassium','Magnesium','Calcium'}; 
else
    error('Date type not regonized. Use 1 for anion, 2 for cation.\n')
end

% print filter ID, IC species
fileID = fopen(ic_zero_fname, 'a');
if fileID == -1
    error('File could not be opened for appending');
end

for ii = 1:length(Master_IDs) % number of filters
    for jj = 1:length(species)
        if data_type == 1 && jj == 1 % Fluoride is always nan for now. No need to print
            continue
        else
            if data_samples(ii,jj) <= 0 || isnan(data_samples(ii,jj))
                fprintf(fileID, '%s   %s   %.1f  \n', Master_IDs{ii},species{jj}, data_samples(ii,jj));
            end
        end
    end
end
fclose(fileID);

