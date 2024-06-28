
function find_ic_zeros(data_samples, Master_IDs, dates, Master_masstype, ic_zero_fname, data_type)
% Sulfate, Nitrate, Chloride, Sodium, Ammonium
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
        if Master_masstype(ii) == 1  ||  Master_masstype(ii) == 2 % only report valid filters (PM25 or PM10)
            if data_samples(ii,jj) <= 0 || isnan(data_samples(ii,jj))
                if data_type == 1 && jj == 7 % sulfate
                    fprintf(fileID, '%s   %s   %s  %.1f  \n', Master_IDs{ii}, dates{jj}, species{jj}, data_samples(ii,jj));

                elseif  data_type == 1 && jj == 5 % nitrate
                    fprintf(fileID, '%s   %s   %s  %.1f  \n', Master_IDs{ii}, dates{jj}, species{jj}, data_samples(ii,jj));

                elseif  data_type == 1 && jj == 2 % chloride
                    fprintf(fileID, '%s   %s   %s  %.1f  \n', Master_IDs{ii}, dates{jj}, species{jj}, data_samples(ii,jj));

                elseif  data_type == 2 && jj == 2 % sodium
                    fprintf(fileID, '%s   %s   %s  %.1f  \n', Master_IDs{ii}, dates{jj}, species{jj}, data_samples(ii,jj));
                    
                elseif  data_type == 2 && jj == 3 % Ammonium
                    fprintf(fileID, '%s   %s   %s  %.1f  \n', Master_IDs{ii}, dates{jj}, species{jj}, data_samples(ii,jj));
                end
            end
        end
    end
end
fclose(fileID);

