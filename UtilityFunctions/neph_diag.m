% This function checks files for various issues and creates plots when an
% issue has been detected for further investigation

function [data_all, FileTracker] = neph_diag(data_all, data_title, FileTracker, Scatmax, filename, fig_dir)
    
warning('off','backtrace') % don't displays a hyperlink with a line number

filename = filename(1:end - 4);
IN_model = FileTracker.Model;

% data_title = {'DateNum','Hour','RF_PMT', 'RF_REF', 'GF_PMT', 'GF_REF', 'BF_PMT', 'BF_REF', 'RB_PMT', 'RB_REF', 'GB_PMT', 'GB_REF', 'BB_PMT', 'BB_REF',... 1-14
%     'Temp' ,'Amb_Pressure', 'RH', ... 15-17
%     'PMT_DARK', 'Forw_Dark_Ref', 'Back_Dark_Ref', ... 18:20
%     'Back_scatt_Red', 'Back_scatt_Green', 'Back_scatt_Blue', 'Total_scatt_Red', 'Total_scatt_Green', 'Total_scatt_Blue', ... 21-26
%     'Fan_RPM', 'Flow_Pressure', 'Flow', 'DAC'}; % 27:30

PMT = data_all(:, contains(data_title, {'F_PMT', 'B_PMT'}, 'IgnoreCase', true)); % RGB_F_Ref; RGB_B_REf
Ref = data_all(:, contains(data_title, {'F_Ref', 'B_Ref'}, 'IgnoreCase', true)); % RGB_F; RGB_B
Dark = data_all(:, contains(data_title, '_Dark', 'IgnoreCase', true)); % PMT_Dark = 1, Fowrd_Dark_Ref = 2; Bkwd_Dark_Ref = 3

TOT = data_all(:, contains(data_title, {'Back_scatt', 'Total_scatt'}, 'IgnoreCase', true));

%% ---------- Check reference sensor signals --------------------
% There is a natural variation in the reference sensor signal, however the variance should not be more that 10 %. 
% If the variance is found to be greater than 10 % there may be an issue with the reference sensor or its power supply.
ref_signals = {'Red_fwd','Green_fwd','Blue_fwd','Red_bkwd','Green_bkwd','Blue_bkwd'};

scat_cols =[21 24; % back_scatter_red + Total_scatt_Red
            22 25; % back_scatter_green + Total_scatt_green
            23 26; % back_scatter_blue + Total_scatt_blue

            21 24; % back_scatter_red + Total_scatt_Red
            22 25; % back_scatter_green + Total_scatt_green
            23 26];% back_scatter_blue + Total_scatt_blue

for i = 1:size(Ref,2)
    % check for more than 10% variability
    avg = mean(Ref(:,i));
    stdev = std(Ref(:,i));
    if stdev/avg > 0.10
        FileTracker.UnstableRef = FileTracker.UnstableRef + 1;
        warning('File %s  %s ref sensor signal has variability greater than 10 %%, scatter in %s rejected', filename, ref_signals{i},ref_signals{i});
        disp('Continuing with diagnostic tests to check for other issues')
        data_all(:,scat_cols(i,:)) = NaN; % set scatter with poor ref signal to NaN
    end

    if sum(sum(~isnan(data_all(:, scat_cols)))) == 0 % all 3 colors are nan out
        FileTracker.Reject = 1;
        warning('URGENT: File %s shows sign of ref sensor death', filename);
    end

    % check how much signal dropped below 10% of average
    idx_badref = find(Ref(:,i) < avg - (0.1*avg));
    if ~isempty(idx_badref) 
        percent_bad = length(idx_badref)/size(Ref,1);
        data_all(idx_badref,scat_cols(i,:)) = NaN; % set scatter with poor ref signal signal to NaN
        
        if percent_bad < 0.2
            FileTracker.LowRef = FileTracker.LowRef + 1;
            warning('File %s flagged; %s reference signal has points below acceptabe limit from mean', filename, ref_signals{i});
        end
        if percent_bad > 0.20
            warning('File %s rejected; %s reference signal has more than 20 percent of points below acceptabe limit from mean\n', filename,ref_signals{i});
            disp('Continuing with diagnostic tests to check for other issues')
            FileTracker.LowRef = FileTracker.LowRef + 1;
            FileTracker.Reject = 1;
        end
    end

end

% ----- plot reference signals to visualize variability issue -----
fig_title{1} = sprintf('File %s Forward Ref Signals', filename);
fig_title{2} = sprintf('File %s Backward Ref Signals', filename);

if FileTracker.UnstableRef > 0
    fig_title{1} = sprintf('%s\nRef signal variance > 10 %%', fig_title{1});
    fig_title{2} = sprintf('%s\nRef signal variance > 10 %%', fig_title{2});
end

if FileTracker.LowRef > 0
    fig_title{1} = sprintf('%s\n%4.2f%% ref signal drops below 10%% of avg', fig_title{1}, percent_bad * 100);
    fig_title{2} = sprintf('%s\n%4.2f%% ref signal drops below 10%% of avg', fig_title{2}, percent_bad * 100);
end

% ---------------------------------------------------------------
clear   percent_bad idx_badref avg stdev

%% ---------- Check dark signals ----------------- 
% Dark reference and dark PMT values are significantly lower than reference measurements at all wavelengths and have little variation (15% based on experience).
dark_signals = {'PMT_dark','Fwd_dark_ref','Bkwd_dark_ref'};
darklim = 0.8; % dark signal should be lower than this faction times mean PMT or mean Ref
for i = 1:size(Dark, 2)
    % compareing dark to ref
    avg = mean(Dark(:, i));
    if avg > darklim * mean(mean(Ref)) 
        warning('File %s rejected; %s signal greater than %d %% of Ref signal', filename, dark_signals{i}, darklim * 100);
        disp('Continuing with diagnostic tests to check for other issues')
        FileTracker.DarkTooHigh = 1;
        FileTracker.Reject = 1;
    end
    % check for more than 15% variability
    stdev = std(Dark(:, i));
    if stdev / avg > 0.15 % allow for 15 % variation
        warning('File %s rejected; %s signal has variability greater than 15 %%', filename, dark_signals{i});
        disp('Continuing with diagnostic tests to check for other issues')
        FileTracker.UnstableDark = 1;
        FileTracker.Reject = 1;
    end
    clear avg stdev
end

% ----- plot reference signals to visualize variability issue -----
if FileTracker.DarkTooHigh ==1 
    fig_title{1} = sprintf('%s\nDark Signal too high', fig_title{1});
    fig_title{2} = sprintf('%s\nDark Signal too high', fig_title{2});
end
if FileTracker.UnstableDark == 1
    fig_title{1} = sprintf('%s\nDark Signal Unstable', fig_title{1});
    fig_title{2} = sprintf('%s\nDark Signal Unstable', fig_title{2});
end

%% ---------- Check PMT signals  ----------------- 
% PMT signal with LED on should be at least 10 % higher than PMT_dark signal
avg_PMTdark = mean(Dark(:,1));
PMTdark_range = avg_PMTdark + (0.10*avg_PMTdark); 
PMT_signals = {'Red_fwd','Green_fwd','Blue_fwd','Red_bkwd','Green_bkwd','Blue_bkwd'};

for i = 1:size(PMT,2)
    avg = mean(PMT(:,i));
    if avg < PMTdark_range
        FileTracker.LowPMT = FileTracker.LowPMT + 1;
        warning('File %s  %s PMT signal is within 10 percent of PMT_dark signal, scatter in %s rejected', filename, PMT_signals{i},PMT_signals{i});
        disp('Continuing with diagnostic tests to check for other issues')
        data_all(:,scat_cols(i,:)) = NaN; % set scatter with poor PMT signal to NaN
        sigid(i) = 1; 
    end
    clear avg
end
if FileTracker.LowPMT == 6 
    FileTracker.Reject = 1;
    warning('URGENT: File %s shows sign of PMT death', filename);
end
if sum(sum(~isnan(data_all(:, scat_cols)))) == 0 
    FileTracker.Reject = 1; % not need to add any warning since this condition can be met by previous checks
end

% ----- plot reference signals to visualize variability issue -----
fig_title{3} = sprintf('File %s Forward PMT signals', filename);
fig_title{4} = sprintf('File %s Backward PMT signals', filename);
if FileTracker.LowPMT > 0
    fig_title{3} = sprintf('%s\nPMT signal too close to Dark',fig_title{3});
    fig_title{4} = sprintf('%s\nPMT signal too close to Dark',fig_title{4});
end


%% ---------- Temp Press and RH at dark ----------
% Time series plots of temperature, relative humidity, and pressure should 
% be compared to that of the dark reference sensor measurements. 
% Dark reference measurements should not show a similar signal/pattern to 
% ambient conditions. 
% A good correlation between dark and ambient conditions indicates a sensor issue.

variables = {'Temp' ,'Amb_Pressure', 'RH'};

for i = 1:length(variables)
    X = data_all(:,ismember(data_title,variables{i}));
    Y = Dark(:,2); % dark forward
    R = corrcoef(X,Y,'rows','complete');
    R = abs(R(1,2)); % R2
    if R > 0.5 % if there is a significant correlation
        if i == 1
            FileTracker.TempCorrDark = FileTracker.TempCorrDark + 1;
        elseif i == 2
            FileTracker.PressCorrDark = FileTracker.PressCorrDark + 1;
        else
            FileTracker.RHCorrDark = FileTracker.RHCorrDark + 1;
        end
        warning('File %s  %s correlate with dark forward reference', filename, variables{i});
        disp('Continuing with diagnostic tests to check for other issues')
    end
    
    Y = Dark(:,3); % dark backward
    R = corrcoef(X,Y,'rows','complete');
    R = abs(R(1, 2)); % R2
    if R > 0.5 % if there is a significant correlation
        if i == 1
            FileTracker.TempCorrDark = FileTracker.TempCorrDark + 1;
        elseif i == 2
            FileTracker.PressCorrDark = FileTracker.PressCorrDark + 1;
        else
            FileTracker.RHCorrDark = FileTracker.RHCorrDark + 1;
        end
        warning('File %s  %s correlate with dark backward reference', filename, variables{i});
        disp('Continuing with diagnostic tests to check for other issues')
    end

    if sum(table2array(FileTracker(1,end-2:end))) == 6
        FileTracker.Reject = 1;
        warning('URGENT: File %s shows similar pattern to ambient conditions', filename);
    end
    clear X Y R
end

% ----- plot ambient conditions to visualize variability issue -----
legends = {'Temp', 'Pressure', 'RH'};
fig_title{5} = sprintf('%s Ambient conditions\n', filename);
if sum(table2array(FileTracker(1, end - 2:end))) > 0
    for ii = 1:length(legends)
        if table2array(FileTracker(1, end - 3 + ii)) > 0
            fig_title{5} = sprintf('%s%s ', fig_title{5}, legends{ii});
        end
    end
    fig_title{5} = sprintf('%scorrelate with dark ref', fig_title{5});
end


%% ---------- Check for poor flow in IN102 models ---------------
if IN_model == 2
    fig_title{6} = sprintf('DAC signal with acceptable flow\n');
    
    idx_noflow = find(data_all(:, ismember(data_title, 'Flow')) < 1.1 | data_all(:, ismember(data_title, 'Flow')) > 11);

    if ~isempty(idx_noflow) % if there are flow rate lower than 1.1
        percent_noflow = length(idx_noflow) / size(data_all, 1);

        if percent_noflow > 0.8
            warning('File %s rejected; over 80 %% of file has flow below 1.1 lpm or above 11 lpm', filename);
            disp('Continuing with diagnostic tests to check for other issues')

            FileTracker.PoorFlow = 1;
            FileTracker.Reject = 1;
            % ----- Plot DAC to visualize the flow issue -------
            fig_title{6} = sprintf('Over 80 %% of %s has flow < 1.1 lpm or > 11 lpm\nChecking for DAC sweeping', filename);
        
        end
    end

end

%% ---------- Sensitivity Check ------------------
% Check that there is variation in the scatter measurements; scatter should 
% not remain at a given number with no variation for more than a few minutes. 
% Lack of sensitivity in scattering is indicative of a larger issue with 
% the nephelometer and requires immediate investigation by the site manager.

if length(TOT) > 1000 % only check the most recent data points
    rows = length(TOT) - 999:length(TOT);
else
    rows = 1:length(TOT);
end
for i = 1:size(TOT, 2)
    rel_change = abs(diff(TOT(rows, i)) ./ TOT(rows(1:end - 1), i));
    ind = find(rel_change < 0.001); % Relative change less than 0.1%
    cont_ind = find(diff(ind)==1);
    if length(cont_ind) > 10*60/15 % if it shows no variation for 10 min (This algorithm is flawed since there could be scatter stagnations for 40 times)
        FileTracker.Insensitive = FileTracker.Insensitive + 1;
        warning('File %s  %s scatter signal show no variation for more than 10 min, scatter in %s rejected', filename, data_title{i + 20}, data_title{i + 20});
        disp('Continuing with diagnostic tests to check for other issues')
        data_all(:,scat_cols(i,:)) = NaN; % set scatter with poor variation to NaN
    end
    if FileTracker.Insensitive == 6
        FileTracker.Reject = 1;
        warning('URGENT: File %s shows sign of low PMT sensitivity', filename);
    end
end


% ----- plot PMT signals to visualize variability issue -----
fig_title{7} = 'Total backward scatter';
fig_title{8} = 'Total forward scatter';
if FileTracker.Insensitive > 0 
    fig_title{7} = sprintf('%s\nSignal show stagnation for more than 10 min', fig_title{7});
    fig_title{8} = sprintf('%s\nSignal show stagnation for more than 10 min', fig_title{8}); 
end


%% ---- check if total scatter higher than maxi --
% upper linear range of green (532nm) scatter. Signal non-linear above 2500 Mm^-1.
if FileTracker.Reject == 0
    signals = {'Red bkwd', 'Green bkwd', 'Blue bkwd', 'Red total', 'Green total', 'Blue total'};
    j = 0;

    for i = 21:26
        j = j + 1;
        over_max = find(data_all(:, i) > Scatmax);
        data_all(over_max, i) = NaN;

        if length(over_max) / size(data_all, 1) > 0.1
            warning('File %s flagged; over 10%% of data in %s over scatmax \n', filename, signals{j});
            fig_title{7} = sprintf('%s\nOver 10%% of data in %s over scatmax', fig_title{7}, signals{j});
            fig_title{8} = sprintf('%s\nOver 10%% of data in %s over scatmax', fig_title{8}, signals{j});

        end

        clear over_max
    end

end

%% now pass the figure titiles to the function making the diagnosis figures.

% figure 1 - PMT, Ref, and Dark
figure('position',[10 10 1600,1600])
% sub plot 1 & 2
y_label = 'Signal';
variable_to_plot = {'RF_REF', 'GF_REF', 'BF_REF', 'Forw_Dark_Ref'};
legends = {'Red Fwd Ref', 'Green Fwd Ref', 'Blue Fwd Ref', 'Dark Fwd Ref'};
ps = 1;
figure_diagosis(ps,data_all, data_title, variable_to_plot, y_label, legends, fig_title{1})

variable_to_plot = {'RB_REF', 'GB_REF', 'BB_REF', 'Back_Dark_Ref'};
legends = {'Red Bkwd Ref', 'Green Bkwd Ref', 'Blue Bkwd Ref', 'Dark Bkwd Ref'};
ps = 2;
figure_diagosis(ps,data_all, data_title, variable_to_plot, y_label, legends, fig_title{2})

% sub plot 3 & 4
variable_to_plot = {'RF_PMT', 'GF_PMT', 'BF_PMT', 'PMT_DARK'};
legends = {'Red Fwd PMT', 'Green Fwd PMT', 'Blue Fwd PMT', 'Dark PMT'};
ps = 3;
figure_diagosis(ps,data_all, data_title, variable_to_plot, y_label, legends, fig_title{3})

variable_to_plot = {'RB_PMT', 'GB_PMT', 'BB_PMT', 'PMT_DARK'};
legends = {'Red Bkwd PMT', 'Green Bkwd PMT', 'Blue Bkwd PMT', 'Dark PMT'};
ps = 4;
figure_diagosis(ps,data_all, data_title, variable_to_plot, y_label, legends, fig_title{4})

% sub plot 5
variable_to_plot = {'Temp', 'Amb_Pressure', 'RH', 'Forw_Dark_Ref', 'Back_Dark_Ref'};
y_label = 'Signal/Conditions';
legends = {'Temp', 'Pressure', 'RH', 'Forw Dark Ref', 'Back Dark Ref'};
ps = 5;
figure_diagosis(ps,data_all, data_title, variable_to_plot, y_label, legends, fig_title{5})

% sub plot 6   
if IN_model == 2
variable_to_plot = {'DAC'};
y_label = 'DAC Signal';
legends = 'DAC Signal';
ps = 6;
figure_diagosis(ps, data_all, data_title, variable_to_plot, y_label, legends, fig_title{6})
end

% save figure
sfname = sprintf('%s/%s.png', fig_dir, filename);
saveas(gcf,sfname)
close 

% figure 2 - total scatter 
if sum(sum(~isnan(data_all(:, scat_cols)))) > 0
    figure('position', [10 10 1000, 1000])
    variable_to_plot = {'Back_scatt_Red', 'Back_scatt_Green', 'Back_scatt_Blue'};
    y_label = 'Signal';
    legends = {'Red Bwd', 'Green Bwd', 'Blue Bwd'};
    ps = 7;
    figure_diagosis(ps, data_all, data_title, variable_to_plot, y_label, legends, fig_title{7})

    variable_to_plot = {'Total_scatt_Red', 'Total_scatt_Green', 'Total_scatt_Blue'};
    legends = {'Red Fwd', 'Green Fwd', 'Blue Fwd'};
    ps = 8;
    figure_diagosis(ps, data_all, data_title, variable_to_plot, y_label, legends, fig_title{8})

    % save figure
    saveas(gcf, sprintf('%s/%s_totalscatter.png', fig_dir, filename));
    close 
end
end    



function figure_diagosis(ps,data_all, data_title, variable_to_plot, y_label, legends, fig_title)
fz = 6;
Colors = {'r','g','b','k','c'};
if ps < 7
    subplot(3,2,ps)
else
    subplot(2,1,ps-6)
end
if contains(y_label, 'conditions', 'IgnoreCase', true) % specifically for plotting ambient conditions
    yyaxis left % ploting ambient conditions
    for ii = 1:3
        x = 1:size(data_all,1);
        plot(x, data_all(:, ismember(data_title, variable_to_plot{ii})),Colors{ii});
        hold on
    end
    ylabel('Temp (C) Pressure (kPa) RH (%)', 'FontWeight', 'bold', 'FontSize', fz);

    yyaxis right % ploting dark scatter
    for ii = 4:5
        x = 1:size(data_all,1);
        plot(x, data_all(:, ismember(data_title, variable_to_plot{ii})), Colors{ii});
        hold on
    end
    ylabel('Signal', 'FontWeight', 'bold', 'FontSize', fz);
else
    for ii = 1:length(variable_to_plot)
        x = 1:size(data_all,1);
        plot(x, data_all(:, ismember(data_title, variable_to_plot{ii})), Colors{ii});
        hold on
    end
    ylabel(y_label, 'FontWeight', 'bold', 'FontSize', fz);
end
xlabel('Data point # ','FontWeight','bold','FontSize',fz);

set(gcf,'outerposition',[520 260 1000 550]);
ax=gca;
ax.FontSize = fz;
ax.FontWeight = 'bold';
xtickangle(45)

legend(legends,'Location','eastoutside','FontWeight','bold','FontSize',fz)

title(fig_title,'FontSize',fz+2)

end

