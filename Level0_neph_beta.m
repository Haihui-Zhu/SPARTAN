%% This script provides an initial, basic screening of raw (from the site operator) nephelometer files
% THe best way to use this script is to look at the generated plots. No
% built-in programming has been provided, it is the responsibility of the
% user to determine validity of nephelometer files

% Written by: Crystal Weagle
% Date updated: February 22, 2022

close all; clear; clc

%% ----------USER SWITCHES ----------
pointsperday = 5760; % approximate number of points generated per day assuming one data point every 15 seconds
pause('on');
% DIRECTORIES
% for use at WashU
direc_input = '//storage1.ris.wustl.edu/rvmartin/Active/SPARTAN-shared/Site_Sampling/Neph_Lvl0/for_testing/';
direc_plot = '//storage1.ris.wustl.edu/rvmartin/Active/SPARTAN-shared/Site_Sampling/Neph_Lvl0/Level0_plots/';
direc_output = '//storage1.ris.wustl.edu/rvmartin/Active/SPARTAN-shared/Site_Sampling/';


site_details = readtable('//storage1.ris.wustl.edu/rvmartin/Active/SPARTAN-shared/Site_Sampling/Site_details.xlsx','PreserveVariableNames',true); 
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));

% - Used by Crystal for development and testing
% direc_input = '/Users/crystalweagle/Documents/SPARTAN/Neph_project_2022/Testing/test_files'; % this is the location of the raw files
% direc_output = '/Users/crystalweagle/Documents/SPARTAN/Neph_project_2022/Testing/Valid'; % this is where files determined to be valid are moved to, this is where they ultimately undergo further processing
% direc_plot = '/Users/crystalweagle/Documents/SPARTAN/Neph_project_2022/Testing/Level0_plots'; % this is where generated plots are saved

%% Get list of files for testing and load


% Get list of files in "Data 0 validation" directory
names = dir(sprintf('%s/IN*',direc_input));
files = {names.name}; clear names
files = files';

for f = 1:length(files)
    
    % find the names and sizes of files in the direc_raw
    fid=fopen(sprintf('%s/%s',direc_input,files{f}));
    filename = sprintf('%s/%s',direc_input,files{f});
    
    test_file=dir(filename);
    file_size=test_file.bytes; clear test_file
    
    if file_size < 5000 % skips over files that are weird computer files
        continue
    end
    % get the name of the file without the ".csv" at the end
    nameB = char(files{f});
    name = nameB(1:end-4); clear nameB

% ------- Read in file --------
    % this file format has a date format with dashes
            A = textscan(fid,'%s  %n/%n/%n  %n:%n:%*n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1);
            if isempty(A{1,8}) % did not read first value after date/time column
                clear A
                % this file format has a date format with slashes
                A = textscan(fid,'%s  %n-%n-%n  %n:%n:%*n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1);
            end
            if isempty(A{1,8}) % did not read first value after date/time column, no seconds in file time
                clear A
                % this file format has a date format with slashes
                A = textscan(fid,'%s  %n-%n-%n  %n:%n  %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','delimiter',',','headerlines',1);
            end
            if isempty(A{1,8}) % has one 'DATE_TIME' column and 6 addtional dates and time columns 
                clear A
                A = textscan(fid,'%s %s %n %n %n %n %n %n %n %n  %n %n %n %n %n %n %n %n %n %n  %n %n %n %n %n %n %n %n %n %n  %n %n %n %n %n %n %n %n %n %n  %n %n %n %n %n %n %n %n %s %n  %n','delimiter',',','headerlines',1);
            end
            % To delete incomplete last raw due to (maybe) SD card memory issue [Haihui 2021-6-6]
            for i = 1:size(A,2)
                if length(A{i}) ~= length(A{end})
                A{i}(end) = []; end
            end
            if length(A) > 40 % raw data with additional dates and time columns
                A([2 8 24 25 26 46:end]) = []; % clear the additional columns 
                scatter = cell2mat(A(25:30));
                met_vars = cell2mat(A(19:21));
                flow_vars = cell2mat(A(31:34));
            else
                scatter = cell2mat(A(25:30));
                met_vars = cell2mat(A(19:21));
                flow_vars = cell2mat(A(31:34));
            end
 % ------ column descriptions -----
    % scatter matrix
    % 1 = backscatter red (1/Mm)
    % 2 = backscatter green (1/Mm)
    % 3 = backscatter blue (1/Mm)
    % 4 = total scatter red (1/Mm)
    % 5 = total scatter green (1/Mm)
    % 6 = total scatter blue (1/Mm)
    % met_vars matrix
    % 1 = Temp (degC)
    % 2 = Atmospheric Pressure (kPa)
    % 3 = RH%
    % flow_vars matrix
    % 1 = fan RPM
    % 2 = Flow pressure
    % 3 = Flow rate (LPM)
    % 4 = DAC
    % -------------------------------   
    
    if isempty(A{1,8}) % did not read first value after date/time column
        prompt = sprintf('Whoops! It looks like file %s was not read properly, move on to next file to be validated?', name);
        definput = {'y/n'};
        dlgTitle = 'User Input';
        dims = [1 50];
        Uinput = inputdlg(prompt,dlgTitle,dims,definput);
        if Uinput{1} == 'y'
            continue
        elseif Uinput{1} == 'n'
            disp('Level0_neph_beta terminated by user due to invalid file')
            return
        end
    end
% ---------    
    
    % estimates how many days are in the data file based on the number of
    % points and asks for user input to determine how many days per time series will be plotted.
    
    file_days = length(scatter)/pointsperday;
    
    prompt = {sprintf(' File %s contains approximately %0.1f days of data \n Enter the number of days do you want to display:',name, file_days)};
    definput = {'999 = all data, 1 = 1 day, 2 = 2 days, etc.'};
    dlgTitle = 'User Input';
    dims = [1 50];
    Uinput = inputdlg(prompt,dlgTitle,dims,definput);
    DaysPerPlot = str2double(Uinput{1});
    clear prompt definput Uinput
    
    %% ---------- Create plots ---------
    
    if DaysPerPlot == 999 % plot all the data on one graph
        x = 1:length(scatter);
        data_index = x;
        file = files{f};
        p = 1;
        fig(p) = neph_plots(x,data_index, scatter, met_vars, flow_vars, p)
        
    else % days per plot is not 999, therefore will divide data up into the set number of points per plot and generate the amount of plots necessary to get all data visiblet
        integ = floor(length(scatter)/(pointsperday*DaysPerPlot)); % approximate number of points per plot based on size of file and user input
        frac = (length(scatter)/(pointsperday*DaysPerPlot)) - integ;
        if frac < 0.2 % Then number of points in the file will be evenly distributed to the other plots, making each plot slightly longer than the user inputclear PPP
            PPP = length(scatter)/integ; % PPP = points per plot
            Xrange = 1:PPP:length(scatter)+PPP;
            Xrange(end) = length(scatter);
            num_plots = integ;
            clear integ frac
        elseif frac > 0.2 % create an additional figure for the remaining points
            num_plots = integ + 1;
            PPP = pointsperday*DaysPerPlot;
            Xrange = 1:PPP:(PPP*integ)+PPP;
            Xrange(end+1) = length(scatter);
            clear integ frac
        end
        for p = 1:num_plots
            x = Xrange(p):Xrange(p+1);
            data_index = x;
            file = files{f};
            neph_plots(x,data_index, scatter, met_vars, flow_vars, p)
            clear x data_index
        end
    end
    
    pause 
    
    prompt = {'Is this file cleared to move onto Level 1 Data Validation?',sprintf('Do you want to save the figure(s) generated from file %s ?',name),'If yes, enter the site code for this file:'};
    definput = {'y/n', 'y/n','(e.g., ZAPR)'};
    dlgTitle = 'User Input';
    dims = [1 50];
    Uinput = inputdlg(prompt,dlgTitle,dims,definput);
    
    if Uinput{1} == 'y' % file is valid, move to site raw folder for further processing
        site_index = find(strcmp(Site_codes, Uinput(3))); 
        movefile(sprintf('%s%s',direc_input,files{f}), sprintf('%s%s_%s/Nephelometer_Data/Raw',direc_output,Uinput{3}, Site_cities{site_index}),'f');
    elseif Uinput{1} == 'n'  % file is not valid, move to site rejected folder?
        site_index = find(strcmp(Site_codes, Uinput(3))); 
        movefile(sprintf('%s%s',direc_input,files{f}), sprintf('%s%s_%s/Nephelometer_Data/Rejected',direc_output,Uinput{3}, Site_cities{site_index}),'f');
    end
    if Uinput{2} == 'y'
        FigList = findobj(allchild(0),'flat','type','figure');
        if length(FigList) == 1
        saveas(FigList(1),sprintf('%s%s/%s_allpoints.png',direc_plot,Uinput{3},name));
        else
            for i = 1:length(FigList)
               saveas(FigList(i),sprintf('%s%s/%s_plot%d.png',direc_plot,Uinput{3},name,i)) ;
            end
        end
    elseif Uinput{2} == 'n'
        close all
    end
    
    fclose(fid)
    close all
    clear prompt definput Uinput dims dlgTitle 
    clear A DaysPerPlot file file_days file_size filename flow_vars met_vars name num_plots p PPP scatter Xrange
end


%% %%%%%%% PLOTTING FUNCTION TO KEEP CODE MORE READABLE %%%%%%%%%%%%
function fig = neph_plots(x,data_index, scatter, met_vars, flow_vars, p);

fig(p) = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% ----- plot of total scatter -----
ax(1) = subplot(6,1,1);
plot(x,scatter(data_index,4),'r')
hold on
plot(x,scatter(data_index,5),'g')
plot(x,scatter(data_index,6),'b')
% set(gcf,'position',[100 267 1400 600]);
ax(1).FontSize = 10;
ax(1).FontWeight = 'bold';
% xtickangle(45)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel('Total Scatter','FontWeight','bold','FontSize',10);
% xlabel('Data point','FontWeight','bold','FontSize',10);
% T(1) = title(sprintf('%s Total Scatter',file), 'FontWeight','bold','FontSize',14);
xlim([x(1) x(end)]);
xtickformat('%.0f')

% ----- plot backscatter -----
ax(2) = subplot(6,1,2);
plot(x,scatter(data_index,1),'r')
hold on
plot(x,scatter(data_index,2),'g')
plot(x,scatter(data_index,3),'b')
% set(gcf,'position',[100 267 1400 600]);
ax(2).FontSize = 10;
ax(2).FontWeight = 'bold';
% xtickangle(45)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel('Backscatter','FontWeight','bold','FontSize',10);
% xlabel('Data point','FontWeight','bold','FontSize',10);
% T(2) = title(sprintf('%s BackScatter',file), 'FontWeight','bold','FontSize',14);
xlim([x(1) x(end)]);
xtickformat('%.0f')

% plot flow rates
ax(3) = subplot(6,1,3);
plot(x,flow_vars(data_index,3),'k')
% set(gcf,'position',[100 267 1400 600]);
ax(3).FontSize = 10;
ax(3).FontWeight = 'bold';
% xtickangle(45)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel('Flow rate (LPM)','FontWeight','bold','FontSize',10);
% xlabel('Data point','FontWeight','bold','FontSize',10);
% T(3) = title(sprintf('%s Flow Rate',file), 'FontWeight','bold','FontSize',14);
xlim([x(1) x(end)]);
ylim([-0.5 5.5])
xtickformat('%.0f')

% plot fan RPM
ax(4) = subplot(6,1,4);
plot(x,flow_vars(data_index,1),'k')
% set(gcf,'position',[100 267 1400 600]);
ax(4).FontSize = 10;
ax(4).FontWeight = 'bold';
% xtickangle(45)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel('Fan RPM','FontWeight','bold','FontSize',10);
% xlabel('Data point','FontWeight','bold','FontSize',10);
% T(4) = title(sprintf('%s Fan RPM',file), 'FontWeight','bold','FontSize',14);
xlim([x(1) x(end)]);
xtickformat('%.0f')

% plot DAC signal
ax(5) = subplot(6,1,5);
plot(x,flow_vars(data_index,4),'k')
% set(gcf,'position',[100 267 1400 600]);
ax(5).FontSize = 10;
ax(5).FontWeight = 'bold';
% xtickangle(45)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel('DAC signal','FontWeight','bold','FontSize',10);
% xlabel('Data point','FontWeight','bold','FontSize',10);
% T(5) = title(sprintf('%s DAC Signal',file), 'FontWeight','bold','FontSize',14);
xlim([x(1) x(end)]);
xtickformat('%.0f')

% plot temp. and %RH
ax(6) = subplot(6,1,6);
yyaxis left
plot(x,met_vars(data_index,1))
ylabel('Temp. (degC)','FontWeight','bold','FontSize',10);
yyaxis right
plot(x,met_vars(data_index,3))
ylabel('RH (%)','FontWeight','bold','FontSize',10);
xlabel('Data point','FontWeight','bold','FontSize',12);
% set(gcf,'position',[100 267 1400 600]);
ax(6).FontSize = 10;
ax(6).FontWeight = 'bold';
xtickangle(45)
xlim([x(1) x(end)]);
% T(6) = title(sprintf('%s Temp. & %%RH',file), 'FontWeight','bold','FontSize',14);
xtickformat('%.0f')

ax(1).Position = [0.03 0.825 0.94 0.15];
ax(2).Position = [0.03 0.67 0.94 0.15];
ax(3).Position = [0.03 0.515 0.94 0.15];
ax(4).Position = [0.03 0.36 0.94 0.15];
ax(5).Position = [0.03 0.205 0.94 0.15];
ax(6).Position = [0.03 0.05 0.94 0.15];

end





