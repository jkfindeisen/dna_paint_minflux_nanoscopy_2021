function results2figures_SI_Fig123()
% Script to get Supplementary Note Figures 1, 2 and 3
%
% Gives Figures of measurement series with varyig laser power,
% pinhole diameter or imager concentration for
% variables: t_btw (time elapsed between events),
% f_bg (estimated background frequency), CFR(center-frequency-ratio),
% sigma_r (radial localization precision)
%
% Requirements: Need .xls-Sheet with results of analysis of MINFLUX
% measurements series (given by test_collection.m)
% --> in results-folder: e.g. code_for_paper\results\laserpower.xls

close all

% determine file locations
root_folder = get_root_folder();

%define format for saving figures
format = {'.png'};

%% User interaction: Selection of analysis parameter
% choose which of the analyzed results would like to plot
list = {'Laser power', 'Pinhole diameter', 'Imager concentration'};

[indx,tf] = listdlg('PromptString',{'Choose a measurement series.',...
    'Only one parameter can be analyzed at a time.',''},...
    'SelectionMode','single','ListString',list);

if isempty(indx)
    error('No measurement series was selected.');
else
    disp(['Starting analysis of: ', list{indx}]);
    if indx == 1
        variable = 'laserpower';
        var_label = ['laser power (',char(181),'W)']; % xlabel
        
        % plotting parameters: l & p are position variables of legend
        l2 = ['NorthEast'];
        l3 = l2;
        p1 = [0.62,0.72,0.15,0];
        p3 = p1;
        % laser power step in plot
        step_x =10;
        
    elseif indx == 2
        variable = 'pinhole';
        var_label = 'pinhole diameter (AU)'; % xlabel
        % plotting parameters: l & p are position variables of legend
        l2 = ['NorthWest'];
        l3 = ['NorthEast'];
        p1 =[0.315,0.72,0.15,0];
        p3 = [0.62,0.72,0.15,0];
        % pinhole diameter step in plot
        step_x =0.15;
        
    elseif indx == 3
        variable = 'concentration';
        var_label = 'imager concentration (nM)'; % xlabel
        % plotting parameters: l & p are position variables of legend
        l2 = ['NorthWest'];
        l3 = ['NorthEast'];
        p1 =[0.315,0.72,0.15,0];
        p3 = [0.62,0.72,0.15,0];
        % concentration step in plot
        step_x =[1,0.2];
        
    else
        error('Wrong measurement series was selected.');
    end
end

%% saving and plot preparations

% folder for saving figures
mkdir(fullfile(root_folder,'figures',variable));

% use fixed axis ranges
ax_t = [0 16]; % tbtw
ax_f = [0 12*1e4]; %fbg
ax_c = [0.4 0.9]; %cfr
ax_s = [2.0 6]; %sigma r

% groupname of series (for legend)
groupn = ' series ';

% colors
cp = [0 0 0];
c1 = [0.619607843137255,0,0];
c2 = [1,0.400000000000000,0];
c3 = [1,0.827450980392157,0.141176470588235];

sm = 30; % size of marker (asteriks)
lw = 1.5; % linewidth
fsize=16; % fontsize
fsize_ax = 12; % fontsize of axis

%% read in results table

% read in results file of previously analyzed measurement series
if isfile(join([fullfile(root_folder, 'results'),'\',variable,'.xls']))
    T_all = readtable(join([fullfile(root_folder, 'results'),'\',variable,'.xls']));
else
    error('Results file not found. Check directories, or rerun analysis script');
end

if indx == 1
    var_x = T_all.MinfluxLaserpower;
    % change percentage to uW
    % within the MINFLUX software the unit of the excitation laser power is
    % a percentage value. We measured the corresponding uW value in the
    % sample as the following: (all values are for the first iteration.
    % within each iteration power is increased)
    var_x(var_x == 4) = 17; % 4% == 17 uW
    var_x(var_x == 6) = 26; % 6% == 26 uW
    var_x(var_x == 8) = 35; % 8% == 35 uW
    var_x(var_x == 10) = 44; % 10% == 44 uW
    var_x(var_x == 12) = 54; % 12% == 54 uW
    var_x(var_x == 14) = 62; % 14% == 62 uW
    var_x(var_x == 16) = 71; % 16% == 71 uW
    
elseif indx == 2
    var_x = T_all.MinfluxPinholesizeAU;
    % according to the manufacturer, the actual pinhole size was half of
    % the value displayed in the MINFLUX software (will be corrected after
    % update). To consider correct values had to divide pinhole size by 2.
    disp('Use half of Pinhole');
    var_x = var_x/2;
elseif indx == 3
    var_x = T_all.ConcentrationnM;
else
    error('Wrong variable for combined analysis given');
end

var_date = T_all.Date;  % date of series (all measurements of same series
% were done on same day)
parameter = unique(var_x); % get all parameters examined in series

% convert var_date to string array
for u = 1:length(var_date)
    array_date(u) = string(var_date(u));
end

% get dates of series to assign data accordingly
series_date = unique(array_date);

% define name of series of series group for legend
groups_names="";
for k = 1:1:numel(series_date)
    groups_names(k) = join([groupn, num2str(k)]);
end

% define new variable, to assign data to corresponding series
for u=1:1:numel(array_date)
    [~,idxx] = find(strcmp(array_date(u), series_date));
    var_series(u) = groups_names(idxx);
end

%% define structure to save data
structPaper = struct();
structPaper.var = parameter;

fig = figure();
fig.Position = [100, 100, 1700, 1000];

%% plot t btw (first 100)
% in previous analysis only first 100 molecules were considered, as all
% series had at least 100 molecules
ax = subplot(2,3,1);
hold on;

% plot mean of measurement series as errorbar
[me,st]=giveMean(var_x, parameter,T_all.MedianOffTimeBtwEv_First2);
gs1 = errorbar(parameter, me, st,'.-','MarkerSize',sm,'Color',cp,'Linewidth',lw, 'DisplayName', 'mean');

%define axis and plot properties and colors
xbounds = xlim;
if indx == 1 % for laserpower series redefine range
    xbounds = [10 80];
end

ybounds = ax_t;
xticks(xbounds(1):step_x(1):xbounds(2));
yticks(ybounds(1):2:ybounds(2));
ylim(ax_t);
ax.Color = 'none';
colors = [c1;c2;c3];
% scatter plot of each individual series
gs2 = gscatter(var_x, T_all.MedianOffTimeBtwEv_First2,var_series(:),colors,'*',sm/5);

for i = 1 : numel(gs2)
    gs2(i).LineWidth = lw;
end

% x & y label
xlabel(var_label,'fontsize',fsize);
ylabel('median{\it t}_{btw} (s)','fontsize',fsize);

ax.FontSize = fsize;
pbaspect(ax,[1 1 1])

% plot legend
leg = legend(ax, gs1, 'Position',p3,'box','off', 'Fontsize',10);
leg.Visible = 'off'; % to plot legends --> comment this line

% save calculated mean and std in structure
structPaper.mean_tbtw = me;
structPaper.std_tbtw = st;

%% FBG estimated background frequency
ax = subplot(2,3,2);
hold on;

% plot mean of measurement series as errorbar
[me,st]=giveMean(var_x, parameter,T_all.MedianFBG);
gs1 = errorbar(parameter, me, st,'.-','MarkerSize',sm,'Color',cp,'Linewidth',lw, 'DisplayName', 'mean');

%define axis and plot properties and colors
xbounds = xlim;
if indx == 1
    xbounds = [10 80];
end

ybounds = ax_f;
xticks(xbounds(1):step_x(1):xbounds(2));
yticks(ybounds(1):2*1e4:ybounds(2));
ylim( ax_f);
xlim(xbounds);
ax.Color = 'none';
colors = [c1;c2;c3];
% scatter plot of each individual series
gs2 = gscatter(var_x, T_all.MedianFBG,var_series(:),colors,'*',sm/5);

for i = 1 : numel(gs2)
    gs2(i).LineWidth = lw;
end

% x & y label
xlabel(var_label,'fontsize',fsize);
ylabel('median {\it f}_{bg} (1/s)','fontsize',fsize);

ax.FontSize = fsize;
pbaspect(ax,[1 1 1])

% legend
leg = legend(ax, gs1, 'Position', p1,'box','off', 'Fontsize',10);
leg.Visible = 'off'; % to plot legends --> comment this line

% save calculated mean and std in structure
structPaper.mean_fbg = me;
structPaper.std_fbg = st;

%% CFR center-frequency-ratio
ax = subplot(2,3,3);
hold on;

% plot mean of measurement series as errorbar
[me,st]=giveMean(var_x, parameter,T_all.MedianCFR);
gs1 = errorbar(parameter, me, st,'.-','MarkerSize',sm,'Color',cp,'Linewidth',lw, 'DisplayName', 'mean');

%define axis and plot properties and colors
xbounds = xlim;
if indx == 1
    xbounds = [10 80];
end

ybounds = ax_c;
xticks(xbounds(1):step_x(1):xbounds(2));
yticks(ybounds(1):0.1:ybounds(2));
ylim( ax_c);
xlim(xbounds);
ax.Color = 'none';
colors = [c1;c2;c3];
gs2 = gscatter(var_x, T_all.MedianCFR,var_series(:),colors,'*',sm/5);

for i = 1 : numel(gs2)
    gs2(i).LineWidth = lw;
end

% x & y label
xlabel(ax,var_label,'fontsize',fsize);
ylabel(ax,'median CFR','fontsize',fsize);

ax.FontSize = fsize;
pbaspect(ax,[1 1 1])

% legend
leg = legend(ax, gs1, 'Position', p1,'box','off','Fontsize',10);
set(leg,'visible','off') % to plot legends --> comment this line

% save calculated mean and std in structure
structPaper.mean_cfr = me;
structPaper.std_cfr = st;

%% SIGMA R (sqrt)
ax = subplot(2,3,4);
hold on;

% plot mean of measurement series as errorbar
[me,st]=giveMean(var_x, parameter,T_all.MedianSigmaR);
gs1 = errorbar(parameter, me, st,'.-','MarkerSize',sm,'Color',cp,'Linewidth',lw, 'DisplayName', 'mean');

%define axis and plot properties and colors
xbounds = xlim;
if indx == 1
    xbounds = [10 80];
end

ybounds = ax_s;
xticks(xbounds(1):step_x(1):xbounds(2));
yticks(ybounds(1):0.5:ybounds(2));
ylim( ax_s);
xlim(xbounds);
ax.Color = 'none';
colors = [c1;c2;c3];
gs2 = gscatter(var_x, T_all.MedianSigmaR,var_series(:),colors,'*',sm/5);

for i = 1 : numel(gs2)
    gs2(i).LineWidth = lw;
end

% x & y label
xlabel(ax,var_label,'fontsize',fsize);
ylabel(ax,'median {\it \sigma}_r (nm)','fontsize',fsize);

ax.FontSize = fsize;
pbaspect(ax,[1 1 1])

% legend
leg = legend(ax, gs1, 'Position', p1,'box','off','Fontsize',10);
set(leg,'visible','off') % to plot legends --> comment this line

% save calculated mean and std in structure
structPaper.mean_sigmar = me;
structPaper.std_sigmar = st;


%% save infos (mean & std) in table
tablePaper = struct2table(structPaper);
disp(tablePaper);

disp('Write all infos into table file');
writetable(tablePaper, fullfile(root_folder,'figures',join([variable,'_mean_std','.xls'])));

end

function [meanList,stdList] = giveMean(var_x, variable_list, originalData)
% calculates mean and std of values belonging to same parameter out of
% different measurement series

el_mean = zeros(size(variable_list));
el_std = el_mean;
get_mean_val = originalData;

% get mean and std to according variable parameter
for g = 1:1:numel(variable_list)
    el_mean(g) = mean(get_mean_val(var_x == variable_list(g)));
    el_std(g) = std(get_mean_val(var_x == variable_list(g)));
end
meanList = el_mean;
stdList = el_std;
end
