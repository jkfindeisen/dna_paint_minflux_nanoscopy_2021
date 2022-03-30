function create_supplemental_figures_minflux_performance()
% Takes output of the analysed Minflux performance measurement series (see
% script analyse_measurement_series.m) and aggregates the data to display
% performance measures.
%
% Recreates figures of Suppl. Note (not exactly).
%
% This file is part of the supplementary software for "DNA-PAINT MINFLUX
% nanoscopy", 2021 by Lynn M. Ostersehlt, Daniel C. Jans, Anna Wittek,
% Jan Keller-Findeisen, Steffen J. Sahl, Stefan W. Hell, and Stefan Jakobs

fprintf('display aggregate Minflux performance measurement series data\n');

close all;

if ~exist('get_root_folder.m', 'file')
    initialize();
end

% determine file locations
root_folder = get_root_folder();
base_folder = [root_folder, filesep, 'results', filesep, 'performance'];

% use fixed axis ranges
ax_t = [0 4]; % tbtw
ax_t2 = [0, 1]; % tloc
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
colors = [c1;c2;c3];

sm = 30; % size of marker (asteriks)
lw = 1.5; % linewidth

names = {'concentration', 'laserpower', 'pinhole'};
for j = 1 : numel(names)
    name = names{j};
    fprintf('work on condition %s\n', name);
    
    % read results file
    file = [base_folder, filesep, name, '.xlsx'];
    assert(isfile(file), 'results file missing, run analyse_measurement_series.m');
    info = readtable(file);
    
    if j == 2 % laser
        var_label = ['laser power (',char(181),'W)']; % xlabel
        
        % plotting parameters: l & p are position variables of legend
        l2 = ['NorthEast'];
        l3 = l2;
        p1 = [0.62,0.72,0.15,0];
        p3 = p1;
        % laser power step in plot
        step_x =10;
        
        var_x = info.MinfluxLaserpower;
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
        
        other_conditions = 'LP, c = 2nM, ph = 0.425AU';
        
    elseif j == 3 % pinhole
        var_label = 'pinhole diameter (AU)'; % xlabel
        % plotting parameters: l & p are position variables of legend
        l2 = ['NorthWest'];
        l3 = ['NorthEast'];
        p1 =[0.315,0.72,0.15,0];
        p3 = [0.62,0.72,0.15,0];
        % pinhole diameter step in plot
        step_x =0.15;
        
        var_x = info.MinfluxPinholesizeAU;
        % according to the manufacturer, the actual pinhole size was half of
        % the value displayed in the MINFLUX software (will be corrected after
        % update). To consider correct values had to divide pinhole size by 2.
        disp('Correction of Pinhole');
        var_x = var_x*360/680;
        
        other_conditions = 'PH, c = 2nM, lp = 16%';
        
    elseif j == 1 % concentration
        var_label = 'imager concentration (nM)'; % xlabel
        % plotting parameters: l & p are position variables of legend
        l2 = ['NorthWest'];
        l3 = ['NorthEast'];
        p1 =[0.315,0.72,0.15,0];
        p3 = [0.62,0.72,0.15,0];
        % concentration step in plot
        step_x =[1,0.2];
        
        var_x = info.ConcentrationnM;
        
        other_conditions = 'C, ph = 0.425AU, lp = 16%';
    end
    
    var_date = info.Date;  % date of series
    % all measurements of same series were done on same day
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
    ax = subplot(3,3,1);
    hold on;
    
    % plot mean of measurement series as errorbar
    [me,st]=giveMean(var_x, parameter,info.MedianScanTime);
    gs1 = errorbar(parameter, me, st,'.-','MarkerSize',sm,'Color',cp,'Linewidth',lw, 'DisplayName', 'mean');
    
    %define axis and plot properties and colors
    xbounds = xlim;
    if j == 2 % for laserpower series redefine range
        xbounds = [10 80];
    end
    
    ybounds = ax_t;
    xticks(xbounds(1):step_x(1):xbounds(2));
    yticks(ybounds(1):0.5:ybounds(2));
    ylim(ax_t);
    ax.Color = 'none';
    % scatter plot of each individual series
    gs2 = gscatter(var_x, info.MedianScanTime,var_series(:),colors,'*',sm/5);
    legend('off');
    
    for i = 1 : numel(gs2)
        gs2(i).LineWidth = lw;
    end
    
    % labels
    decorate(var_label, 'median{\it t}_{btw} (s)');
    
    % save calculated mean and std in structure
    structPaper.mean_tbtw = me;
    structPaper.std_tbtw = st;
    
    %% FBG estimated background frequency
    ax = subplot(3,3,2);
    hold on;
    
    % plot mean of measurement series as errorbar
    [me,st]=giveMean(var_x, parameter,info.MedianFBG);
    gs1 = errorbar(parameter, me, st,'.-','MarkerSize',sm,'Color',cp,'Linewidth',lw, 'DisplayName', 'mean');
    
    %define axis and plot properties and colors
    xbounds = xlim;
    if i == 1
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
    gs2 = gscatter(var_x, info.MedianFBG,var_series(:),colors,'*',sm/5);
    legend('off');
    
    for i = 1 : numel(gs2)
        gs2(i).LineWidth = lw;
    end
    
    % labels
    decorate(var_label, 'median {\it f}_{bg} (1/s)');
    
    % save calculated mean and std in structure
    structPaper.mean_fbg = me;
    structPaper.std_fbg = st;
    
    %% CFR center-frequency-ratio
    ax = subplot(3,3,3);
    hold on;
    
    % plot mean of measurement series as errorbar
    [me,st]=giveMean(var_x, parameter,info.MedianCFR);
    gs1 = errorbar(parameter, me, st,'.-','MarkerSize',sm,'Color',cp,'Linewidth',lw, 'DisplayName', 'mean');
    
    %define axis and plot properties and colors
    xbounds = xlim;
    if i == 1
        xbounds = [10 80];
    end
    
    ybounds = ax_c;
    xticks(xbounds(1):step_x(1):xbounds(2));
    yticks(ybounds(1):0.1:ybounds(2));
    ylim( ax_c);
    xlim(xbounds);
    ax.Color = 'none';
    gs2 = gscatter(var_x, info.MedianCFR,var_series(:),colors,'*',sm/5);
    legend('off');
    
    for i = 1 : numel(gs2)
        gs2(i).LineWidth = lw;
    end
    
    % labels
    decorate(var_label, 'median CFR');
    
    % save calculated mean and std in structure
    structPaper.mean_cfr = me;
    structPaper.std_cfr = st;
    
    %% SIGMA R (sqrt)
    ax = subplot(3,3,4);
    hold on;
    
    % plot mean of measurement series as errorbar
    [me,st]=giveMean(var_x, parameter,info.MedianSigmaR);
    gs1 = errorbar(parameter, me, st,'.-','MarkerSize',sm,'Color',cp,'Linewidth',lw, 'DisplayName', 'mean');
    
    %define axis and plot properties and colors
    xbounds = xlim;
    if i == 1
        xbounds = [10 80];
    end
    
    ybounds = ax_s;
    xticks(xbounds(1):step_x(1):xbounds(2));
    yticks(ybounds(1):0.5:ybounds(2));
    ylim( ax_s);
    xlim(xbounds);
    ax.Color = 'none';
    gs2 = gscatter(var_x, info.MedianSigmaR,var_series(:),colors,'*',sm/5);
    legend('off');
    
    for i = 1 : numel(gs2)
        gs2(i).LineWidth = lw;
    end
    
    % labels
    decorate(var_label, 'median {\it \sigma}_r (nm)');
    
    % save calculated mean and std in structure
    structPaper.mean_sigmar = me;
    structPaper.std_sigmar = st;
    
    %% t_loc
    ax = subplot(3,3,5);
    hold on;
    
    % plot mean of measurement series as errorbar
    [me,st]=giveMean(var_x, parameter,info.MedianLocTime);
    gs1 = errorbar(parameter, me, st,'.-','MarkerSize',sm,'Color',cp,'Linewidth',lw, 'DisplayName', 'mean');
    
    %define axis and plot properties and colors
    xbounds = xlim;
    if i == 1 % for laserpower series redefine range
        xbounds = [10 80];
    end
    
    ybounds = ax_t2;
    xticks(xbounds(1):step_x(1):xbounds(2));
    yticks(ybounds(1):0.1:ybounds(2));
    ylim(ax_t2);
    ax.Color = 'none';
    % scatter plot of each individual series
    gs2 = gscatter(var_x, info.MedianLocTime,var_series(:),colors,'*',sm/5);
    legend('off');
    
    for i = 1 : numel(gs2)
        gs2(i).LineWidth = lw;
    end
    
    % labels
    decorate(var_label, 'median{\it t}_{loc} (s)');
    
    %% info box with legend and everything
    ax = subplot(3,3,6);
    ax.Visible = 'off';
    h = [gs2;gs1];
    legend(ax, h, 'box', 'off', 'Fontsize', 12);
    
    text(0.5, 0.4, other_conditions, 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
    
    %% FRC
    ax = subplot(3,3,7);
    hold on;
    
    % plot mean of measurement series as errorbar
    [me,st]=giveMean(var_x, parameter,info.FRCcombined);
    gs1 = errorbar(parameter, me, st,'.-','MarkerSize',sm,'Color',cp,'Linewidth',lw, 'DisplayName', 'mean');
    
    %define axis and plot properties and colors
    xbounds = xlim;
    if i == 1 % for laserpower series redefine range
        xbounds = [10 80];
    end
    
    xticks(xbounds(1):step_x(1):xbounds(2));
    ylim([0, 14]);
    yticks(0:2:14);
    ax.Color = 'none';
    % scatter plot of each individual series
    gs2 = gscatter(var_x, info.FRCcombined,var_series(:),colors,'*',sm/5);
    legend('off');
    
    for i = 1 : numel(gs2)
        gs2(i).LineWidth = lw;
    end
    
    % labels
    decorate(var_label, 'FRC combined (nm)');
    
    %% Loc per Events
    ax = subplot(3,3,8);
    hold on;
    
    % plot mean of measurement series as errorbar
    [me,st]=giveMean(var_x, parameter,info.MeanLocPerEvent);
    gs1 = errorbar(parameter, me, st,'.-','MarkerSize',sm,'Color',cp,'Linewidth',lw, 'DisplayName', 'mean');
    
    %define axis and plot properties and colors
    xbounds = xlim;
    if i == 1 % for laserpower series redefine range
        xbounds = [10 80];
    end
    
    xticks(xbounds(1):step_x(1):xbounds(2));
    ylim([0, 35]);
    yticks(0:5:35);
    ax.Color = 'none';
    % scatter plot of each individual series
    gs2 = gscatter(var_x, info.MeanLocPerEvent,var_series(:),colors,'*',sm/5);
    legend('off');
    
    for i = 1 : numel(gs2)
        gs2(i).LineWidth = lw;
    end
    
    % labels
    decorate(var_label, 'Mean loc./event');
    
    
    %% save figure
    exportgraphics(fig, [file(1:end-5), '.pdf']);
end

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

function decorate(labelx, labely)
%

xlabel(labelx);
ylabel(labely);

grid on;
box on;

ax = gca;
ax.FontSize = 16;
pbaspect(ax,[1 1 1])
end