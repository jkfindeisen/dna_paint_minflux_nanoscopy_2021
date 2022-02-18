function series2results_SI_Fig123()
%% Script to analyze measurement series of MINFLUX data
% Script imports all measurements belonging to one measurement series and
% analyses each individually. Data for three repetitions of each
% measurement series is provided.
% Provides .xls -sheet with resulting information and plots of steps
% in-between
% Requirements: Files need to be in ordered file structure in a "data" folder.
% All measurements of the same series must be in an according subfolder.
% For each measurement an exported .json file (exported from Abberior
% Imspector Software) and a .txt file with additional information in
% tabular format must be provided. The corresponding .json and .txt file
% must be saved with the same name.
% Example: data\laserpower\...

close all

% determine file locations
root_folder = get_root_folder();

% formats with which figures should be saved
format = {'.png'};

% for measurement series only data recorded for 1h after starting the
% measurement is considered -> filter each file accordingly
timefilter = 1 ; % h

% plot information
bw = 0.25; % binwidth (nm) for histogram of sigma_r
nb = 50; % number of bins for a histogram
lw = 1; % linewidth of median line

% Colors for plot
cf = [0.3010 0.7450 0.9330]; % histogram face color
ce = [0 0.4470 0.7410]; % histogram edge color
cm = [0.6350 0.0780 0.1840]; % plot median line color
cp = [0.6350 0.0780 0.1840];

%% User interaction: Selection of analysis parameter
% select the measurement series which should be analyzed
list = {'Laser power', 'Pinhole diameter', 'Imager concentration'};

indx = listdlg('PromptString',{'Choose a measurement series.',...
    'Only one parameter can be analyzed at a time.',''},...
    'SelectionMode','single','ListString',list);

if isempty(indx)
    error('No measurement series was selected.');
else
    disp(['Starting analysis of: ', list{indx}]);
    if indx == 1
        variable = 'laserpower';
    elseif indx == 2
        variable = 'pinhole';
    elseif indx == 3
        variable = 'concentration';
    else
        error('Wrong measurement series was selected.');
    end
end

%% analysis of each measurement file

% get all text files of the corresponding folder
files = dir(fullfile(root_folder, 'data',variable, '*.mat')); % filenames with .mat
filenames = {files(:).name};

mkdir(fullfile(root_folder,'results',variable));

for j=1:1:numel(filenames)
    close all
    
    % display filename
    fname = filenames{j};
    disp(['### Start with: ',fname, '. (', num2str(j), '/', num2str(numel(filenames)), ')' ]);
    
    %% read MINFLUX measurement file
    disp('### Loading file... ###')
    disp(fname)
    minflux = load(fullfile(root_folder, 'data',variable, fname));
    minflux = minflux.minflux; % positions are already in nm
    minflux.cfr = minflux.efc ./ minflux.efo;
    
    % filter on time passed after starting measurement
    if timefilter > 0
        filterval = timefilter*60*60; % in s
        disp(['Filter on time at ', num2str(timefilter), ' h = ', num2str(filterval), ' s.']);
        minflux = filter_time(minflux, filterval);
    end
    
    %% calculate some interesting properties
    
    % time between events (t_btw)
    consider_ids = 100; % all measurements had at least 100 ids, for t_btw filter on first 100 molecules
    
    % get index of localization in unfiltered data set
    orig_idx = 1:length(minflux.id);
    id_act = 0; % current id
    
    on_list = zeros(numel(minflux.id),1);
    off_list2 = zeros(numel(minflux.id),1); % for considering the end as the start of the next loc
    
    % run through all ids of valid localizations
    for i=1:1:numel(minflux.id)
        % if trace id is differs from the previous one...
        if minflux.id(i) > id_act
            id_act = minflux.id(i); % update current id
            if i>1 %not consider very first localization
                idx_end = orig_idx(i-1); % get index of last valid localization of previous molecule
                off_list2(i-1) = minflux.t(idx_end+1); % get starting point of first non-valid localization after previous molecule
                % settings of MINFLUX microscope ensure that at least 2 non-valid localizations follow a valid one
            end
            on_list(i) = minflux.t(i);
        end
    end
    
    idx_end = orig_idx(length(minflux.t));
    off_list2(end) = minflux.t(idx_end); %get last time point
    % remove unrelevant data points (just determine as many time differences as there are molecules)
    on_list = on_list(on_list>0);
    off_list2 = off_list2(off_list2>0);
    
    off_time2=zeros(size(off_list2,1),1);
    for i =1:1:numel(on_list)
        if i ==1
            off_time2(i) = on_list(i);
        else
            off_time2(i) = on_list(i)-off_list2(i-1);
        end
    end
    
    % only considering first 100 localizations for off time
    off_time_first2 = off_time2(1:consider_ids);
    
    % localization precision x,y,z,r
    uid = unique(minflux.id);
    
    prec_calc.mean_x=zeros(length(uid),1);
    prec_calc.mean_y=zeros(length(uid),1);
    prec_calc.mean_z=zeros(length(uid),1);
    
    prec_calc.std_x=zeros(length(uid),1);
    prec_calc.std_y=zeros(length(uid),1);
    prec_calc.std_z=zeros(length(uid),1);
    
    prec_calc.nbID = zeros(length(uid),1);
    
    % run through all trace ids
    for ttid = 1:length(uid)
        act_locs = minflux.pos(minflux.id == uid(ttid),:);
        prec_calc.nbID(ttid) = length(act_locs);
        %prec_calc.mean_x(ttid) = mean(act_locs(:,1));
        prec_calc.std_x(ttid) = std(act_locs(:,1));
        %prec_calc.mean_y(ttid) = mean(act_locs(:,2));
        prec_calc.std_y(ttid) = std(act_locs(:,2));
        %prec_calc.mean_z(ttid) = mean(act_locs(:,3));
        prec_calc.std_z(ttid) = std(act_locs(:,3));
    end
    min_sample = 4; % need at least 5 localizations
    
    % filter: need at least 5 localizations
    prec_calc.filter_std_x = prec_calc.std_x(prec_calc.std_x>0 & prec_calc.nbID > min_sample);
    prec_calc.filter_std_y = prec_calc.std_y(prec_calc.std_y>0 & prec_calc.nbID > min_sample);
    prec_calc.filter_std_z = prec_calc.std_z(prec_calc.std_z>0 & prec_calc.nbID > min_sample);
    prec_calc.filter_std_r = sqrt((prec_calc.filter_std_x.^2 + prec_calc.filter_std_y.^2)./2);
    
    median_x = median(median(prec_calc.filter_std_x));
    median_y = median(median(prec_calc.filter_std_y));
    median_z = median(median(prec_calc.filter_std_z));
    median_r = median(median(prec_calc.filter_std_r));
    
    %% background frequency - plot
    fig = figure();
    fig.Position = [100, 100, 1600, 1100];
    
    subplot(3, 4, 1);
    hold on;
    m = median(minflux.fbg);
    
    g = 0:1e3:3e4;
    histogram(minflux.fbg, g, 'FaceColor', cf, 'EdgeColor',ce)
    xlim(g([1,end]));
    
    plot(m*[1,1],ylim(),'LineWidth',3,'Color',cm)
    title(sprintf('median: %.0f',m));
    xlabel('background signal (Hz)');
    ylabel('occurence');
    grid on;
    box on;
    
    %% CFR center-frequency-ratio - plot
    subplot(3, 4, 2);
    hold on
    m = median(minflux.cfr);
    
    g = 0:0.05:1.5;
    histogram(minflux.cfr,g,'FaceColor',	cf, 'EdgeColor',ce)
    xlim(g([1,end]));
    
    plot(m*[1,1],ylim(),'LineWidth',3,'Color',cm)
    title(sprintf('median: %.2f', m));
    xlabel('center-frequency-ratio');
    ylabel('occurence');
    grid on;
    box on;
    
    %% time between events (t_btw)
    subplot(3, 4, 3);
    hold on;
    m = median(off_time_first2);
    
    g = 0:0.2:10;
    histogram(off_time_first2,g,'FaceColor',	cf, 'EdgeColor',ce)
    xlim(g([1,end]));
    
    plot(m*[1, 1],ylim(),'LineWidth',3,'Color',cm)
    title(sprintf('median {\\it t}_{btw}= %.2f s', m));
    xlabel('time elapsed between events {\it t}_{btw} (s)');
    ylabel('occurence');
    grid on;
    box on;
    
    %% plot sigma xyz
    g = 0:0.5:20; % limit of x-range for histogram of sigma_r
    subplot(3, 4, 4);
    histogram(prec_calc.filter_std_x,g,'FaceColor',	cf, 'EdgeColor',ce);
    hold on
    plot(median_x*[1,1],ylim(),'LineWidth',lw,'Color',cm)
    title(sprintf('median \\sigma_x = %.1f nm ', median_x));
    xlim(g([1,end]));
    
    xlabel('\sigma_x (nm)');
    ylabel('occurence')
    pbaspect([1 1 1])
    grid on;
    box on;
    
    subplot(3, 4, 5);
    histogram(prec_calc.filter_std_y,g,'FaceColor',	cf, 'EdgeColor',ce);
    hold on
    plot(median_y*[1,1],ylim(),'LineWidth',lw,'Color',cm)
    title(sprintf('median \\sigma_y = %.1f nm ', median_y));
    xlim(g([1,end]));
    xlabel('\sigma_y (nm)');
    ylabel('occurence')
    pbaspect([1 1 1])
    grid on;
    box on;
    
    subplot(3, 4, 6);
    histogram(prec_calc.filter_std_z,g,'FaceColor',	cf, 'EdgeColor',ce);
    hold on
    plot(median_z*[1,1],ylim(),'LineWidth',lw,'Color',cm)
    title(sprintf('median \\sigma_z = %.1f nm ', median_z));
    xlim(g([1,end]));
    xlabel('\sigma_z (nm)');
    ylabel('occurence')
    pbaspect([1 1 1])
    grid on;
    box on;
    
    %% plot sigma r
    subplot(3, 4, 7);
    histogram(prec_calc.filter_std_r,g,'FaceColor',	cf, 'EdgeColor',ce);
    hold on
    plot(median_r*[1,1],ylim(),'LineWidth',lw,'Color',cm)
    title(sprintf('median \\sigma_r = %.1f nm ', median_r));
    xlim(g([1,end]));
    xlabel('\sigma_r (nm)');
    ylabel('occurence')
    pbaspect([1 1 1])
    grid on;
    box on;
    
    %% Save parameters in structure
    fields = fieldnames(minflux.meta);
    for g = 1:1:numel(fields)
        field = fields{g};
        structInfo(j).(field)= minflux.meta.(field); % structure with info of all files of measurement series combined
    end
    structInfo(j).MedianFBG = median(minflux.fbg);
    structInfo(j).MedianCFR = median(minflux.cfr);
    structInfo(j).MedianOffTimeBtwEv_First2 = median(off_time_first2);
    structInfo(j).considerIDs = consider_ids;
    structInfo(j).MedianSigmaR = median_r;
    structInfo(j).MedianSigmaZ = median_z;
    structInfo(j).MedianSigmaX = median_x;
    structInfo(j).MedianSigmaY = median_y;
    structInfo(j).TimeFilter = timefilter;
    
    %% last message
    disp(['### DONE WITH :',fname]);
    
end

%% Save analysis results in combinded table file
T_all = struct2table(structInfo);
disp(T_all)
disp('Write results in combined file');

table_path = fullfile(root_folder,'results',join([variable,'.xls']));
writetable(T_all, table_path);

end

function m = filter_time(m, t)

ix = m.t < t;
fields = {'t', 'id', 'pos', 'rpos', 'dt', 'efo', 'efc', 'cfr', 'fbg'};
for i = 1 : numel(fields)
    field = fields{i};
    assert(isfield(m, field));
    x = m.(field);
    m.(field) = x(ix, :);
end

end