function create_supplemental_figure_localization_precisions_fig1()
% Computes histograms of the localization precision for all subpanels in
% Fig1
%
% Recreates Suppl. Figure 2
%
% This file is part of the supplementary software for "DNA-PAINT MINFLUX
% nanoscopy", 2021 by Lynn M. Ostersehlt, Daniel C. Jans, Anna Wittek,
% Jan Keller-Findeisen, Steffen J. Sahl, Stefan W. Hell, and Stefan Jakobs

fprintf('calculates localization precision histograms for Minflux data in Fig 1&2\n');

if ~exist('get_root_folder.m', 'file')
    initialize();
end

% determine file locations
root_folder = get_root_folder();
figure_folder = [root_folder, filesep, 'results', filesep, 'precisions'];


% create figure
fig = figure();
fig.Position = [100, 100, 1200, 1400];

% 2D data from Fig1
data_folder = [root_folder, filesep, 'data', filesep, '2D', filesep];
T = 5;

% iterate over all subpanels
panels = {'a', 'b', 'c', 'd', 'e', 'f'};
for i = 1 : numel(panels)
    panel = panels{i};
    
    % load Minflux file
    files = dir([data_folder, sprintf('Fig1%s*_Minflux.mat', panel)]);
    assert(numel(files) == 1);
    file = [data_folder, files(1).name];
    minflux = load(file);
    minflux = minflux.minflux;
    
    % calculate extended statistics
    minflux = calculate_extended_statistics(minflux);
    
    % sigma_r (see Supplement for definitions)
    s = minflux.combined.std_xyz(minflux.combined.n >= T, :);
    sr = sqrt((s(:, 1).^2 + s(:, 2).^2)/2);
    med = median(sr);
    
    % display histogram
    ax = subplot(4,3,i);
    histogram(sr, 0:0.4:12);
    hold on;
    plot(med*[1,1], ylim(), 'r-', 'LineWidth', 1.3);
    
    xlim([0, 12]);
    xticks(0:4:12);
    xlabel('\sigma_r (nm)');
    ylabel('occurence');
    title(sprintf('median \\sigma_r=%.1f nm', med));
    
    grid on;
    box on;
    ax.FontSize = 12;
end

% 3D data from Fig2
data_folder = [root_folder, filesep, 'data', filesep, '3D', filesep];
names = {'GFP', 'Mic60', 'ATP5B'};
for i = 1 : numel(names)
    
    % load Minflux file
    file_name = sprintf('Fig2_U2OS_Tom70-Dreiklang_%s_AB_Minflux3D.mat', names{i});
    file = [data_folder, file_name];
    minflux = load(file);
    minflux = minflux.minflux;    
    
    % calculate extended statistics
    minflux = calculate_extended_statistics(minflux);    
    
    % sigma_r, sigma_z (see Supplement for definitions)
    s = minflux.combined.std_xyz(minflux.combined.n >= T, :);
    sr = sqrt((s(:, 1).^2 + s(:, 2).^2)/2);    
    med = median(sr);
    sz = s(:, 3);
    medz = median(sz);
    
    % display sr histogram
    ax = subplot(4,3,6+i);
    histogram(sr, 0:0.4:12);
    hold on;
    plot(med*[1,1], ylim(), 'r-', 'LineWidth', 1.3);
    
    xlim([0, 12]);
    xticks(0:4:12);
    xlabel('\sigma_r (nm)');
    ylabel('occurence');
    title(sprintf('median \\sigma_r=%.1f nm', med));
    
    grid on;
    box on;
    ax.FontSize = 12;
    
    % display sz histogram
    ax = subplot(4,3,9+i);
    histogram(sz, 0:0.4:12);
    hold on;
    plot(medz*[1,1], ylim(), 'r-', 'LineWidth', 1.3);
    
    xlim([0, 12]);
    xticks(0:4:12);
    xlabel('\sigma_z (nm)');
    ylabel('occurence');
    title(sprintf('median \\sigma_z=%.1f nm', medz));
    
    grid on;
    box on;
    ax.FontSize = 12;    
end

exportgraphics(fig, [figure_folder, filesep, 'Localization_precision_histograms_Fig12.png']);

end