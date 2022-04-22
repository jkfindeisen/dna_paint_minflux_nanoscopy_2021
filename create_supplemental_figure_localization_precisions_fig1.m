function create_supplemental_figure_localization_precisions_fig1()
% Computes histograms of the localization precision for all subpanels in
% Fig1
%
% Recreates Suppl. Figure 2
%
% This file is part of the supplementary software for "DNA-PAINT MINFLUX
% nanoscopy", 2021 by Lynn M. Ostersehlt, Daniel C. Jans, Anna Wittek,
% Jan Keller-Findeisen, Steffen J. Sahl, Stefan W. Hell, and Stefan Jakobs

fprintf('calculates localization precision histograms for Minflux data in Fig 1\n');

if ~exist('get_root_folder.m', 'file')
    initialize();
end

% determine file locations
root_folder = get_root_folder();
data_folder = [root_folder, filesep, 'data', filesep, '2D', filesep];
figure_folder = [root_folder, filesep, 'results', filesep, 'precisions'];

% create figure
fig = figure();
fig.Position = [100, 100, 1200, 700];

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
    s = minflux.combined.std_xyz;
    sr = sqrt((s(:, 1).^2 + s(:, 2).^2)/2);
    sr = sr(sr > 0); % remove those where no localization precision could be calculated
    msr = median(sr);
    
    % display histogram
    ax = subplot(2,3,i);
    histogram(sr, 0:0.4:12);
    hold on;
    plot(msr*[1,1], ylim(), 'r-', 'LineWidth', 1.3);
    
    xlim([0, 12]);
    xticks(0:4:12);
    xlabel('\sigma_r (nm)');
    ylabel('occurence');
    title(sprintf('median \\sigma_r=%.2f nm', msr));
    
    grid on;
    box on;
    ax.FontSize = 12;
end

exportgraphics(fig, [figure_folder, filesep, 'Localization_precision_histograms_Fig1.png']);

end