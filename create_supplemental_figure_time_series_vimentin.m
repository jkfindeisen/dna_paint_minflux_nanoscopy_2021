function create_supplemental_figure_time_series_vimentin()
% Takes the 2D Vimentin Minflux measurement and cuts the data to various
% time segments. Then computes the FRC and displays it.
%
% Recreates part of Suppl. Figure 3
%
% This file is part of the supplementary software for "DNA-PAINT MINFLUX
% nanoscopy", 2021 by Lynn M. Ostersehlt, Daniel C. Jans, Anna Wittek,
% Jan Keller-Findeisen, Steffen J. Sahl, Stefan W. Hell, and Stefan Jakobs

fprintf('calculates the FRC (Fourier ring correlation) vs recording time for one data set.\n');

if ~exist('get_root_folder.m', 'file')
    initialize();
end

% determine file locations
root_folder = get_root_folder();
data_file = [root_folder, filesep, 'data', filesep, '2D', filesep, '04_U2OS_Vim-rsEGFP2_het1_NB_A655_2nM_PH0.85_LP14_16.03.2021_Minflux.mat'];
figure_folder = [root_folder, filesep, 'results', filesep, 'frc'];

% load Minflux file
minflux = load(data_file);
minflux = minflux.minflux;

% calculate extended statistics
minflux = calculate_extended_statistics(minflux);

% reasonable boundaries
x = minflux.pos(:, 1); % in nm
y = minflux.pos(:, 2);
a = 0.02;
Rx = quantile(x, [a, 1-a]);
Ry = quantile(y, [a, 1-a]);

% for different times, calculate FRC of combined
ti = (0.5:0.5:8.5).'*3600;
Ni = numel(ti);

frc = zeros(Ni, 1);
sxy = 1; % pixelation of histogram
for i = 1 : Ni
    fprintf('%d/%d\n', i, Ni);
    
    % get all localizations until that time
    ix = minflux.combined.t < ti(i);
    cpos = minflux.combined.pos(ix, :);
    x = cpos(:, 1);
    y = cpos(:, 2);
    
    % split dataset in two random populations (means the FRC can look
    % different every time, maybe we should repeat and average here)
    ix = rand(size(x)) < 0.5;
    h1 = render_xy(x(ix), y(ix), sxy, sxy, Rx, Ry);
    h2 = render_xy(x(~ix), y(~ix), sxy, sxy, Rx, Ry);

    % estimat the FRC
    estimated_resolution = img_fourier_ring_correlation(h1, h2, size(h1)*sxy*1e-9);
    frc(i) = estimated_resolution / 1e-9; % in nm
end

x = ti / 3600;
y = frc;
idx = x < 8.6;
x = x(idx);
y = y(idx);

% display FRC over time
fig = figure();
fig.Position = [100, 100, 400, 400];
plot(x, y, 'bo-', 'LineWidth', 1.5);
xlabel('time (h)');
ylabel('FRC (nm)');
box on;
grid on;
xlim([0, 9]);
ylim([0, 60]);
ax = gca;
ax.FontSize = 12;
exportgraphics(fig, [figure_folder, filesep, 'FRC_vs_recording_time.png']);

end