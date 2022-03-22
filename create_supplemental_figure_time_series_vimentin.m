function create_supplemental_figure_time_series_vimentin()
% Analyses the measurement series (laser power, pinhole, ...) localization
% data and outputs an Excel file containing a table of information about
% the results for the various conditions. Can be used to create Supl.
% Figures.
%
%
%
% This file is part of the supplementary software for
% "DNA Paint Minflux nanoscopy"  by Lynn M. Ostersehlt, Daniel C. Jans, Anna Wittek,
% Jan Keller-Findeisen, Steffen J. Sahl, Stefan W. Hell, and Stefan Jakobs

fprintf('Analyse measurement series\n');

if ~exist('get_root_folder.m', 'file')
    initialize();
end

% determine file locations
root_folder = get_root_folder();
file = [root_folder, filesep, 'data', filesep, '2D', filesep, '04_U2OS_Vim-rsEGFP2_het1_NB_A655_2nM_PH0.85_LP14_16.03.2021_Minflux.mat'];

% load Minflux file
minflux = load(file);
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
n = zeros(Ni, 1);
% clusters
nc = zeros(Ni, 1);

sxy = 1;
for i = 1 : Ni
% for i = Ni
    fprintf('%d/%d\n', i, Ni);
    ix = minflux.combined.t < ti(i);
    cpos = minflux.combined.pos(ix, :);
    
    x = cpos(:, 1);
    y = cpos(:, 2);
    t = minflux.combined.t(ix);
    n(i) = numel(x);
    
    % do some clustering
%     eps = 1:0.5:10;
%     neps = zeros(numel(eps), 1);
%     for j = 1 : numel(eps)
%         idx = dbscan([x,y], eps(j), 1);
%         neps(j) = idx(end);
%     end
    
    idx = dbscan([x,y], 4, 1);
    nc(i) = numel(unique(idx));
    
    ix = rand(size(x)) < 0.5;
%     ix = t < t(end) / 2;
% ix = (1:n(i)) < n(i)/2;
    h1 = render_xy(x(ix), y(ix), sxy, sxy, Rx, Ry);
    h2 = render_xy(x(~ix), y(~ix), sxy, sxy, Rx, Ry);
    
    estimated_resolution = img_fourier_ring_correlation(h1, h2, size(h1)*sxy*1e-9);
    frc(i) = estimated_resolution / 1e-9; % in nm
end


x = ti / 3600;
y = frc;
idx = x < 8.6;
x = x(idx);
y = y(idx);

% for supplement
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
exportgraphics(fig, 'P:\Shared\Daniel Jans\2D Vimentin MINFLUX\suppl time evolution figure\FRC.pdf');

fig = figure();
fig.Position = [100, 100, 400, 400];
plot(x, nc, 'bo-', 'LineWidth', 1.5);
xlabel('time (h)');
ylabel('number of clusters');
box on;
grid on;
xlim([0, 9]);
ax = gca;
ax.FontSize = 12;
exportgraphics(fig, 'P:\Shared\Daniel Jans\2D Vimentin MINFLUX\suppl time evolution figure\NC.pdf');


fig = figure();
fig.Position = [100, 100, 400, 400];
plot(x, y, 'bo');
xlabel('time (h)');
ylabel('FRC (nm)');
box on;
grid on;
% fit to a*exp(-x/b)+c
model = @(p) p(1)*exp(-x/p(2))+p(3);
minimizer = @(p) mean((y-model(p)).^2);
p0 = [50, 4, 10];
p = fmincon(minimizer, p0);
m = model(p);
hold on;
plot(x, m, 'r-', 'LineWidth', 1.5);
x2 = 8.5:0.1:16;
model = @(p) p(1)*exp(-x2/p(2))+p(3);
m = model(p);
plot(x2, m, 'r:', 'LineWidth', 1.5);
ylim([0, 60]);
tx = p(2)*(log(p(1)) - log(5 - p(3)));
title(sprintf('FRC(t→∞) = %.1f nm, t(FRC=5) =  %.1f h', p(3), tx));
legend('data', 'fit', 'extrapol.');
% exportgraphics(fig, 'P:\FRC_vs_t.png');

fig = figure();
fig.Position = [200, 200, 400, 400];
plot(x, n, 'bo');
xlabel('time (h)');
ylabel('#binding events measured');
box on;
grid on;
% fit with exponential model
model = @(p) p(1)*(1-exp(-x/p(2)));
minimizer = @(p) mean((n-model(p)).^2);
p0 = [50000, 5];
lb = [n(end), 1];
ub = [1e6, 100];
p = fmincon(minimizer, p0, [], [], [], [], lb, ub);
m = model(p);
hold on;
plot(x, m, 'r-', 'LineWidth', 1.5);
model = @(p) p(1)*(1-exp(-x2/p(2)));
m = model(p);
plot(x2, m, 'r:', 'LineWidth', 1.5);
% plot([0, 16], p(1)*[1,1], 'k-');
tx = -log(1-0.9)*p(2);
title(sprintf('speed %.1f Hz, 8.5h %.1f%%, %.1fh 90%%', numel(minflux.combined.n) / minflux.t(end), n(end)/p(1)*100, tx));
legend('data', 'fit', 'extrapol.', 'Location', 'NorthEast');
% exportgraphics(fig, 'P:\N_vs_t.png');

end