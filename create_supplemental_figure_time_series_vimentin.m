function create_supplemental_figure_time_series_vimentin()
%

profile off;
profile on;

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
a = 0.01;
Rx = quantile(x, [a, 1-a]);
Ry = quantile(y, [a, 1-a]);

% for different times, calculate FRC of combined
ti = (0.5:0.5:10).'*3600;
Ni = numel(ti);

frc = zeros(Ni, 1);
n = zeros(Ni, 1);

sxy = 0.5;
for i = 1 : Ni
    ix = minflux.ct < ti(i);
    cpos = minflux.cpos(ix, :);
    
    x = cpos(:, 1);
    y = cpos(:, 2);
    n(i) = numel(x);
    ix = rand(size(x)) < 0.5;
    h1 = render_xy(x(ix), y(ix), sxy, sxy, Rx, Ry);
    h2 = render_xy(x(~ix), y(~ix), sxy, sxy, Rx, Ry);
    
    estimated_resolution = img_fourier_ring_correlation(h1, h2, size(h1)*sxy*1e-9);
    frc(i) = estimated_resolution / 1e-9; % in nm
    
end


x = ti / 3600;
y = frc;
fig = figure();
fig.Position = [100, 100, 600, 600];
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
title(sprintf('lb = %.1f nm', p(3)));
end