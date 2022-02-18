function minflux = calculate_extended_statistics(minflux)
% Calculates additional statistics to a gather Minflux experiment

assert(nargin == 1);

%% center frequency ratio
minflux.cfr = minflux.efc ./ minflux.efo;

%% ranges
x = minflux.pos(:, 1); % in nm
y = minflux.pos(:, 2);
z = minflux.pos(:, 3);
t = minflux.t;

% reasonable boundaries (will also be used further below)
a = 0.01;
Rx = quantile(x, [a, 1-a]);
Ry = quantile(y, [a, 1-a]);
Rz = quantile(z, [a, 1-a]);
Rt = t([1,end]);

%% drift correction
T = numel(unique(minflux.id))*diff(Rx)*diff(Ry)/3e6; % heuristic for optimal length fo time window
T = min([T, diff(Rt)/2, 3600]); % need at least two time windows
T = max([T, 600]); % but at least 10 minutes long
sxy = 2;
sxyz = 5;
use_gpu = true;

% switch for 2D/3D
is_3D = diff(Rz) > 1e-6;

if is_3D
    % 3D drift correction
    [dx, dy, dz, dxt, dyt, dzt, ti, fig] = drift_correction_time_windows_3D(x, y, z, t, Rx, Ry, Rz, T, sxyz, use_gpu);
    close(fig);
    minflux.dpos = minflux.pos - [dx, dy, dz];
    minflux.drift = [ti(:), dxt(:), dyt(:), dzt(:)];
else
    % 2D drift correction
    [dx, dy, dxt, dyt, ti, fig] = drift_correction_time_windows_2D(x, y, t, Rx, Ry, T, sxy, use_gpu);
    close(fig);
    minflux.dpos = minflux.pos - [dx, dy, zeros(size(dx))];
    minflux.drift = [ti(:), dxt(:), dyt(:)];
end

% will work on minflux.dpos from now one

%% combine localization of each event, compute std deviation in x,y,z
uid = unique(minflux.id);
N = numel(uid);

minflux.std_xyz = zeros(N, 4);
minflux.cpos = zeros(N, 4);
minflux.ct = zeros(N, 1);
for i = 1 : N
    ix = minflux.id == uid(i);
    minflux.std_xyz(i, :) = [std(minflux.dpos(ix, :), [], 1), sum(ix)];
    minflux.cpos(i, :) = [mean(minflux.dpos(ix, :), 1), sum(ix)];
    minflux.ct(i) = mean(minflux.t(ix));
end

%% start and end times of each binding event (first and last time)
minflux.t_event = zeros(N, 2);
for i = 1 : N
    ti = minflux.t(minflux.id == uid(i));
    minflux.t_event(i, :) = [min(ti), max(ti)];
end
minflux.t_scan = minflux.t_event(2:end, 1) - minflux.t_event(1:end-1, 2); % start of next - end of previous
minflux.t_loc = minflux.t_event(:, 2) - minflux.t_event(:, 1); % end - start

%% fourier ring correlation (only in x,y)

% x, y coordinates (in nm)
x = minflux.dpos(:, 1);
y = minflux.dpos(:, 2);
ix = rand(size(x)) < 0.5;

% 2D render of x,y (histogram, 1nm pixel size)
sxy = 0.5;
h1 = render_xy(x(ix), y(ix), sxy, sxy, Rx, Ry);
h2 = render_xy(x(~ix), y(~ix), sxy, sxy, Rx, Ry);

[estimated_resolution, ~, qi, ci] = img_fourier_ring_correlation(h1, h2, size(h1)*sxy*1e-9);
minflux.frc.resolution = estimated_resolution / 1e-9; % in nm
minflux.frc.qi = qi;
minflux.frc.ci = ci;

% now on combined
x = minflux.cpos(:, 1);
y = minflux.cpos(:, 2);
ix = rand(size(x)) < 0.5;
h1 = render_xy(x(ix), y(ix), sxy, sxy, Rx, Ry);
h2 = render_xy(x(~ix), y(~ix), sxy, sxy, Rx, Ry);

[estimated_resolution, ~, qi, ci] = img_fourier_ring_correlation(h1, h2, size(h1)*sxy*1e-9);
minflux.frc_combined.resolution = estimated_resolution / 1e-9; % in nm
minflux.frc_combined.qi = qi;
minflux.frc_combined.ci = ci;
    
end