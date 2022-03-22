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
[~, ~, uid] = unique(minflux.id);
c.n = accumarray(uid, 1);
N = numel(c.n);
c.t = accumarray(uid, minflux.t) ./ c.n;
c.pos = zeros(N, 3);
c.std_xyz = zeros(N, 3);
for i = 1 : 3
    c.pos(:, i) = accumarray(uid, minflux.dpos(:, i)) ./ c.n;
    c.std_xyz(:, i) = accumarray(uid, minflux.dpos(:, i), [], @std);
end
minflux.combined = c;

%% start and end times of each binding event (first and last time)
minflux.t_start = accumarray(uid, minflux.t, [], @min);
minflux.t_end = accumarray(uid, minflux.t, [], @max);
minflux.t_scan = minflux.t_start(2:end, 1) - minflux.t_end(1:end-1); % start of next - end of previous
minflux.t_loc = minflux.t_end - minflux.t_start; % end - start

%% fourier ring correlation (only in x,y)

% x, y coordinates (in nm)
x = minflux.dpos(:, 1);
y = minflux.dpos(:, 2);
ix = rand(size(x)) < 0.5;

% 2D render of x,y (histogram, 1nm pixel size)
sxy = 1;
h1 = render_xy(x(ix), y(ix), sxy, sxy, Rx, Ry);
h2 = render_xy(x(~ix), y(~ix), sxy, sxy, Rx, Ry);

[estimated_resolution, ~, qi, ci] = img_fourier_ring_correlation(h1, h2, size(h1)*sxy*1e-9);
minflux.frc.resolution = estimated_resolution / 1e-9; % in nm
minflux.frc.qi = qi;
minflux.frc.ci = ci;

% now on combined
x = minflux.combined.pos(:, 1);
y = minflux.combined.pos(:, 2);
ix = rand(size(x)) < 0.5;
h1 = render_xy(x(ix), y(ix), sxy, sxy, Rx, Ry);
h2 = render_xy(x(~ix), y(~ix), sxy, sxy, Rx, Ry);

[estimated_resolution, ~, qi, ci] = img_fourier_ring_correlation(h1, h2, size(h1)*sxy*1e-9);
minflux.frc_combined.resolution = estimated_resolution / 1e-9; % in nm
minflux.frc_combined.qi = qi;
minflux.frc_combined.ci = ci;
  
end