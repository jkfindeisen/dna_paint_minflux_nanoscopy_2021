function minflux = calculate_extended_statistics(minflux)
% Calculates additional statistics on a gather Minflux experiment

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
% minimal boundaries 200nm (because of CR in drift_correction)
if diff(Rx) < 200
    Rx = Rx + (200 - diff(Rx))/2*[-1,1];
end
if diff(Ry) < 200
    Ry = Ry + (200 - diff(Ry))/2*[-1,1];
end
if diff(Rz) < 200 && diff(Rz) > 1e-6 % only if 3D
    Rz = Rz + (200 - diff(Rz))/2*[-1,1];
end
Rt = t([1,end]);

%% drift correction
T = numel(unique(minflux.id))*diff(Rx)*diff(Ry)/3e6; % heuristic for optimal length of time window
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

% compute std-error in x,y and in z
s = c.std_xyz(c.n > 1, :); % ignore only seen once (likely false positives anyway)
n = c.n(c.n > 1);
sr = sqrt(s(:, 1).^2+s(:, 2).^2);
sre = sr ./ sqrt(n); % standard error = standard deviation / sqrt(N)
sz = s(:, 3);
sze = sz ./ sqrt(n);
c.ste_rz = [sre, sze];

minflux.combined = c;

%% start and end times of each binding event (first and last time)
minflux.t_start = accumarray(uid, minflux.t, [], @min);
minflux.t_end = accumarray(uid, minflux.t, [], @max);
minflux.t_scan = minflux.t_start(2:end, 1) - minflux.t_end(1:end-1); % start of next - end of previous
minflux.t_loc = minflux.t_end - minflux.t_start; % end - start

%% fourier ring correlation (only in x,y)
R = 8; % five repetitions to get more stable results

% x, y coordinates (in nm)
x = minflux.dpos(:, 1);
y = minflux.dpos(:, 2);

for r = 1 : R
    % new split
    ix = rand(size(x)) < 0.5;
    
    % 2D render of x,y (histogram, 1nm pixel size)
    sxy = 1;
    h1 = render_xy(x(ix), y(ix), sxy, sxy, Rx, Ry);
    h2 = render_xy(x(~ix), y(~ix), sxy, sxy, Rx, Ry);
    
    [estimated_resolution, ~, qi, ci] = img_fourier_ring_correlation(h1, h2, size(h1)*sxy*1e-9);
    
    minflux.frc.resolution(r) = estimated_resolution / 1e-9; % in nm
    minflux.frc.qi{r} = qi;
    minflux.frc.ci{r} = ci;
    
    % now on combined
    x = minflux.combined.pos(:, 1);
    y = minflux.combined.pos(:, 2);
    ix = rand(size(x)) < 0.5;
    h1 = render_xy(x(ix), y(ix), sxy, sxy, Rx, Ry);
    h2 = render_xy(x(~ix), y(~ix), sxy, sxy, Rx, Ry);
    
    [estimated_resolution, ~, qi, ci] = img_fourier_ring_correlation(h1, h2, size(h1)*sxy*1e-9);
    minflux.frc_combined.resolution(r) = estimated_resolution / 1e-9; % in nm
    minflux.frc_combined.qi{r} = qi;
    minflux.frc_combined.ci{r} = ci;
end

end