function figs = display_single_measurement(minflux)
%

% TODO efo, efc, loc/event, combined precision, FRC image

assert(nargin == 1);
figs = [];
labels = {'x', 'y', 'z', 'r'};

%% parameters

% plot information
lw = 2; % linewidth of median line

% Colors for plot
cf = [0.3010 0.7450 0.9330]; % histogram face color
ce = [0 0.4470 0.7410]; % histogram edge color
cm = [0.6350 0.0780 0.1840]; % plot median line color

%% 2D xy histogram rendering

% x, y coordinates (in nm)
x = minflux.dpos(:, 1);
y = minflux.dpos(:, 2);

% reasonable boundaries (rounded to next 100nm)
a = 0.02;
R = 100;
Rx = quantile(x, [a, 1-a]);
Rx = [floor(Rx(1)/R)*R, ceil(Rx(2)/R)*R];
Ry = quantile(y, [a, 1-a]);
Ry = [floor(Ry(1)/R)*R, ceil(Ry(2)/R)*R];

% 2D render of x,y (histogram, 1nm pixel size)
sxy = 1;
[h, xh, yh] = render_xy(x, y, sxy, sxy, Rx, Ry);

% display
fig = figure(345);
clf('reset');
fig.WindowState = 'maximized';
im = imagesc(xh, yh, h.');
im.Parent.YDir = 'normal';
axis image;
colormap(hot);
caxis([0, max(h(:))* 0.5]);
hold on;
xlabel('x (nm)');
xlim(Rx);
ylabel('y (nm)');
ylim(Ry);
title(sprintf('%d events recorded in %.1f min', numel(unique(minflux.id)), minflux.t(end)/60));

figs = [figs; fig];

%% interesting histograms
fig = figure(346);
fig.Position = [100, 100, 1500, 1000];
clf('reset');

% background frequency - plot

subplot(3, 4, 1);
hold on;
m = median(minflux.fbg);

g = 0:1e3:3e4;
histogram(minflux.fbg, g, 'FaceColor', cf, 'EdgeColor',ce)
xlim(g([1,end]));

plot(m*[1,1],ylim(),'LineWidth',3,'Color',cm)

decorate('background signal (Hz)', 'occurence', sprintf('median: %.0f',m));


%% CFR center-frequency-ratio - plot
subplot(3, 4, 2);
hold on
m = median(minflux.cfr);

g = 0:0.05:1.5;
histogram(minflux.cfr,g,'FaceColor',	cf, 'EdgeColor',ce)
xlim(g([1,end]));

plot(m*[1,1],ylim(),'LineWidth',lw,'Color',cm)

decorate('center-frequency-ratio', 'occurence', sprintf('median: %.2f', m));

%% time between events (t_scan)
subplot(3, 4, 3);
hold on;

m = median(minflux.t_scan);

g = 0:0.2:10;
histogram(minflux.t_scan,g,'FaceColor',	cf, 'EdgeColor', ce)
xlim(g([1,end]));

plot(m*[1, 1],ylim(),'LineWidth',lw,'Color',cm)
decorate('time elapsed between events t_{scan} (s)', 'occurence', sprintf('median t_{scan}= %.2f s', m));

%% time within events (t_loc)
subplot(3, 4, 4);
hold on;

m = median(minflux.t_loc);

g = 0:0.3:8;
histogram(minflux.t_loc,g,'FaceColor', cf, 'EdgeColor', ce)
xlim(g([1,end]));

plot(m*[1, 1],ylim(),'LineWidth',lw,'Color',cm)
decorate('time elapsed within events t_{loc} (s)', 'occurence', sprintf('median t_{loc}= %.2f s', m));


%% plot sigma x, y, z, r
T = 5;
sigmas = minflux.combined.std_xyz(minflux.combined.n >= T, :); % all events with at least T localizations
sigmas = [sigmas, sqrt(sum(sigmas(:, 1:2).^2, 2)/2)]; % adds sigma_r
med = median(sigmas, 1);

g = 0:0.5:15; % limit of x-range for histogram of sigma_r

for i = 1 : 4
    subplot(3, 4, 4+i);
    histogram(sigmas(:, i),g,'FaceColor',	cf, 'EdgeColor',ce);
    hold on
    plot(med(i)*[1,1],ylim(),'LineWidth',lw,'Color',cm)
    xlim(g([1,end]));
    decorate('\sigma_x (nm)', 'occurence', sprintf('median \\sigma_%s = %.1f nm ', labels{i}, med(i)));
end

%% FRC
subplot(3,4,9);
plot(minflux.frc.qi, minflux.frc.ci);
decorate('k', 'FRC', sprintf('est. res. %.1f nm', minflux.frc.resolution));

subplot(3,4,10);
plot(minflux.frc_combined.qi, minflux.frc_combined.ci);
decorate('k', 'FRC of combined', sprintf('est. res. %.1f nm', minflux.frc_combined.resolution));

figs = [figs; fig];

%% drift
d = minflux.drift;
subplot(3,4,11);
hold on;
for i = 2 : size(d, 2)
    plot(d(:, 1), d(:, i), 'DisplayName', labels{i-1});
end
decorate('time (s)', 'Est. drift (nm)');

end

function decorate(labelx, labely, plot_title)
% some often used  functionality together

if nargin > 2
    title(plot_title);
end

xlabel(labelx);
ylabel(labely);
grid on;
box on;
pbaspect([1 1 1]);
end