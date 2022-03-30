function create_supplemental_figures_cfr_simulations()
% Using precomputed excitation and detection PSFs, can estimate the CFR
% (center frequency ratio) for various conditions including pinhole size
% and image concentration.
%
% Recreates Suppl. Note figure III.
%
% This file is part of the supplementary software for "DNA-PAINT MINFLUX
% nanoscopy", 2021 by Lynn M. Ostersehlt, Daniel C. Jans, Anna Wittek,
% Jan Keller-Findeisen, Steffen J. Sahl, Stefan W. Hell, and Stefan Jakobs

close all;

if ~exist('get_root_folder.m', 'file')
    initialize();
end

disp('simulation/estimation of the CFR (center frequency ratio) vs. laser power, pinhole, image concentration');
disp('reproduces Suppl. Note Figure III');

% determine file locations
root_folder = get_root_folder();
psf_folder = [root_folder, filesep, 'data', filesep, 'psf'];
figure_folder = [root_folder, filesep, 'results', filesep, 'simulation'];

% import excitation & detection psf collections
load([psf_folder, filesep, 'excitation.mat']);
load([psf_folder, filesep, 'detection.mat']);

% define ranges (x-variable) for simulation and L
lp = [exc.lp];
ph_au = [det.ph_au]; % au
c = (1:10)*1e-9; % mol
L = [300 150 75 40]; % nm

% define figure information
fontsize = 14; %fontsize
colors(1,:) = [0.2 0.2 0.2];
colors(2,:) = [0.07 0.62 1];
colors(3,:) = [0 1 0];
colors(4,:) = [0.8 0.8 0.8];

%% get default value
default_c = 2*1e-9; % mol
idx_lp = find(lp==max(lp)); % use maximal laser power
default_exc = exc(idx_lp).h;
idx_ph = find(ph_au == 0.85*360/680);
default_det = det(idx_ph).h;

%% calculate CFR laser power
disp('calculate CFR for varying laser powers');
CFR_lp = zeros(length(L),length(lp));

for l=1:1:length(L)
    for po = 1:1:length(lp)
        var.L = L(l);
        var.exc = exc(po).h;
        cfr_lp_L = CFR(var.exc, default_det, default_c, var.L);
        CFR_lp(l,po)=cfr_lp_L;
    end
end

%% plot CFR laser power for different L
fig = figure();
fig.Position = [100, 100, 1400, 350];
ax = subplot(1,3,1);
hold on

for k = 1:size(CFR_lp,1)
    plot(1e6*lp, CFR_lp(k,:),'linewidth',3,'color', colors(k,:))
end
leg=legend(string(L), 'location','eastoutside');
title(leg,'{\it L} (nm)')
leg.Title.Visible = 'on';
xlabel(' laser power (µW)');
ylabel('CFR');
ylim([0 0.8]);
ax.FontSize = fontsize;

%% calculate CFR pinhole size
disp('calculate CFR for varying pinhole sizes');
CFR_ph = zeros(length(L),length(ph_au));
for l=1:1:length(L)
    for p = 1:1:length(ph_au)
        var.L = L(l);
        var.det = det(p).h;
        cfr_ph_L = CFR(default_exc,var.det,default_c,var.L);
        CFR_ph(l,p)=cfr_ph_L;
    end
end

%% plot CFR pinhole for different L
ax = subplot(1,3,2);
hold on

for k = 1:size(CFR_ph,1)
    plot(ph_au, CFR_ph(k,:),'linewidth',3,'color', colors(k,:))
end
leg=legend(string(L),'location','eastoutside');
title(leg,'{\it L} (nm)')
leg.Title.Visible = 'on';
xlabel('{\it d}_{pinhole} (AU)');
xlim([0.2, 0.8]);
xticks(0.2:0.2:0.8);
ylabel('CFR');
ax.FontSize = fontsize;

%% calculate CFR concentration
disp('calculate CFR for varying imager concentrations');
CFR_conc = zeros(length(L),length(c));

for cc = 1:1:length(c)
    for g=1:1:length(L)
        var.L = L(g);
        var.c = c(cc);
        cfr_L = CFR(default_exc,default_det,var.c,var.L);
        CFR_conc(g,cc)=cfr_L;
    end
end

%% plot CFR concentration for different L
ax = subplot(1,3,3);
hold on

for k = 1:size(CFR_conc,1)
    plot(c*1e9, CFR_conc(k,:),'linewidth',3,'color', colors(k,:))
end
leg=legend(string(L),'location','eastoutside');
title(leg,'{\it L} (nm)')
leg.Title.Visible = 'on';
xlabel('imager concentration (nM)');
ylabel('CFR');
xlim([0 11]);
ax.FontSize = fontsize;

exportgraphics(fig, [figure_folder, filesep, 'CFR_simulations.png']);

% end message
disp('all calculations done');

end

function cfr = CFR(h_exc, h_det, c, L)
%% CFR calculation (with pre-calculated excitation and detection PSFs)
% h_exc: PSF of excitation
% h_det: PSF of confocal detection
% c: imager strand concentration [mol]
% L: Targeted coordinate pattern diameter [nm]

h_exc = single(h_exc);
% h_exc = h_exc(41:121, 41:121, 71:211);
h_det = single(h_det);
% h_det = h_det(41:121, 41:121, 71:211);

%% sets default information
n_a = 6.022140857e23; %1/mol

%general
p.L = L;
p.step = 10; % according to step-size of functions [nm]
p.thresh = 0.1; % need to be optimized
p.power = 1;
p.bgc = 1; % percentage of background considered
p.B = 150; % brightness (photons) of the dye (here the CFR is independent
           % of the brightness because both signal and bg scale with it

% region
p.xlim = 400;
p.ylim = p.xlim;
p.zlim = 700;

%% define volume of one pixel and # imager strands
voxel_nm = p.step^3; % nm^3
voxel_l = voxel_nm * 1e-24; % liter
n_bg_voxel = c .* voxel_l .* n_a ; % number background molecules in voxel

%% get targeted coordinate pattern positions

% get exposure positions of TCP for given L and voxel size
L_2 = p.L/2;
angle = deg2rad(60);
x_hex = L_2 * cos(angle);
y_hex = L_2 * sin(angle);
center_position(1,:) = [0,0,0]; % central position
center_position(2,:) = [L_2,0,0];
center_position(3,:) = [-L_2,0,0];
center_position(4,:) = [x_hex,y_hex,0];
center_position(5,:) = [-x_hex,y_hex,0];
center_position(6,:) = [x_hex,-y_hex,0];
center_position(7,:) = [-x_hex,-y_hex,0];
center_position = round(center_position ./ p.step) * p.step;

x = -p.xlim:p.step:p.xlim;
y = -p.ylim:p.step:p.ylim;
z = -p.zlim:p.step:p.zlim;

idxs = zeros(size(center_position));
% get indexes for positions
for i =1:1:length(center_position)
    [~,idx]= find(x == center_position(i,1));
    idxs(i,1)= idx;
    [~,idx]= find(y == center_position(i,2));
    idxs(i,2)= idx;
    [~,idx]= find(z == center_position(i,3));
    idxs(i,3)= idx;
end

if all(center_position(1,:) == [0, 0 , 0])
    idx_det = idxs(1,:);
else
    error('confocal position not defined!');
end

%% for each exposure (of TCP) get background and intensity
room = 200;
sx = -p.xlim-room:p.step:p.xlim+room;
sy = -p.ylim-room:p.step:p.ylim+room;
sz = -p.zlim-room:p.step:p.zlim+room;
space = zeros(length(sx),length(sy),length(sz));

I_bg_shift = zeros(size(center_position,1),1);
I_sg_shift = zeros(size(center_position,1),1);

% shift confocal detection - to make same size
h_det_shift = shift2([0,0,0], h_det, space, sx, sy, sz);

% get background for each exposure position
for i_cp=1:1:size(center_position,1)
    % current center position
    cp = center_position(i_cp,:);
    
    % shift the exc beam
    h_exc_shift = shift2(cp, h_exc, space, sx, sy, sz);
    
    % get h_eff
    h_eff_shift = h_exc_shift.*h_det_shift;
    
    % calculate background as the sum of the effective PSF
    % assuming homogeneous background
    h_bg_shift = h_eff_shift .* (p.B * n_bg_voxel);
    I_bg_shift(i_cp) = sum(h_bg_shift, 'all');
end


sg_det = h_det(idx_det(1), idx_det(2), idx_det(3));

% get signal for each exposure position
for i_cp=1:1:size(center_position,1)
    % get exc signal at position
    sg_exc = h_exc(idxs(i_cp,1),idxs(i_cp,2),idxs(i_cp,3));
    
    % to get h_eff multiply with central h_det
    sg = sg_exc*sg_det;
    
    I_sg_shift(i_cp) = sg * p.B;
end

%% combine all obtained signal (signal + background)

I_all = I_bg_shift + I_sg_shift;

%% calculate CFR result (center-frequency-ratio)
cfr = I_all(1) ./ mean(I_all(2:end));

end

function shifted = shift2(center,beam, space, x,y,z)
% shifts psf to a defined position in a larger space

sI = size(beam);
%find the right position in the space
ind_space = [find(x == center(1)), find(y == center(2)), find(z == center(3))];
space(ind_space(1)-(sI(1)-1)/2:ind_space(1)+(sI(1)-1)/2,ind_space(2)-(sI(2)-1)/2:ind_space(2)+(sI(2)-1)/2,ind_space(3)-(sI(3)-1)/2:ind_space(3)+(sI(3)-1)/2)=beam;
shifted=space;
end