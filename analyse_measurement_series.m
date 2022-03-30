function analyse_measurement_series()
% Analyses the Minflux performance measurement series
% (varying laser power, pinhole, ...). Uses localization
% data and outputs an Excel file containing a table of information about
% the results for the various conditions.
%
% Recreates Figures of Suppl. Note III (not exactly).
%
% This file is part of the supplementary software for "DNA-PAINT MINFLUX
% nanoscopy", 2021 by Lynn M. Ostersehlt, Daniel C. Jans, Anna Wittek,
% Jan Keller-Findeisen, Steffen J. Sahl, Stefan W. Hell, and Stefan Jakobs

fprintf('analyse Minflux performance measurement series\n');

if ~exist('get_root_folder.m', 'file')
    initialize();
end

% determine file locations
root_folder = get_root_folder();
data_folder = [root_folder, filesep, 'data'];
output_folder = [root_folder, filesep, 'results', filesep, 'performance'];
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% parameter
Tmax = 1 * 60 * 60; % time filter [s]
Nmax = 100; % consider only that many molecules for estiamating of scan/loc times

folder_names = {'concentration', 'laserpower', 'pinhole'};
for i = 1 : numel(folder_names)
    folder_name = folder_names{i};
    fprintf('work on condition: %s\n', folder_name);
    folder = [data_folder, filesep, folder_name];
    
    files = dir([folder, filesep, '*.mat']);
    Nf = numel(files);
    
    info = [];
    for j = 1 : Nf
        file_name = files(j).name;
        fprintf(' work on file %s (%d/%d)\n', file_name, j, Nf);
        file = [folder, filesep, file_name];
        
        % load Minflux file
        minflux = load(file);
        minflux = minflux.minflux;
        
        % filter in time
        minflux = minflux_time_filter(minflux, Tmax);
        
        % calculate extended statistics
        minflux = calculate_extended_statistics(minflux);
        
        % sigma computation
        T = 5;
        sigmas = minflux.combined.std_xyz(minflux.combined.n >= T, :); % all events with at least T localizations
        sigmas = [sigmas, sqrt(sum(sigmas(:, 1:2).^2, 2)/2)]; % adds sigma_r
        
        % save properties to structure
        fields = fieldnames(minflux.meta);
        for g = 1:1:numel(fields)
            field = fields{g};
            info(j).(field)= minflux.meta.(field); % structure with info of all files of measurement series combined
        end
        info(j).MedianFBG = median(minflux.fbg);
        info(j).MedianCFR = median(minflux.cfr);
        info(j).MedianScanTime = median(minflux.t_scan(1:min(end, Nmax)));
        info(j).MedianLocTime = median(minflux.t_loc);
        info(j).MedianSigmaX = median(sigmas(:, 1));
        info(j).MedianSigmaY = median(sigmas(:, 2));
        info(j).MedianSigmaZ = median(sigmas(:, 3));
        info(j).MedianSigmaR = median(sigmas(:, 4));
        info(j).FRCcombined = minflux.frc_combined.resolution;
        info(j).MeanLocPerEvent = mean(minflux.combined.n);
        info(j).TimeFilter = Tmax;
    end
    
    % save analysis results in combined table Excel file
    info = struct2table(info);
    disp(info)
    
    table_path = [output_folder, filesep, folder_name, '.xlsx'];
    writetable(info, table_path);
    
end

end

function m = minflux_time_filter(m, tmax)
% t - maximal time in s

ix = m.t < tmax;
fields = {'t', 'id', 'pos', 'rpos', 'dt', 'efo', 'efc', 'fbg'};
for i = 1 : numel(fields)
    field = fields{i};
    assert(isfield(m, field));
    x = m.(field);
    m.(field) = x(ix, :);
end

end