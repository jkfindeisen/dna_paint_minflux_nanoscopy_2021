function analyse_single_measurements()
% Analyses all contained single measurements (data for main text figures
% and data for measurement series for supplementary figures) and outputs a
% standardized analysis (localization precision, crude visualization) for
% every dataset.
%
% This file is part of the supplementary software for "DNA-PAINT MINFLUX
% nanoscopy", 2021 by Lynn M. Ostersehlt, Daniel C. Jans, Anna Wittek,
% Jan Keller-Findeisen, Steffen J. Sahl, Stefan W. Hell, and Stefan Jakobs

fprintf('analyse all single measurements\n');

if ~exist('get_root_folder.m', 'file')
    initialize();
end

% determine file locations
root_folder = get_root_folder();
data_folder = [root_folder, filesep, 'data'];
output_folder = [root_folder, filesep, 'results', filesep, 'single-dataset'];
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% get all subfolders of data_folder
folders = dir(data_folder);
for i = 1 : numel(folders)
    folder = folders(i);
    if ~folder.isdir || strcmp(folder.name, '.') || strcmp(folder.name, '..') || strcmp(folder.name, 'psf')
        continue; % need to ignore the PSF folder
    end
    
    folder_name = folder.name;
    fprintf('work on folder %s\n', folder_name);
    folder = [data_folder, filesep, folder_name];
    
    files = dir([folder, filesep, '*.mat']);
    Nf = numel(files);
    
    for j = 1 : Nf
        file_name = files(j).name;
        fprintf(' work on file %s (%d/%d)\n', file_name, j, Nf);
        file = [folder, filesep, file_name];
        
        % load Minflux file
        minflux = load(file);
        minflux = minflux.minflux;
        
        % calculate extended statistics
        minflux = calculate_extended_statistics(minflux);
        
        % show as figure
        figs = display_single_measurement(minflux);
        
        % store images
        out_folder = [output_folder, filesep, folder_name];
        if ~exist(out_folder, 'dir')
            mkdir(out_folder);
        end
        exportgraphics(figs(1), [out_folder, filesep, file_name(1:end-4), '.render-xy.png']);
        exportgraphics(figs(2), [out_folder, filesep, file_name(1:end-4), '.properties.png']);
        
    end
end

end