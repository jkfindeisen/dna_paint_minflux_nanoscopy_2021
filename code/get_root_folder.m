function folder = get_root_folder()
% Returns the folder, one folder up (which is the root folder of the
% project)
%
% This file is part of the supplementary software for
% Optimal Precision and Accuracy in 4Pi-STORM using Dynamic Spline PSF Models
% by Mark Bates, Jan Keller-Findeisen, Adrian Przybylski et al.

folder = fileparts(mfilename('fullpath'));
parts = strsplit(folder, filesep);
parts = parts(1:end-1);
folder = strjoin(parts, filesep);

end