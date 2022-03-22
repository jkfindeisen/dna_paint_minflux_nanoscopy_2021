# Software package and datasets for "DNA Paint Minflux nanoscopy"

by Lynn M. Ostersehlt1, Daniel C. Jans1,, Anna Wittek, Jan Keller-Findeisen, Steffen J. Sahl, 
Stefan W. Hell, and Stefan Jakobs

to be published

This is the software package containing raw localization data and analysis scripts in Matlab (https://www.mathworks.com/products/matlab.html).

1. Data description

Files *.mat in folder data (subfolder 2D for 2D data and 3D for 3D data). The data files are Matlab "mat" files (based on and compatible with HDF5)
and contain the following fields

t    - time of each event [s]
id   - event id (equal value indicates consecutive localization within a single binding event)
pos  - localization positions X,Y,Z [nm] (for 2D Z = 0)
rpos - position of localization relative to beam center in last iteration X,Y,Z [nm]
dt   - duration of illumination in last iteration [s]
efo  - average frequency of counts in periphery of illumination pattern [Hz]
efc  - average frequency of counts in center of illumination pattern [Hz]
fbg  - est. frequency of counts of the background [Hz]

for each of the localization (so Nx1 or Nx3 vectors)

"meta" struct with meta information (recording data, laser power, pinhole size, used Minflux sequence, used dye and image concentration)

2. Code description

run initialize.m at the beginning

analyse_single_measurements.m 
  Creates a rendering and a figure containing key data (localization precision, brightness, FRC) for each single measurement used in the main text figures.

analyse_measurement_series.m
  Creates an Excel file that contains a table of relevant information for all measurement series (laser power, pinhole, ...).

create_supplemental_figure_time_series_vimentin.m
  Recreates Supl. Fig. X, the time series of the Vimentin measurement.
  
create_supplemental_figures.m
  Recreates Supl. Figs. YY, the measurement series  (laser power, pinhole, ...).

3. License

The code is licensed under the MIT license (see file LICENSE). The data is all rights reserved and only presented for visualization purposes.
