# Software package and datasets for "DNA-PAINT MINFLUX nanoscopy"

by Lynn M. Ostersehlt (1), Daniel C. Jans (1,2), Anna Wittek (1,2), Jan Keller-Findeisen (1), Steffen J. Sahl (1), 
Stefan W. Hell (1,3,4), and Stefan Jakobs (1,2,4,5)

1 Max Planck Institute for Multidisciplinary Sciences, Department of NanoBiophotonics, Göttingen, Germany
2 University Medical Center Göttingen, Department of Neurology, Göttingen, Germany
3 Max Planck Institute for Medical Research, Department of Optical Nanoscopy, Heidelberg, Germany
4 Cluster of Excellence "Multiscale Bioimaging: from Molecular Machines to Networks of Excitable Cells" (MBExC),
  University of Göttingen, Germany
5 Fraunhofer Institute for Translational Medicine and Pharmacology ITMP, Translational Neuroinflammation and
  Automated Microscopy, Göttingen, Germany 


1. System requirements and licenses

This is the software package containing raw localization data and analysis scripts in Matlab
(https://www.mathworks.com/products/matlab.html). It has been tested with Matlab R2020b.

The source code of the software is licensed under the MIT license (see file LICENSE).

The data is provided for academic and visualization purposes only. Commercial usage of the provided DNA-Paint Minflux
nanoscopy data (reproduction outside of the publication, etc. ) is forbidden. Please contact the authors for further inquiries.


2. Installation guide

No specific installation is necessary. Execute "initialize.m" at least once before running other scripts to include all paths.


3. Demonstration and usage instructions

Make sure to run initialize.m at least once. 

Call "analyse_single_measurements.m"
  Creates a rendering and a figure containing key data (localization precision, brightness, FRC) for each single measurement
  used in the main text figures and for all measurements used in the performance analysis in the supplement of the publication.

Call "analyse_measurement_series.m"
  Creates an Excel file that contains a table of relevant information for all measurement series (laser power, pinhole, ...).
  
Call "create_supplemental_figures_minflux_performance.m"
  Recreates Supl. Note figues, performance measures for aggregated measurement series  (laser power, pinhole, ...).
  
Call "create_supplemental_figure_time_series_vimentin.m"
  Recreates FRC measurements in Supl. Fig. 3, the time series of the Vimentin measurement.
  
Call "create_supplemental_figures_cfr_simulations.m"
  Recreates the CFR simulation/estimation from Suppl. Note figure III



4. Data description

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

