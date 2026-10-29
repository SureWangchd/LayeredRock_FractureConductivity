% Replot the figure from the saved inputs in this bundle; no COMSOL needed.
bundle_dir=fileparts(mfilename('fullpath'));
addpath(bundle_dir);
if ~isfolder(fullfile(bundle_dir,'source_notes')),mkdir(fullfile(bundle_dir,'source_notes'));end
figure08_finite3d_matlab(bundle_dir,bundle_dir);
