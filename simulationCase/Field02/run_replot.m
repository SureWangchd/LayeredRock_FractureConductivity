% Replot saved figure data only; no COMSOL connection or solve.
bundle_dir=fileparts(mfilename('fullpath'));
addpath(bundle_dir);
figure_field3d_refined_matlab(bundle_dir,bundle_dir);
