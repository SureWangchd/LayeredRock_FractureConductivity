% Replot saved figure data only; no COMSOL connection or solve.
bundle_dir=fileparts(mfilename('fullpath'));
addpath(fullfile(bundle_dir,'scripts'));
figure02_intervention_acd_matlab(bundle_dir,bundle_dir);
