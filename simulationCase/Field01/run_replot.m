% Replot the bundled, existing COMSOL fields. No COMSOL solve is launched.
bundle_dir=fileparts(mfilename('fullpath'));
addpath(bundle_dir);
figure_field_sequence_matlab(bundle_dir,bundle_dir);
