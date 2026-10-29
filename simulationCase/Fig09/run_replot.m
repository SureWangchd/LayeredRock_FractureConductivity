% Replot Fig. 9 from the compact, existing inputs in this directory.
% No COMSOL connection, new solution, or source CT lookup is required.
% MATLAB R2025b was used for the delivered render.
% For the high-resolution mixed vector/raster PDF, Python needs PyMuPDF.
% Set its interpreter if required:
% setenv('FIG09_PYTHON','C:/path/to/python.exe');
bundle_dir=fileparts(mfilename('fullpath'));
addpath(bundle_dir);
if ~isfolder(fullfile(bundle_dir,'source_notes'))
    mkdir(fullfile(bundle_dir,'source_notes'));
end
figure09_ct_matlab(bundle_dir,bundle_dir);
