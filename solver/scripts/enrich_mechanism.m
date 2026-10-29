function report = enrich_mechanism(runDirectory, caseIDs)
% ENRICH_MECHANISM Audit saved solutions without rebuilding or solving.
% Usage: enrich_mechanism
%        enrich_mechanism(runDirectory, {'mechanism_cut_mesh3'})
% Writes only derived JSON/CSV files in the selected existing case folders.
% Geometry and mphinterp coordinates are mm; COMSOL components x,y mean x,z.
% No study.run, sol.run, geometry/mesh build, mphsave, or ModelUtil.clear.
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(runDirectory)
    runDirectory = fullfile(root, 'runs', '2026-10-21');
end
if nargin < 2 || isempty(caseIDs)
    caseIDs = {};
    for level = 1:3
        for kind = {'full', 'cut', 'side'}
            caseIDs{end+1} = sprintf('mechanism_%s_mesh%d', kind{1}, level); %#ok<AGROW>
        end
    end
end
if ischar(caseIDs) || isstring(caseIDs), caseIDs = cellstr(caseIDs); end
base = 'D:/Program Files/COMSOL/COMSOL63/Multiphysics';
addpath(fullfile(base, 'mli'));
import com.comsol.model.util.*
ownsConnection = false;
try
    mphstart('127.0.0.1', 2037, base);
    ownsConnection = true;
catch connectError
    % An already-connected MATLAB process may also call this function.
    try
        ModelUtil.tags;
        fprintf('Using the existing LiveLink connection.\n');
    catch
        rethrow(connectError);
    end
end
connectionCleanup = onCleanup(@() disconnect_if_owned(ownsConnection)); %#ok<NASGU>
report = repmat(struct('case_id', '', 'status', '', 'message', ''), numel(caseIDs), 1);
for k = 1:numel(caseIDs)
    id = char(caseIDs{k});
    report(k).case_id = id;
    folder = fullfile(runDirectory, id);
    modelFile = fullfile(folder, [id '.mph']);
    if ~isfile(fullfile(folder, 'SUCCESS.json')) || ~isfile(modelFile)
        report(k).status = 'skipped';
        report(k).message = 'A completed SUCCESS.json and solved MPH are required.';
        fprintf('ENRICH_SKIP %s: %s\n', id, report(k).message);
        continue;
    end
    tag = ['Enrich_' strrep(char(java.util.UUID.randomUUID), '-', '')];
    modelCleanup = onCleanup(@() remove_owned_model(tag));
    fprintf('ENRICH_START %s\n', id);
    started = tic;
    try
        p = jsondecode(fileread(fullfile(folder, 'parameters.json')));
        sel = jsondecode(fileread(fullfile(folder, 'selections.json')));
        model = mphload(modelFile, tag);
        metadata = solver_metadata(model, id, modelFile);
        write_json(fullfile(folder, 'solver_metadata.json'), metadata);
        [metrics, boundaryTable] = geometry_metrics(model, p, sel);
        femResult = jsondecode(fileread(fullfile(folder, 'SUCCESS.json')));
        metrics.cad_areas = annotate_geometry_measurement_quality(metrics.cad_areas, femResult);
        writetable(boundaryTable, fullfile(folder, 'interface_segments.csv'));
        metrics.case_id = id;
        metrics.source_mph = modelFile;
        metrics.postprocessing_only = true;
        metrics.status = 'success';
        metrics.coordinates = 'COMSOL x,y = physical x,z; geometry/input coordinates mm';
        metrics.current_units = 'A/m of out-of-plane thickness';
        try
            metrics.section = partitioned_section_current(model, p, sel);
        catch sectionError
            metrics.section = struct('status', 'failed', 'error', sectionError.message);
            metrics.status = 'partial';
            warning('enrich_mechanism:SectionFailed', '%s: %s', id, sectionError.message);
        end
        metrics.postprocessing_seconds = toc(started);
        write_json(fullfile(folder, 'mechanism_geometry_metrics.json'), metrics);
        report(k).status = metrics.status;
        report(k).message = sprintf('Measured %d fracture components; design check %d.', ...
            metrics.fracture_component_count, metrics.expected_component_count_for_design);
        fprintf('ENRICH_%s %s: %s %.1fs\n', upper(metrics.status), id, report(k).message, toc(started));
    catch caseError
        report(k).status = 'failed';
        report(k).message = getReport(caseError, 'extended', 'hyperlinks', 'off');
        warning('enrich_mechanism:CaseFailed', '%s: %s', id, report(k).message);
    end
    clear modelCleanup model;
end
end

function [m, T] = geometry_metrics(model, p, sel)
geomTag = 'geom1';
allIDs = row_ids(sel.all_domains);
fracIDs = row_ids(sel.fracture_domains);
lowIDs = row_ids(sel.low_domains);
highIDs = row_ids(sel.high_domains);
compIDs = [];
if isfield(sel, 'compensation_domains'), compIDs = row_ids(sel.compensation_domains); end
mainIDs = setdiff(fracIDs, compIDs);
assert(numel(unique([lowIDs highIDs fracIDs])) == numel([lowIDs highIDs fracIDs]), ...
    'Material domain selections overlap.');
assert(isequal(sort([lowIDs highIDs fracIDs]), sort(allIDs)), 'Material domain coverage is incomplete.');
assert(all(ismember(compIDs, fracIDs)), 'Compensation domains must be fracture material.');
% Union the adjacency of each domain; this also retains artificial partition
% boundaries, which are reported but never counted as material interfaces.
boundaryIDs = [];
for d = allIDs
    boundaryIDs = union(boundaryIDs, row_ids(mphgetadj(model, geomTag, 'boundary', 'domain', d)));
end
boundaryIDs = sort(boundaryIDs(:));
nb = numel(boundaryIDs);
length_mm = zeros(nb, 1);
adjacent_domain_count = zeros(nb, 1);
adjacent_domain_ids = cell(nb, 1);
adjacent_materials = cell(nb, 1);
fracture_low_interface = false(nb, 1);
fracture_high_interface = false(nb, 1);
fracture_internal_partition = false(nb, 1);
fracture_on_exterior = false(nb, 1);
main_ellipse_low_interface = false(nb, 1);
main_ellipse_high_interface = false(nb, 1);
compensation_low_interface = false(nb, 1);
compensation_high_interface = false(nb, 1);
adjacency = false(numel(fracIDs));
unexpected = [];
for n = 1:nb
    bnd = boundaryIDs(n);
    adj = row_ids(mphgetadj(model, geomTag, 'domain', 'boundary', bnd));
    assert(all(ismember(adj, allIDs)), 'Boundary references an unclassified domain.');
    L = mphint2(model, '1', 1, 'selection', bnd, 'intorder', 8, 'unit', 'mm', 'solnum', 'end');
    assert(isscalar(L) && isfinite(L) && isreal(L) && L >= 0, 'Invalid boundary length.');
    length_mm(n) = double(L);
    adjacent_domain_count(n) = numel(adj);
    adjacent_domain_ids{n} = strtrim(sprintf('%d ', adj));
    labels = cell(1, numel(adj));
    for q = 1:numel(adj)
        if ismember(adj(q), fracIDs), labels{q} = 'fracture';
        elseif ismember(adj(q), lowIDs), labels{q} = 'low';
        else, labels{q} = 'high'; end
    end
    adjacent_materials{n} = strjoin(labels, ' ');
    af = intersect(adj, fracIDs);
    hasLow = any(ismember(adj, lowIDs));
    hasHigh = any(ismember(adj, highIDs));
    fracture_low_interface(n) = ~isempty(af) && hasLow;
    fracture_high_interface(n) = ~isempty(af) && hasHigh;
    fracture_internal_partition(n) = numel(af) >= 2 && numel(af) == numel(adj);
    fracture_on_exterior(n) = ~isempty(af) && numel(adj) == 1;
    main_ellipse_low_interface(n) = any(ismember(adj, mainIDs)) && hasLow;
    main_ellipse_high_interface(n) = any(ismember(adj, mainIDs)) && hasHigh;
    compensation_low_interface(n) = any(ismember(adj, compIDs)) && hasLow;
    compensation_high_interface(n) = any(ismember(adj, compIDs)) && hasHigh;
    % Only positive-length shared boundaries connect two fracture domains;
    % touching at a point is not a finite-width conducting connection.
    if L > 1e-10 && numel(af) >= 2
        [~, indices] = ismember(af, fracIDs);
        adjacency(indices, indices) = true;
    end
    if ~ismember(numel(adj), [1 2]) || L <= 1e-10
        unexpected(end+1) = bnd; %#ok<AGROW>
    end
end
adjacency(1:size(adjacency, 1)+1:end) = false;
G = graph(adjacency);
if isempty(fracIDs), bins = []; componentCount = 0;
else, bins = conncomp(G); componentCount = max(bins); end
components = cell(1, componentCount);
for n = 1:componentCount, components{n} = fracIDs(bins == n); end
switch char(p.intervention)
    case {'full', 'none'}, expected = 1;
    case 'cut', expected = 3;
    case 'side', expected = 2;
    otherwise, expected = NaN;
end
m = struct;
m.boundary_count = nb;
m.fracture_domain_ids = fracIDs;
m.fracture_component_count = componentCount;
m.fracture_components_domain_ids = components;
m.fracture_domain_component_labels = bins;
m.expected_component_count_for_design = expected;
m.component_count_matches_design = componentCount == expected;
m.connectivity_definition = 'Fracture domains connected by a shared boundary of length >1e-10 mm; point contact excluded.';
m.fracture_low_interface_length_mm = sum(length_mm(fracture_low_interface));
m.fracture_high_interface_length_mm = sum(length_mm(fracture_high_interface));
m.main_ellipse_low_interface_length_mm = sum(length_mm(main_ellipse_low_interface));
m.main_ellipse_high_interface_length_mm = sum(length_mm(main_ellipse_high_interface));
m.compensation_low_interface_length_mm = sum(length_mm(compensation_low_interface));
m.compensation_high_interface_length_mm = sum(length_mm(compensation_high_interface));
m.fracture_internal_partition_length_mm = sum(length_mm(fracture_internal_partition));
m.fracture_exterior_boundary_length_mm = sum(length_mm(fracture_on_exterior));
m.fracture_exterior_boundary_ids = boundaryIDs(fracture_on_exterior)';
m.unexpected_or_negligible_boundary_ids = unexpected;
m.interface_length_definition = 'Each geometric boundary counted once; same-material partitions excluded from material-interface totals.';
m.interface_length_measurement = 'mphint2 on the saved finite-element solution/mesh, requested unit mm; these are not CAD measure lengths.';
m.cad_areas = cad_area_metrics(model, p, sel);
T = table(boundaryIDs, adjacent_domain_ids, adjacent_materials, adjacent_domain_count, length_mm, ...
    fracture_low_interface, fracture_high_interface, fracture_internal_partition, fracture_on_exterior, ...
    main_ellipse_low_interface, main_ellipse_high_interface, compensation_low_interface, compensation_high_interface, ...
    'VariableNames', {'boundary_id', 'adjacent_domain_ids', 'adjacent_materials', 'adjacent_domain_count', 'length_mm', ...
    'fracture_low_interface', 'fracture_high_interface', 'fracture_internal_partition', 'fracture_on_exterior', ...
    'main_ellipse_low_interface', 'main_ellipse_high_interface', 'compensation_low_interface', 'compensation_high_interface'});
end

function out = cad_area_metrics(model, p, sel)
% Every mphmeasure call explicitly disables geometry building. Unsupported
% two-dimensional measures or unverified units remain unavailable, not data.
out = struct('status', 'unavailable', ...
    'method', 'mphmeasure on finalized geometry: entity=domain, build=off, usefinal=on', ...
    'distinction', 'Approximate rendering-mesh geometry measures, separate from finite-element mphint2 areas; not higher-accuracy CAD ground truth.', ...
    'errors', {{}});
try
    out.geometry_length_unit = char(model.component('comp1').geom('geom1').lengthUnit);
catch ME
    out.errors{end+1} = ['Geometry length unit: ' ME.message];
end
try
    allIDs = row_ids(sel.all_domains);
    rawTotal = raw_cad_measure(model, allIDs);
    out.raw_total_domain_measure = rawTotal;
    expectedTotal = p.L_mm^2;
    % Validate against the known rectangular domain before assigning mm^2.
    % This also detects a 3D-style volume/perimeter return on a 2D model.
    if abs(rawTotal-expectedTotal) <= expectedTotal*1e-8
        factor = 1;
        out.raw_measure_unit_verified_as = 'mm^2';
    elseif abs(rawTotal-expectedTotal*1e-6) <= expectedTotal*1e-14
        factor = 1e6;
        out.raw_measure_unit_verified_as = 'm^2';
    else
        error('enrich_mechanism:CADUnitUnknown', ...
            '2D total-domain measure %.16g is neither the known %.16g mm^2 nor its SI equivalent.', rawTotal, expectedTotal);
    end
    out.raw_to_mm2_multiplier = factor;
    out.unit_validation = 'Verified units against the known 120 mm by 120 mm rectangle; this does not validate curved-boundary area accuracy. No geometry rebuild.';
    out.area_total_mm2 = rawTotal*factor;
    fracIDs = row_ids(sel.fracture_domains);
    compIDs = []; removedIDs = [];
    if isfield(sel, 'compensation_domains'), compIDs = row_ids(sel.compensation_domains); end
    if isfield(sel, 'removed_domains'), removedIDs = row_ids(sel.removed_domains); end
    out.area_low_mm2 = factor*raw_cad_measure(model, row_ids(sel.low_domains));
    out.area_high_mm2 = factor*raw_cad_measure(model, row_ids(sel.high_domains));
    out.area_fracture_mm2 = factor*raw_cad_measure(model, fracIDs);
    out.area_remaining_main_ellipse_mm2 = factor*raw_cad_measure(model, setdiff(fracIDs, compIDs));
    out.removed_area_mm2 = factor*raw_cad_measure(model, removedIDs);
    out.compensation_area_mm2 = factor*raw_cad_measure(model, compIDs);
    originalIDs = double(sel.all_domains(:)');
    backgroundLayers = double(sel.domain_background_layer(:)');
    [present, locations] = ismember(fracIDs, originalIDs);
    assert(all(present) && numel(originalIDs) == numel(backgroundLayers), 'Invalid background-layer mapping for CAD areas.');
    replacedLayers = backgroundLayers(locations);
    out.replaced_low_mm2 = factor*raw_cad_measure(model, fracIDs(mod(replacedLayers, 2) == 1));
    out.replaced_high_mm2 = factor*raw_cad_measure(model, fracIDs(mod(replacedLayers, 2) == 0));
    out.phase_area_closure_error_mm2 = out.area_low_mm2+out.area_high_mm2+out.area_fracture_mm2-out.area_total_mm2;
    out.compensation_minus_removed_mm2 = out.compensation_area_mm2-out.removed_area_mm2;
    out.analytic_reference_ellipse_area_mm2 = pi*p.a_mm*p.b_mm;
    out.reassembled_ellipse_minus_analytic_mm2 = out.area_remaining_main_ellipse_mm2+out.removed_area_mm2-out.analytic_reference_ellipse_area_mm2;
    if ~isempty(compIDs)
        out.analytic_compensation_disk_area_mm2 = pi*p.comp_radius_mm^2;
        out.compensation_minus_analytic_disk_mm2 = out.compensation_area_mm2-out.analytic_compensation_disk_area_mm2;
    end
    out.measurement_status = 'completed';
    out.status = 'approximate_not_used_for_matching';
catch ME
    out.errors{end+1} = ME.message;
    if isfield(out, 'area_total_mm2'), out.status = 'partial'; end
end
end

function out = annotate_geometry_measurement_quality(out, fem)
% Keep every originally measured number. The geometry measurement API uses
% a rendering mesh and must not be promoted to a CAD area reference.
out.quality_status = 'approximate_not_used_for_matching';
out.eligible_for_phase_matching = false;
out.preserved_for = 'Low-precision geometry-measurement diagnostic only.';
out.distinction = 'Approximate rendering-mesh geometry measures, separate from finite-element mphint2 areas; not higher-accuracy CAD ground truth.';
out.documentation_source = 'D:/Program Files/COMSOL/COMSOL63/Multiphysics/doc/help/wtpwebapps/ROOT/doc/com.comsol.help.comsol/comsol_api_geom.48.046.html';
out.documentation_section = 'Measurement Methods';
out.documentation_quote_excerpt = 'The volumes, areas, and lengths are approximate because they are based on approximating the geometry with a triangular or polygonal mesh';
out.documentation_precision_summary = 'The documentation identifies the approximation as a rendering mesh with rendering detail set to fine.';
out.accuracy_control_audit = 'No measurement-accuracy/tolerance argument is documented for mphmeasure or GeomMeasureBase/GeomMeasureFinal in the installed COMSOL 6.3 API. No repair tolerance, rendering setting, or geometry rebuild was changed.';
out.mphmeasure_reference = 'D:/Program Files/COMSOL/COMSOL63/Multiphysics/doc/help/wtpwebapps/ROOT/doc/com.comsol.help.llmatlab/llmatlab_ug_ref.9.43.html';
out.geometry_measure_api_reference = 'D:/Program Files/COMSOL/COMSOL63/Multiphysics/doc/help/wtpwebapps/ROOT/doc/com.comsol.help.comsol/api/com/comsol/model/GeomMeasureBase.html';
out.matching_evidence_policy = 'Use separately labeled finite-element integrals with analytic geometry checks and mesh convergence. Passing total area closure establishes neither curved-area accuracy nor exact phase matching.';
if ismember(out.status, {'success', 'approximate_not_used_for_matching'})
    out.measurement_status = 'completed';
    out.status = 'approximate_not_used_for_matching';
end
if isfield(out, 'analytic_reference_ellipse_area_mm2') && isfield(out, 'reassembled_ellipse_minus_analytic_mm2')
    reference = out.analytic_reference_ellipse_area_mm2;
    out.ellipse_analytic_check_absolute_error_mm2 = abs(out.reassembled_ellipse_minus_analytic_mm2);
    out.ellipse_analytic_check_relative_error = out.ellipse_analytic_check_absolute_error_mm2/reference;
    out.ellipse_analytic_check_percent_error = 100*out.ellipse_analytic_check_relative_error;
    femReassembled = fem.area_fracture_mm2;
    if isfield(fem, 'compensation_area_mm2'), femReassembled = femReassembled-fem.compensation_area_mm2; end
    if isfield(fem, 'removed_area_mm2'), femReassembled = femReassembled+fem.removed_area_mm2; end
    out.fem_reassembled_ellipse_area_mm2 = femReassembled;
    out.fem_ellipse_analytic_check_absolute_error_mm2 = abs(femReassembled-reference);
    out.fem_ellipse_analytic_check_relative_error = out.fem_ellipse_analytic_check_absolute_error_mm2/reference;
    out.ellipse_measurement_error_exceeds_fem = out.ellipse_analytic_check_absolute_error_mm2 > out.fem_ellipse_analytic_check_absolute_error_mm2;
    out.accuracy_rejection_reason = sprintf(['Rendering-mesh ellipse area check error is %.9g percent; ' ...
        'the separately stored FEM integral check error is %.9g percent. Geometry-measure values are diagnostic only and excluded from phase-matching evidence.'], ...
        out.ellipse_analytic_check_percent_error, 100*out.fem_ellipse_analytic_check_relative_error);
else
    out.accuracy_rejection_reason = 'Approximate rendering-mesh measurement; analytic validation unavailable or incomplete. Excluded from phase matching.';
end
end

function value = raw_cad_measure(model, ids)
if isempty(ids), value = 0; return; end
value = mphmeasure(model, 'geom1', 'domain', 'selection', ids, 'build', 'off', 'usefinal', 'on');
assert(isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value) && value >= 0, ...
    'mphmeasure did not return a nonnegative scalar domain measure.');
value = double(value);
end

function out = partitioned_section_current(model, p, sel)
s = 6;
assert(abs(s) < p.a_mm, 'Section does not intersect the reference ellipse.');
halfWidth = p.b_mm * sqrt(1 - s^2/p.a_mm^2);
breaks = [-6 -3 -halfWidth halfWidth 3 6];
if strcmp(p.intervention, 'side') && p.side_width_mm > -halfWidth && p.side_width_mm < halfWidth
    breaks(end+1) = p.side_width_mm;
end
breaks = unique(sort(breaks));
% Verify that the coordinate convention actually used by interpolation is mm.
xy = local_coordinates(s, 0, p);
[qx, qz] = mphinterp(model, {'x', 'y'}, 'coord', xy, 'edim', 2, 'ext', 0, ...
    'solnum', 'end', 'unit', {'mm', 'mm'});
coordinateError = max(abs([qx(:); qz(:)] - xy(:)));
assert(isfinite(coordinateError) && coordinateError < 1e-6, 'mphinterp coordinate-unit validation failed.');
steps = [0.01 0.005];
passes = cell(1, numel(steps));
for pass = 1:numel(steps)
    segments = repmat(struct('t_start_mm', 0, 't_end_mm', 0, 'category', '', ...
        'material_selection', '', 'sample_count', 0, 'sample_step_mm', 0, ...
        'current_signed_A_per_m', 0, 'mean_Js_A_m2', 0, 'min_Js_A_m2', 0, 'max_Js_A_m2', 0), numel(breaks)-1, 1);
    Ifrac = 0; Irestored = 0; Ioutside = 0;
    for j = 1:numel(breaks)-1
        lo = breaks(j); hi = breaks(j+1); mid = (lo+hi)/2;
        inEllipse = abs(mid) < halfWidth;
        removed = false;
        if strcmp(p.intervention, 'cut')
            removed = inEllipse && abs(s-p.cut_u_mm) < p.cut_width_mm/2;
        elseif strcmp(p.intervention, 'side')
            removed = inEllipse && s > 4 && s < 8 && mid > p.side_width_mm;
        end
        if inEllipse && ~removed
            category = 'fracture'; selection = row_ids(sel.fracture_domains); material = 'fracture';
        elseif removed
            category = 'restored_low'; selection = row_ids(sel.low_domains); material = 'low';
        else
            category = 'surrounding_rock_outside_original_ellipse';
            selection = unique([row_ids(sel.low_domains) row_ids(sel.high_domains)]); material = 'low_or_high';
        end
        n = max(16, ceil((hi-lo)/steps(pass)));
        dt = (hi-lo)/n;
        % Midpoints partition the entire segment without sampling either
        % material boundary. Do not trapezoid across a conductivity jump.
        t = lo + ((1:n)-0.5)*dt;
        xy = local_coordinates(s*ones(size(t)), t, p);
        [jx, jz] = mphinterp(model, {'ec.Jx', 'ec.Jy'}, 'coord', xy, ...
            'edim', 2, 'selection', selection, 'ext', 0, 'recover', 'off', 'solnum', 'end', 'unit', {'A/m^2', 'A/m^2'});
        assert(numel(jx) == n && numel(jz) == n, 'Unexpected solution count in interpolation.');
        Js = jx(:)*cosd(p.theta_deg) + jz(:)*sind(p.theta_deg);
        assert(all(isfinite(Js)) && isreal(Js), 'Interpolation contains nonfinite/complex values; material classification may be wrong.');
        current = sum(Js)*dt*1e-3;
        switch category
            case 'fracture', Ifrac = Ifrac+current;
            case 'restored_low', Irestored = Irestored+current;
            otherwise, Ioutside = Ioutside+current;
        end
        segments(j) = struct('t_start_mm', lo, 't_end_mm', hi, 'category', category, ...
            'material_selection', material, 'sample_count', n, 'sample_step_mm', dt, ...
            'current_signed_A_per_m', current, 'mean_Js_A_m2', mean(Js), ...
            'min_Js_A_m2', min(Js), 'max_Js_A_m2', max(Js));
    end
    passes{pass} = struct('maximum_sample_step_mm', steps(pass), ...
        'fracture_current_signed_A_per_m', Ifrac, 'restored_low_current_signed_A_per_m', Irestored, ...
        'surrounding_rock_outside_original_ellipse_current_signed_A_per_m', Ioutside, ...
        'all_rock_current_signed_A_per_m', Irestored+Ioutside, ...
        'total_current_signed_A_per_m', Ifrac+Irestored+Ioutside, 'segments', segments);
end
fine = passes{end}; coarse = passes{1};
out = fine;
out.status = 'success';
out.s_mm = s;
out.t_limits_mm = [-6 6];
out.reference_ellipse_halfwidth_mm = halfWidth;
out.breakpoints_mm = breaks;
out.coordinate_unit_check_max_error_mm = coordinateError;
out.direction = 'Positive along e_s=(cos(theta),sin(theta)); Js=ec.Jx*cos(theta)+ec.Jy*sin(theta), where COMSOL y is physical z.';
out.unit_conversion = 'Integral of Js [A/m^2] times dt [mm]*1e-3 -> A/m; multiply by recorded out-of-plane thickness for A.';
out.out_of_plane_thickness_m = p.out_of_plane_thickness_m;
out.total_current_signed_A = fine.total_current_signed_A_per_m*p.out_of_plane_thickness_m;
out.integration_method = 'Composite midpoint, separate material segments, domain selection restricted, ext=0, recover=off; no boundary endpoint samples.';
out.coarse_integration = coarse;
out.quadrature_total_change_A_per_m = fine.total_current_signed_A_per_m-coarse.total_current_signed_A_per_m;
out.quadrature_component_absolute_change_sum_A_per_m = ...
    abs(fine.fracture_current_signed_A_per_m-coarse.fracture_current_signed_A_per_m) + ...
    abs(fine.restored_low_current_signed_A_per_m-coarse.restored_low_current_signed_A_per_m) + ...
    abs(fine.surrounding_rock_outside_original_ellipse_current_signed_A_per_m-coarse.surrounding_rock_outside_original_ellipse_current_signed_A_per_m);
out.quadrature_change_relative_to_total = out.quadrature_component_absolute_change_sum_A_per_m / max(abs(fine.total_current_signed_A_per_m), eps);
end

function meta = solver_metadata(model, id, modelFile)
meta = struct('case_id', id, 'source_mph', modelFile, 'postprocessing_only', true, ...
    'solver_sequences', {{}}, 'errors', {{}});
try
    meta.comsol_version = char(com.comsol.model.util.ModelUtil.getComsolVersion);
catch ME
    meta.errors{end+1} = ['COMSOL version: ' ME.message];
end
try
    solTags = text_list(model.sol.tags);
    sequences = cell(1, numel(solTags));
    for k = 1:numel(solTags)
        seq = model.sol(solTags{k});
        info = struct('tag', solTags{k}, 'nodes', {{}}, 'errors', {{}});
        try, info.label = char(seq.label); catch ME, info.errors{end+1} = ME.message; end
        try, info.type = char(seq.getType); catch ME, info.errors{end+1} = ['Sequence type: ' ME.message]; end
        try, info.solution_size = double(seq.getSize(1)); catch ME, info.errors{end+1} = ['Solution size: ' ME.message]; end
        info.nodes = solver_nodes(seq, solTags{k}, 0);
        sequences{k} = info;
    end
    meta.solver_sequences = sequences;
catch ME
    meta.errors{end+1} = ['Solver tree: ' ME.message];
end
meta.physics = struct('component_tag', 'comp1', 'physics_tag', 'ec', ...
    'coordinate_mapping', 'COMSOL x,y are physical x,z');
try
    ec = model.component('comp1').physics('ec');
    try, meta.physics.type = char(ec.getType); catch ME, meta.physics.type_error = ME.message; end
    try
        shape = ec.prop('ShapeProperty');
        meta.physics.ShapeProperty = property_snapshot(shape, {'order_electricpotential', 'order'}, true);
    catch ME
        meta.physics.ShapeProperty = struct('error', ME.message);
    end
    try
        meta.physics.out_of_plane_thickness_expression = char(ec.prop('d').getString('d'));
    catch ME
        meta.physics.out_of_plane_thickness_error = ME.message;
    end
catch ME
    meta.errors{end+1} = ['Physics properties: ' ME.message];
end
try, meta.physics.d_out_parameter_expression = char(model.param.get('d_out'));
catch ME, meta.physics.d_out_parameter_expression_error = ME.message; end
try, meta.physics.d_out_parameter_evaluated_SI_m = double(model.param.evaluate('d_out'));
catch ME, meta.physics.d_out_parameter_evaluation_error = ME.message; end
meta.interpretation = 'Stored node properties, including linsolver and active state where supported; unavailable API values remain explicit errors. No solver settings are inferred from defaults.';
end

function nodes = solver_nodes(parent, prefix, depth)
nodes = {};
if depth > 12
    nodes = {struct('path', prefix, 'error', 'Solver-tree depth limit reached.')}; return;
end
try
    tags = text_list(parent.feature.tags);
catch ME
    nodes = {struct('path', prefix, 'error', ['Feature enumeration: ' ME.message])}; return;
end
for k = 1:numel(tags)
    path = [prefix '/' tags{k}];
    info = struct('tag', tags{k}, 'path', path, 'errors', {{}});
    try
        node = parent.feature(tags{k});
        try, info.type = char(node.getType); catch ME, info.errors{end+1} = ['Type: ' ME.message]; end
        try, info.label = char(node.label); catch ME, info.errors{end+1} = ['Label: ' ME.message]; end
        try, info.active = logical(node.isActive); catch ME, info.errors{end+1} = ['Active: ' ME.message]; end
        keys = {'linsolver', 'precond', 'rtol', 'atol', 'maxiter', 'iter', 'method', 'relin', 'dtech', 'btech', 'rowperm', 'colperm'};
        info.properties = property_snapshot(node, keys, false);
        info.children = solver_nodes(node, path, depth+1);
    catch ME
        info.errors{end+1} = ME.message;
    end
    nodes{end+1} = info; %#ok<AGROW>
end
end

function snap = property_snapshot(node, requested, readAll)
snap = struct('property_names', {{}}, 'values', struct, 'unavailable', struct, 'errors', {{}});
namesKnown = false;
try
    names = text_list(node.properties);
    snap.property_names = names;
    namesKnown = true;
catch ME
    names = {};
    snap.errors{end+1} = ['Property enumeration: ' ME.message];
end
if readAll && namesKnown
    keys = unique([names requested], 'stable');
elseif namesKnown
    keys = intersect(requested, names, 'stable');
else
    keys = requested;
end
for k = 1:numel(keys)
    key = keys{k}; field = matlab.lang.makeValidName(key);
    if namesKnown && ~ismember(key, names)
        snap.unavailable.(field) = 'Not listed by this node.'; continue;
    end
    try
        snap.values.(field) = char(node.getString(key));
        continue;
    catch
    end
    try
        snap.values.(field) = text_list(node.getStringArray(key));
        continue;
    catch
    end
    try
        snap.values.(field) = double(node.getDouble(key));
        continue;
    catch
    end
    try
        snap.values.(field) = double(node.getInt(key));
        continue;
    catch
    end
    try
        snap.values.(field) = logical(node.getBoolean(key));
    catch ME
        snap.unavailable.(field) = ME.message;
    end
end
end

function xy = local_coordinates(s, t, p)
xy = [p.cx_mm+s*cosd(p.theta_deg)-t*sind(p.theta_deg); ...
      p.cz_mm+s*sind(p.theta_deg)+t*cosd(p.theta_deg)];
end

function ids = row_ids(values)
ids = unique(double(values(:)'));
ids = ids(ids > 0);
end

function result = text_list(values)
if isempty(values), result = {}; return; end
if ischar(values), result = cellstr(values)'; return; end
if isstring(values), result = cellstr(values(:))'; return; end
result = cell(1, numel(values));
if iscell(values)
    for k = 1:numel(values), result{k} = char(values{k}); end
else
    for k = 1:numel(values), result{k} = char(values(k)); end
end
end

function write_json(path, value)
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid >= 0, 'Cannot open derived JSON output: %s', path);
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', jsonencode(value, 'PrettyPrint', true));
end

function remove_owned_model(tag)
try, com.comsol.model.util.ModelUtil.remove(tag); catch, end
end

function disconnect_if_owned(owned)
if owned
    try, com.comsol.model.util.ModelUtil.disconnect; catch, end
end
end
