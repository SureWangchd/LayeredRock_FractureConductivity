function report=enrich_round2_exchange(caseIDs)
% Read-only postprocessing of saved Round 2 solutions. No solve or MPH write.
% Usage: enrich_round2_exchange({'case1','case2'}), or no args for all successes.
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'runs','2026-10-21_round2');
if nargin<1 || isempty(caseIDs)
 files=dir(fullfile(out,'*','SUCCESS.json'));caseIDs=cell(1,numel(files));
 for k=1:numel(files),[~,caseIDs{k}]=fileparts(files(k).folder);end
end
if ischar(caseIDs)||isstring(caseIDs),caseIDs=cellstr(caseIDs);end
base='D:/Program Files/COMSOL/COMSOL63/Multiphysics';addpath(fullfile(base,'mli'));
import com.comsol.model.util.*
owned=false;
try,mphstart('127.0.0.1',2037,base);owned=true;
catch err,try,ModelUtil.tags;catch,rethrow(err);end,end
connectionCleanup=onCleanup(@()disconnect_if_owned(owned)); %#ok<NASGU>
report=repmat(struct('case_id','','status','','error',''),numel(caseIDs),1);
for k=1:numel(caseIDs)
 id=char(caseIDs{k});report(k).case_id=id;folder=fullfile(out,id);
 mphFile=fullfile(folder,[id '.mph']);success=fullfile(folder,'SUCCESS.json');
 stamp=fullfile(folder,'enrichment_complete.json');
 if ~isfile(success)||~isfile(mphFile),report(k).status='missing_success';continue;end
 if isfile(stamp)
  previous=jsondecode(fileread(stamp));mi=dir(mphFile);
  if strcmp(previous.status,'success')&&previous.source_mph_datenum==mi.datenum
   report(k).status='already_complete';fprintf('ENRICH_ROUND2_SKIP %s\n',id);continue;
  end
 end
 tag=['Round2Extra_' strrep(char(java.util.UUID.randomUUID),'-','')];
 modelCleanup=onCleanup(@()remove_owned_model(tag));
 fprintf('ENRICH_ROUND2_START %s\n',id);started=tic;
 try
  p=jsondecode(fileread(fullfile(folder,'parameters.json')));
  sel=jsondecode(fileread(fullfile(folder,'selections.json')));
  model=mphload(mphFile,tag);
  endpoint=finite_cap_potentials(model,sel);
  endpoint.case_id=id;endpoint.passive_ellipse=p.passive_ellipse;
  endpoint.source_mph=mphFile;endpoint.postprocessing_only=true;
  write_json(fullfile(folder,'endpoint_potential_averages.json'),endpoint);
  metadata=solver_metadata(model,id,mphFile);
  metadata=complete_metadata(metadata,p,sel);
  write_json(fullfile(folder,'metadata.json'),metadata);
  exchangeStatus='success';
  try
   [exchange,T]=ellipse_arc_exchange(model,p,sel);
   exchange.case_id=id;exchange.source_mph=mphFile;exchange.postprocessing_only=true;
   write_json(fullfile(folder,'ellipse_boundary_exchange.json'),exchange);
   writetable(T,fullfile(folder,'ellipse_boundary_exchange_segments.csv'));
  catch err
   exchangeStatus='failed';
   write_json(fullfile(folder,'ellipse_boundary_exchange.json'),struct('case_id',id,'status','failed','error',getReport(err,'extended','hyperlinks','off')));
   warning('round2:ExchangeFailed','%s: %s',id,err.message);
  end
  mi=dir(mphFile);completion=struct('case_id',id,'status','success','source_mph_datenum',mi.datenum,...
    'source_mph_bytes',mi.bytes,'postprocessing_seconds',toc(started),'exchange_status',exchangeStatus);
  if strcmp(exchangeStatus,'failed'),completion.status='partial';end
  write_json(stamp,completion);report(k).status=completion.status;
  fprintf('ENRICH_ROUND2_%s %s %.1fs\n',upper(completion.status),id,toc(started));
 catch err
  report(k).status='failed';report(k).error=getReport(err,'extended','hyperlinks','off');
  warning('round2:EnrichmentFailed','%s: %s',id,report(k).error);
 end
 clear modelCleanup model;
end
end

function e=finite_cap_potentials(model,sel)
minus=row_ids(sel.capminus_domains);plus=row_ids(sel.capplus_domains);
assert(~isempty(minus)&&~isempty(plus),'Finite cap selections are empty');
am=mphint2(model,'1',2,'selection',minus,'intorder',8,'unit','m^2');
ap=mphint2(model,'1',2,'selection',plus,'intorder',8,'unit','m^2');
vm=mphint2(model,'V',2,'selection',minus,'intorder',8,'unit','V*m^2')/am;
vp=mphint2(model,'V',2,'selection',plus,'intorder',8,'unit','V*m^2')/ap;
assert(all(isfinite([am ap vm vp]))&&am>0&&ap>0,'Invalid cap integrals');
e=struct('status','success','capminus_V',vm,'capplus_V',vp,'dV_V',vp-vm,'dV2_V2',(vp-vm)^2,...
 'capminus_area_mm2',am*1e6,'capplus_area_mm2',ap*1e6,...
 'definition','Area average of V on actual ellipse caps of axial length 3 mm. dV=plus-minus, not a mathematical point-tip voltage.',...
 'integration_order',8);
end

function meta=complete_metadata(meta,p,sel)
interfaces=cell(1,numel(sel.z_interfaces_mm));
for j=1:numel(interfaces)
 interfaces{j}=struct('index',j,'expression_geometry_mm',sprintf('%.17g + %.17g*sin(2*pi*x/120)',...
  sel.z_interfaces_mm(j),sel.interface_amplitudes_mm(j)),...
  'phase_below',sel.layer_materials(j),'phase_above',sel.layer_materials(j+1));
end
meta.interfaces=interfaces;meta.parameters_recorded=p;
meta.direct_solver_nodes={};
for k=1:numel(meta.solver_sequences)
 sequence=meta.solver_sequences{k};
 meta.direct_solver_nodes=[meta.direct_solver_nodes collect_direct_nodes(sequence.nodes)]; %#ok<AGROW>
end
meta.units=struct('geometry_input','mm','COMSOL_internal_length','m','fields_csv_coordinates','m',...
 'voltage','V','conductivity','S/m','electric_field','V/m','current_density','A/m^2',...
 'boundary_current_2D','A/m of out-of-plane thickness','reaction_current','A, including d_out',...
 'area','mm^2 unless explicitly suffixed m2','volume_dissipation','W/m of out-of-plane thickness',...
 'local_dissipation_density','W/m^3','out_of_plane_thickness','m','theta','degrees counterclockwise from +x toward +z');
meta.electrodes=struct('load',p.load,'positive_V',1,'ground_V',0,'positive_boundaries',sel.positive_electrode,...
 'ground_boundaries',sel.ground_electrode,'other_exterior_boundaries','Electric insulation');
meta.linear_solver_evidence='Read solver_sequences[*].nodes recursively for actual Direct node properties.linsolver; do not infer from defaults.';
meta.geometry_coordinate_convention='COMSOL x,y correspond to physical x,z. Rotation has s=(x-cx)cos(theta)+(z-cz)sin(theta).';
end

function found=collect_direct_nodes(nodes)
found={};
for k=1:numel(nodes)
 n=nodes{k};
 if isfield(n,'type')&&strcmp(n.type,'Direct')
  f=struct('path',n.path,'active',false,'linsolver','unavailable');
  if isfield(n,'active'),f.active=n.active;end
  if isfield(n,'properties')&&isfield(n.properties,'values')&&isfield(n.properties.values,'linsolver')
   f.linsolver=n.properties.values.linsolver;
  end
  found{end+1}=f; %#ok<AGROW>
 end
 if isfield(n,'children'),found=[found collect_direct_nodes(n.children)];end %#ok<AGROW>
end
end

function [e,T]=ellipse_arc_exchange(model,p,sel)
names={'boundary_id','matrix_phase','length_mm','ellipse_residual_mean','net_out_A_per_m','absolute_exchange_A_per_m','trace_net_difference_A_per_m','trace_absolute_difference_A_per_m'};
rows=zeros(0,numel(names));
e=struct('status','success','normal_definition','Analytic gradient of (s/a)^2+(t/b)^2, pointing from ellipse into matrix.',...
 'scope','Only surviving original-ellipse arcs adjoining matrix. Newly created cut/side faces and compensation-circle boundaries are excluded.',...
 'units','A/m of out-of-plane thickness','trace_method','Mean of up and down current traces; their difference is exported as an interface discretization diagnostic.',...
 'integration_order',8,'low_net_out_A_per_m',0,'high_net_out_A_per_m',0,...
 'low_absolute_exchange_A_per_m',0,'high_absolute_exchange_A_per_m',0,...
 'excluded_new_face_count',0,'selected_arc_count',0,'max_mean_ellipse_residual',0);
main=setdiff(row_ids(sel.fracture_domains),row_ids(sel.compensation_domains));
if isempty(main)
 e.status='not_applicable_passive_ellipse';e.scope='No physical fracture/matrix conductivity interface in this passive background model.';
 T=array2table(rows,'VariableNames',names);return;
end
low=row_ids(sel.low_domains);high=row_ids(sel.high_domains);
bnds=[];for d=main,bnds=union(bnds,row_ids(mphgetadj(model,'geom1','boundary','domain',d)));end
% union([],row) may return a column vector; MATLAB for iterates columns.
% Normalize before scalar iteration, otherwise an entire boundary set is
% passed to mphgetadj and legitimate interfaces appear to have many domains.
bnds=row_ids(bnds);
s='((x-cx)*cos(theta)+(y-cz)*sin(theta))';
t='(-(x-cx)*sin(theta)+(y-cz)*cos(theta))';
f=['((' s '/a)^2+(' t '/b)^2-1)'];
gx=['(' s '*cos(theta)/a^2-' t '*sin(theta)/b^2)'];
gz=['(' s '*sin(theta)/a^2+' t '*cos(theta)/b^2)'];
den=['sqrt((' s '/a^2)^2+(' t '/b^2)^2)'];
n1=['(' gx '/' den ')'];n2=['(' gz '/' den ')'];
ju=['(up(ec.Jx)*' n1 '+up(ec.Jy)*' n2 ')'];
jd=['(down(ec.Jx)*' n1 '+down(ec.Jy)*' n2 ')'];
jmean=['0.5*(' ju '+' jd ')'];
for boundaryIndex=1:numel(bnds)
 b=bnds(boundaryIndex);
 adj=row_ids(mphgetadj(model,'geom1','domain','boundary',b));
 if ~any(ismember(adj,main)),continue;end
 matrix=intersect(adj,[low high]);if isempty(matrix),continue;end
 assert(numel(adj)==2&&numel(matrix)==1,'Unexpected matrix/fracture boundary adjacency: boundary=%d adjacent=[%s] matrix=[%s]',b,num2str(adj),num2str(matrix));
 length=mphint2(model,'1',1,'selection',b,'intorder',8,'unit','m');
 residual=mphint2(model,['abs(' f ')'],1,'selection',b,'intorder',8,'unit','m')/length;
 if residual>1e-6,e.excluded_new_face_count=e.excluded_new_face_count+1;continue;end
 phase=double(ismember(matrix,high));
 net=mphint2(model,jmean,1,'selection',b,'intorder',8,'unit','A/m');
 absolute=mphint2(model,['abs(' jmean ')'],1,'selection',b,'intorder',8,'unit','A/m');
 traceNet=mphint2(model,['(' ju '-' jd ')'],1,'selection',b,'intorder',8,'unit','A/m');
 traceAbs=mphint2(model,['abs(' ju '-' jd ')'],1,'selection',b,'intorder',8,'unit','A/m');
 assert(all(isfinite([length residual net absolute traceNet traceAbs])),'Nonfinite interface integral');
 rows(end+1,:)=[b phase length*1000 residual net absolute traceNet traceAbs]; %#ok<AGROW>
 if phase==0,e.low_net_out_A_per_m=e.low_net_out_A_per_m+net;e.low_absolute_exchange_A_per_m=e.low_absolute_exchange_A_per_m+absolute;
 else,e.high_net_out_A_per_m=e.high_net_out_A_per_m+net;e.high_absolute_exchange_A_per_m=e.high_absolute_exchange_A_per_m+absolute;end
end
T=array2table(rows,'VariableNames',names);e.selected_arc_count=height(T);
assert(e.selected_arc_count>0,'No original ellipse boundary arcs were identified');
e.max_mean_ellipse_residual=max(T.ellipse_residual_mean);
e.total_net_out_A_per_m=e.low_net_out_A_per_m+e.high_net_out_A_per_m;
e.total_absolute_exchange_A_per_m=e.low_absolute_exchange_A_per_m+e.high_absolute_exchange_A_per_m;
e.total_trace_absolute_difference_A_per_m=sum(T.trace_absolute_difference_A_per_m);
e.trace_difference_relative_to_absolute_exchange=e.total_trace_absolute_difference_A_per_m/max(e.total_absolute_exchange_A_per_m,eps);
e.net_to_absolute_exchange=abs(e.total_net_out_A_per_m)/max(e.total_absolute_exchange_A_per_m,eps);
e.has_new_cut_faces=ismember(p.intervention,{'cut','side'});
e.net_interpretation='For full ellipses the net should tend to zero with refinement. For cut/side cases these arcs are an open boundary subset; net need not vanish.';
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
