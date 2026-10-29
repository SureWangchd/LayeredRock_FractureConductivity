function r6_truth_run(root,action)
% Round 6 full-target validation only. Never constructs prediction inputs.
% Reuses the existing geometry/solver and R3 truth/refinement routines.
% Caller owns the LiveLink connection. This entry never starts a server.
root=char(root);action=char(action);out=fullfile(root,'runs','round6_applicability');
design=jsondecode(fileread(fullfile(out,'design.json')));
cases=as_cells(design.cases);
assert(isfile(fullfile(out,'method_frozen.json')),'Freeze Round 6 settings before target validation');
assert(isfile(fullfile(out,'criteria_ready.json')),'All 32 prospective applicability decisions and intervals must precede target validation');
jsondecode(fileread(fullfile(out,'method_frozen.json')));
jsondecode(fileread(fullfile(out,'criteria_ready.json')));
assert(numel(cases)==16,'Round 6 contains exactly 16 preregistered geometries');
if strcmp(action,'initial')
 readyFile=fullfile(out,'predictions_ready.json');
 assert(isfile(readyFile),'All 32 predictions must be saved before any new full target');
 ready=jsondecode(fileread(readyFile));states=as_cells(design.new_states);
 assert(numel(states)==32,'The preregistered initial batch must contain exactly 32 states');
 preflight_predictions(out,states,ready);
 run_states(out,cases,states,'initial_32_prediction_batch');
 write_completion(out,states,fullfile(out,'truth_initial_complete.json'));
elseif strcmp(action,'refine')
 assert(isfile(fullfile(out,'truth_initial_complete.json')),'Finish the 32 main states before target refinement');
 request=jsondecode(fileread(fullfile(out,'refinement_targets.json')));
 if isstruct(request)&&isscalar(request)&&isfield(request,'targets'),entries=as_cells(request.targets);else,entries=as_cells(request);end
 assert(numel(entries)<=6,'Round 6 permits at most six target/background refinement pairs');
 planned=as_cells(design.new_states);seen={};
 for j=1:numel(entries)
  e=entries{j};
  assert(isfield(e,'case_id')&&isfield(e,'sigma_f_S_m'),'Each refinement must identify a prescribed state');
  key=[char(e.case_id) '/' sf_token(e.sigma_f_S_m)];
  assert(~ismember(key,seen),'Duplicate refinement state');seen{end+1}=key; %#ok<AGROW>
  assert(any(cellfun(@(s)strcmp(char(s.case_id),char(e.case_id))&&abs(s.sigma_f_S_m-e.sigma_f_S_m)<1e-12,planned)), ...
   'Refinements must belong to the 32 prescribed states');
  if isfield(e,'source'),source=char(e.source);else,source=fullfile(out,'truth',char(e.case_id),sf_token(e.sigma_f_S_m));end
  if isfield(e,'output'),destination=char(e.output);else,destination=fullfile(out,'precision',[char(e.case_id) '_' sf_token(e.sigma_f_S_m)]);end
  source=absolute_path(root,source);destination=absolute_path(root,destination);
  assert(inside_folder(destination,out)&&inside_folder(source,fullfile(out,'truth')), ...
   'Round 6 refinement paths must use only this round''s new truth and output directories');
  r3_refine_case(source,destination);
 end
 r3_json(fullfile(out,'truth_refinement_complete.json'),struct('status','success','count',numel(entries), ...
  'maximum_allowed_count',6,'completed_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'))));
else,error('Unknown Round 6 truth action: %s',action);end
end

function run_states(out,cases,states,phase)
ids=cellfun(@(s)char(s.case_id),states,'UniformOutput',false);
for caseId=reshape(unique(ids,'stable'),1,[])
 id=caseId{1};found=find(cellfun(@(c)strcmp(case_id(c),id),cases));
 assert(numel(found)==1,'Target case ID must identify one preregistered geometry');
 entry=cases{found};if isfield(entry,'c'),c=entry.c;else,c=entry;end
 assert(strcmp(c.id,id),'Case ID mismatch');
 c.passive_ellipse=true;assert(c.fracture&&strcmp(c.load,'z')&&c.L_mm==120,'Expected the 120 mm passive ellipse z background');
 assert(ismember(2*c.a_mm,[40 64])&&c.b_mm==1.5&&c.theta_deg==60&&c.cx_mm==60&&c.cz_mm==60, ...
  'Round 6 source geometry differs from the prescribed lengths/width/angle/center');
 assert(c.out_of_plane_thickness_m==1,'Existing truth energy convention requires t=1 m');
 folder=fullfile(out,'truth',id);if ~isfolder(folder),mkdir(folder);end
 selected=states(strcmp(ids,id));pending=cellfun(@(s)~isfile(fullfile(folder,sf_token(s.sigma_f_S_m),'SUCCESS.json')),selected);
 if ~any(pending),fprintf('R6_TRUTH_CASE_ALREADY_COMPLETE %s\n',id);continue;end
 [m,d,tag]=prepare_light_background(c,folder);
 modelCleanup=onCleanup(@()remove_model(tag));
 for j=1:numel(selected)
  s=selected{j};sf=s.sigma_f_S_m;target=fullfile(folder,sf_token(sf));
  if isfile(fullfile(target,'SUCCESS.json')),continue;end
  result=r3_truth(m,d,sf,target,'');
  result.round6_phase=phase;result.round6_target_only=true;
  if isfield(c,'layout'),result.round6_layout=c.layout;end
  if isfield(c,'configuration_role'),result.round6_role=c.configuration_role;end
  result.round6_full_length_mm=2*c.a_mm;
  result.round6_source_key=c.source_key;
  result.prediction_gate='All 32 predictions, applicability criteria and frozen settings verified saved before any new full target';
  r3_json(fullfile(target,'SUCCESS.json'),result);
  marginal=struct('sigma_f_S_m',sf,'Sf_m',mphint2(m,'ec.Ex^2+ec.Ey^2',2, ...
   'selection',d.sel.inclusion_geometry_domains,'intorder',8)*c.out_of_plane_thickness_m, ...
   'use','Posthoc dG/dsigma_f only; never used as a prediction input');
  r3_json(fullfile(target,'marginal_scalar.json'),marginal);
 end
 clear m d modelCleanup;
end
end

function [m,d,tag]=prepare_light_background(c,folder)
% No auxiliary basis responses. Only the matrices needed by r3_truth for
% posthoc same-mesh energy accounting are assembled here.
tag=['R6Truth_' strrep(char(java.util.UUID.randomUUID),'-','')];
cache=fullfile(folder,'background_responses.mat');
if isfile(cache)
 d=load(cache);assert(strcmp(d.c.id,c.id),'Background cache belongs to another case');
 m=mphload(fullfile(folder,'background.mph'),tag);return;
end
source=fullfile(folder,'source_background');if ~isfolder(source),mkdir(source);end
if ~isfile(fullfile(source,'SUCCESS.json'))
 r3_json(fullfile(source,'parameters.json'),c);
 assert(~any(strcmp(cell(com.comsol.model.util.ModelUtil.tags),'Rock2D')), ...
  'Rock2D tag is already owned by another operation; finish it before truth background construction');
 sourceCleanup=onCleanup(@()remove_model('Rock2D'));
 started=tic;r=solve_round2_case(c,source);r.status='success';r.elapsed_seconds=toc(started);
 r3_json(fullfile(source,'SUCCESS.json'),r);clear sourceCleanup;
end
sel=jsondecode(fileread(fullfile(source,'selections.json')));
m=mphload(fullfile(source,[c.id '.mph']),tag);
try
 r3_setup_materials(m,sel);
 started=tic;m.sol('sol1').runAll;backgroundSolveSeconds=toc(started);
 u0=mphgetu(m);u0=u0(:);xm=mphxmeshinfo(m);
 G0=mphint2(m,'reacf(V)',1,'selection',sel.positive_electrode,'method','summation');
 ground=mphint2(m,'reacf(V)',1,'selection',sel.ground_electrode,'method','summation');
 power=mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',2,'intorder',8)*c.out_of_plane_thickness_m;
 assert(abs(G0+ground)/G0<1e-6&&abs(G0-power)/G0<1e-6,'Background balance failed');
 S=mphmatrix(m,'sol1','out',{'K','N','uscale'},'rowscale','off', ...
  'initmethod','sol','initsol','zero','extractafter','v1');
 K0=S.K;assert(all(abs(S.uscale(:)-1)<1e-12),'Expected unscaled physical V');
 assert(all(sum(spones(S.N),2)<=1),'Expected direct voltage constraints');
 fixed=find(any(S.N,1));free=setdiff(1:numel(u0),fixed);free=free(:);
 assert(norm(K0-K0','fro')/norm(K0,'fro')<1e-10,'Background K symmetry');
 assert(abs(u0'*(K0*u0)-G0)/G0<1e-6,'Background energy scaling');
 m.param.set('r3bgscale','0');m.param.set('r3dlow','1[S/m]');DL=r3_matrix(m);
 m.param.set('r3dlow','0[S/m]');m.param.set('r3dhigh','1[S/m]');DH=r3_matrix(m);
 m.param.set('r3dhigh','0[S/m]');m.param.set('r3bgscale','1');m.param.set('r3unit','0[S/m]');
 % mphmatrix(...,'initsol','zero') can clear stored reaction information.
 % Restore and solve the passive background before saving any solved MPH or
 % fields. The mesh/DOF coordinates have not changed during matrix extraction.
 restoredTimer=tic;m.sol('sol1').runAll;restoredBackgroundSolveSeconds=toc(restoredTimer);
 restoredU=mphgetu(m);restoredX=mphxmeshinfo(m);
 restoredU=r3_map_dofs(restoredU,restoredX,xm);
 assert(norm(restoredU-u0)/max(norm(u0),1e-30)<1e-9,'Restored background differs from cached matrix reference');
 restoredG=mphint2(m,'reacf(V)',1,'selection',sel.positive_electrode,'method','summation');
 restoredGround=mphint2(m,'reacf(V)',1,'selection',sel.ground_electrode,'method','summation');
 restoredPower=mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',2,'intorder',8)*c.out_of_plane_thickness_m;
 assert(abs(restoredG-G0)/G0<1e-9&&abs(restoredG+restoredGround)/restoredG<1e-6&& ...
  abs(restoredG-restoredPower)/restoredG<1e-6,'Restored background balance failed');
 capV=zeros(1,2);caps={sel.capminus_domains,sel.capplus_domains};
 for k=1:2,capV(k)=mphint2(m,'V',2,'selection',caps{k},'intorder',8)/mphint2(m,'1',2,'selection',caps{k},'intorder',8);end
 provenance=struct('purpose','Truth-only same-mesh energy accounting after the whole prediction batch is frozen', ...
  'prediction_inputs_from_this_cache',false,'auxiliary_response_count',0,'background_solve_seconds',backgroundSolveSeconds, ...
  'restored_background_solve_seconds',restoredBackgroundSolveSeconds,'background_resolved_after_matrix_extraction',true, ...
  'source_folder',source,'units','t=1 m; DV=1 V; G in S; conductivity in S/m; geometry inputs in mm');
 mphsave(m,fullfile(folder,'background.mph'));
 bgfields=mpheval(m,{'V','ec.Ex','ec.Ey','ec.Jx','ec.Jy'},'edim',2,'refine',1,'smooth','none');
 save(fullfile(folder,'background_fields.mat'),'bgfields','-v7.3');
 r3_json(fullfile(folder,'parameters.json'),c);r3_json(fullfile(folder,'selections.json'),sel);
 summary=struct('case_id',c.id,'G0_S',G0,'cap_dV_V',diff(capV),'dofs',numel(u0), ...
  'reaction_flux_relative_error',abs(G0+ground)/G0,'energy_relative_error',abs(G0-power)/G0,'provenance',provenance);
 r3_json(fullfile(folder,'background_summary.json'),summary);
 areas=jsondecode(fileread(fullfile(source,'SUCCESS.json')));areas.round6_passive_background=true;
 assert(abs(areas.background_high_area_mm2-7200)<.001&&abs(areas.background_low_area_mm2-7200)<.001,'Actual background phase amounts changed');
 ellipseArea=pi*c.a_mm*c.b_mm;
 assert(abs(areas.original_ellipse_area_mm2-ellipseArea)<.001,'Actual ellipse area differs from its new source length');
 assert(abs(areas.original_ellipse_replaced_high_mm2-c.expected_replaced_high_area_mm2)<.001&& ...
  abs(areas.original_ellipse_replaced_low_mm2-c.expected_replaced_low_area_mm2)<.001,'Actual replacement areas differ from this geometry''s analytical values');
 assert(abs(areas.original_ellipse_replaced_high_mm2+areas.original_ellipse_replaced_low_mm2-ellipseArea)<.001, ...
  'Replacement phase areas do not cover the actual ellipse');
 % Finite cap contacts are measured for this length, never inherited from 60 mm.
 % A/B endcap equality is explicitly NOT required (40 mm changes B contacts).
 assert(abs(areas.capminus_high_fraction-c.expected_endcap_high_fraction_lower)<1e-4&& ...
  abs(areas.capplus_high_fraction-c.expected_endcap_high_fraction_upper)<1e-4, ...
  'Measured finite endcap contacts differ from the new geometry calculation');
 areas.round6_geometry_scope='Actual COMSOL areas for this source length and role; no A/B endcap matching assumption';
 areas.round6_full_length_mm=2*c.a_mm;
 areas.round6_source_key=c.source_key;
 r3_json(fullfile(folder,'actual_geometry_and_background.json'),areas);
 save(cache,'K0','DL','DH','u0','xm','free','fixed','c','sel','G0','capV','provenance','-v7.3');
 d=load(cache);
catch err
 remove_model(tag);rethrow(err);
end
end

function preflight_predictions(out,states,ready)
seen={};
for j=1:numel(states)
 s=states{j};assert(isfield(s,'case_id')&&isfield(s,'sigma_f_S_m'),'State needs case_id and sigma_f_S_m');
 key=[char(s.case_id) '/' sf_token(s.sigma_f_S_m)];assert(~ismember(key,seen),'Duplicate state in prediction gate');seen{end+1}=key; %#ok<AGROW>
 path=prediction_path(out,s,ready);assert(isfile(path),'Missing pre-target prediction: %s',path);
 p=jsondecode(fileread(path));
 assert(isstruct(p)&&isscalar(p)&&isfield(p,'case_id')&&isfield(p,'sigma_f_S_m')&&isfield(p,'methods'), ...
  'Prediction file must contain the completed state and method predictions');
 assert(strcmp(char(p.case_id),char(s.case_id))&&abs(p.sigma_f_S_m-s.sigma_f_S_m)<1e-12, ...
  'Pre-target prediction belongs to a different prescribed state');
end
fprintf('R6_PREDICTION_GATE_PASSED %d states\n',numel(states));
end

function path=prediction_path(out,s,ready)
if isfield(s,'prediction_file'),path=char(s.prediction_file);
else
 path=fullfile(out,'predictions',char(s.case_id),sf_token(s.sigma_f_S_m),'prediction.json');
 for field={'entries','predictions'}
  if ~isfield(ready,field{1}),continue;end
  entries=as_cells(ready.(field{1}));
  for k=1:numel(entries)
   e=entries{k};
   if isfield(e,'case_id')&&isfield(e,'sigma_f_S_m')&&strcmp(char(e.case_id),char(s.case_id))&&abs(e.sigma_f_S_m-s.sigma_f_S_m)<1e-12*max(1,s.sigma_f_S_m)
    if isfield(e,'prediction_file'),path=char(e.prediction_file);end
   end
  end
 end
end
path=absolute_path(out,path);
end

function write_completion(out,states,file)
paths=cell(size(states));
for j=1:numel(states)
 s=states{j};paths{j}=fullfile(out,'truth',char(s.case_id),sf_token(s.sigma_f_S_m),'SUCCESS.json');
 assert(isfile(paths{j}),'Incomplete target batch');
end
r3_json(file,struct('status','success','count',numel(states),'success_files',{paths}, ...
 'completed_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'))));
end

function cells=as_cells(value)
if iscell(value),cells=value;else,cells=num2cell(value);end
end
function id=case_id(entry)
if isfield(entry,'c'),id=char(entry.c.id);else,id=char(entry.id);end
end
function path=absolute_path(base,path)
if ~java.io.File(path).isAbsolute,path=fullfile(base,path);end
end
function inside=inside_folder(path,parent)
resolved=lower(char(java.io.File(path).getCanonicalPath));
base=lower(char(java.io.File(parent).getCanonicalPath));
inside=strcmp(resolved,base)||startsWith(resolved,[base filesep]);
end
function token=sf_token(sf)
token=['sf_' strrep(sprintf('%.10g',sf),'.','p')];
end
function remove_model(tag)
try,com.comsol.model.util.ModelUtil.remove(tag);catch,end
end
