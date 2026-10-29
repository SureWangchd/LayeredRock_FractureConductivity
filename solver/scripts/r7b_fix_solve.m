function r7b_fix_solve(root,action)
%R7B_SOLVE Prescribed native-ellipsoid CT T/H background/target closure.
% Uses the already connected single LiveLink service. Starts no service.
% Geometry is loaded from inputs/<sample>_<M>/geometry_ready.mph; this file
% never builds, remeshes, smooths, reclassifies, or modifies the geometry.
% Actions: pilot, main, refinements, all, S1_M1_0p18, S5_M1_10, etc.
% Fixed scientific target budget: 8 main + 8 representative refined states.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,action='pilot';end
out=fullfile(root,'runs','round7b_localfit');
assert(isfile(fullfile(out,'uniform_S1','PILOT_ELECTRICAL_SUCCESS.json')), ...
 'New cut-mesh uniform S1 electrical comparison must pass before formal CT targets');
jobs=prescribed_jobs(action);
for ji=1:size(jobs,1)
 sample=jobs{ji,1};variant=jobs{ji,2};water=jobs{ji,3};
 input=fullfile(out,'inputs',[sample '_' variant]);
 folder=fullfile(out,'states',[sample '_' variant]);
 if ~isfolder(folder),mkdir(folder);end
 wf=fullfile(folder,water_token(water));
 assert(isfile(fullfile(input,'GEOMETRY_SUCCESS.json')),'Missing actual geometry acceptance for %s %s',sample,variant);
 assert(isfile(fullfile(input,'geometry_ready.mph')),'Missing accepted physical mesh MPH');
 p=jsondecode(fileread(fullfile(input,'parameters.json')));
 sel=jsondecode(fileread(fullfile(input,'selections.json')));
 gate=jsondecode(fileread(fullfile(input,'GEOMETRY_SUCCESS.json')));
 % Every explicit-cut or fixed-midpoint mesh passes its own shared-face
 % check before solving. Native Refine's rejected meshes are never used.
 sourceFaces=fullfile(input,'all_faces_verification.json');
 assert(isfile(sourceFaces),'Missing full-volume shared-face check for this mesh');
 faceCheck=jsondecode(fileread(sourceFaces));
 assert(faceCheck.all_volume_faces_checked && faceCheck.all_face_topology_checks_passed,...
  'Current mesh has not passed the full-volume shared-face check');
 assert(sel.min_mesh_quality>0,'Nonpositive actual mesh quality cannot enter electrical solving');
 if isfield(gate,'status'),assert(~contains(lower(char(gate.status)),'fail'),'Geometry status is failure');end
 if isfield(gate,'geometry_passed'),assert(gate.geometry_passed,'Actual geometric checks have not passed');end
 assert(strcmp(p.sample,sample),'Wrong sample in geometry input');
 p.input_variant=variant;p.model_output_folder=folder;
 p.geometry_input_folder=input;p.geometry_model_path=fullfile(input,'geometry_ready.mph');
 p.geometry_acceptance_file=fullfile(input,'GEOMETRY_SUCCESS.json');
 p.id=sprintf('r7bfix_CT_%s_%s_native_ellipsoid',sample,variant);
 if ~isfield(p,'ellipsoid_representation') && isfield(p,'geometry_representation'),p.ellipsoid_representation=p.geometry_representation;end
 assert(isfield(p,'ellipsoid_representation'),'Geometry representation must be explicitly recorded');
 assert(~contains(lower(p.ellipsoid_representation),'centroid staircase'),'Staircase representation is disallowed');
 if ~isfield(p,'mesh_tag'),p.mesh_tag='mesh1';end
 if ~isfield(p,'geometry_tag'),p.geometry_tag='geom1';end
 assert(p.skeleton_sigma_S_m==1e-5,'Round 7b keeps the original skeleton conductivity');
 assert(p.potential_order==1,'Round 7b uses existing P1 potential convention');
 if isfile(fullfile(wf,'T_target','SUCCESS.json')) && isfile(fullfile(wf,'H_target','SUCCESS.json'))
  tr=jsondecode(fileread(fullfile(wf,'T_target','SUCCESS.json')));
  hr=jsondecode(fileread(fullfile(wf,'H_target','SUCCESS.json')));
  assert(isfile(tr.model_path)&&isfile(hr.model_path),'Successful target checkpoint MPH is missing');
  if ~isfile(fullfile(wf,'PAIR_SUCCESS.json'))
   tbg=jsondecode(fileread(fullfile(wf,'T_background','SUCCESS.json')));
   frozen=jsondecode(fileread(fullfile(wf,'kappa_frozen.json')));
   pair_summary(p,water,wf,tbg,frozen,tr,hr);
  end
  fprintf('R7B_REUSE %s %s water %.8g requested pair complete\n',sample,variant,water);continue;
 end
 tag='R7BCTSolve';cleanup=onCleanup(@()remove_model(tag)); %#ok<NASGU>
 try
  m=mphload(p.geometry_model_path,tag);
  com.comsol.model.util.ModelUtil.showProgress(fullfile(folder,'solver_progress.log'));
  sel=install_electrical(m,p,sel,folder);
  if ~isfolder(wf),mkdir(wf);end
  numericalTol=1e-10;if water>=10,numericalTol=1e-8;end
  % An optional stricter native tolerance is allowed, never a looser one.
  requestedTol=str2double(getenv('R7B_CT_TOL'));
  if isfinite(requestedTol)&&requestedTol>0&&requestedTol<numericalTol,numericalTol=requestedTol;end
  solver=r7_ct_native_solver(m,numericalTol);
  r3_json(fullfile(wf,'solver_configuration.json'),solver);
  tbg=solve_state(m,p,sel,water,0,0,p.skeleton_sigma_S_m,fullfile(wf,'T_background'),'T_background',NaN,true);
  frozen=calibrate(m,p,sel,water,wf,tbg);
  assert(isfile(fullfile(wf,'kappa_frozen.json')),'Background-only kappa must be committed before targets');
  tr=solve_state(m,p,sel,water,1,0,frozen.kappa_eq_S_m,fullfile(wf,'T_target'),'T_target',tbg.G_S,true);
  hr=solve_state(m,p,sel,water,1,1,frozen.kappa_eq_S_m,fullfile(wf,'H_target'),'H_target',frozen.G0_H_S,true);
  pair_summary(p,water,wf,tbg,frozen,tr,hr);
 catch ME
  fail=struct('status','FAILED','sample',sample,'input_variant',variant,'water_sigma_S_m',water,'action',action,...
   'message',ME.message,'identifier',ME.identifier,'report',getReport(ME,'extended','hyperlinks','off'),...
   'time_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
  r3_json(fullfile(folder,['FAILURE_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '.json']),fail);rethrow(ME);
 end
 clear m cleanup;
end
end

function jobs=prescribed_jobs(action)
main={'S1','M1',.18;'S1','M1',10;'S5','M1',.18;'S5','M1',10};
fine={'S1','M2',.18;'S1','M3',.18;'S5','M2',10;'S5','M3',10};
allJobs=[main;fine];
switch action
 case 'pilot',jobs=main(1,:);
 case 'main',jobs=main;
 case 'refinements',jobs=fine;
 case 'all',jobs=allJobs;
 otherwise
  jobs=cell(0,3);
  for k=1:size(allJobs,1)
   token=sprintf('%s_%s_%s',allJobs{k,1},allJobs{k,2},strrep(sprintf('%.10g',allJobs{k,3}),'.','p'));
   if strcmp(action,token),jobs=allJobs(k,:);break;end
  end
  assert(~isempty(jobs),'Action is outside fixed Round 7b scientific target budget: %s',action);
end
end

function sel=install_electrical(m,p,sel,folder)
co=m.component('comp1');co.sorder('linear');
me=co.mesh(p.mesh_tag);
assert(strcmp(char(me.geometricModel()),''),'Accepted geometry mesh must define its own topology');
co.geometricModel(['mesh/' p.mesh_tag]);
assert(strcmp(char(co.geometricModel()),['mesh/' p.mesh_tag]),'Component does not use accepted material mesh');
ids=reshape(sel.all_domains,1,[]);phase=reshape(sel.original_phase,1,[]);
region=reshape(sel.region,1,[]);ell=logical(reshape(sel.is_ellipse,1,[]));
assert(numel(ids)==numel(phase)&&numel(ids)==numel(region)&&numel(ids)==numel(ell),'Domain metadata lengths disagree');
assert(numel(unique(ids))==numel(ids)&&all(ismember(phase,[0 1]))&&all(ismember(region,[0 1 2])),'Invalid original CT labels or regions');
assert(all(region(ell)==0),'Actual ellipsoid escaped the protected core');
assert(isequal(sort(reshape(sel.inclusion_domains,1,[])),sort(ids(ell))),'Inclusion selection differs from accepted domain metadata');
assert(~isempty(sel.top)&&~isempty(sel.bottom)&&isempty(intersect(sel.top,sel.bottom)),'Invalid or overlapping electrodes');
assert(double(me.getNumElem('tet'))==sel.tetrahedra&&double(me.getNumVertex())==sel.vertices,'Accepted geometry mesh counts changed');
sel.domain_code=1+phase+2*region+6*double(ell);
sel.pore_domains=ids(phase==0);sel.skeleton_domains=ids(phase==1);
sel.remote_lower_domains=ids(region==1);sel.remote_upper_domains=ids(region==2);
m.param.set('sigma_water','0.18[S/m]');m.param.set('sigma_solid',sprintf('%.17g[S/m]',p.skeleton_sigma_S_m));
m.param.set('kappa',sprintf('%.17g[S/m]',p.skeleton_sigma_S_m));
m.param.set('frac_on','0');m.param.set('H_on','0');m.param.set('DV','1[V]');
assert(isempty(cell(co.material.tags)),'geometry_ready.mph must not contain preassigned target materials');
assert(isempty(cell(co.physics.tags)),'geometry_ready.mph must not contain preassigned target physics');
for code=reshape(unique(sel.domain_code),1,[])
 members=ids(sel.domain_code==code);pore=mod(code-1,2)==0;reg=floor(mod(code-1,6)/2);inc=code>=7;
 native='sigma_solid';if pore,native='sigma_water';end
 expr=native;
 if reg>0,expr=[native '+H_on*(kappa-' native ')'];end
 if inc,expr=[native '+frac_on*(sigma_water-' native ')'];end
 addmat(co,sprintf('r7b_code%d',code),members,expr);
end
% COMSOL physics create takes the component's geometry tag, even when its
% analyzed topology is supplied by a mesh. Assert the resulting physics
% selection includes every accepted material domain, not merely source CAD.
co.physics.create('ec','ConductiveMedia',p.geometry_tag);ec=co.physics('ec');
ec.prop('ShapeProperty').set('order_electricpotential',1);
actual=reshape(double(ec.selection.entities(3)),1,[]);
assert(isequal(sort(actual),sort(ids)),'Electric Currents is not attached to all accepted mesh domains');
ec.create('pot1','ElectricPotential',2);ec.feature('pot1').selection.set(sel.top);ec.feature('pot1').set('V0','DV');
ec.create('gnd1','Ground',2);ec.feature('gnd1').selection.set(sel.bottom);
assert(isempty(cell(m.study.tags)),'Geometry input unexpectedly contains studies');
m.study.create('std1');m.study('std1').create('stat','Stationary');
m.study('std1').feature('stat').set('usestol',true);m.study('std1').feature('stat').set('stol',1e-10);
% Explicitly select the accepted mesh for this physical geometry.
m.study('std1').feature('stat').mesh(p.geometry_tag,p.mesh_tag);
r3_json(fullfile(folder,'selections.json'),sel);r3_json(fullfile(folder,'parameters.json'),p);
r3_json(fullfile(folder,'electrical_setup.json'),struct('component_geometric_model',char(co.geometricModel()),...
 'mesh_geometric_model',char(me.geometricModel()),'physics_geometry_tag',p.geometry_tag,...
 'verified_all_physics_domains',actual,'geometry_changed_by_solver',false,'same_geometry_all_material_states',true));
end

function pair_summary(p,water,wf,tbg,frozen,tr,hr)
delta=hr.g_S-tr.g_S;
nearZero=abs(tr.g_S)<=1e-12*abs(tbg.G_S);
epsilon=NaN;if ~nearZero,epsilon=delta/tr.g_S;end
res=frozen.G0_H_S-tbg.G_S;
rr=struct('schema','r7b_same_mesh_TH_pair_v1','status','success','sample',p.sample,'input_variant',p.input_variant,...
 'water_sigma_S_m',water,'G0_T_S',tbg.G_S,'G0_H_S',frozen.G0_H_S,'Gf_T_S',tr.G_S,'Gf_H_S',hr.G_S,...
 'gT_S',tr.g_S,'gH_S',hr.g_S,'delta_gH_minus_gT_S',delta,'epsilon_delta_over_gT',epsilon,...
 'etaT',tr.g_S/tbg.G_S,'etaH',hr.g_S/frozen.G0_H_S,'whole_target_relative_difference',(hr.G_S-tr.G_S)/tr.G_S,...
 'background_match_residual_S',res,'background_match_residual_over_gT',safe_ratio(abs(res),abs(tr.g_S)),...
 'background_match_residual_over_delta',safe_ratio(abs(res),abs(delta)),...
 'unresolved_small_difference_chasing_disallowed',true,'near_zero_gT',nearZero,...
 'geometry_acceptance_file',p.geometry_acceptance_file,'geometry_model_path',p.geometry_model_path,...
 'T_target_model',tr.model_path,'H_target_model',hr.model_path,...
 'completed_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
r3_json(fullfile(wf,'PAIR_SUCCESS.json'),rr);writetable(struct2table(rr),fullfile(wf,'pair_response.csv'));
end
function y=safe_ratio(a,b),y=NaN;if b>0,y=a/b;end,end


function frozen=calibrate(m,p,sel,water,wf,tbg)
freezePath=fullfile(wf,'kappa_frozen.json');
if isfile(freezePath)
 frozen=jsondecode(fileread(freezePath));assert(strcmp(frozen.input_variant,p.input_variant),'Frozen calibration mesh changed');assert(abs(frozen.G0_T_S/tbg.G_S-1)<1e-13,'Frozen calibration baseline changed');return;
end
cal=fullfile(wf,'calibration');if ~isfolder(cal),mkdir(cal);end
tol=1e-7;maxSolves=40;target=tbg.G_S;
% Commit the initial-bracket rule before any new calibration solve. Existing
% step_001 without a rule always replays the original full-range schedule.
rule=calibration_rule(p,water,wf,cal,tol,maxSolves);
logLo=log(rule.initial_kappa_interval_S_m(1));logHi=log(rule.initial_kappa_interval_S_m(2));
history=struct([]);lo=run_point(1,logLo);hi=run_point(2,logHi);
best=lo;if abs(hi.G_S/target-1)<abs(best.G_S/target-1),best=hi;end
% A reused background kappa is only an initial guess. Establish the bracket
% using THIS fluid/mesh's true no-fracture responses, expanding as needed.
while ~(lo.G_S<=target*(1+tol) && hi.G_S>=target*(1-tol))
 assert(numel(history)<maxSolves,'CT background matching exhausted 40 solves during bracket expansion');
 if lo.G_S>target*(1+tol)
  newK=max(p.skeleton_sigma_S_m,exp(logLo)/2);
  assert(newK<exp(logLo)*(1-1e-14),'H background target below physical lower-bound response');
  hi=lo;logHi=logLo;logLo=log(newK);lo=run_point(numel(history)+1,logLo);trialPoint=lo;
 elseif hi.G_S<target*(1-tol)
  newK=min(water,exp(logHi)*2);
  assert(newK>exp(logHi)*(1+1e-14),'H background target above physical upper-bound response');
  lo=hi;logLo=logHi;logHi=log(newK);hi=run_point(numel(history)+1,logHi);trialPoint=hi;
 else,error('CT background bracketing encountered inconsistent endpoint responses');end
 if abs(trialPoint.G_S/target-1)<abs(best.G_S/target-1),best=trialPoint;end
end
for step=numel(history)+1:maxSolves
 if abs(best.G_S/target-1)<=tol,break;end
 % Safeguarded log(kappa)/log(G) secant; endpoints and every trial are
 % genuine no-fracture responses. Clip to the inner 80% of the bracket.
 denom=log(hi.G_S)-log(lo.G_S);fraction=.5;
 if isfinite(denom)&&abs(denom)>1e-14
  trial=(log(target)-log(lo.G_S))/denom;
  if isfinite(trial),fraction=min(.9,max(.1,trial));end
 end
 logMid=logLo+fraction*(logHi-logLo);mid=run_point(step,logMid);
 if abs(mid.G_S/target-1)<abs(best.G_S/target-1),best=mid;end
 if mid.G_S<target,logLo=logMid;lo=mid;else,logHi=logMid;hi=mid;end %#ok<NASGU>
end
residual=best.G_S/target-1;
assert(abs(residual)<=tol,'CT background matching not reached within 40 solves');
frozen=struct('schema','r7b_CT_background_only_kappa_v1','status','FROZEN_BEFORE_TARGETS','sample',p.sample,'input_variant',p.input_variant,...
 'fluid_and_fracture_sigma_S_m',water,'skeleton_sigma_S_m',p.skeleton_sigma_S_m,...
 'kappa_eq_S_m',best.kappa_S_m,'G0_T_S',target,'G0_H_S',best.G_S,'relative_background_mismatch',residual,...
 'absolute_background_mismatch_S',best.G_S-target,'required_relative_tolerance',tol,...
 'rule','Frozen initial interval, background-only bracket verification/expansion by factor2, then safeguarded log-kappa/log-G secant clipped0.1..0.9; bisection if degenerate; max40 total background solves; no target result read',...
 'calibration_rule_json',fullfile(wf,'calibration_rule.json'),'calibration_rule',rule,...
 'final_kappa_bracket_S_m',[exp(logLo) exp(logHi)],...
 'number_background_trials',numel(history),'H_background_solution_mph',best.model_path,...
 'H_background_success_json',best.success_path,'T_background_solution_mph',fullfile(wf,'T_background','solution.mph'),...
 'same_kappa_in_both_remote_regions',true,'full_tensor_matching_claimed',false,...
 'committed_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
frozen.geometry_acceptance_file=p.geometry_acceptance_file;frozen.geometry_model_path=p.geometry_model_path;
r3_json(freezePath,frozen);writetable(struct2table(history),fullfile(wf,'calibration_history.csv'));
fprintf('R7B_CT_KAPPA_FROZEN %s water %.8g kappa %.12g residual %.3g backgrounds %d\n',p.sample,water,best.kappa_S_m,residual,numel(history));
 function rr=run_point(index,logK)
  assert(index==numel(history)+1 && index<=maxSolves,'Calibration checkpoint indices must be consecutive within budget');
  k=exp(logK);path=fullfile(cal,sprintf('step_%03d',index));
  rr=solve_state(m,p,sel,water,0,1,k,path,'H_calibration',NaN,false);
  assert(abs(log(rr.kappa_S_m)-logK)<1e-12,'Stored calibration point does not match deterministic resume');
  row=struct('step',index,'kappa_S_m',k,'G_S',rr.G_S,'relative_background_mismatch',rr.G_S/target-1,...
   'solve_seconds',rr.solve_seconds,'model_path',rr.model_path);
  if isempty(history),history=row;else,history(end+1)=row;end %#ok<AGROW>
  r3_json(fullfile(wf,'calibration_progress.json'),history);
end
end

function rule=calibration_rule(p,water,wf,cal,tol,maxSolves)
% Only prior background calibration records are inspected; no old target
% truth or fields enter this new-mesh/new-fluid background matching.
path=fullfile(wf,'calibration_rule.json');
if isfile(path)
 rule=jsondecode(fileread(path));
 assert(strcmp(rule.sample,p.sample)&&rule.water_sigma_S_m==water&&strcmp(rule.input_variant,p.input_variant),'Saved calibration rule identity changed');
 assert(rule.skeleton_sigma_S_m==p.skeleton_sigma_S_m&&rule.relative_background_tolerance==tol&&rule.max_background_solves==maxSolves,'Saved background numerical rule changed');
 interval=rule.initial_kappa_interval_S_m;
 assert(numel(interval)==2&&all(isfinite(interval))&&interval(1)>=p.skeleton_sigma_S_m&&interval(2)<=water&&interval(1)<interval(2),'Invalid saved bracket');return;
end
legacy=isfolder(fullfile(cal,'step_001'));seed=NaN;source=struct();mode='full_physical_interval';
if legacy
 mode='existing_step001_without_rule_preserves_full_interval';
else
 % Prefer a preceding R7b mesh background, then old R7 same-sample same-
 % fluid background. Only a seed: both endpoints are recomputed on this mesh.
 out=fileparts(fileparts(p.model_output_folder));project=fileparts(fileparts(out));
 candidates={};kinds={};
 if ~strcmp(p.input_variant,'M1')
  candidates{end+1}=fullfile(out,'states',[p.sample '_M1'],water_token(water),'kappa_frozen.json');
  kinds{end+1}='r7b_same_sample_fluid_M1_background';
 end
 candidates{end+1}=fullfile(project,'runs','round7_3d_ct','ct',[p.sample '_base'],water_token(water),'kappa_frozen.json');
 kinds{end+1}='old_round7_same_sample_fluid_background_initial_guess_only';
 others=[.18 10];others=others(abs(others-water)>1e-12);
 for other=others
  candidates{end+1}=fullfile(p.model_output_folder,water_token(other),'kappa_frozen.json');
  kinds{end+1}='r7b_same_mesh_other_fluid_background_scaled_initial_guess';
 end
 for j=1:numel(candidates)
  candidate=candidates{j};if ~isfile(candidate),continue;end
  old=jsondecode(fileread(candidate));
  if ~isfield(old,'status')||~strcmp(old.status,'FROZEN_BEFORE_TARGETS')||~strcmp(old.sample,p.sample),continue;end
  if old.skeleton_sigma_S_m~=p.skeleton_sigma_S_m,continue;end
  oldWater=old.fluid_and_fracture_sigma_S_m;oldK=old.kappa_eq_S_m;
  if ~isscalar(oldWater)||~isfinite(oldWater)||oldWater<=0||~isscalar(oldK)||~isfinite(oldK)||oldK<=0,continue;end
  if ~isfield(old,'H_background_success_json')||~isfile(old.H_background_success_json),continue;end
  bg=jsondecode(fileread(old.H_background_success_json));
  if ~strcmp(bg.status,'success')||bg.fracture_on~=0||bg.H_on~=1,continue;end
  if bg.water_sigma_S_m~=oldWater||abs(log(bg.kappa_S_m/oldK))>1e-12||~isfile(bg.model_path),continue;end
  seed=min(water,max(p.skeleton_sigma_S_m,oldK*water/oldWater));
  source=struct('kind',kinds{j},'frozen_background_json',candidate,'H_background_success_json',old.H_background_success_json,...
   'source_water_sigma_S_m',oldWater,'source_kappa_S_m',oldK,'fluid_ratio_applied',water/oldWater,...
   'clipped_seed_kappa_S_m',seed,'source_contains_target_response',false);
  mode='background_only_initial_interval';break;
 end
end
interval=[p.skeleton_sigma_S_m water];if isfinite(seed),interval=[max(p.skeleton_sigma_S_m,seed/2),min(water,seed*2)];end
assert(interval(1)<interval(2),'Invalid positive calibration bracket');
rule=struct('schema','r7b_background_initial_bracket_v1','sample',p.sample,'input_variant',p.input_variant,...
 'water_sigma_S_m',water,'skeleton_sigma_S_m',p.skeleton_sigma_S_m,'mode',mode,...
 'physical_kappa_bounds_S_m',[p.skeleton_sigma_S_m water],'initial_kappa_interval_S_m',interval,'seed_background_source',source,...
 'existing_step001_without_rule',legacy,'endpoint_rule','Recompute both endpoints on current actual accepted geometry and current fluid, then expand by factor 2 if needed',...
 'iteration_rule','Safeguarded log-kappa/log-G secant clipped to 0.1..0.9; bisection if degenerate',...
 'relative_background_tolerance',tol,'max_background_solves',maxSolves,'target_results_read',false,...
 'geometry_acceptance_file',p.geometry_acceptance_file,'committed_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
r3_json(path,rule);
end

function r=solve_state(m,p,sel,water,fracOn,HOn,kappa,folder,label,G0,saveFields)
done=fullfile(folder,'SUCCESS.json');
if isfile(done)
 r=jsondecode(fileread(done));assert(strcmp(r.input_variant,p.input_variant),'Checkpoint mesh variant mismatch');assert(r.water_sigma_S_m==water && r.fracture_on==fracOn && r.H_on==HOn,'Checkpoint state mismatch');
 assert(abs(log(r.kappa_S_m/kappa))<1e-12,'Checkpoint kappa mismatch');assert(isfile(r.model_path),'Successful checkpoint MPH is missing');return;
end
if ~isfolder(folder),mkdir(folder);end
m.param.set('sigma_water',sprintf('%.17g[S/m]',water));m.param.set('kappa',sprintf('%.17g[S/m]',kappa));
m.param.set('frac_on',num2str(fracOn));m.param.set('H_on',num2str(HOn));
timer=tic;m.sol('sol1').runAll;seconds=toc(timer);
I=mphint2(m,'reacf(V)',2,'selection',sel.top,'method','summation');Ig=mphint2(m,'reacf(V)',2,'selection',sel.bottom,'method','summation');
expr='ec.Jx*ec.Ex+ec.Jy*ec.Ey+ec.Jz*ec.Ez';power=mphint2(m,expr,3,'intorder',2);
area=mphint2(m,'1',2,'selection',sel.top,'intorder',2,'unit','m^2');G=I;
q=mphint2(m,'ec.Jz',2,'selection',sel.top,'intorder',2);qg=-mphint2(m,'ec.Jz',2,'selection',sel.bottom,'intorder',2);
flux=abs(I+Ig)/abs(I);energy=abs(I-power)/abs(I);
assert(G>0 && abs(area/p.L_m^2-1)<1e-9,'CT electrode area/current check failed');
if max(flux,energy)>=1e-7
 mphsave(m,fullfile(folder,'failed_conservation_solution.mph'),'optimize','size');
 r3_json(fullfile(folder,'CONSERVATION_FAILED.json'),struct('status','solved_but_conservation_failed',...
  'state',label,'water_sigma_S_m',water,'kappa_S_m',kappa,'G_S',G,'reaction_ground_A',Ig,...
  'volume_dissipation_W',power,'reaction_flux_relative_error',flux,'energy_relative_error',energy,'solve_seconds',seconds));
 error('r7b:ct:conservation','CT reaction flux/energy balance >1e-7: %.3g %.3g',flux,energy);
end
r=struct('status','success','case_id',p.id,'sample',p.sample,'input_variant',p.input_variant,'state',label,'dimension',3,...
 'formal_CT_target',logical(fracOn),'water_sigma_S_m',water,'fracture_sigma_S_m',water,...
 'skeleton_sigma_S_m',p.skeleton_sigma_S_m,'fracture_on',fracOn,'H_on',HOn,'kappa_S_m',kappa,...
 'G_S',G,'G0_S',G0,'g_S',G-G0,'relative_gain',(G-G0)/G0,'sigma_z_app_S_m',G*p.L_m/area,...
 'L_m',p.L_m,'electrode_area_m2',area,'reaction_positive_A',I,'reaction_ground_A',Ig,...
 'reaction_flux_relative_error',flux,'energy_relative_error',energy,'volume_dissipation_W',power,...
 'gradient_positive_out_A',q,'gradient_ground_out_A',qg,'gradient_flux_relative_error',abs(q+qg)/max(abs([q qg])),...
 'source_CT',p.source,'source_raw',p.source_raw,'source_spacing_m',p.spacing_m,...
 'geometry_representation',p.ellipsoid_representation,'geometry_acceptance_file',p.geometry_acceptance_file,'geometry_model_path',p.geometry_model_path,'potential_order',1,'geometry_order','linear_on_frozen_accepted_ellipsoid_surface',...
 'same_pore_and_fracture_fluid',true,'background_recomputed_for_each_fluid',true,...
 'tetrahedra',sel.tetrahedra,'vertices',sel.vertices,'min_mesh_quality',sel.min_mesh_quality,...
 'solve_seconds',seconds,'stationary_relative_tolerance',double(m.sol('sol1').feature('s1').getDouble('stol')),'comsol_version','6.3.0.290',...
 'model_path',fullfile(folder,'solution.mph'),'success_path',done,'completed_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
% Geometry is identical for every fluid/calibration/target state. Perform
% actual COMSOL phase-volume integrals once, then reuse that explicit record.
volumePath=fullfile(p.model_output_folder,'actual_integrated_volumes.json');
volumeReused=isfile(volumePath);
if volumeReused
 volumeRecord=jsondecode(fileread(volumePath));vc=reshape(volumeRecord.encoded_state,1,[]);vol=reshape(volumeRecord.volume_m3,1,[]);
else
 ids=reshape(sel.all_domains,1,[]);codes=reshape(sel.domain_code,1,[]);vc=reshape(unique(codes),1,[]);vol=zeros(size(vc));
 for j=1:numel(vc),vol(j)=mphint2(m,'1',3,'selection',ids(codes==vc(j)),'intorder',2,'unit','m^3');end
 volumeRecord=struct('encoded_state',vc,'volume_m3',vol,'first_evaluated_state',label,...
  'method','COMSOL mphint2 on the accepted boundary-fitted fixed geometry, reused only on this same M-level mesh');
 r3_json(volumePath,volumeRecord);
end
r.total_volume_m3=sum(vol);r.actual_inclusion_volume_m3=sum(vol(vc>=7));r.original_pore_volume_m3=sum(vol(mod(vc-1,2)==0));
r.phase_volume_source_file=volumePath;r.phase_volume_reused=volumeReused;
r.original_skeleton_volume_m3=sum(vol(mod(vc-1,2)==1));
r.actual_replaced_skeleton_volume_m3=sum(vol(vc>=7 & mod(vc-1,2)==1));
regions=floor(mod(vc-1,6)/2);homogeneous=logical(HOn)&regions>0;
fluid=(mod(vc-1,2)==0 | (logical(fracOn)&vc>=7)) & ~homogeneous;
r.final_fluid_volume_m3=sum(vol(fluid));
r.final_skeleton_volume_m3=sum(vol(~fluid & ~homogeneous));
r.final_uniform_remote_volume_m3=sum(vol(homogeneous));
r.material_volume_note='Original CT pore/skeleton volumes and ellipsoid overlap are bookkeeping; H replaces remote material by kappa. Final phase volumes above account for that replacement without double-counting pore/fracture fluid.';
assert(abs(r.total_volume_m3/p.total_volume_m3-1)<1e-8,'CT physical-volume integration mismatch');
r.original_pore_volume_relative_error=r.original_pore_volume_m3/p.background_pore_volume_m3-1;
assert(abs(r.original_pore_volume_relative_error)<1e-8,'Native CT pore fraction changed');
writetable(table(vc(:),vol(:),'VariableNames',{'encoded_state','actual_volume_m3'}),fullfile(folder,'phase_volumes.csv'));
% Preserve all mesh and solution data; native size compression controls the
% disk cost of the prescribed calibration and fixed-geometry refinement.
mphsave(m,r.model_path,'optimize','size');
if saveFields
 u=mphgetu(m);save(fullfile(folder,'solution_vector.mat'),'u','-v7.3');clear u;
 save_slices(m,p,folder);
 r.inclusion_dissipation_W=mphint2(m,expr,3,'selection',sel.inclusion_domains,'intorder',2);
 r.lower_remote_dissipation_W=mphint2(m,expr,3,'selection',sel.remote_lower_domains,'intorder',2);
 r.upper_remote_dissipation_W=mphint2(m,expr,3,'selection',sel.remote_upper_domains,'intorder',2);
 r.Jy_squared_volume_integral_A2_per_m=mphint2(m,'ec.Jy^2',3,'intorder',2);
end
params=p;params.current_state=label;params.current_fluid_sigma_S_m=water;params.current_kappa_S_m=kappa;params.fracture_on=fracOn;params.H_on=HOn;
r3_json(fullfile(folder,'parameters.json'),params);r3_json(done,r);
fprintf('R7B_CT_SOLVED %s %s water %.8g kappa %.8g G %.12g gain %.12g flux %.2g energy %.2g seconds %.2f\n',p.sample,label,water,kappa,G,G-G0,flux,energy,seconds);
end

function save_slices(m,p,folder)
L=p.L_m;[x,z]=meshgrid(linspace(0,L,129));ys=L*[.375 .5 .625];expr={'V','ec.Ex','ec.Ey','ec.Ez','ec.Jx','ec.Jy','ec.Jz'};slices=struct([]);
for j=1:3
 xyz=[x(:)';repmat(ys(j),1,numel(x));z(:)'];values=cell(1,numel(expr));
 [values{:}]=mphinterp(m,expr,'coord',xyz,'coorderr','off');
 slices(j).y_m=ys(j);slices(j).coordinates_m=xyz;slices(j).expressions=expr;slices(j).values=vertcat(values{:});
 T=array2table([xyz' vertcat(values{:})'],'VariableNames',{'x_m','y_m','z_m','V_V','Ex_V_m','Ey_V_m','Ez_V_m','Jx_A_m2','Jy_A_m2','Jz_A_m2'});
 writetable(T,fullfile(folder,sprintf('slice_y_%03d_permil.csv',round(1000*ys(j)/L))));
end
save(fullfile(folder,'slices.mat'),'slices','-v7.3');
end

function addmat(co,tag,ids,expr)
co.material.create(tag,'Common');ma=co.material(tag);ma.selection.set(reshape(ids,1,[]));ma.propertyGroup('def').set('electricconductivity',{expr '0' '0' '0' expr '0' '0' '0' expr});
end
function token=water_token(w),token=['water_' strrep(sprintf('%.10g',w),'.','p')];end
function remove_model(tag),try,com.comsol.model.util.ModelUtil.remove(tag);catch,end,end
