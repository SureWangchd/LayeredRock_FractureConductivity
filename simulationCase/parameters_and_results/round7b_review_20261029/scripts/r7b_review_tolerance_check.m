function result=r7b_review_tolerance_check(root,successJson,newTol)
%R7B_FIX_TOLERANCE_CHECK One already-solved NEW CT state, fixed mesh/physics.
% Called on the root's existing serial COMSOL LiveLink connection. This
% helper neither connects nor starts a service. It does not invoke the main
% solve/checkpoint loop, recalibrate kappa, or overwrite a saved target.
% Example, after an actual target SUCCESS exists:
% r7b_review_tolerance_check(root,fullfile(root,'runs','round7b_localfit',...
%   'states','S1_M1','water_0p18','T_target','SUCCESS.json'));
% Review adapter: 1e-8 -> 1e-10 or 1e-10 -> 2e-11 is allowed.
% The original 1e-11 attempt remains separately uncompleted in its own folder.
% Also accepts the real T_background and exactly the H calibration state
% named by this pair's kappa_frozen.H_background_success_json. It never
% selects another calibration trial or changes the frozen kappa.
% A single state's change DOES NOT establish the whole T/H gain-difference
% solver uncertainty. The report explicitly leaves that applicability false.
if nargin<1||isempty(root),root=fileparts(fileparts(mfilename('fullpath')));end
assert(nargin>=2 && isfile(successJson),'Supply one actual completed local-fit state SUCCESS.json');
out=fullfile(root,'runs','round7b_localfit');
jroot=javaObject('java.io.File',out);jinput=javaObject('java.io.File',successJson);
canonicalRoot=char(jroot.getCanonicalPath());canonicalInput=char(jinput.getCanonicalPath());
assert(startsWith(lower(canonicalInput),[lower(canonicalRoot) filesep]),'Only new round7b_localfit results may be checked');
old=jsondecode(fileread(successJson));
assert(strcmp(old.status,'success'),'Expected one successful actual CT state');
assert(any(strcmp(old.sample,{'S1','S5'})) && any(strcmp(old.input_variant,{'M1','M2','M3'})),'Unexpected sample or mesh scope');
assert(any(old.water_sigma_S_m==[.18 10]) && old.skeleton_sigma_S_m==1e-5,'Unexpected fluid/skeleton material parameters');
state=char(old.state);role=state;
wf=fullfile(out,'states',[old.sample '_' old.input_variant],['water_' strrep(sprintf('%.10g',old.water_sigma_S_m),'.','p')]);
if strcmp(state,'T_target')
 expectedFrac=1;expectedH=0;expectedJson=fullfile(wf,'T_target','SUCCESS.json');
elseif strcmp(state,'H_target')
 expectedFrac=1;expectedH=1;expectedJson=fullfile(wf,'H_target','SUCCESS.json');
elseif strcmp(state,'T_background')
 expectedFrac=0;expectedH=0;expectedJson=fullfile(wf,'T_background','SUCCESS.json');
elseif strcmp(state,'H_calibration')
 role='H_background';expectedFrac=0;expectedH=1;
 frozenPath=fullfile(wf,'kappa_frozen.json');assert(isfile(frozenPath),'H background must already be selected by an actual frozen calibration');
 frozen=jsondecode(fileread(frozenPath));
 assert(strcmp(frozen.status,'FROZEN_BEFORE_TARGETS') && strcmp(frozen.sample,old.sample) && strcmp(frozen.input_variant,old.input_variant),...
  'Frozen H background sample/mesh/status mismatch');
 assert(frozen.fluid_and_fracture_sigma_S_m==old.water_sigma_S_m && abs(log(old.kappa_S_m/frozen.kappa_eq_S_m))<1e-12,...
  'Only the selected frozen H background may be tolerance checked');
 expectedJson=frozen.H_background_success_json;
 assert(startsWith(lower(canonical_path(expectedJson)),[lower(canonical_path(fullfile(wf,'calibration'))) filesep]),...
  'Frozen H background must be a calibration state in this exact local-fit pair');
 assert(strcmpi(canonical_path(old.model_path),canonical_path(frozen.H_background_solution_mph)),...
  'Selected H background MPH does not match the frozen calibration record');
 assert(abs(old.G_S/frozen.G0_H_S-1)<1e-12,'Frozen H background response differs from saved selected state');
else,error('Only actual T/H targets, T_background and the selected frozen H_background are allowed');end
assert(strcmpi(canonicalInput,canonical_path(expectedJson)),'The supplied state is not the exact current pair state');
assert(old.fracture_on==expectedFrac && old.H_on==expectedH && logical(old.formal_CT_target)==logical(expectedFrac),...
 'Saved state identity and fracture/background flags disagree');
assert(isfile(old.model_path),'State solution MPH missing');
assert(strcmpi(canonical_path(old.model_path),canonical_path(fullfile(fileparts(expectedJson),'solution.mph'))),...
 'State MPH is not its own local-fit solution file');
assert(isfield(old,'stationary_relative_tolerance') && isscalar(old.stationary_relative_tolerance),'Saved stationary tolerance missing');
previousTol=double(old.stationary_relative_tolerance);
if abs(log10(previousTol)+8)<1e-9,expectedTol=1e-10;
elseif abs(log10(previousTol)+10)<1e-9,expectedTol=2e-11;
else,error('This bounded helper only checks saved 1e-8 or 1e-10 targets');end
if nargin<3||isempty(newTol),newTol=expectedTol;end
assert(isscalar(newTol) && isfinite(newTol) && abs(log10(newTol/expectedTol))<1e-12,...
 'Allowed stricter tolerance is %.12g for this saved target',expectedTol);
input=fullfile(out,'inputs',[old.sample '_' old.input_variant]);
assert(isfile(fullfile(input,'GEOMETRY_SUCCESS.json')) && isfile(fullfile(input,'selections.json')),'Accepted geometry evidence missing');
sel=jsondecode(fileread(fullfile(input,'selections.json')));
assert(strcmpi(canonical_path(old.geometry_model_path),canonical_path(fullfile(input,'geometry_ready.mph'))),...
 'State geometry path does not match this accepted local-fit mesh');
token=sprintf('%s_%s_water_%s_%s_fgmres_tol_%s',old.sample,old.input_variant,...
 strrep(sprintf('%.10g',old.water_sigma_S_m),'.','p'),role,strrep(sprintf('%.0e',newTol),'-','m'));
folder=fullfile(root,'runs','round7b_review_20261029','tolerance_checks',token);
done=fullfile(folder,'comparison.json');
if isfile(done)
 result=jsondecode(fileread(done));
 assert(strcmp(result.status,'success') && strcmp(result.source_success_json,canonicalInput),...
  'Existing comparison does not match this target');
 assert(result.strict_stationary_tolerance==newTol && isfile(result.strict_model_path),'Saved stricter comparison/MPH mismatch');
 fprintf('R7BFIX_TOLERANCE_REUSE %s\n',done);return;
end
if ~isfolder(folder),mkdir(folder);end
import com.comsol.model.util.*
tag='R7BfixTolerance';assert(~any(strcmp(cell(ModelUtil.tags),tag)),'Tolerance model tag already in use');
cleanup=onCleanup(@()remove_model(tag)); %#ok<NASGU>
started=tic;
try
 m=mphload(old.model_path,tag);co=m.component('comp1');me=co.mesh('mesh1');
 ModelUtil.showProgress(fullfile(folder,'solver_progress.log'));
 assert(strcmp(char(co.geometricModel()),'mesh/mesh1'),'Target is not attached to the accepted local-fit mesh');
 countsBefore=[double(me.getNumVertex()),double(me.getNumElem('tet'))];
 assert(isequal(countsBefore,[sel.vertices sel.tetrahedra]),'Saved target mesh counts differ from accepted geometry');
 assert(old.vertices==countsBefore(1)&&old.tetrahedra==countsBefore(2),'Saved target metadata differs from actual MPH');
 paramNames={'sigma_water','sigma_solid','kappa','frac_on','H_on','DV'};
 before=struct();for k=1:numel(paramNames),before.(paramNames{k})=double(m.param.evaluate(paramNames{k}));end
 assert(abs(before.sigma_water/old.water_sigma_S_m-1)<1e-12 && before.sigma_solid==old.skeleton_sigma_S_m,...
  'Target conductivity parameters do not match saved result');
 assert(abs(log(before.kappa/old.kappa_S_m))<1e-12 && before.frac_on==expectedFrac && before.H_on==expectedH && before.DV==1,...
  'Saved target state/kappa/loading mismatch');
 assert(abs(log(double(m.sol('sol1').feature('s1').getDouble('stol'))/previousTol))<1e-12,...
  'Actual saved MPH tolerance differs from result record');
 oldI=mphint2(m,'reacf(V)',2,'selection',sel.top,'method','summation');
 assert(abs(oldI/old.G_S-1)<1e-10,'Loaded saved solution conductance does not reproduce its SUCCESS record');
 r3_json(fullfile(folder,'request_before_solve.json'),struct('sample',old.sample,'input_variant',old.input_variant,...
  'state',state,'pair_state_role',role,'source_success_json',canonicalInput,'source_model_path',old.model_path,...
  'old_stationary_tolerance',previousTol,'strict_stationary_tolerance',newTol,...
  'parameters_before',before,'vertices',countsBefore(1),'tetrahedra',countsBefore(2),...
  'mesh_geometry_materials_loading_kappa_fixed',true,'background_recalibrated',false,...
  'main_target_files_overwritten',false,'single_state_only',true,'single_target_only',logical(expectedFrac),'applies_to_pair_difference',false));
 solver=r7_ct_native_solver(m,newTol);
 nativeIter=m.sol('sol1').feature('s1').feature(solver.iterative_feature);
 nativeIter.set('linsolver','fgmres');nativeIter.set('prefuntype','right');
 assert(strcmp(char(nativeIter.getString('linsolver')),'fgmres')&&strcmp(char(nativeIter.getString('prefuntype')),'right'),...
  'Native FGMRES/right-preconditioning settings were not accepted');
 solver.linear_solver='fgmres';solver.actual_iterative_properties.prefuntype='right';
 solver.check_type='fixed-mesh solver crosscheck: GMRES/left 1e-10 to FGMRES/right 2e-11';
 solver.tolerance_only_comparison=false;
 % Print per-linear-iteration convergence; this changes logging only.
 st=m.sol('sol1').feature('s1');children=cell(st.feature.tags);advancedTag='';
 for jj=1:numel(children)
  if strcmp(char(st.feature(children{jj}).getType),'Advanced'),advancedTag=children{jj};break;end
 end
 if isempty(advancedTag),advancedTag='adReview';st.create(advancedTag,'Advanced');end
 st.feature(advancedTag).set('convinfo','detailed');
 solver.convergence_logging='detailed';
 % Reuse the accepted field only as an initial guess. sol1 is cleared and
 % runAll executes a new solve with the tightened tolerance on the same mesh.
 initialTag='solReviewInitial';
 assert(~any(strcmp(cell(m.sol.tags),initialTag)),'Initial guess tag already exists');
 sourceU=double(m.sol('sol1').getU());
 m.sol('sol1').copySolution(initialTag);
 initialU=double(m.sol(initialTag).getU());
 assert(isequal(size(sourceU),size(initialU))&&isequal(sourceU,initialU),...
  'The copied initial solution must equal every accepted degree of freedom');
 assert(max(abs(initialU))>0,'Copied initial potential field is zero');
 sequenceTags=cell(m.sol('sol1').feature.tags);variableTag='';
 for jj=1:numel(sequenceTags)
  if strcmp(char(m.sol('sol1').feature(sequenceTags{jj}).getType),'Variables'),variableTag=sequenceTags{jj};break;end
 end
 assert(~isempty(variableTag),'Variables feature is required for a saved initial guess');
 vv=m.sol('sol1').feature(variableTag);vv.set('control','user');
 variableProperties=cell(vv.properties);
 if any(strcmp(variableProperties,'useinitsol')),vv.set('useinitsol','on');end
 vv.set('initmethod','sol');vv.set('initsol',initialTag);vv.set('initsoluse','current');
 assert(strcmp(char(vv.getString('initmethod')),'sol')...
  &&strcmp(char(vv.getString('initsol')),initialTag),'Saved-field initialization settings were not accepted');
 solver.initial_guess='copied accepted source solution; sol1 cleared and genuinely re-solved';
 solver.initial_guess_solution_tag=initialTag;
 solver.initial_guess_source_model=old.model_path;
 solver.initial_guess_DOF_count=numel(initialU);
 solver.initial_guess_linf=max(abs(initialU));
 solver.initial_guess_difference_from_saved_linf=max(abs(initialU-sourceU));
 solver.initial_guess_variables_feature=variableTag;
 solver.initial_guess_available_properties=variableProperties;
 solver.initial_guess_initmethod=char(vv.getString('initmethod'));
 solver.initial_guess_initsol=char(vv.getString('initsol'));
 fprintf('R7B_REVIEW_INITIAL_GUESS %s: %d DOFs, max potential %.17g, max source difference %.3g\n',...
  role,numel(initialU),solver.initial_guess_linf,solver.initial_guess_difference_from_saved_linf);
 clear sourceU initialU;
 r3_json(fullfile(folder,'strict_solver_configuration.json'),solver);
 % Force a genuine solve; do not mistake a stored converged field for a new
 % tolerance measurement. No geometry/mesh/study rebuild is requested.
 m.sol('sol1').clearSolutionData();
 timer=tic;m.sol('sol1').runAll();seconds=toc(timer);
 after=struct();for k=1:numel(paramNames),after.(paramNames{k})=double(m.param.evaluate(paramNames{k}));end
 assert(isequal(before,after),'A physical parameter changed during tolerance check');
 assert(isequal(countsBefore,[double(me.getNumVertex()),double(me.getNumElem('tet'))]),'Mesh counts changed during tolerance check');
 I=mphint2(m,'reacf(V)',2,'selection',sel.top,'method','summation');
 Ig=mphint2(m,'reacf(V)',2,'selection',sel.bottom,'method','summation');
 power=mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey+ec.Jz*ec.Ez',3,'intorder',2);
 area=mphint2(m,'1',2,'selection',sel.top,'unit','m^2','intorder',2);
 G=I/before.DV;flux=abs(I+Ig)/abs(I);energy=abs(before.DV*I-power)/abs(before.DV*I);
 strictPath=fullfile(folder,'solution_strict.mph');mphsave(m,strictPath,'optimize','size');
 scale=abs(G-old.G_S);
 oldG0=number_or_null(old,'G0_S');oldGain=number_or_null(old,'g_S');
 strictGain=NaN;relativeGainChange=NaN;gainScale=NaN;targetChange=NaN;targetScale=NaN;
 if expectedFrac==1
  if isfinite(oldG0),strictGain=G-oldG0;end
  if isfinite(oldGain)&&oldGain~=0,relativeGainChange=scale/abs(oldGain);end
  gainScale=scale;targetChange=G-old.G_S;targetScale=scale;
 end
 result=struct('schema','r7b_localfit_single_state_tolerance_v2','status','COMPUTED',...
  'sample',old.sample,'input_variant',old.input_variant,'water_sigma_S_m',old.water_sigma_S_m,'state',state,'pair_state_role',role,...
  'source_success_json',canonicalInput,'source_model_path',old.model_path,'strict_model_path',strictPath,...
  'old_stationary_tolerance',previousTol,'strict_stationary_tolerance',newTol,...
  'old_G_S',old.G_S,'strict_G_S',G,'conductance_change_S',G-old.G_S,'target_conductance_change_S',targetChange,...
  'observed_conductance_change_scale_S',scale,'observed_target_conductance_change_scale_S',targetScale,'observed_gain_or_pair_change_scale_S',gainScale,...
  'old_G0_S',oldG0,'old_g_S',oldGain,'strict_g_using_same_saved_background_S',strictGain,...
  'relative_target_conductance_change',scale/abs(old.G_S),...
  'relative_gain_change_with_saved_background',relativeGainChange,...
  'reaction_positive_A',I,'reaction_ground_A',Ig,'dissipation_W',power,'electrode_area_m2',area,...
  'reaction_flux_relative_error',flux,'energy_relative_error',energy,...
  'vertices',countsBefore(1),'tetrahedra',countsBefore(2),'parameters_before',before,'parameters_after',after,...
  'native_solver',solver,'solve_seconds',seconds,'total_seconds',toc(started),...
  'evidence_model_paths',{{old.model_path,strictPath}},'source_geometry_modified',false,...
  'background_recalibrated',false,'background_tolerance_rechecked',logical(expectedFrac==0),'single_state_only',true,'single_target_only',logical(expectedFrac),...
  'fracture_on',expectedFrac,'H_on',expectedH,'kappa_S_m',old.kappa_S_m,...
  'applies_to_pair_difference',false,'observed_scale_is_rigorous_bound',false,...
  'scope_note','One actual state at fixed accepted mesh, materials, loading and its saved kappa. Native GMRES/left at 1e-10 is compared with native FGMRES/right at 2e-11; this is a solver crosscheck, not a pure tolerance-only comparison. A single-state record does not measure the complete T/H gain-difference solver scale.',...
  'completed_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
 r3_json(fullfile(folder,'computed.json'),result);
 assert(isfinite(G)&&G>0&&abs(area/old.L_m^2-1)<1e-9,'Stricter target current/electrode area check failed');
 assert(max(flux,energy)<1e-7,'Stricter target reaction flux/energy balance failed');
 result.status='success';r3_json(done,result);
 fprintf('R7BFIX_TOLERANCE_SOLVED %s %s %s sigma %.8g tol %.3g->%.3g dG %.12g seconds %.2f\n',...
  old.sample,old.input_variant,state,old.water_sigma_S_m,previousTol,newTol,G-old.G_S,seconds);
catch ME
 r3_json(fullfile(folder,['FAILED_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')) '.json']),...
  struct('status','FAILED','source_success_json',canonicalInput,'source_model_path',old.model_path,...
  'requested_strict_tolerance',newTol,'message',ME.message,'identifier',ME.identifier,...
  'report',getReport(ME,'extended','hyperlinks','off'),'main_target_files_overwritten',false));
 rethrow(ME);
end
end

function remove_model(tag),try,com.comsol.model.util.ModelUtil.remove(tag);catch,end,end
function p=canonical_path(path),f=javaObject('java.io.File',path);p=char(f.getCanonicalPath());end
function x=number_or_null(s,name)
x=NaN;if isfield(s,name)&&isnumeric(s.(name))&&isscalar(s.(name))&&isfinite(s.(name)),x=double(s.(name));end
end
