function r4_periphery(root,action)
% Task-authorized equal-area peripheral exchange; frozen R3 prediction only.
import com.comsol.model.util.*
out=fullfile(root,'runs','round4_crossover');manifest=jsondecode(fileread(fullfile(out,'periphery_manifest.json')));cases=manifest.cases;
if ~iscell(cases),cases=num2cell(cases);end
designFile=fullfile(out,'periphery_initial_design.json');
if ~isfile(designFile)
 original=jsondecode(fileread(fullfile(out,'original_brackets.json')));
 sigmas=unique([.1 .3 1 10 reshape(original.brackets_S_m,1,[])]);
 design=struct('sigma_f_S_m',sigmas,'chosen_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')),...
  'selection_basis','Original strict COMSOL crossover bracket plus 0.1,0.3,1,10; before any new periphery target solve',...
  'method_source',fullfile(root,'runs','round3_predict','method_frozen','r3_enrich_predict_v1.m'),...
  'method_version','r3_background_vector4_residual_rankguard_v1','all_initial_predictions_before_all_initial_targets',true);
 r3_json(designFile,design);
else,design=jsondecode(fileread(designFile));sigmas=reshape(design.sigma_f_S_m,1,[]);end
if strcmp(action,'predict')
 ready=fullfile(out,'periphery_prediction_batch_ready.json');
 if isfile(ready),fprintf('PERIPHERY_PREDICTION_BATCH_ALREADY_FROZEN\n');return;end
 assert(isempty(dir(fullfile(out,'periphery','*','sf_*','SUCCESS.json'))),'New target truth already exists; do not rewrite initial predictions');
 count=0;
 for j=1:numel(cases)
  c=cases{j};folder=fullfile(out,'periphery',c.id);source=fullfile(folder,'source_background');
  if ~isfolder(source),mkdir(source);end
  if ~isfile(fullfile(source,'SUCCESS.json'))
   assert(c.passive_ellipse,'Background only before frozen predictions');
   r3_json(fullfile(source,'parameters.json'),c);ts=tic;r=solve_round2_case(c,source);
   r.status='success';r.elapsed_seconds=toc(ts);r3_json(fullfile(source,'SUCCESS.json'),r);ModelUtil.clear;
  end
  [m,d]=r3_prepare(folder,source,false);
  bg=jsondecode(fileread(fullfile(source,'SUCCESS.json')));
  matchfile=fullfile(folder,'background_matching.json');
  if ~isfile(matchfile)
   ec=m.component('comp1').physics('ec');ec.feature('pot1').selection.set(d.sel.right_boundaries);ec.feature('gnd1').selection.set(d.sel.left_boundaries);
   ts=tic;m.sol('sol1').runAll;gxSeconds=toc(ts);
   gx=mphint2(m,'reacf(V)',1,'selection',d.sel.right_boundaries,'method','summation');
   neg=mphint2(m,'reacf(V)',1,'selection',d.sel.left_boundaries,'method','summation');
   pow=mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',2,'intorder',8);
   mphsave(m,fullfile(folder,'background_x.mph'));
   fields=mpheval(m,{'V','ec.Ex','ec.Ey'},'edim',2,'refine',1,'smooth','none');save(fullfile(folder,'background_x_fields.mat'),'fields','-v7.3');
   originalCache=fullfile(out,'original',c.configuration_role,'background_responses.mat');ori=load(originalCache,'G0','capV');
   match=struct('case_id',c.id,'role',c.configuration_role,'layout',c.layout,'G0z_S',d.G0,'G0x_S',gx,...
    'G0z_original_S',ori.G0,'G0x_analytic_S',.055,'G0z_relative_to_original',abs(d.G0/ori.G0-1),...
    'G0x_relative_to_analytic',abs(gx/.055-1),'cap_dV_V',diff(d.capV),'original_cap_dV_V',diff(ori.capV),...
    'background_high_area_mm2',bg.background_high_area_mm2,'background_low_area_mm2',bg.background_low_area_mm2,...
    'ellipse_area_mm2',bg.original_ellipse_area_mm2,'replaced_high_area_mm2',bg.original_ellipse_replaced_high_mm2,...
    'replaced_low_area_mm2',bg.original_ellipse_replaced_low_mm2,'capminus_high_fraction',bg.capminus_high_fraction,...
    'capplus_high_fraction',bg.capplus_high_fraction,'x_flux_relative_error',abs(gx+neg)/gx,...
    'x_energy_relative_error',abs(gx-pow)/gx,'x_solve_seconds',gxSeconds,...
    'background_response_tolerance_relative',1e-6,'phase_area_tolerance_mm2',.001,...
    'areas_source','Actual COMSOL domain quadrature on new mesh, not analytic candidate areas');
   assert(match.G0z_relative_to_original<1e-6 && match.G0x_relative_to_analytic<1e-6,'Background matching failed');
   assert(abs(match.background_high_area_mm2-7200)<.001 && abs(match.background_low_area_mm2-7200)<.001,'Background phase areas changed');
   assert(abs(match.replaced_high_area_mm2-pi*30*1.5/2)<.001 && abs(match.replaced_low_area_mm2-pi*30*1.5/2)<.001,'Ellipse phase areas changed');
   assert(abs(match.capminus_high_fraction-1)<1e-8 && abs(match.capplus_high_fraction)<1e-8,'Endcap phases changed');
   match.strict_matching_pass=true;r3_json(matchfile,match);
  end
  freeze_predictions(d,sigmas,folder,c,'initial_prespecified');count=count+numel(sigmas);
  ModelUtil.clear;clear m d;
 end
 r3_json(ready,struct('status','All initial new-periphery predictions saved before target solves','count',count,...
  'saved_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')),'design_file',designFile));
elseif strcmp(action,'truth')
 assert(isfile(fullfile(out,'periphery_prediction_batch_ready.json')),'Freeze all initial predictions first');
 for j=1:numel(cases),truth_points(out,cases{j},sigmas,'initial_prespecified');end
 pair_summary(out,cases);
elseif strcmp(action,'locate')
 for iteration=1:12
  [~,brackets]=pair_summary(out,cases);pending=struct([]);
  for j=1:numel(brackets)
   b=brackets(j);
   if (b.upper-b.lower)/mean([b.upper b.lower])>.01
    e=struct('layout',b.layout,'sf',mean([b.lower b.upper]));
    if isempty(pending),pending=e;else,pending(end+1)=e;end %#ok<AGROW>
   end
  end
  if isempty(pending),break;end
  for j=1:numel(pending)
   todo=cases(cellfun(@(c)strcmp(c.layout,pending(j).layout),cases));sf=pending(j).sf;
   % Both A/B predictions are recorded before either target at this new point.
   for k=1:numel(todo)
    c=todo{k};folder=fullfile(out,'periphery',c.id);d=load(fullfile(folder,'background_responses.mat'));
    freeze_predictions(d,sf,folder,c,'adaptive_location_after_initial_truth');clear d;
   end
   for k=1:numel(todo),truth_points(out,todo{k},sf,'adaptive_location_after_initial_truth');end
  end
 end
 [~,brackets]=pair_summary(out,cases);r3_json(fullfile(out,'periphery_brackets.json'),brackets);
elseif strcmp(action,'refine')
 [~,brackets]=pair_summary(out,cases);
 for j=1:numel(brackets)
  todo=cases(cellfun(@(c)strcmp(c.layout,brackets(j).layout),cases));
  for k=1:numel(todo)
   c=todo{k};
   for sf=[brackets(j).lower brackets(j).upper]
    source=fullfile(out,'periphery',c.id,sf_token(sf));dest=fullfile(out,'precision',[c.id '_' sf_token(sf)]);
    r3_refine_case(source,dest);ModelUtil.clear;
   end
  end
 end
else,error('Unknown periphery action %s',action);end
end

function freeze_predictions(d,sigmas,folder,c,phase)
pending=[];
for sf=reshape(sigmas,1,[])
 target=fullfile(folder,sf_token(sf));if ~isfolder(target),mkdir(target);end
 if ~isfile(fullfile(target,'prediction.json'))
  assert(~isfile(fullfile(target,'SUCCESS.json')),'Cannot first predict a seen target');pending(end+1)=sf; %#ok<AGROW>
 end
end
if isempty(pending),return;end
[pred,fields]=r3_enrich_predict_v1(d,pending);
for k=1:numel(pred)
 sf=pred(k).sigma_f_S_m;target=fullfile(folder,sf_token(sf));
 pred(k).phase=phase;pred(k).layout=c.layout;pred(k).role=c.configuration_role;
 pred(k).prediction_saved_local=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
 pred(k).independent_geometry=true;pred(k).method_version='r3_background_vector4_residual_rankguard_v1';
 r3_json(fullfile(target,'prediction.json'),pred(k));
 key=['s' strrep(num2str(sf),'.','p')];prediction_fields=fields.(key); %#ok<NASGU>
 save(fullfile(target,'prediction_fields.mat'),'prediction_fields','-v7.3');
end
end

function truth_points(out,c,sigmas,phase)
folder=fullfile(out,'periphery',c.id);d=load(fullfile(folder,'background_responses.mat'));
m=mphload(fullfile(folder,'background.mph'),'R4');
for sf=reshape(sigmas,1,[])
 target=fullfile(folder,sf_token(sf));assert(isfile(fullfile(target,'prediction.json')),'Predict before truth');
 if isfile(fullfile(target,'SUCCESS.json')),continue;end
 r=r3_truth(m,d,sf,target,'');r.round4_phase=phase;r.round4_layout=c.layout;r.round4_role=c.configuration_role;
 r.round4_reused=false;r3_json(fullfile(target,'SUCCESS.json'),r);
 marginal=struct('sigma_f_S_m',sf,'Sf_m',mphint2(m,'ec.Ex^2+ec.Ey^2',2,'selection',d.sel.inclusion_geometry_domains,'intorder',8)*c.out_of_plane_thickness_m,...
  'definition','t/DV^2 integral ellipse |E|^2; DV=1 V; posthoc derivative, not background prediction');
 r3_json(fullfile(target,'marginal_scalar.json'),marginal);
end
com.comsol.model.util.ModelUtil.clear;
end

function [rows,brackets]=pair_summary(out,cases)
rows=struct([]);brackets=struct([]);layouts=unique(cellfun(@(c)c.layout,cases,'UniformOutput',false));
for j=1:numel(layouts)
 layout=layouts{j};a=cases{find(cellfun(@(c)strcmp(c.layout,layout)&&strcmp(c.configuration_role,'A'),cases),1)};
 b=cases{find(cellfun(@(c)strcmp(c.layout,layout)&&strcmp(c.configuration_role,'B'),cases),1)};
 list=dir(fullfile(out,'periphery',a.id,'sf_*'));part=struct([]);
 for k=1:numel(list)
  af=fullfile(list(k).folder,list(k).name,'SUCCESS.json');bf=fullfile(out,'periphery',b.id,list(k).name,'SUCCESS.json');
  if ~isfile(af)||~isfile(bf),continue;end
  ra=jsondecode(fileread(af));rb=jsondecode(fileread(bf));
  r=struct('layout',layout,'sigma_f_S_m',ra.sigma_f_S_m,'G0_A_S',ra.G0_S,'G0_B_S',rb.G0_S,...
   'deltaG_A_S',ra.deltaG_S,'deltaG_B_S',rb.deltaG_S,'D_B_minus_A_S',rb.deltaG_S-ra.deltaG_S,'source_A',af,'source_B',bf);
  if isempty(part),part=r;else,part(end+1)=r;end %#ok<AGROW>
 end
 if isempty(part),continue;end
 [~,ix]=sort([part.sigma_f_S_m]);part=part(ix);
 if isempty(rows),rows=part;else,rows=[rows part];end %#ok<AGROW>
 for k=1:numel(part)-1
  if part(k).D_B_minus_A_S*part(k+1).D_B_minus_A_S<0
   v=struct('layout',layout,'lower',part(k).sigma_f_S_m,'upper',part(k+1).sigma_f_S_m,...
    'D_lower_S',part(k).D_B_minus_A_S,'D_upper_S',part(k+1).D_B_minus_A_S);
   if isempty(brackets),brackets=v;else,brackets(end+1)=v;end %#ok<AGROW>
  end
 end
end
if ~isempty(rows),writetable(struct2table(rows),fullfile(out,'periphery_paired_response.csv'));end
fprintf('PERIPHERY_COMSOL_BRACKETS\n');disp(brackets);
end

function s=sf_token(sf)
s=['sf_' strrep(sprintf('%.10g',sf),'.','p')];
end
