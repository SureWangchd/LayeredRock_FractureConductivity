function r6_predict(root,phase)
% R6 orchestration only: unchanged r5_engine and newly computed core sources.
out=fullfile(root,'runs','round6_applicability');design=jsondecode(fileread(fullfile(out,'design.json')));
mesh='base';baseFolder='predictions';
if strcmp(phase,'precision')
 mesh='fine';baseFolder='prediction_precision';req=jsondecode(fileread(fullfile(out,'refinement_targets.json')));selected=req.targets;
else
 selected=design.new_states;
 if isfile(fullfile(out,'predictions_ready.json')),fprintf('R6_PREDICTIONS_ALREADY_SEALED\n');return;end
end
names={'P11','P77','P17','P71'};Ks=[0 8 32 64];layouts=load_layouts(fullfile(out,'port_coefficients.csv'),names);
for ig=1:numel(design.source_groups)
 group=design.source_groups(ig);d=load(fullfile(out,'reference',[group.source_key '_' mesh],'reference.mat'));
 [backgrounds,bgFields]=r5_engine(d,NaN,layouts,Ks);
 bgFolder=fullfile(out,'predicted_backgrounds',[group.source_key '_' mesh]);if ~isfolder(bgFolder),mkdir(bgFolder);end
 r3_json(fullfile(bgFolder,'background.json'),backgrounds);save(fullfile(bgFolder,'background_fields.mat'),'bgFields','-v7.3');clear bgFields;
 for sf=[.18 10]
  take=false(1,4);expected=cell(1,4);
  for il=1:4
   id=['r6_' group.source_key '_' names{il}];expected{il}=fullfile(out,baseFolder,id,sf_token(sf));
   take(il)=any(arrayfun(@(s)strcmp(s.case_id,id)&&s.sigma_f_S_m==sf,selected));
  end
  if ~any(take)||all(cellfun(@(f)isfile(fullfile(f,'prediction.json')),expected(take))),continue;end
  [results,fields]=r5_engine(d,sf,layouts,Ks);
  diagnostics=r6_domination_check(d,sf,layouts(1));
  r3_json(fullfile(bgFolder,[sf_token(sf) '_conditions.json']),diagnostics);
  assert(diagnostics.numerical_domination_pass,'Finite-FE domination check failed');
  for il=1:4
   if ~take(il),continue;end
   folder=expected{il};if isfile(fullfile(folder,'prediction.json')),continue;end
   if ~isfolder(folder),mkdir(folder);end
   pick=(il-1)*4+(1:4);methods=results(pick);predictedFields=fields(pick);
   for j=1:4
    methods(j).G0_S=backgrounds(pick(j)).G_S;methods(j).deltaG_S=methods(j).G_S-methods(j).G0_S;
    methods(j).method=sprintf('M%d',methods(j).K);
   end
   prediction=struct('case_id',['r6_' group.source_key '_' names{il}],'role',group.role,'layout',names{il},...
    'length_mm',group.full_length_mm,'full_length_mm',group.full_length_mm,'sigma_f_S_m',sf,...
    'phase',phase,'mesh',mesh,'methods',methods,'prediction_saved_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')),...
    'method_version','r5_selfconsistent_modal_v1_unchanged','target_truth_used',false,...
    'field_description','Physical core, common M0 outer realization; nonlocal modal correction on core interfaces.',...
    'finite_model_conditions',diagnostics);
   save(fullfile(folder,'predicted_fields.mat'),'predictedFields','Ks','-v7.3');r3_json(fullfile(folder,'prediction.json'),prediction);
   fprintf('R6_PREDICTED %s %s sf %.8g M0 %.12g M64 %.12g\n',group.source_key,names{il},sf,methods(1).deltaG_S,methods(4).deltaG_S);
  end
  clear results fields;
 end
 clear d;
end
if strcmp(phase,'new')
 entries=design.new_states;
 for j=1:numel(entries)
  entries(j).prediction_file=fullfile(out,'predictions',entries(j).case_id,sf_token(entries(j).sigma_f_S_m),'prediction.json');
  assert(isfile(entries(j).prediction_file));
  assert(~isfile(fullfile(out,'truth',entries(j).case_id,sf_token(entries(j).sigma_f_S_m),'SUCCESS.json')),'Target before prediction seal');
 end
 r3_json(fullfile(out,'prediction_fields_ready.json'),struct('count',32,'entries',entries,'saved_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'))));
end
end
function layouts=load_layouts(path,names)
t=readtable(path,'TextType','string');layouts=struct([]);
for i=1:numel(names)
 actual=zeros(2,65);reference=actual;
 for s=1:2
  side={'lower','upper'};rows=t(t.layout==string(names{i}) & t.side==string(side{s}),:);rows=sortrows(rows,'n');
  assert(height(rows)==65&&isequal(rows.n(:),(0:64)'));actual(s,:)=rows.Y_actual_S_m2;reference(s,:)=rows.Y_reference_S_m2;
 end
 layouts(i).layout=names{i};layouts(i).Y_actual=actual;layouts(i).Y_reference=reference;
end
end
function s=sf_token(sf),s=['sf_' strrep(sprintf('%.10g',sf),'.','p')];end
