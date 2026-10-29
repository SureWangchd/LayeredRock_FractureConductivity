function r5_predict(root,phase)
% All prediction inputs come from the common reference and analytic layers.
out=fullfile(root,'runs','round5_modal_closure');
design=jsondecode(fileread(fullfile(out,'design.json')));
if strcmp(phase,'new') && isfile(fullfile(out,'predictions_ready.json'))
 ready=jsondecode(fileread(fullfile(out,'predictions_ready.json')));
 assert(ready.count==48,'Unexpected saved prediction batch');
 for j=1:numel(ready.entries),assert(isfile(ready.entries(j).prediction_file),'A sealed prediction file is missing');end
 fprintf('R5_INITIAL_PREDICTION_BATCH_ALREADY_SEALED\n');return;
end
if strcmp(phase,'development')
 names={'NN','NF','FN','FF'};sfs=[.1 .159375 .3];Ks=[0 1 2 4 8 16 32 64];mesh='base';
elseif strcmp(phase,'new')
 frozen=jsondecode(fileread(fullfile(out,'method_frozen.json')));
 names={'P11','P77','P17','P71','P41','P14'};sfs=[.1 .18 .3 10];Ks=unique([0 frozen.K_star frozen.K_reference]);mesh='base';
elseif strcmp(phase,'precision')
 frozen=jsondecode(fileread(fullfile(out,'method_frozen.json')));
 names={'NN','NF','FN','FF','P11','P77'};sfs=[.159375 10];Ks=unique([0 frozen.K_star 32 64]);mesh='fine';
else,error('Unknown prediction phase');end
layouts=load_layouts(fullfile(out,'port_coefficients.csv'),names);
for role={'A','B'}
 r=role{1};d=load(fullfile(out,'reference',[r '_' mesh],'reference.mat'));
 [backgrounds,bgFields]=r5_engine(d,NaN,layouts,Ks);
 bgFolder=fullfile(out,'predicted_backgrounds',[r '_' mesh]);if ~isfolder(bgFolder),mkdir(bgFolder);end
 r3_json(fullfile(bgFolder,[phase '.json']),backgrounds);
 save(fullfile(bgFolder,[phase '_fields.mat']),'bgFields','-v7.3');clear bgFields;
 for sf=sfs
  expected=cellfun(@(name) output_folder(out,phase,r,name,sf),names,'UniformOutput',false);
  selected=true(size(names));
  if strcmp(phase,'precision'),selected=startsWith(names,'P')==(sf==10);end
  if all(cellfun(@(p)isfile(fullfile(p,'prediction.json')),expected(selected))),fprintf('R5_PREDICTION_REUSE %s %s %.8g\n',phase,r,sf);continue;end
  [results,fields]=r5_engine(d,sf,layouts,Ks);
  for il=1:numel(names)
   if ~selected(il),continue;end
   pick=(il-1)*numel(Ks)+(1:numel(Ks));methods=results(pick);predictedFields=fields(pick);
   for j=1:numel(methods)
    methods(j).G0_S=backgrounds(pick(j)).G_S;methods(j).deltaG_S=methods(j).G_S-methods(j).G0_S;
    methods(j).method=sprintf('M%d',methods(j).K);
   end
   folder=expected{il};if ~isfolder(folder),mkdir(folder);end
   if isfile(fullfile(folder,'prediction.json')),continue;end
   prediction=struct('case_id',['r5_strict_' r '_' names{il}], 'role',r,'layout',names{il}, ...
    'sigma_f_S_m',sf,'phase',phase,'mesh',mesh,'methods',methods, ...
    'prediction_saved_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'method_version','r5_selfconsistent_modal_v1','target_truth_used',false, ...
    'field_description','Core field is physical; outer FE field is the common M0 realization. Analytic nonlocal correction lives on the two core interfaces.');
   save(fullfile(folder,'predicted_fields.mat'),'predictedFields','Ks','-v7.3');
   r3_json(fullfile(folder,'prediction.json'),prediction);
   fprintf('R5_PREDICTED %s %s %s sf=%.8g g(M0)=%.12g g(Kmax)=%.12g\n',phase,r,names{il},sf,methods(1).deltaG_S,methods(end).deltaG_S);
  end
  clear fields results;
 end
 clear d;
end
if strcmp(phase,'new')
 entries=design.new_states;
 for j=1:numel(entries)
  entries(j).prediction_file=fullfile(output_folder(out,phase,entries(j).role,entries(j).layout,entries(j).sigma_f_S_m),'prediction.json');
  assert(isfile(entries(j).prediction_file),'Incomplete initial prediction batch');
  assert(~isfile(fullfile(out,'truth',entries(j).case_id,sf_token(entries(j).sigma_f_S_m),'SUCCESS.json')),'Target existed before batch prediction was sealed');
 end
 r3_json(fullfile(out,'predictions_ready.json'),struct('status','all_48_predictions_saved_before_targets', ...
  'count',48,'saved_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')),'entries',entries, ...
  'method_file',fullfile(out,'method_frozen.json')));
end
end

function layouts=load_layouts(path,names)
t=readtable(path,'TextType','string');layouts=struct([]);
for i=1:numel(names)
 actual=zeros(2,65);reference=actual;
 for s=1:2
  side={'lower','upper'};rows=t(t.layout==string(names{i}) & t.side==string(side{s}),:);rows=sortrows(rows,'n');
  assert(height(rows)==65&&isequal(rows.n(:),(0:64)'),'Missing analytic coefficients');
  actual(s,:)=rows.Y_actual_S_m2;reference(s,:)=rows.Y_reference_S_m2;
 end
 item=struct('layout',names{i},'Y_actual',actual,'Y_reference',reference);
 if isempty(layouts),layouts=item;else,layouts(end+1)=item;end %#ok<AGROW>
end
end
function path=output_folder(out,phase,role,name,sf)
if strcmp(phase,'development'),base=fullfile(out,'development',[role '_' name]);
elseif strcmp(phase,'new'),base=fullfile(out,'predictions',['r5_strict_' role '_' name]);
else,base=fullfile(out,'prediction_precision',[role '_' name]);end
path=fullfile(base,sf_token(sf));
end
function s=sf_token(sf),s=['sf_' strrep(sprintf('%.10g',sf),'.','p')];end
