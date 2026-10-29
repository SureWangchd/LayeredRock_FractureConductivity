function r5_crossover_predict(root)
% Limited crossover prediction from frozen R5 matrices and analytic layers.
% No full-target files or target fields are read; no COMSOL solve is called.
% Caller owns any existing LiveLink connection. This function starts none.
root=char(root);out=fullfile(root,'runs','round5_modal_closure');
assert(isfile(fullfile(out,'predictions_ready.json')),'Seal all 48 initial predictions first');
frozen=jsondecode(fileread(fullfile(out,'method_frozen.json')));
Ks=unique([double(frozen.K_star) double(frozen.K_reference)]);
assert(isequal(Ks,[8 32]),'This limited design expects frozen K8 and K32');
names={'P11','P77','P17'};layouts=load_layouts(fullfile(out,'port_coefficients.csv'),names);
refs=cell(1,2);backgrounds=cell(1,2);roles={'A','B'};
for ir=1:2
 refs{ir}=load(fullfile(out,'reference',[roles{ir} '_base'],'reference.mat'));
 allbg=jsondecode(fileread(fullfile(out,'predicted_backgrounds',[roles{ir} '_base'],'new.json')));
 backgrounds{ir}=zeros(numel(names),numel(Ks));
 for il=1:numel(names)
  for ik=1:numel(Ks)
   matches=find(arrayfun(@(q)strcmp(q.layout,names{il})&&q.K==Ks(ik),allbg));
   assert(numel(matches)==1,'Missing true background for crossover prediction');
   backgrounds{ir}(il,ik)=allbg(matches).G_S;
  end
 end
end
initialPoints=[.1 .14 .16 .18 .2 .3];points=initialPoints;
for sf=points,ensure_point(out,layouts,Ks,refs,backgrounds,sf);end
targetWidth=.003;maxAdaptivePoints=12;adaptivePoints=[];
for iteration=1:maxAdaptivePoints
 [brackets,~]=collect_curve(out,names,Ks,points);
 widths=arrayfun(@(q)q.interval_upper_S_m-q.interval_lower_S_m,brackets);
 eligible=find(arrayfun(@(q)q.found,brackets)&widths>targetWidth);
 if isempty(eligible),break;end
 [~,where]=max(widths(eligible));chosen=brackets(eligible(where));
 sf=(chosen.interval_lower_S_m+chosen.interval_upper_S_m)/2;
 assert(all(abs(points-sf)>1e-12),'Adaptive midpoint already exists');
 ensure_point(out,layouts,Ks,refs,backgrounds,sf);
 points=sort([points sf]);adaptivePoints(end+1)=sf; %#ok<AGROW>
end
[brackets,~]=collect_curve(out,names,Ks,points);
% Endpoints surround both frozen-method intervals, without targeting an
% exact root. All endpoint predictions are saved before the target manifest.
entries=struct([]);targets=struct([]);margin=.003;minGap=3e-7;
for il=1:numel(names)
 name=names{il};selected=brackets(strcmp({brackets.layout},name));
 info=struct('layout',name,'prediction_brackets',selected,'status','no_sign_change_in_prespecified_range', ...
  'target_lower_sigma_f_S_m',NaN,'target_upper_sigma_f_S_m',NaN, ...
  'endpoint_margin_S_m',0,'minimum_requested_abs_D_S',minGap, ...
  'endpoint_point_signs_clear',false,'endpoint_predictions',struct([]));
 if all([selected.found])
  lo=max(.1,min([selected.interval_lower_S_m])-margin);
  hi=min(.3,max([selected.interval_upper_S_m])+margin);
  ensure_point(out,layouts,Ks,refs,backgrounds,lo);
  ensure_point(out,layouts,Ks,refs,backgrounds,hi);
  points=sort(unique([points lo hi]));
  endpoint=endpoint_values(out,name,Ks,lo,hi);clearSigns=all([endpoint.clear_signs]);
  if ~clearSigns
   % One deterministic widening, still inside the originally prescribed range.
   lo=max(.1,min([selected.interval_lower_S_m])-.005);
   hi=min(.3,max([selected.interval_upper_S_m])+.005);
   ensure_point(out,layouts,Ks,refs,backgrounds,lo);
   ensure_point(out,layouts,Ks,refs,backgrounds,hi);
   points=sort(unique([points lo hi]));
   endpoint=endpoint_values(out,name,Ks,lo,hi);clearSigns=all([endpoint.clear_signs]);
   info.endpoint_margin_S_m=.005;
  else,info.endpoint_margin_S_m=margin;end
  info.status='predicted_crossing_bracket_selected_for_validation';
  info.target_lower_sigma_f_S_m=lo;info.target_upper_sigma_f_S_m=hi;
  info.endpoint_point_signs_clear=clearSigns;info.endpoint_predictions=endpoint;
  for sf=[lo hi]
   for ir=1:2
    role=roles{ir};id=['r5_strict_' role '_' name];
    entry=struct('case_id',id,'role',role,'layout',name,'sigma_f_S_m',sf, ...
     'prediction_file',fullfile(out,'crossover_predictions',id,sf_token(sf),'prediction.json'), ...
     'selection_basis','Frozen K8/K32 prediction intervals plus fixed exterior margin; no target field or result used', ...
     'phase','limited_crossover_prediction_selected');
    assert(isfile(entry.prediction_file),'Endpoint prediction missing');
    if isempty(targets),targets=entry;else,targets(end+1)=entry;end %#ok<AGROW>
   end
  end
 end
 if isempty(entries),entries=info;else,entries(end+1)=info;end %#ok<AGROW>
end
assert(numel(targets)<=12,'At most 4 targets per three selected layouts');
[finalBrackets,curve]=collect_curve(out,names,Ks,points);
stamp=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
summary=struct('status','all_limited_crossover_predictions_saved_before_targets', ...
 'saved_local',stamp,'method_version',frozen.version,'K_values',Ks, ...
 'sigma_search_interval_S_m',[.1 .3],'initial_prediction_points_S_m',initialPoints, ...
 'adaptive_prediction_points_S_m',adaptivePoints,'all_prediction_points_S_m',points, ...
 'adaptive_point_limit',maxAdaptivePoints,'desired_prediction_bracket_width_S_m',targetWidth, ...
 'target_count',numel(targets),'layouts',entries,'final_sampled_prediction_brackets',finalBrackets, ...
 'target_truth_used',false,'initial_48_predictions_reused',true, ...
 'interpretation','Prediction-side brackets only; full COMSOL validation has not been read. They are not certified roots or proof of uniqueness.');
r3_json(fullfile(out,'crossover_predictions_summary.json'),summary);
writetable(struct2table(curve),fullfile(out,'crossover_prediction_curve.csv'));
r3_json(fullfile(out,'crossover_targets.json'),struct('status','prediction_only_targets_ready', ...
 'saved_local',stamp,'targets',targets,'target_count',numel(targets), ...
 'prediction_summary',fullfile(out,'crossover_predictions_summary.json'),'no_target_truth_read',true));
fprintf('R5_CROSSOVER_PREDICTIONS_READY points=%d targets=%d\n',numel(points),numel(targets));
end

function ensure_point(out,layouts,Ks,refs,backgrounds,sf)
roles={'A','B'};
for ir=1:2
 role=roles{ir};folders=cell(1,numel(layouts));pending=false(size(folders));
 for il=1:numel(layouts)
  name=layouts(il).layout;id=['r5_strict_' role '_' name];
  folders{il}=fullfile(out,'crossover_predictions',id,sf_token(sf));
  saved=fullfile(folders{il},'prediction.json');
  if isfile(saved),continue;end
  old=fullfile(out,'predictions',id,sf_token(sf),'prediction.json');
  if isfile(old)
   p=jsondecode(fileread(old));assert(~p.target_truth_used,'Only independently saved predictions can be reused');
   assert(abs(p.sigma_f_S_m-sf)<1e-12*max(1,sf),'Rounded prediction path collision');
   keep=ismember([p.methods.K],Ks);assert(sum(keep)==numel(Ks),'Frozen K missing in initial prediction');
   p.methods=p.methods(keep);p.phase='limited_crossover_prediction_reuse';p.source_prediction_file=old;
   p.prediction_copied_local=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
   if ~isfolder(folders{il}),mkdir(folders{il});end
   r3_json(saved,p);
  else,pending(il)=true;end
 end
 if ~any(pending),continue;end
 [results,fields]=r5_engine(refs{ir},sf,layouts,Ks);
 for il=find(pending)
  name=layouts(il).layout;pick=(il-1)*numel(Ks)+(1:numel(Ks));methods=results(pick);
  for ik=1:numel(Ks)
   methods(ik).G0_S=backgrounds{ir}(il,ik);methods(ik).deltaG_S=methods(ik).G_S-methods(ik).G0_S;
   methods(ik).method=sprintf('M%d',methods(ik).K);
  end
  prediction=struct('case_id',['r5_strict_' role '_' name],'role',role,'layout',name, ...
   'sigma_f_S_m',sf,'phase','limited_crossover_prediction_only','mesh','base','methods',methods, ...
   'prediction_saved_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
   'method_version','r5_selfconsistent_modal_v1','target_truth_used',false, ...
   'source','Frozen COMSOL-FE M0 matrices and analytical peripheral layers only');
  if ~isfolder(folders{il}),mkdir(folders{il});end
  predictedFields=fields(pick); %#ok<NASGU>
  save(fullfile(folders{il},'predicted_fields.mat'),'predictedFields','Ks','-v7.3');
  r3_json(fullfile(folders{il},'prediction.json'),prediction);
  fprintf('R5_CROSSOVER_PREDICTED %s %s sf=%.10g gK8=%.12g gK32=%.12g\n',role,name,sf,methods(1).deltaG_S,methods(end).deltaG_S);
 end
 clear fields results;
end
end

function [brackets,curve]=collect_curve(out,names,Ks,points)
brackets=struct([]);curve=struct([]);points=sort(unique(points));
for il=1:numel(names)
 name=names{il};
 for ik=1:numel(Ks)
  K=Ks(ik);values=zeros(size(points));
  for j=1:numel(points)
   [values(j),gA,gB]=pair_value(out,name,K,points(j));
   row=struct('layout',name,'K',K,'sigma_f_S_m',points(j),'gA_prediction_S',gA,'gB_prediction_S',gB, ...
    'D_prediction_S',values(j),'target_truth_used',false);
   if isempty(curve),curve=row;else,curve(end+1)=row;end %#ok<AGROW>
  end
  changes=find(values(1:end-1).*values(2:end)<=0);
  record=struct('layout',name,'K',K,'found',~isempty(changes),'interval_lower_S_m',NaN, ...
   'interval_upper_S_m',NaN,'D_lower_S',NaN,'D_upper_S',NaN,'number_sampled_sign_change_intervals',numel(changes), ...
   'selection_note','First sign change in increasing sigma; all sampled D values retained, no uniqueness assertion');
  if ~isempty(changes)
   index=changes(1);record.interval_lower_S_m=points(index);record.interval_upper_S_m=points(index+1);
   record.D_lower_S=values(index);record.D_upper_S=values(index+1);
  end
  if isempty(brackets),brackets=record;else,brackets(end+1)=record;end %#ok<AGROW>
 end
end
end

function records=endpoint_values(out,name,Ks,lo,hi)
records=struct([]);
for K=Ks
 low=pair_value(out,name,K,lo);high=pair_value(out,name,K,hi);
 item=struct('K',K,'lower_sigma_f_S_m',lo,'upper_sigma_f_S_m',hi,'D_lower_S',low,'D_upper_S',high, ...
  'opposite_point_signs',low*high<0,'clear_signs',low*high<0&&min(abs([low high]))>=3e-7);
 if isempty(records),records=item;else,records(end+1)=item;end %#ok<AGROW>
end
end

function [D,gA,gB]=pair_value(out,name,K,sf)
values=zeros(1,2);roles={'A','B'};
for ir=1:2
 file=fullfile(out,'crossover_predictions',['r5_strict_' roles{ir} '_' name],sf_token(sf),'prediction.json');
 p=jsondecode(fileread(file));assert(~p.target_truth_used,'Target field used in prediction');
 found=find([p.methods.K]==K);assert(numel(found)==1,'Missing requested frozen order');
 values(ir)=p.methods(found).deltaG_S;
end
gA=values(1);gB=values(2);D=gB-gA;
end

function layouts=load_layouts(path,names)
t=readtable(path,'TextType','string');layouts=struct([]);
for i=1:numel(names)
 actual=zeros(2,65);reference=actual;
 for s=1:2
  side={'lower','upper'};rows=t(t.layout==string(names{i})&t.side==string(side{s}),:);rows=sortrows(rows,'n');
  assert(height(rows)==65&&isequal(rows.n(:),(0:64)'),'Missing analytic coefficients');
  actual(s,:)=rows.Y_actual_S_m2;reference(s,:)=rows.Y_reference_S_m2;
 end
 item=struct('layout',names{i},'Y_actual',actual,'Y_reference',reference);
 if isempty(layouts),layouts=item;else,layouts(end+1)=item;end %#ok<AGROW>
end
end
function value=sf_token(sf),value=['sf_' strrep(sprintf('%.10g',sf),'.','p')];end
