function report=enrich_round2_topology(caseIDs)
% Read existing Round 2 MPH files; export connectivity and all interface lengths.
% No solve, geometry rebuild, model save, or COMSOL process launch is performed.
% The caller must run this serially after other clients finish with port 2037.
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
report=repmat(struct('case_id','','status','','error','','seconds',0),numel(caseIDs),1);
for k=1:numel(caseIDs)
 id=char(caseIDs{k});report(k).case_id=id;folder=fullfile(out,id);
 mphFile=fullfile(folder,[id '.mph']);jsonFile=fullfile(folder,'topology_interfaces.json');
 csvFile=fullfile(folder,'topology_interface_boundaries.csv');
 if ~isfile(fullfile(folder,'SUCCESS.json'))||~isfile(mphFile)
  report(k).status='missing_success';continue;
 end
 mi=dir(mphFile);
 if isfile(jsonFile)&&isfile(csvFile)
  previous=jsondecode(fileread(jsonFile));
  if isfield(previous,'status')&&strcmp(previous.status,'success')&& ...
     previous.source_mph_datenum==mi.datenum&&previous.source_mph_bytes==mi.bytes
   report(k).status='already_complete';fprintf('TOPOLOGY_ROUND2_SKIP %s\n',id);continue;
  end
 end
 started=tic;tag=['Round2Topo_' strrep(char(java.util.UUID.randomUUID),'-','')];
 modelCleanup=onCleanup(@()remove_owned_model(tag));
 fprintf('TOPOLOGY_ROUND2_START %s\n',id);
 try
  p=jsondecode(fileread(fullfile(folder,'parameters.json')));
  sel=jsondecode(fileread(fullfile(folder,'selections.json')));
  model=mphload(mphFile,tag);
  [data,T]=topology_interfaces(model,p,sel,folder);
  data.case_id=id;data.status='success';data.source_mph=mphFile;
  data.source_mph_datenum=mi.datenum;data.source_mph_bytes=mi.bytes;
  data.postprocessing_only=true;data.postprocessing_seconds=toc(started);
  writetable(T,csvFile);
  write_json(jsonFile,data);
  report(k).status='success';report(k).seconds=data.postprocessing_seconds;
  fprintf('TOPOLOGY_ROUND2_DONE %s %.1fs main=%d total=%d\n',id, ...
   report(k).seconds,data.main_fracture_component_count,data.all_fracture_component_count);
 catch err
  report(k).status='failed';report(k).error=getReport(err,'extended','hyperlinks','off');
  report(k).seconds=toc(started);
  write_json(jsonFile,struct('case_id',id,'status','failed','error',report(k).error, ...
   'source_mph',mphFile,'postprocessing_only',true));
  warning('round2:TopologyFailed','%s: %s',id,err.message);
 end
 clear modelCleanup model;
end
end

function [data,T]=topology_interfaces(model,p,sel,folder)
allfrac=row_ids(sel.fracture_domains);comp=row_ids(sel.compensation_domains);
assert(all(ismember(comp,allfrac)),'Compensation domain is not assigned fracture material');
main=setdiff(allfrac,comp);low=row_ids(sel.low_domains);high=row_ids(sel.high_domains);
removed=row_ids(sel.removed_domains);
assert(isempty(intersect(allfrac,[low high])),'Fracture/matrix selections overlap');
assert(isempty(intersect(low,high)),'Low/high matrix selections overlap');
bnds=[];
for k=1:numel(allfrac)
 bnds=union(bnds,row_ids(mphgetadj(model,'geom1','boundary','domain',allfrac(k))));
end
bnds=row_ids(bnds);graph=false(numel(allfrac));
% Columns: boundary, fracture domain, matrix domain, matrix phase, owner,
% source. Owner: 0 main / 1 compensation. Source: 0 original ellipse arc /
% 1 new face adjoining the removed-and-restored matrix / 2 compensation.
rows=zeros(0,6);outside=[];internal=[];
for k=1:numel(bnds)
 boundary=bnds(k);
 adj=row_ids(mphgetadj(model,'geom1','domain','boundary',boundary));
 frac=intersect(adj,allfrac);matrix=intersect(adj,[low high]);
 if numel(frac)>1
  assert(numel(adj)==2&&numel(frac)==2,'Unexpected fracture/fracture boundary adjacency');
  [~,indices]=ismember(frac,allfrac);graph(indices(1),indices(2))=true;graph(indices(2),indices(1))=true;
  internal(end+1)=boundary; %#ok<AGROW>
 elseif numel(frac)==1&&numel(matrix)==1
  assert(numel(adj)==2,'Unexpected fracture/matrix boundary adjacency');
  owner=double(ismember(frac,comp));phase=double(ismember(matrix,high));
  if owner==1,source=2;
  elseif ismember(matrix,removed),source=1;
  else,source=0;
  end
  rows(end+1,:)=[boundary frac matrix phase owner source]; %#ok<AGROW>
 elseif numel(frac)==1&&numel(adj)==1
  outside(end+1)=boundary; %#ok<AGROW>
 else
  error('round2:TopologyAdjacency','Unclassified boundary %d adjacent=[%s]',boundary,num2str(adj));
 end
end
assert(isempty(outside),'Fracture material unexpectedly touches the specimen exterior');
T=array2table(rows,'VariableNames',{'boundary_id','fracture_domain_id','matrix_domain_id','matrix_phase','owner_code','source_code'});
data=struct('main_fracture_domains',main,'compensation_domains',comp,'all_fracture_domains',allfrac, ...
 'physical_fracture_present',~isempty(allfrac),'intervention',p.intervention, ...
 'connectivity_definition','Material domains are connected only across a shared one-dimensional geometric boundary; contact at a vertex does not create a conducting bridge.', ...
 'matrix_phase_definition','0=low; 1=high, from actual final matrix domain selections. Restored cut domains retain their original matrix phase.', ...
 'owner_code_definition','0=main finite ellipse remnant; 1=compensation inclusion.', ...
 'source_code_definition','0=surviving original ellipse/matrix interface; 1=new main-fracture face adjoining a removed/restored matrix domain; 2=compensation-circle/matrix interface.', ...
 'new_face_classifier','Exact domain adjacency against selections.removed_domains; not an ellipse-residual tolerance.', ...
 'units','Interface lengths in mm; domain areas remain in the original SUCCESS.json.', ...
 'fracture_internal_partition_boundaries',internal,'exterior_contact_boundaries',outside);
data.main_fracture_components=components(main,allfrac,graph);
data.compensation_components=components(comp,allfrac,graph);
data.all_fracture_components=components(allfrac,allfrac,graph);
data.main_fracture_component_count=numel(data.main_fracture_components);
data.compensation_component_count=numel(data.compensation_components);
data.all_fracture_component_count=numel(data.all_fracture_components);
data.components_include_only_sigma_frac_material=true;
data.background_high_conductivity_layers_are_not_fracture_components=true;
data.interface_boundary_count=height(T);
data.main_original_arc_boundaries=row_ids(rows(rows(:,6)==0,1));
data.main_new_face_boundaries=row_ids(rows(rows(:,6)==1,1));
data.compensation_matrix_boundaries=row_ids(rows(rows(:,6)==2,1));
data.main_original_arc_count=numel(data.main_original_arc_boundaries);
data.main_new_face_count=numel(data.main_new_face_boundaries);
data.compensation_matrix_boundary_count=numel(data.compensation_matrix_boundaries);
data.lengths_mm=struct;
lengthGroups=cell(1,6);position=0;
for source=0:2
 for phase=0:1
  position=position+1;ids=row_ids(rows(rows(:,6)==source&rows(:,4)==phase,1));
  if isempty(ids),len=0;else,len=mphint2(model,'1',1,'selection',ids,'intorder',8,'unit','mm');end
  assert(isfinite(len)&&len>=0,'Invalid interface length');
  ownerNames={'main_original','main_new','compensation'};phaseNames={'low','high'};
  field=[ownerNames{source+1} '_' phaseNames{phase+1}];data.lengths_mm.(field)=len;
  lengthGroups{position}=struct('source_code',source,'matrix_phase',phase,'boundary_ids',ids,'length_mm',len);
 end
end
data.interface_groups=lengthGroups;
v=data.lengths_mm;
data.lengths_mm.main_low=v.main_original_low+v.main_new_low;
data.lengths_mm.main_high=v.main_original_high+v.main_new_high;
data.lengths_mm.all_low=data.lengths_mm.main_low+v.compensation_low;
data.lengths_mm.all_high=data.lengths_mm.main_high+v.compensation_high;
data.lengths_mm.main_total=data.lengths_mm.main_low+data.lengths_mm.main_high;
data.lengths_mm.compensation_total=v.compensation_low+v.compensation_high;
data.lengths_mm.all_total=data.lengths_mm.all_low+data.lengths_mm.all_high;
data.compensation_present=~isempty(comp);
if ~isempty(comp)
 data.compensation_analytic_circle_circumference_mm=2*pi*p.comp_radius_mm;
 data.compensation_circumference_error_mm=data.lengths_mm.compensation_total-2*pi*p.comp_radius_mm;
else
 data.compensation_analytic_circle_circumference_mm=0;
 data.compensation_circumference_error_mm=0;
end
if isempty(allfrac),expectedMain=0;elseif strcmp(p.intervention,'cut'),expectedMain=2;else,expectedMain=1;end
data.designed_main_component_count=expectedMain;
data.main_component_count_matches_design=data.main_fracture_component_count==expectedMain;
data.compensation_is_disconnected_from_main=data.all_fracture_component_count==data.main_fracture_component_count+data.compensation_component_count;
data.interpretation='Connectivity and interface changes are measured geometric controls. They do not establish a purely topological electrical effect independent of shape or restored material.';
% Keep the older flux selection visible as a diagnostic, never overwrite it.
fluxCsv=fullfile(folder,'ellipse_boundary_exchange_segments.csv');
fluxJson=fullfile(folder,'ellipse_boundary_exchange.json');
data.flux_arc_comparison=struct('status','unavailable');
if isfile(fluxCsv)&&isfile(fluxJson)
 exchange=jsondecode(fileread(fluxJson));
 if strcmp(exchange.status,'success')
  fluxT=readtable(fluxCsv);fluxIDs=row_ids(fluxT.boundary_id);
  data.flux_arc_comparison=struct('status','available','exchange_arc_ids',fluxIDs, ...
   'topology_original_arcs_missing_in_exchange',setdiff(data.main_original_arc_boundaries,fluxIDs), ...
   'exchange_arcs_not_topology_original',setdiff(fluxIDs,data.main_original_arc_boundaries), ...
   'all_ids_match',isequal(data.main_original_arc_boundaries,fluxIDs));
 end
end
end

function groups=components(subset,allIDs,graph)
groups={};if isempty(subset),return;end
[~,indices]=ismember(subset,allIDs);g=graph(indices,indices);seen=false(1,numel(subset));
for k=1:numel(subset)
 if seen(k),continue;end
 queue=k;seen(k)=true;members=[];
 while ~isempty(queue)
  current=queue(1);queue(1)=[];members(end+1)=current; %#ok<AGROW>
  neighbors=find(g(current,:)&~seen);seen(neighbors)=true;queue=[queue neighbors]; %#ok<AGROW>
 end
 groups{end+1}=row_ids(subset(members)); %#ok<AGROW>
end
end

function ids=row_ids(values)
ids=unique(double(values(:)'));ids=ids(ids>0);
end

function write_json(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0,'Cannot open %s',path);
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end

function remove_owned_model(tag)
try,com.comsol.model.util.ModelUtil.remove(tag);catch,end
end

function disconnect_if_owned(owned)
if owned,try,com.comsol.model.util.ModelUtil.disconnect;catch,end,end
end
