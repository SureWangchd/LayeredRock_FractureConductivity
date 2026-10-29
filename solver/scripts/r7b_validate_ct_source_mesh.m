function [result,sel] = r7b_validate_ct_source_mesh(model,meshTag,src,outputFolder,mapping)
%R7B_VALIDATE_CT_SOURCE_MESH Verify original CT after real mesh Union/FreeTet.
% No fracture classification is performed here. In particular, this helper
% cannot manufacture a smooth fracture using centroids. It independently
% checks that each ACTUAL fitted tet retains one original CT phase and one
% protected-region class, using the existing complete tet/voxel SAT routine
% whenever its AABB reaches an opposite-phase voxel. Same-phase voxel
% partitions may disappear without changing the categorical CT geometry.
%
% src is returned by r7b_import_ct_source. meshTag must be the actual final
% volume mesh. result.domainCode is indexed by real COMSOL domain number;
% mapping.domainCode (indexed by actual domain ID) or .domainIDsByCode{1..6}
% provides independently inherited original phase/region classes.
% mapping.ellipseIDs is the ACTUAL mesh-Boolean inside-Ellipsoid selection.
% Alternative parallel mapping arrays: all_domains, original_phase, region,
% is_ellipse. No ellipse centre test is used anywhere in this function.
% Optional mapping absent: original phase verification only, sel is empty.
if nargin<4,outputFolder='';end
if nargin<5,mapping=[];end
if ~isempty(outputFolder) && ~isfolder(outputFolder),mkdir(outputFolder);end
timer=tic;me=model.component('comp1').mesh(meshTag);
vertices=double(me.getVertex);nt=double(me.getNumElem('tet'));
L=src.L_m;h=src.spacing_m;labels=src.labels_xyz;
assert(nt>0,'r7b:original:empty','No final tetrahedra');
assert(size(vertices,1)==3,'r7b:original:shape','Expected 3D physical vertices');
assert(all(abs(min(vertices,[],2))<L*1e-9) && all(abs(max(vertices,[],2)-L)<L*1e-9),'r7b:original:bbox','Final mesh outer sample bounds changed');
chunk=100000;domainCode=zeros(0,1);domainVolume=zeros(0,1);domainCount=zeros(0,1);
stats=[];planeCrossings=0;hasRegions=src.p.protection_regions_present;
hasMapping=~isempty(mapping);sel=[];
inheritedCode=[];inheritedEllipse=[];
if hasMapping,[inheritedCode,inheritedEllipse]=read_mapping(mapping);end
ellipseTetBlocks={};ellipsePhaseBlocks={};ellipseRegionBlocks={};
volumeByCodeEllipse=zeros(6,2);
for first=0:chunk:nt-1
 count=min(chunk,nt-first);tet=double(me.getElem('tet',int32(first),int32(count)))+1;
 ent=reshape(double(me.getElemEntity('tet',int32(first),int32(count))),[],1);
 q=reshape(vertices(:,tet),3,4,count);
 [phase,checked,rejected]=r7_ct_validate_native_phase(q,h,labels);
 if any(rejected)
  failure=struct('status','FAILED_ORIGINAL_CT_PHASE_PRESERVATION','mesh_tag',char(meshTag),...
   'chunk_first_tetrahedron_0based',first,'tetrahedra_in_chunk',count,'phase_diagnostics',checked,...
   'rejected_indices_0based',first+find(rejected)-1,'no_labels_modified',true);
  if ~isempty(outputFolder)
   r3_json(fullfile(outputFolder,'ORIGINAL_CT_PHASE_FAILED.json'),failure);
   actual_bad_tetrahedra_m=q(:,:,rejected); %#ok<NASGU>
   save(fullfile(outputFolder,'actual_phase_crossing_tets.mat'),'actual_bad_tetrahedra_m','-v7');
  end
  error('r7b:original:phase','%d actual final tets cross original pore/skeleton interfaces in chunk %d',nnz(rejected),first);
 end
 if isempty(stats),stats=checked;
 else
  fields=fieldnames(checked);
  for j=1:numel(fields)
   name=fields{j};if strcmp(name,'coordinate_tolerance_native_voxels'),continue;end
   stats.(name)=stats.(name)+checked.(name);
  end
 end
 center=reshape(mean(q,2),3,count);
 if hasRegions
  z=center(3,:);region=double(z>.1*L & z<.25*L)+2*double(z>.75*L & z<.9*L);
  lo=reshape(min(q(3,:,:),[],2),1,count);hi=reshape(max(q(3,:,:),[],2),1,count);
  for plane=[.1 .25 .75 .9]*L
   crossing=lo<plane-L*1e-10 & hi>plane+L*1e-10;
   planeCrossings=planeCrossings+nnz(crossing);
   if any(crossing)
    if ~isempty(outputFolder)
     r3_json(fullfile(outputFolder,'PROTECTED_PLANE_FAILED.json'),struct('mesh_tag',char(meshTag),...
      'plane_z_m',plane,'chunk_first_tetrahedron_0based',first,'crossing_indices_0based',first+find(crossing)-1));
    end
    error('r7b:original:plane','Actual final mesh crosses protected plane z=%.17g m',plane);
   end
  end
 else,region=zeros(1,count);end
 code=1+phase+2*region;
 if hasMapping
  assert(max(ent)<=numel(inheritedCode),'r7b:original:mapping','Inherited selection misses actual final domain IDs');
  assert(all(reshape(inheritedCode(ent),1,[])==code),'r7b:original:mapping','Inherited source selection disagrees with actual original CT phase/region');
  isEllipse=reshape(inheritedEllipse(ent),1,[]);
  assert(all(region(isEllipse)==0),'r7b:original:ellipse','Inherited ellipse selection enters simplifiable remote region');
  if any(isEllipse)
   ellipseTetBlocks{end+1}=int32(tet(:,isEllipse)); %#ok<AGROW>
   ellipsePhaseBlocks{end+1}=uint8(phase(isEllipse)); %#ok<AGROW>
   ellipseRegionBlocks{end+1}=uint8(region(isEllipse)); %#ok<AGROW>
  end
 end
 mn=accumarray(ent,code(:),[],@min,NaN);mx=accumarray(ent,code(:),[],@max,NaN);valid=find(~isnan(mn));
 assert(all(mn(valid)==mx(valid)),'r7b:original:domain','Actual final domain merges different original phase/region classes');
 if numel(domainCode)<numel(mn)
  domainCode(end+1:numel(mn),1)=0;domainVolume(end+1:numel(mn),1)=0;domainCount(end+1:numel(mn),1)=0;
 end
 assert(all(domainCode(valid)==0 | domainCode(valid)==mn(valid)),'r7b:original:domain','Actual final domain code changes between chunks');
 domainCode(valid)=mn(valid);
 a=reshape(q(:,2,:)-q(:,1,:),3,count);b=reshape(q(:,3,:)-q(:,1,:),3,count);c=reshape(q(:,4,:)-q(:,1,:),3,count);
 volume=abs(sum(a.*cross(b,c,1),1))/6;
 if any(volume<=0)
  bad=volume<=0;actual_bad_tetrahedra_m=q(:,:,bad);
  if ~isempty(outputFolder)
   r3_json(fullfile(outputFolder,'DEGENERATE_TET_FAILED.json'),struct('status','FAILED_ZERO_VOLUME_TET','chunk_first_tetrahedron_0based',first,'indices_0based',first+find(bad)-1,'actual_domain_ids',ent(bad),'count',nnz(bad),'minimum_tet_volume_m3',min(volume),'no_degenerate_tets_removed',true));
   save(fullfile(outputFolder,'actual_degenerate_tets.mat'),'actual_bad_tetrahedra_m','-v7');
  end
  error('r7b:original:degenerate','%d zero-volume actual tets in chunk %d; no failed cells removed',nnz(bad),first);
 end
 domainVolume=domainVolume+accumarray(ent,volume(:),[numel(domainCode) 1],@sum,0);
 domainCount=domainCount+accumarray(ent,ones(count,1),[numel(domainCode) 1],@sum,0);
 if hasMapping
  volumeByCodeEllipse=volumeByCodeEllipse+accumarray([code(:),double(isEllipse(:))+1],volume(:),[6 2],@sum,0);
 end
end
ids=find(domainCode>0);codes=domainCode(ids);phase=mod(codes-1,2);region=floor((codes-1)/2);
actual=accumarray(codes,domainVolume(ids),[6 1],@sum,0);
expected=reshape(src.verification.expected_phase_region_volume_m3,6,1);
err=max(abs(actual-expected))/L^3;
% Original categorical CT is piecewise planar. This threshold verifies
% preservation, not permission to smooth/change its geometry.
assert(err<1e-8,'r7b:original:volume','Original CT phase/region volume changed by %.3g L^3',err);
byCode=cell(1,6);for code=1:6,byCode{code}=reshape(ids(codes==code),1,[]);end
result=struct('schema','r7b_original_CT_after_fitted_remeshing_v1','status','ORIGINAL_PHASE_AND_REGION_VERIFIED',...
 'mesh_tag',char(meshTag),'vertices',size(vertices,2),'tetrahedra',nt,'domainCode',domainCode,...
 'all_domains',ids','domainIDsByCode',{byCode},'domain_volume_m3',domainVolume,...
 'domain_tetrahedron_count',domainCount,'actual_phase_region_volume_m3',actual',...
 'expected_phase_region_volume_m3',expected','maximum_volume_error_over_L3',err,...
 'original_phase_all_tetrahedra',stats,'protected_plane_crossings',planeCrossings,...
 'protection_regions_present',hasRegions,'fracture_geometry_not_classified_by_this_helper',true,...
 'all_tetrahedra_checked_against_original_labels',true,'seconds',toc(timer));
if hasMapping
 isEllipse=logical(inheritedEllipse(ids));
 assert(any(isEllipse),'r7b:original:ellipse','No actual ellipse domain after Boolean/intersection');
 ms=mphmeshstats(model,char(meshTag));
 sel=struct('all_domains',ids','domain_code',codes','original_phase',phase','region',region',...
  'is_ellipse',reshape(isEllipse,1,[]),'inclusion_domains',reshape(ids(isEllipse),1,[]),...
  'pore_domains',reshape(ids(phase==0),1,[]),'skeleton_domains',reshape(ids(phase==1),1,[]),...
  'remote_lower_domains',reshape(ids(region==1),1,[]),'remote_upper_domains',reshape(ids(region==2),1,[]),...
  'top',[],'bottom',[],'electrode_selections_pending_caller',true,...
  'tetrahedra',nt,'vertices',size(vertices,2),'min_mesh_quality',ms.minquality,...
  'mean_mesh_quality',ms.meanquality,'verified_phase_region_mask_all_tetrahedra',true,...
  'ellipse_selection_source','Inherited actual Boolean/mesh domain selection; no centroid classification',...
  'actual_phase_region_ellipse_volume_m3',volumeByCodeEllipse,'mesh_seconds',toc(timer));
 globalTet=cat(2,ellipseTetBlocks{:});
 [used,~,localIndex]=unique(globalTet(:));
 vertices_m=vertices(:,double(used)).'; %#ok<NASGU>
 tetrahedra_inside=int32(reshape(localIndex,size(globalTet)).'-1); %#ok<NASGU>
 original_phase_inside=cat(2,ellipsePhaseBlocks{:}).'; %#ok<NASGU>
 region_inside=cat(2,ellipseRegionBlocks{:}).'; %#ok<NASGU>
 index_base=0; %#ok<NASGU>
 result.inherited_phase_region_selection_verified=true;
 result.actual_phase_region_ellipse_volume_m3=volumeByCodeEllipse;
 result.ellipse_tetrahedra_exported=size(globalTet,2);
 result.ellipse_local_vertices_exported=numel(used);
 result.ellipse_label_from_actual_boolean_selection=true;
 result.ellipse_geometry_measurement_still_required=true;
 if ~isempty(outputFolder)
  exportFile=fullfile(outputFolder,'actual_boundary.mat');
  save(exportFile,'vertices_m','tetrahedra_inside','original_phase_inside','region_inside','index_base','-v7');
  result.actual_boundary_export_file=exportFile;
  r3_json(fullfile(outputFolder,'verified_selections_pending_electrodes.json'),sel);
  rows=[];
  for code=1:6
   for ell=0:1
    rows(end+1,:)=[code mod(code-1,2) floor((code-1)/2) ell volumeByCodeEllipse(code,ell+1)]; %#ok<AGROW>
   end
  end
  writetable(array2table(rows,'VariableNames',{'source_code','native_phase','region','inside_ellipsoid','actual_volume_m3'}),...
   fullfile(outputFolder,'actual_phase_region_ellipse_volumes.csv'));
 end
end
if ~isempty(outputFolder)
 r3_json(fullfile(outputFolder,'ORIGINAL_CT_PHASE_PRESERVATION.json'),result);
 writetable(table(ids,codes,phase,region,domainCount(ids),domainVolume(ids),...
  'VariableNames',{'domain_id','original_source_code','native_phase','region','tetrahedra','volume_m3'}),...
  fullfile(outputFolder,'actual_domains_original_phase.csv'));
end
fprintf('R7B_ORIGINAL_CT_VERIFIED mesh %s tets %d domains %d volumeError/L3 %.3g seconds %.2f\n',char(meshTag),nt,numel(ids),err,toc(timer));
end

function [code,ellipse]=read_mapping(mapping)
% Accept either complete domain-indexed inherited selections or parallel
% domain arrays. Domain IDs are 1-based COMSOL entities, not tet indices.
if isfield(mapping,'domainCode')
 code=reshape(double(mapping.domainCode),[],1);
elseif isfield(mapping,'domainIDsByCode')
 code=zeros(0,1);
 for c=1:6
  ids=reshape(double(mapping.domainIDsByCode{c}),[],1);
  if isempty(ids),continue;end
  if numel(code)<max(ids),code(end+1:max(ids),1)=0;end
  assert(all(code(ids)==0),'r7b:original:mapping','Inherited original code selections overlap');code(ids)=c;
 end
elseif all(isfield(mapping,{'all_domains','original_phase','region'}))
 ids=reshape(double(mapping.all_domains),[],1);code=zeros(max(ids),1);
 code(ids)=1+reshape(double(mapping.original_phase),[],1)+2*reshape(double(mapping.region),[],1);
else,error('r7b:original:mapping','Need domainCode, domainIDsByCode, or all_domains/original_phase/region');end
assert(all(code==0 | ismember(code,1:6)),'r7b:original:mapping','Invalid source code; expected 1..6');
ellipse=false(size(code));
if isfield(mapping,'ellipseIDs'),eids=reshape(double(mapping.ellipseIDs),[],1);
elseif isfield(mapping,'inclusion_domains'),eids=reshape(double(mapping.inclusion_domains),[],1);
elseif isfield(mapping,'all_domains') && isfield(mapping,'is_ellipse')
 ids=reshape(double(mapping.all_domains),[],1);eids=ids(logical(mapping.is_ellipse(:)));
else,error('r7b:original:mapping','Actual inherited ellipseIDs selection is required');end
assert(~isempty(eids) && all(eids>=1 & eids==floor(eids) & eids<=numel(code)),'r7b:original:mapping','Invalid or empty ellipse IDs');
assert(all(code(eids)>0),'r7b:original:mapping','Ellipse domains are missing original phase labels');
ellipse(eids)=true;
end
