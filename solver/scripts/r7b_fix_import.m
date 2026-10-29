function r7b_fix_import(root,sample,mode)
% Import explicitly cut volume elements; NO surface Union and NO FreeTet.
% Called on an existing LiveLink connection. Does not start/stop a service.
if nargin<3,mode='formal';end
import com.comsol.model.util.*
out=fullfile(root,'runs','round7b_localfit');
isRefined=any(strcmp(mode,{'M2','M3'}));
if strcmp(mode,'technical32'),folder=fullfile(out,'technical32');
elseif isRefined,folder=fullfile(out,'inputs',[sample '_' mode]);
else,folder=fullfile(out,'inputs',[sample '_M1']);end
assert(isfile(fullfile(folder,'PREPROCESS_COMPLETE.json')),'Run r7b_fix_prepare.py first');
if isfile(fullfile(folder,'GEOMETRY_SUCCESS.json')),fprintf('Accepted input exists; do not rebuild\n');return;end
if isfile(fullfile(folder,'geometry_ready.mph')),error('Existing imported MPH awaits validation. Do not overwrite; resume geometry checking.');end
p=jsondecode(fileread(fullfile(folder,'parameters.json')));
rec=jsondecode(fileread(fullfile(folder,'PREPROCESS_COMPLETE.json')));
check_binary_size(fullfile(folder,'vertices_m.bin'),rec.vertices*3*8);
check_binary_size(fullfile(folder,'tets0.bin'),rec.tetrahedra*4*4);
check_binary_size(fullfile(folder,'entity_code.bin'),rec.tetrahedra);
labels=load(fullfile(folder,'original_labels_and_parameters.mat'),'labels_xyz','expected_phase_region_volume_m3');
L=p.L_m;src=struct('p',p,'L_m',L,'spacing_m',p.spacing_m,'labels_xyz',labels.labels_xyz,...
 'verification',struct('expected_phase_region_volume_m3',reshape(labels.expected_phase_region_volume_m3,1,6)));
save(fullfile(folder,'original_CT_source_meta.mat'),'src','-v7');
if isRefined
 assert(p.refinement_preserves_frozen_surface,'Refined source lacks fixed-interface proof');
 anchorFolder=fullfile(out,'inputs',[sample '_M1']);
 anchorP=jsondecode(fileread(fullfile(anchorFolder,'parameters.json')));
 assert(strcmp(p.fixed_surface_identifier,anchorP.fixed_surface_identifier),'Frozen M1 interface identity changed');
 faceProof=jsondecode(fileread(fullfile(folder,'all_faces_verification.json')));
 assert(faceProof.all_face_topology_checks_passed&&faceProof.all_volume_faces_checked,'All-volume shared-face verification is required');
 copyfile(fullfile(anchorFolder,'native_ellipsoid_definition.mph'),fullfile(folder,'native_ellipsoid_definition.mph'));
else
% Retain the native CAD primitive as a separate definition record, not a
% second overlapping material body in the electrical component.
cadtag='R7BfixCAD';if any(strcmp(cell(ModelUtil.tags),cadtag)),ModelUtil.remove(cadtag);end
cm=ModelUtil.create(cadtag);cm.component.create('comp1',true);cm.component('comp1').geom.create('geom1',3);
g=cm.component('comp1').geom('geom1');g.lengthUnit('m');g.create('ell_ct','Ellipsoid');
e=g.feature('ell_ct');e.set('type','solid');e.set('semiaxes',reshape(p.ellipsoid_semiaxes_m,1,[]));e.set('pos',reshape(p.ellipsoid_center_m,1,[]));
e.set('axistype','cartesian');e.set('axis',[-sqrt(3)/2 0 .5]);e.set('rot',0);g.run();
mphsave(cm,fullfile(folder,'native_ellipsoid_definition.mph'));ModelUtil.remove(cadtag);
end

tag='R7BfixImport';if any(strcmp(cell(ModelUtil.tags),tag)),ModelUtil.remove(tag);end
cl=onCleanup(@()remove_model(tag)); %#ok<NASGU>
m=ModelUtil.create(tag);m.component.create('comp1',true);co=m.component('comp1');co.sorder('linear');
co.geom.create('geom1',3);co.geom('geom1').lengthUnit('m');co.mesh.create('cutsource','geom1');me=co.mesh('cutsource');
nv=double(rec.vertices);nt=double(rec.tetrahedra);chunk=100000;
fv=fopen(fullfile(folder,'vertices_m.bin'),'r','ieee-le');ft=fopen(fullfile(folder,'tets0.bin'),'r','ieee-le');fc=fopen(fullfile(folder,'entity_code.bin'),'r');
assert(fv>0&&ft>0&&fc>0);fileclose=onCleanup(@()close_files([fv ft fc])); %#ok<NASGU>
for first=floor((nv-1)/chunk)*chunk:-chunk:0
 count=min(chunk,nv-first);fseek(fv,first*24,'bof');v=fread(fv,[3 count],'double=>double');assert(size(v,2)==count);
 me.data.setVertex(int32(first),v);
end
for first=floor((nt-1)/chunk)*chunk:-chunk:0
 count=min(chunk,nt-first);fseek(ft,first*16,'bof');t=fread(ft,[4 count],'int32=>int32');fseek(fc,first,'bof');c=fread(fc,[1 count],'uint8=>int32');
 assert(size(t,2)==count&&numel(c)==count);me.data.setElem('tet',int32(first),t);me.data.setElemEntity('tet',int32(first),c);
end
me.data.createMesh;me.data.clearData;
assert(me.getNumElem('tet')==nt&&me.getNumVertex()==nv,'MeshData changed counts');
% A geometry-defining editable import sequence is retained for fixed-surface
% COMSOL Refine. No boolean or free remeshing is involved.
co.mesh.create('mesh1','geom1');ms=co.mesh('mesh1');ms.geometricModel('');ms.create('cutimp','Import');
im=ms.feature('cutimp');im.set('source','sequence');im.set('sequence','cutsource');im.set('buildsource','off');im.set('domelemsequence','on');im.set('unmesheddom','on');im.importData;ms.run('cutimp');
assert(ms.getNumElem('tet')==nt&&ms.getNumVertex()==nv,'Editable import changed counts');
domainFullCode=zeros(0,1);
for first=0:chunk:nt-1
 count=min(chunk,nt-first);fseek(ft,first*16,'bof');expected=fread(ft,[4 count],'int32=>int32');actual=ms.getElem('tet',int32(first),int32(count));
 assert(isequal(sort(actual,1),sort(expected,1)),'COMSOL changed element ordering; stop rather than guess material mapping');
 fseek(fc,first,'bof');codes=double(fread(fc,[1 count],'uint8=>uint8'));ent=double(ms.getElemEntity('tet',int32(first),int32(count)));ent=ent(:);
 lo=accumarray(ent,codes(:),[],@min,NaN);hi=accumarray(ent,codes(:),[],@max,NaN);ids=find(~isnan(lo));assert(all(lo(ids)==hi(ids)),'Imported domain merged material identities');
 if numel(domainFullCode)<numel(lo),domainFullCode(end+1:numel(lo),1)=0;end
 assert(all(domainFullCode(ids)==0|domainFullCode(ids)==lo(ids)));domainFullCode(ids)=lo(ids);
end
for first=0:chunk:nv-1
 count=min(chunk,nv-first);fseek(fv,first*24,'bof');expected=fread(fv,[3 count],'double=>double');actual=double(ms.getVertex(int32(first),int32(count)));
 assert(max(abs(actual(:)-expected(:)))<=L*1e-12,'COMSOL changed physical source vertices');
end
mapping=struct('domainIDsByCode',{cell(1,6)},'ellipseIDs',find(domainFullCode>6)');
for code=1:6
 ids=find(domainFullCode>0 & mod(domainFullCode-1,6)+1==code)';mapping.domainIDsByCode{code}=ids;
 if isempty(ids),continue;end
 name=sprintf('orig%d',code);ms.create(name,'ExplicitSelection');sel=ms.feature(name).selection;sel.geom(3);sel.set(ids);ms.feature(name).set('selshow','on');ms.run(name);
end
ms.create('ellvol','ExplicitSelection');sel=ms.feature('ellvol').selection;sel.geom(3);sel.set(mapping.ellipseIDs);ms.feature('ellvol').set('selshow','on');ms.run('ellvol');
co.geometricModel('mesh/mesh1');ms.run();
% Finalize may canonize domain IDs. Obtain the propagated source selections
% from the actual component topology, never reclassify material by centroid.
preFinalizeMapping=mapping;mapping=struct('domainIDsByCode',{cell(1,6)},'ellipseIDs',[]);
for code=1:6
 name=sprintf('orig%d',code);if ~any(strcmp(cell(ms.feature.tags),name)),continue;end
 tags=cell(ms.feature(name).componentSelection());assert(~isempty(tags),'Original source selection did not propagate through Finalize');
 ids=[];for j=1:numel(tags),ids=union(ids,double(co.selection(tags{j}).entities(3)));end
 mapping.domainIDsByCode{code}=reshape(ids,1,[]);
end
tags=cell(ms.feature('ellvol').componentSelection());assert(~isempty(tags),'Ellipsoid cut selection did not propagate through Finalize');
for j=1:numel(tags),mapping.ellipseIDs=union(mapping.ellipseIDs,double(co.selection(tags{j}).entities(3)));end
r3_json(fullfile(folder,'domain_mapping_after_finalize.json'),struct('before_finalize',preFinalizeMapping,'after_finalize',mapping,'method','Native propagated ExplicitSelection; no geometry or material classification changed'));
mphsave(m,fullfile(folder,'imported_before_validation.mph'),'optimize','size');
actualStats=mphmeshstats(m,'mesh1');
r3_json(fullfile(folder,'actual_import_mesh_quality.json'),struct('minquality',actualStats.minquality,...
 'measured_after_COMSOL_import',true,'tetrahedra',double(ms.getNumElem('tet'))));
assert(actualStats.minquality>0,'Actual imported mesh has nonpositive quality; never solve or delete failed cells');
[verified,sel]=r7b_validate_ct_source_mesh(m,'mesh1',src,folder,mapping);
% Correct the provenance wording in the reused validator. Checks unchanged.
verified.ellipse_label_from_actual_boolean_selection=false;
verified.ellipse_label_from_explicit_conforming_cut=true;
sel.ellipse_selection_source='Explicit cut-cell domain inheritance, not a centroid indicator and not COMSOL mesh Union';
tol=L*1e-7;sel.top=reshape(mphselectbox(m,'geom1',[-tol L+tol;-tol L+tol;L-tol L+tol],'boundary'),1,[]);
sel.bottom=reshape(mphselectbox(m,'geom1',[-tol L+tol;-tol L+tol;-tol tol],'boundary'),1,[]);assert(~isempty(sel.top)&&~isempty(sel.bottom));
sel.electrode_selections_pending_caller=false;sel.sample=sample;sel.formal_CT_target=strcmp(mode,'formal')||isRefined;
if isRefined
 anchor=jsondecode(fileread(fullfile(anchorFolder,'native_phase_verification.json')));
 before=anchor.actual_phase_region_ellipse_volume_m3;after=verified.actual_phase_region_ellipse_volume_m3;
 delta=max(abs(before(:)-after(:)));
 assert(delta/L^3<1e-9,'Explicit midpoint refinement changed phase/region/ellipse volumes');
 assert(abs(sum(after(:,2))/sum(before(:,2))-1)<1e-8,'Frozen M1 ellipsoid volume changed');
 proofPath=fullfile(folder,'fixed_geometry_refinement.json');proof=jsondecode(fileread(proofPath));
 proof.actual_COMSOL_import_verified=true;proof.actual_COMSOL_min_mesh_quality=actualStats.minquality;
 proof.maximum_phase_region_ellipse_volume_change_m3=delta;
 proof.relative_ellipse_volume_change=sum(after(:,2))/sum(before(:,2))-1;
 proof.external_geometry_verification_pending=true;proof.PDE_solved=false;
 r3_json(proofPath,proof);
end
r3_json(fullfile(folder,'selections.json'),sel);r3_json(fullfile(folder,'native_phase_verification.json'),verified);
r3_json(fullfile(folder,'ORIGINAL_CT_PHASE_PRESERVATION.json'),verified);
mphsave(m,fullfile(folder,'geometry_ready.mph'),'optimize','size');
r7b_mesh_sections(m,'mesh1',sel,folder,L);
r3_json(fullfile(folder,'IMPORT_COMPLETE_NOT_GEOMETRY_ACCEPTED.json'),struct('status','IMPORTED_AND_PHASE_CHECKED','sample',sample,...
 'mode',mode,'new_ellipsoid_surface_version',true,'old_union_skipped',true,'old_freetet_skipped',true,'PDE_solved',false));
fprintf('Run original geometry checker on actual_boundary.mat, then accept_geometry. No target PDE run yet.\n');
end
function check_binary_size(p,n),s=dir(p);assert(~isempty(s)&&s.bytes==n,'Binary file size mismatch: %s',p);end
function close_files(ids),for x=ids,try,fclose(x);catch,end;end;end
function remove_model(tag),try,com.comsol.model.util.ModelUtil.remove(tag);catch,end;end
