function summary = r7b_mesh_sections(model,meshTag,sel,outputFolder,L)
%R7B_MESH_SECTIONS Actual tetrahedral/CT phase sections for internal viewing.
% No solve, no geometry modification, no ellipse indicator classification.
% Reads tet connectivity in chunks and exports only intersections with the
% three physical planes x/y/z=L/2. Material labels come from actual domain
% selections supplied by sel. Coordinates are SI metres.
%
% At a material face exactly coincident with a section, use its positive
% coordinate side: negative-side touching tets are omitted, preventing
% duplicate coincident polygons with ambiguous material labels. This is an
% explicitly recorded one-sided trace, not an averaged measured CT slice.
if ~isfolder(outputFolder),mkdir(outputFolder);end
required={'all_domains','original_phase','region','is_ellipse'};
assert(all(isfield(sel,required)),'r7b:section:labels','Actual domain selections required');
assert(isscalar(L) && isfinite(L) && L>0,'r7b:section:scale','Positive sample L in metres required');
ids=reshape(double(sel.all_domains),[],1);assert(~isempty(ids) && all(ids>=1 & ids==floor(ids)));
codePhase=-ones(max(ids),1);codeRegion=-ones(max(ids),1);codeEllipse=false(max(ids),1);
codePhase(ids)=double(sel.original_phase(:));codeRegion(ids)=double(sel.region(:));codeEllipse(ids)=logical(sel.is_ellipse(:));
assert(all(ismember(codePhase(ids),[0 1])) && all(ismember(codeRegion(ids),[0 1 2])),'r7b:section:labels','Invalid inherited phase/region');
me=model.component('comp1').mesh(meshTag);vertices=double(me.getVertex);nt=double(me.getNumElem('tet'));
assert(size(vertices,1)==3 && nt>0,'r7b:section:mesh','Actual 3D tet mesh required');
timer=tic;chunk=100000;plane=L/2;tol=L*1e-11;
polys=cell(3,1);counts=cell(3,1);domains=cell(3,1);tetids=cell(3,1);
for axis=1:3,polys{axis}={};counts{axis}={};domains{axis}={};tetids{axis}={};end
for first=0:chunk:nt-1
 count=min(chunk,nt-first);t=double(me.getElem('tet',int32(first),int32(count)))+1;
 ent=reshape(double(me.getElemEntity('tet',int32(first),int32(count))),[],1);
 assert(max(ent)<=numel(codePhase) && all(codePhase(ent)>=0),'r7b:section:mapping','Unlabelled actual domain');
 q=reshape(vertices(:,t),3,4,count);
 for axis=1:3
  d=reshape(q(axis,:,:),4,count)-plane;
  candidates=find(min(d,[],1)<=tol & max(d,[],1)>tol);
  if isempty(candidates),continue;end
  block=nan(3,4,numel(candidates));nv=zeros(numel(candidates),1,'uint8');keep=false(numel(candidates),1);
  for k=1:numel(candidates)
   j=candidates(k);[p,n]=cut_tet(q(:,:,j),axis,plane,tol);
   if n<3,continue;end
   other=setdiff(1:3,axis);xy=p(other,1:n);area=abs(sum(xy(1,:).*xy(2,[2:n 1])-xy(2,:).*xy(1,[2:n 1])))/2;
   if area<=L^2*1e-25,continue;end
   block(:,1:n,k)=p(:,1:n);nv(k)=uint8(n);keep(k)=true;
  end
  if any(keep)
   polys{axis}{end+1}=block(:,:,keep); %#ok<AGROW>
   counts{axis}{end+1}=nv(keep); %#ok<AGROW>
   domains{axis}{end+1}=int32(ent(candidates(keep))); %#ok<AGROW>
   tetids{axis}{end+1}=int32(first+candidates(keep)-1); %#ok<AGROW>
  end
 end
end
[sample,scope]=scope_info(sel,L);
metadata=struct('schema','r7b_actual_tetrahedral_sections_v1','sample',sample,'geometry_scope',scope,...
 'mesh_tag',char(meshTag),'sample_L_m',L,'actual_source_tetrahedra',nt,'actual_source_vertices',size(vertices,2),...
 'planes','x=L/2; y=L/2; z=L/2','purpose','INTERNAL_GEOMETRY_AND_MESH_CHECK_NOT_PUBLICATION_FIGURE',...
 'coordinates','Actual plane/tetrahedron intersections, SI metres; no analytic ellipse substituted',...
 'phase_and_ellipse_source','Actual COMSOL domain sel; no point-based ellipse relabeling',...
 'on_interface_material_convention','positive-coordinate-side limit when plane coincides with an actual material/mesh face',...
 'intersection_coordinate_tolerance_m',tol,'no_COMSOL_solve',true);
records=cell(3,1);
for axis=1:3
 assert(~isempty(polys{axis}),'r7b:section:empty','No actual polygons for section axis %d',axis);
 q=cat(3,polys{axis}{:});polygon_vertices_m=permute(q,[3 2 1]); %#ok<NASGU>
 polygon_vertex_count=cat(1,counts{axis}{:});domain_id=cat(1,domains{axis}{:});tetrahedron_index0=cat(2,tetids{axis}{:});tetrahedron_index0=tetrahedron_index0(:); %#ok<NASGU>
 original_phase=uint8(codePhase(double(domain_id)));region=uint8(codeRegion(double(domain_id)));inside_ellipsoid=logical(codeEllipse(double(domain_id))); %#ok<NASGU>
 axis_index0=int32(axis-1);plane_coordinate_m=plane;sample_L_m=L;coordinate_tolerance_m=tol; %#ok<NASGU>
 metadata_json=jsonencode(metadata); %#ok<NASGU>
 names='xyz';filename=fullfile(outputFolder,['actual_mesh_section_' names(axis) '.mat']);
 save(filename,'polygon_vertices_m','polygon_vertex_count','domain_id','tetrahedron_index0',...
  'original_phase','region','inside_ellipsoid','axis_index0','plane_coordinate_m',...
  'sample_L_m','coordinate_tolerance_m','metadata_json','-v7');
 other=setdiff(1:3,axis);areas=zeros(numel(polygon_vertex_count),1);
 for n=3:4
  take=find(polygon_vertex_count==n);xy=q(other,1:n,take);
  areas(take)=reshape(abs(sum(xy(1,:,:).*xy(2,[2:n 1],:)-xy(2,:,:).*xy(1,[2:n 1],:),2))/2,[],1);
 end
 records{axis}=struct('axis',names(axis),'file',filename,'polygons',numel(polygon_vertex_count),...
  'triangles',nnz(polygon_vertex_count==3),'quadrilaterals',nnz(polygon_vertex_count==4),...
  'actual_total_polygon_area_m2',sum(areas),'expected_whole_section_area_m2',L^2,...
  'section_area_relative_error',sum(areas)/L^2-1,...
  'original_pore_section_area_m2',sum(areas(original_phase==0)),...
  'original_skeleton_section_area_m2',sum(areas(original_phase==1)),...
  'inside_ellipsoid_section_area_m2',sum(areas(inside_ellipsoid)),...
  'internal_check_only',true);
end
summary=metadata;summary.sections=vertcat(records{:});summary.seconds=toc(timer);
r3_json(fullfile(outputFolder,'actual_mesh_sections.json'),summary);
fprintf('R7B_MESH_SECTIONS %s polygons %d/%d/%d seconds %.2f\n',scope,summary.sections.polygons,toc(timer));
end

function [p,n]=cut_tet(v,axis,plane,tol)
d=v(axis,:)-plane;p=zeros(3,0);
for j=1:4
 if abs(d(j))<=tol,x=v(:,j);x(axis)=plane;p(:,end+1)=x;end %#ok<AGROW>
end
edges=[1 2;1 3;1 4;2 3;2 4;3 4];
for k=1:6
 i=edges(k,1);j=edges(k,2);
 if (d(i)<-tol && d(j)>tol)||(d(j)<-tol && d(i)>tol)
  x=v(:,i)+d(i)/(d(i)-d(j))*(v(:,j)-v(:,i));x(axis)=plane;p(:,end+1)=x; %#ok<AGROW>
 end
end
uniquePoints=zeros(3,0);
for k=1:size(p,2)
 if isempty(uniquePoints)||all(vecnorm(uniquePoints-p(:,k),2,1)>tol),uniquePoints(:,end+1)=p(:,k);end %#ok<AGROW>
end
p=uniquePoints;n=size(p,2);assert(n<=4,'r7b:section:polygon','Tet-plane section has more than four distinct vertices');
if n<3,return;end
other=setdiff(1:3,axis);mid=mean(p(other,:),2);angles=atan2(p(other(2),:)-mid(2),p(other(1),:)-mid(1));
[~,order]=sort(angles);p=p(:,order);
end

function [sample,scope]=scope_info(sel,L)
sample='UNKNOWN';scope='ACTUAL_SCALE_INTERNAL_CHECK_NOT_FORMAL_CT';
if isfield(sel,'sample'),sample=char(sel.sample);end
if isfield(sel,'geometry_scope'),scope=char(sel.geometry_scope);return;end
names={'S1','S5'};spacing=[8.683e-6 3.997e-6];
for k=1:2
 n=L/spacing(k);
 if abs(n-128)<1e-8,sample=names{k};scope='prescribed_native_centre128_actual_mesh';return;
 elseif abs(n-32)<1e-8,sample=names{k};scope='technical_native_centre32_NOT_FORMAL_CT';return;end
end
end
