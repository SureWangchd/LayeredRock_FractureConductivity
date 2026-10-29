function r=r7_thinmesh_check(md,ms,c,ellipseDomainIDs,folder)
%R7_THINMESH_CHECK Read-only central-thickness test on mphmeshstats tetrahedra.
% r=r7_thinmesh_check(md,ms,c,ellipseDomainIDs,folder)
% Optional domain IDs restrict tetrahedra to the original ellipsoid; folder
% writes diagnostics only. No model, geometry, material or mesh is modified.
% c uses cx_mm/cy_mm/cz_mm, basis_n, c_mm and L_mm. For nonstandard meshes,
% set c.mesh_coordinate_to_mm explicitly. A mesh-statistics tetrahedron uses
% its four corner vertices: this is a linear-chord resolution check, not an
% exact intersection with a curved quadratic geometric element.
if nargin<4,ellipseDomainIDs=[];end
if nargin<5,folder='';end
it=find(strcmp(ms.types,'tet'));
assert(numel(it)==1,'r7:thinmesh:tet','Expected one tetrahedron element type');
vertices=double(md.vertex);
assert(size(vertices,1)==3,'r7:thinmesh:dimension','Expected a 3D mesh');
tet=double(md.elem{it}(1:4,:))+1;
entities=reshape(double(md.elementity{it}),1,[]);
assert(size(tet,2)==numel(entities),'r7:thinmesh:entities','Element/entity lengths disagree');
center=[getvalue(c,'cx_mm',60);getvalue(c,'cy_mm',60);getvalue(c,'cz_mm',60)];
if isfield(c,'basis_n'),normal=double(c.basis_n(:));
else,theta=getvalue(c,'theta_deg',60);normal=[-sind(theta);0;cosd(theta)];end
assert(numel(normal)==3 && norm(normal)>0,'r7:thinmesh:normal','Invalid thin-axis direction');
normal=normal/norm(normal);
halfThickness=getvalue(c,'c_mm',1.5);
assert(isscalar(halfThickness) && halfThickness>0,'r7:thinmesh:thickness','Positive c_mm is required');
if isfield(c,'mesh_coordinate_to_mm')
 scale=double(c.mesh_coordinate_to_mm);
else
 expectedL=getvalue(c,'L_mm',120);
 span=max(max(vertices,[],2)-min(vertices,[],2));
 if abs(span-expectedL)<max(1e-7,1e-7*expectedL),scale=1;
 elseif abs(span*1000-expectedL)<max(1e-7,1e-7*expectedL),scale=1000;
 else,error('r7:thinmesh:units','Mesh span %.16g does not identify mm or m; set c.mesh_coordinate_to_mm',span);end
end
vertices=vertices*scale;
assert(isfinite(scale) && scale>0,'r7:thinmesh:units','Invalid mesh-coordinate scale');
tol=max(1e-8,halfThickness*1e-8); % mm; also clusters numerically shared faces
p1=center-normal*halfThickness;p2=center+normal*halfThickness;
lineMin=min(p1,p2)-tol;lineMax=max(p1,p2)+tol;
if isempty(ellipseDomainIDs),eligible=1:size(tet,2);
else,eligible=find(ismember(entities,double(ellipseDomainIDs(:))));end
raw=zeros(0,5);candidateCount=0;degenerateCount=0;
for offset=1:20000:numel(eligible)
 ids=eligible(offset:min(offset+19999,numel(eligible)));corners=tet(:,ids);
 inside=true(1,numel(ids));
 for ax=1:3
  vals=reshape(vertices(ax,corners(:)),4,[]);
  inside=inside & max(vals,[],1)>=lineMin(ax) & min(vals,[],1)<=lineMax(ax);
 end
 ids=ids(inside);candidateCount=candidateCount+numel(ids);
 for k=ids
  v=vertices(:,tet(:,k));edge=v(:,2:4)-v(:,1);
  if rcond(edge)<1e-14,degenerateCount=degenerateCount+1;continue;end
  aa=edge\(center-v(:,1));bb=edge\normal;
  alpha=[1-sum(aa);aa];beta=[-sum(bb);bb];
  lo=-halfThickness;hi=halfThickness;valid=true;
  betaTol=1e-13/max(max(vecnorm(edge)),eps);
  for q=1:4
   if beta(q)>betaTol,lo=max(lo,-alpha(q)/beta(q));
   elseif beta(q)<-betaTol,hi=min(hi,-alpha(q)/beta(q));
   elseif alpha(q)<-1e-10,valid=false;break;
   end
  end
  if valid && hi-lo>tol
   raw(end+1,:)=[k entities(k) lo hi hi-lo]; %#ok<AGROW>
  end
 end
end
if ~isempty(raw),raw=sortrows(raw,[3 4 1]);end

% All endpoints are clustered, then each physical segment is counted once.
% This removes duplicate coverage if the line lies on a shared face/edge.
breaks=sort([-halfThickness;halfThickness;reshape(raw(:,3:4),[],1)]);
clustered=zeros(0,1);i=1;
while i<=numel(breaks)
 j=i;while j<numel(breaks) && breaks(j+1)-breaks(i)<=tol,j=j+1;end
 clustered(end+1,1)=mean(breaks(i:j)); %#ok<AGROW>
 i=j+1;
end
segments=zeros(0,5);gaps=zeros(0,3);
for i=1:numel(clustered)-1
 lo=clustered(i);hi=clustered(i+1);mid=(lo+hi)/2;
 if hi-lo<=tol,continue;end
 covering=find(raw(:,3)<mid & raw(:,4)>mid);
 if isempty(covering),gaps(end+1,:)=[lo hi hi-lo]; %#ok<AGROW>
 else,segments(end+1,:)=[lo hi hi-lo numel(covering) raw(covering(1),1)]; %#ok<AGROW>
 end
end
covered=sum(segments(:,3));gap=sum(gaps(:,3));expected=2*halfThickness;
if isempty(segments),stats=[NaN NaN NaN NaN];
else,stats=[min(segments(:,3)) median(segments(:,3)) mean(segments(:,3)) max(segments(:,3))];end
if isempty(raw),uniqueIntervals=0;
else
 % For a conforming straight-sided mesh shared-face duplicates have equal
 % interval endpoints. Distinct endpoints stay visible in atomic segments.
 canonical=round(raw(:,3:4)/tol);uniqueIntervals=size(unique(canonical,'rows'),1);
end
r=struct('schema','round7_central_thinmesh_v1','case_id',getvalue(c,'id',''),...
 'read_only_mesh_diagnostic',true,'method','Linear-tetrahedron barycentric line clipping; union coverage with shared-face duplicate removal',...
 'center_mm',center','normal_unit_vector',normal','s_range_mm',[-halfThickness halfThickness],...
 'mesh_coordinate_to_mm',scale,'expected_central_thickness_mm',expected,...
 'eligible_tetrahedra',numel(eligible),'bbox_candidate_tetrahedra',candidateCount,...
 'degenerate_candidates_skipped',degenerateCount,'positive_chord_tetrahedra',size(raw,1),...
 'unique_chord_intervals',uniqueIntervals,'resolved_nonoverlapping_thickness_segments',size(segments,1),...
 'covered_length_mm',covered,'uncovered_length_mm',gap,'relative_uncovered_length',gap/expected,...
 'raw_sum_chord_lengths_mm',sum(raw(:,5)),'duplicate_overlap_length_mm',max(0,sum(raw(:,5))-covered),...
 'min_segment_length_mm',stats(1),'median_segment_length_mm',stats(2),...
 'mean_segment_length_mm',stats(3),'max_segment_length_mm',stats(4),...
 'endpoint_clustering_tolerance_mm',tol,'entire_thickness_covered_within_1e_5',gap/expected<1e-5,...
 'six_or_more_resolved_segments',size(segments,1)>=6,...
 'interpretation','This checks the central thin-axis chord only. Curved-element boundary geometry and tip regions require separate geometry/mesh convergence checks; it does not establish solution convergence.');
if ~isempty(folder)
 if ~isfolder(folder),mkdir(folder);end
 writetable(array2table(raw,'VariableNames',{'tetrahedron_index_1based','domain_id','s_lower_mm','s_upper_mm','chord_length_mm'}),fullfile(folder,'thinmesh_raw_tetrahedron_intervals.csv'));
 writetable(array2table(segments,'VariableNames',{'s_lower_mm','s_upper_mm','length_mm','covering_tetrahedron_count','representative_tetrahedron_index_1based'}),fullfile(folder,'thinmesh_unique_intervals.csv'));
 writetable(array2table(gaps,'VariableNames',{'s_lower_mm','s_upper_mm','length_mm'}),fullfile(folder,'thinmesh_uncovered_intervals.csv'));
 fid=fopen(fullfile(folder,'thinmesh_check.json'),'w','n','UTF-8');assert(fid>=0,'Cannot write thinmesh JSON');
 cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
 fprintf(fid,'%s\n',jsonencode(r,'PrettyPrint',true));
end
end

function x=getvalue(s,name,fallback)
if isfield(s,name),x=s.(name);else,x=fallback;end
end
