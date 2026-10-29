function [mapped,info,targetToSource]=r3_map_dofs(u,sourceXmesh,targetXmesh)
% Explicit one-to-one transfer between unchanged-mesh DOF orderings.
% Match the physical field name AND every DOF coordinate. No interpolation.
% Returns targetToSource(j): source row used for target row j.
sc=double(sourceXmesh.dofs.coords);tc=double(targetXmesh.dofs.coords);
assert(size(sc,1)==size(tc,1),'Spatial coordinate dimensions differ');
assert(size(sc,2)==size(tc,2),'Different DOF counts: not a same-mesh transfer');
n=size(sc,2);assert(size(u,1)==n,'Solution row count does not match source DOFs');
snames=string(sourceXmesh.dofs.dofnames);tnames=string(targetXmesh.dofs.dofnames);
si=double(sourceXmesh.dofs.nameinds(:))+1;ti=double(targetXmesh.dofs.nameinds(:))+1;
assert(numel(si)==n&&numel(ti)==n,'Incomplete DOF field indices');
assert(all(si>=1&si<=numel(snames))&&all(ti>=1&ti<=numel(tnames)),'Unexpected field-name index convention');
sf=reshape(snames(si),[],1);tf=reshape(tnames(ti),[],1);
orderedFieldsSame=isequal(sf,tf);
orderedCoordsExact=isequal(sc,tc);
orderedCoordsEqualNaN=isequaln(sc,tc);
info=struct('method','physical_field_name_and_coordinate_bijection_no_interpolation', ...
 'source_ndofs',n,'target_ndofs',n,'ordered_field_names_equal',orderedFieldsSame, ...
 'ordered_coordinates_isequal',orderedCoordsExact,'ordered_coordinates_isequaln',orderedCoordsEqualNaN, ...
 'source_nan_coordinate_count',sum(isnan(sc(:))),'target_nan_coordinate_count',sum(isnan(tc(:))));
if orderedFieldsSame&&orderedCoordsEqualNaN
 targetToSource=(1:n)';
 info.mapping_kind='identity_exact_coordinates';
 if ~orderedCoordsExact,info.mapping_kind='identity_equaln_only';end
 info.absolute_coordinate_tolerance=0;info.max_matched_coordinate_difference=0;
else
 assert(all(isfinite(sc(:)))&&all(isfinite(tc(:))), ...
  'Nonfinite coordinates with changed ordering cannot establish a coordinate-based bijection');
 extent=max([max(sc,[],2)-min(sc,[],2);max(tc,[],2)-min(tc,[],2)]);
 scale=max([1;abs(sc(:));abs(tc(:))]);
 % Raw xmesh coordinate units. This is far below any physical mesh length;
 % matching must still have EXACTLY ONE candidate and a global bijection.
 atol=max(1e-11*max(1,extent),128*eps(scale));
 targetToSource=zeros(n,1);fields=unique(tf);
 assert(isequal(sort(unique(sf)),sort(fields)),'Physical field sets differ');
 for k=1:numel(fields)
  sids=find(sf==fields(k));tids=find(tf==fields(k));
  assert(numel(sids)==numel(tids),'Field DOF counts differ for %s',char(fields(k)));
  [matched,candidates]=ismembertol(tc(:,tids)',sc(:,sids)',atol, ...
   'ByRows',true,'DataScale',1,'OutputAllIndices',true);
  counts=cellfun(@numel,candidates);
  assert(all(matched)&&all(counts==1), ...
   'Coordinate matching is incomplete or ambiguous for %s: missing=%d ambiguous=%d', ...
   char(fields(k)),sum(~matched),sum(counts>1));
  local=cellfun(@(v)v(1),candidates);
  targetToSource(tids)=sids(local);
 end
 assert(all(targetToSource>=1)&&numel(unique(targetToSource))==n, ...
  'Coordinate match is not a bijection');
 assert(isequal(sf(targetToSource),tf),'Mapped field names differ');
 err=max(abs(sc(:,targetToSource)-tc),[],'all');
 assert(err<=atol,'Matched coordinates exceed mapping tolerance');
 info.absolute_coordinate_tolerance=atol;info.max_matched_coordinate_difference=err;
 if all(targetToSource==(1:n)'),info.mapping_kind='identity_with_coordinate_roundoff';
 else,info.mapping_kind='explicit_coordinate_field_permutation';end
end
info.permuted_dof_count=sum(targetToSource~=(1:n)');
info.bijective=numel(unique(targetToSource))==n;
mapped=u(targetToSource,:);
end
