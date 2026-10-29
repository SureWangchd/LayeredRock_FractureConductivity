function [phase,stats,rejected] = r7_ct_validate_native_phase(q,h,labels_xyz)
% Validate the ORIGINAL piecewise-constant CT material, not artificial faces
% between neighbouring voxels of identical phase. No labels are modified.
% q: 3 x 4 x N physical tet coordinates. SAT is used only for opposite-phase
% AABB candidate voxels; centre/vertex sampling is never used to accept them.
n=size(q,3);shape=size(labels_xyz);assert(numel(shape)==3 && all(shape>1));
tol=1e-10;centers=reshape(mean(q,2),3,n)/h;
idx=min(shape(:)-1,max(0,floor(centers)))+1;
phase=reshape(double(labels_xyz(sub2ind(shape,idx(1,:),idx(2,:),idx(3,:)))),1,n);
lo=reshape(min(q,[],2),3,n)/h;hi=reshape(max(q,[],2),3,n)/h;
ilow=max(0,floor(lo+tol));ihigh=min(shape(:)-1,floor(hi-tol));
crossing=find(any(ilow~=ihigh,1));rejected=false(1,n);
stats=struct('tetrahedra_checked',n,'single_voxel_tetrahedra',n-numel(crossing),...
    'cross_voxel_tetrahedra',numel(crossing),'same_phase_aabb_tetrahedra',0,...
    'SAT_tested_tetrahedra',0,'SAT_opposite_phase_candidate_voxels',0,...
    'SAT_positive_opposite_phase_intersections',0,'rejected_tetrahedra',0,...
    'same_phase_crossing_tetrahedra',0,'coordinate_tolerance_native_voxels',tol);
for j=reshape(crossing,1,[])
    [xx,yy,zz]=ndgrid(ilow(1,j):ihigh(1,j),ilow(2,j):ihigh(2,j),ilow(3,j):ihigh(3,j));
    vox=[xx(:) yy(:) zz(:)]';
    cl=double(labels_xyz(sub2ind(shape,vox(1,:)+1,vox(2,:)+1,vox(3,:)+1)));
    other=vox(:,cl~=phase(j));
    if isempty(other)
        stats.same_phase_aabb_tetrahedra=stats.same_phase_aabb_tetrahedra+1;continue;
    end
    stats.SAT_tested_tetrahedra=stats.SAT_tested_tetrahedra+1;
    x=q(:,:,j)/h;axes=sat_axes(x);projection=axes'*x;
    tlo=min(projection,[],2);thi=max(projection,[],2);radius=.5*sum(abs(axes),1)';
    for k=1:size(other,2)
        mid=axes'*(other(:,k)+.5);
        overlap=min(thi,mid+radius)-max(tlo,mid-radius);
        stats.SAT_opposite_phase_candidate_voxels=stats.SAT_opposite_phase_candidate_voxels+1;
        if all(overlap>tol)
            % A strict positive intersection with another phase is forbidden.
            % Merely touching a face/edge/vertex has zero overlap on some axis.
            rejected(j)=true;
            stats.SAT_positive_opposite_phase_intersections=stats.SAT_positive_opposite_phase_intersections+1;
        end
    end
end
stats.rejected_tetrahedra=nnz(rejected);
stats.same_phase_crossing_tetrahedra=numel(crossing)-nnz(rejected);
end

function axes=sat_axes(q)
% Complete separating axes for a convex tet and an axis-aligned cube.
axes=eye(3);faces=[1 2 3;1 2 4;1 3 4;2 3 4];edges=nchoosek(1:4,2);
for j=1:size(faces,1)
    f=faces(j,:);axes(:,end+1)=cross(q(:,f(2))-q(:,f(1)),q(:,f(3))-q(:,f(1))); %#ok<AGROW>
end
for j=1:size(edges,1)
    e=q(:,edges(j,2))-q(:,edges(j,1));
    axes=[axes cross(repmat(e,1,3),eye(3),1)]; %#ok<AGROW>
end
norms=sqrt(sum(axes.^2,1));axes=axes(:,norms>1e-13)./norms(norms>1e-13);
end
