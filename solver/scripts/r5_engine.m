function [results,fields]=r5_engine(d,sf,layouts,Ks)
% Self-consistent COMSOL-FE / analytic-periphery combination, no COMSOL calls.
% d.C integrates physical V edge shapes against cos(n*pi*x/.12), n=0..64,
% lower side first. This function never reads a complete target solution.
started=tic;
assert(isscalar(sf) && (isnan(sf) || (isfinite(sf)&&sf>0)), 'Invalid fracture conductivity');
assert(isstruct(layouts)&&~isempty(layouts),'At least one layout is required');
Ks=double(Ks(:)');assert(~isempty(Ks)&&all(Ks==floor(Ks)&Ks>=0&Ks<=64),'K must be 0..64');
assert(numel(unique(Ks))==numel(Ks),'Duplicate K values');
W=.12;t=1;DV=1;nmax=64;norms=[W repmat(W/2,1,nmax)];
n=numel(d.u0);free=double(d.free(:));fixed=double(d.fixed(:));
top=double(d.top(:));bottom=double(d.bottom(:));
assert(size(d.C,1)==n && size(d.C,2)==130,'Expected two interfaces, n=0..64');
assert(numel(free)+numel(fixed)==n && numel(unique([free;fixed]))==n,'Invalid DOF partition');
assert(isequal(sort(fixed),sort([top;bottom])),'Fixed DOFs must be exactly the electrodes');
assert(isempty(intersect(top,bottom)),'Electrodes overlap');
assert(norm(d.C(fixed,:),'fro')==0,'Port integrals must not include electrode DOFs');
if isnan(sf)
 K=d.K0;
else
 K=d.K0+(sf-.01)*d.DL+(sf-.1)*d.DH;
end
assert(issparse(K),'Physical operator must remain sparse');
symmetry=norm(K-K','fro')/max(norm(K,'fro'),realmin);
assert(symmetry<1e-10,'Physical operator is not symmetric');
K=(K+K')/2;
uFixed=zeros(n,1);uFixed(top)=DV; % Real external affine drive, not port truth.
assert(max(abs(d.u0(fixed)-uFixed(fixed)))<1e-8,'Background electrode values differ');
b=-K(free,fixed)*uFixed(fixed);
affineConstant=uFixed(fixed)'*(K(fixed,fixed)*uFixed(fixed));
A0=K(free,free);ts=tic;F=decomposition(A0,'chol');factorSeconds=toc(ts);
ts=tic;uf0=F\b;referenceSolveSeconds=toc(ts);
uReference=uFixed;uReference(free)=uf0;
referencePower=uReference'*(K*uReference);
referenceAffinePower=affineConstant-b'*uf0;
assert(abs(referenceAffinePower-referencePower)<2e-7*max(abs(referencePower),realmin),...
 'Reference affine energy does not close');
modeColumns=[2:65 67:130];B=d.C(free,modeColumns);
ts=tic;Q=F\B;R=full(B'*Q);R=(R+R')/2;p=full(B'*uf0);
responseSeconds=toc(ts);
[Vr,Er]=eig(R);er=real(diag(Er));
assert(min(er)>-1e-10*max(max(er),realmin),'Port response Gram matrix is not PSD');
Rhalf=Vr*diag(sqrt(max(er,0)))*Vr';Rhalf=(Rhalf+Rhalf')/2;
assemblySeconds=toc(started);
results=struct([]);fields=cell(numel(layouts)*numel(Ks),1);count=0;
for il=1:numel(layouts)
 item=layouts(il);actual=double(item.Y_actual);reference=double(item.Y_reference);
 assert(isequal(size(actual),[2 65])&&isequal(size(reference),[2 65]),'Expected Y arrays 2 by 65');
 assert(all(isfinite(actual(:))&actual(:)>0)&&all(isfinite(reference(:))&reference(:)>0),'All complete admittances must be positive');
 assert(max(abs(actual(:,1)-reference(:,1)))<1e-11,'The common zero mode must match');
 expectedY0=[.013157894736842105/.03;.029411764705882353/.03];
 assert(max(abs(reference(:,1)-expectedY0))<1e-10,'Wrong fixed-buffer zero-mode reference');
 if il>1
  assert(norm(reference-layouts(1).Y_reference,'fro')<1e-11*norm(reference,'fro'),'Layouts must share one M0 exterior');
 end
 for ik=1:numel(Ks)
  solveClock=tic;order=Ks(ik);coefficient=zeros(128,1);
  if order>0
   coefficient(1:order)=t*(actual(1,2:order+1)-reference(1,2:order+1))'/(W/2);
   coefficient(65:64+order)=t*(actual(2,2:order+1)-reference(2,2:order+1))'/(W/2);
  end
  % Signed Woodbury update. A negative response correction is retained.
  active=find(coefficient~=0);q=zeros(128,1);smallRcond=1;
  if ~isempty(active)
   scales=sqrt(abs(coefficient(active)));signs=sign(coefficient(active));
   small=diag(signs)+(scales*scales').*R(active,active);small=(small+small')/2;
   smallRcond=rcond(small);
   assert(smallRcond>1e-14,'Ill-conditioned modal correction solve');
   alpha=small\(scales.*p(active));q(active)=scales.*alpha;
  end
  uf=uf0-Q*q;u=uFixed;u(free)=uf;
  trace=full(d.C'*u);v=reshape(trace,65,2)'./norms;
  projected=full(B'*uf);hFree=B*(coefficient.*projected);
  kForce=K*u;reaction=kForce;reaction(free)=reaction(free)+hFree;
  referenceVolumePower=u'*kForce;
  correctionPower=sum(coefficient.*projected.^2);
  power=referenceVolumePower+correctionPower;
  affinePower=affineConstant-b'*uf;
  topCurrent=sum(reaction(top));bottomCurrent=-sum(reaction(bottom));
  topMeanCurrent=t*W*actual(2,1)*(DV-v(2,1));
  bottomMeanCurrent=t*W*actual(1,1)*v(1,1);
  freeResidual=reaction(free);
  residualScale=max(norm(A0*uf)+norm(b)+norm(hFree),realmin);
  relativeResidual=norm(freeResidual)/residualScale;
  currentScale=max([abs(topCurrent),abs(bottomCurrent),abs(power/DV),realmin]);
  fluxBalance=abs(topCurrent-bottomCurrent)/currentScale;
  reactionPowerBalance=abs(power-DV*topCurrent)/max(abs(power),realmin);
  affinePowerBalance=abs(power-affinePower)/max(abs(power),realmin);
  meanCurrentBalance=max(abs([topCurrent-topMeanCurrent,bottomCurrent-bottomMeanCurrent]))/currentScale;
  passive=eye(128)+Rhalf*(coefficient.*Rhalf);passive=(passive+passive')/2;
  [~,cholFlag]=chol(passive);passiveMin=min(eig(passive));
  assert(cholFlag==0&&passiveMin>0,'Combined physical operator lost positive definiteness');
  assert(relativeResidual<1e-8,'Self-consistent free-DOF residual is too large');
  assert(max([fluxBalance,reactionPowerBalance,affinePowerBalance,meanCurrentBalance])<2e-6,...
   'Combined model current, affine drive or complete power does not close');
  count=count+1;
  row=struct('layout',char(string(item.layout)),'K',order,'sigma_f_S_m',sf,...
   'is_true_background',isnan(sf),'G_S',power/DV^2,'G0_S',double(d.G0),...
   'deltaG_S',power/DV^2-double(d.G0),'M0_G_S',referencePower/DV^2,...
   'reference_volume_power_W',referenceVolumePower,'modal_correction_power_W',correctionPower,...
   'total_power_W',power,'affine_power_W',affinePower,'affine_constant_W',affineConstant,...
   'top_current_in_A',topCurrent,'bottom_current_out_A',bottomCurrent,...
   'top_current_from_mean_A',topMeanCurrent,'bottom_current_from_mean_A',bottomMeanCurrent,...
   'flux_balance_relative',fluxBalance,'reaction_power_balance_relative',reactionPowerBalance,...
   'affine_power_balance_relative',affinePowerBalance,'mean_current_balance_relative',meanCurrentBalance,...
   'minimum_passivity_eigenvalue',passiveMin,'passivity_cholesky_flag',cholFlag,...
   'free_residual_relative',relativeResidual,'free_residual_norm_A',norm(freeResidual),...
   'free_residual_max_abs_A',max(abs(freeResidual)),'modal_solve_rcond',smallRcond,...
   'port_v',v,'mode_numbers',0:64,'active_correction_rank',numel(active),...
   'negative_correction_count',sum(coefficient<0),'physical_matrix_symmetry_relative',symmetry,...
   'factor_seconds_shared',factorSeconds,'reference_solve_seconds_shared',referenceSolveSeconds,...
   'response_backsolve_seconds_shared',responseSeconds,'preparation_seconds_shared',assemblySeconds,...
   'layout_solve_seconds',toc(solveClock),'ndofs',n,'nfree',numel(free),...
   'thickness_m',t,'DV_V',DV,'target_truth_used',false,...
   'method','COMSOL physical FE M0 plus signed analytic modal response correction');
  if isempty(results),results=row;else,results(end+1)=row;end %#ok<AGROW>
  fields{count}=u;
 end
end
end
