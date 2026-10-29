function [rows,fields]=r3_enrich_predict_v1(d,sigmas)
% Background-only Galerkin prediction with background-inverse residual bases.
% No COMSOL call, target solution, target factorization, or truth-file access.
% Initial basis: 4 axial groups x original phase x {phi0,xSI,zSI} sources.
% Version 1: normalized initial Gram cutoff 1e-10; residual energy cutoff 1e-14.
% Fixed development rule: <=40 extra vectors; E/(prediction-E)<=0.03.
% K0,DL,DH contain the physical t=1 m factor. DV=1 V in this project;
% every quadratic power divided by DV^2 is conductance in S, not S/m.
assert(all(sigmas>=d.c.sigma_high_S_m),'Residual bound requires nonnegative material replacement');
assert(mod(d.nseg,4)==0,'Background segment count must aggregate exactly into four groups');
assert(size(d.responses,2)==3*size(d.blockinfo,1),'Response/block layout mismatch');
free=d.free(:);nfull=numel(d.u0);u0=d.u0(:);DV=1;
Kf=d.K0(free,free);Lf=d.DL(free,free);Hf=d.DH(free,free);
assert(norm(Kf-Kf','fro')/norm(Kf,'fro')<1e-10,'Background operator is not physical symmetric');
Kf=(Kf+Kf')/2;
started=tic;factor=decomposition(Kf,'chol');factorSeconds=toc(started);
% The preparation cache contains response vectors, not a serialized factor.
% This function makes exactly one additional background factorization and
% explicitly records its shared cost; it never factors a target operator.
started=tic;P=zeros(numel(free),0);nonemptySourceGroups=0;
for group=1:4
 for phase=0:1
  blocks=find(ceil(d.blockinfo(:,1)/(d.nseg/4))==group & d.blockinfo(:,2)==phase);
  if isempty(blocks),continue;end
  nonemptySourceGroups=nonemptySourceGroups+1;
  for direction=1:3
   P(:,end+1)=sum(d.responses(:,3*blocks-3+direction),2); %#ok<AGROW>
  end
 end
end
% Version 1 rank guard: normalize column energies before a Gram eigentest.
% In a flat phase the phi0 response is proportional to its z-coordinate
% response; solver/assembly noise must not become an extra physical mode.
candidateColumns=size(P,2);diagE=sum(P.*(Kf*P),1);
positiveColumns=diagE>max(diagE)*1e-16;
assert(any(positiveColumns),'Initial response columns have zero energy');
Pn=P(:,positiveColumns)./sqrt(diagE(positiveColumns));
C=Pn'*(Kf*Pn);C=(C+C')/2;
[Q,eigenvalues]=eig(C,'vector');gramThreshold=1e-10*max(eigenvalues);
keep=eigenvalues>gramThreshold;
assert(any(keep),'Initial response space has zero rank');
preBasis=Pn*(Q(:,keep)./sqrt(eigenvalues(keep)'));
Z0=zeros(numel(free),0);
for j=1:size(preBasis,2)
 [Z0,kept]=append_energy_basis(Z0,preBasis(:,j),Kf);
 assert(kept,'Retained Gram direction became dependent on reorthogonalization');
end
discardedInitial=candidateColumns-size(Z0,2);
initialMinimumKeptEigenvalue=min(eigenvalues(keep));
if any(~keep),initialMaximumDiscardedEigenvalue=max(eigenvalues(~keep));else,initialMaximumDiscardedEigenvalue=0;end
assert(~isempty(Z0),'Initial response space has zero rank');
initialRank=size(Z0,2);initialOrthogonality=norm(Z0'*(Kf*Z0)-eye(initialRank),inf);
assert(initialOrthogonality<1e-7,'Background energy orthogonalization failed');
initialBasisSeconds=toc(started);
UL=u0'*(d.DL*u0);UH=u0'*(d.DH*u0);
loadL=d.DL(free,:)*u0;loadH=d.DH(free,:)*u0;
backgroundResidual=d.K0(free,:)*u0;
backgroundResidualNorm=norm(backgroundResidual);
rows=struct([]);fields=struct;
relativeTarget=0.03;maxExtra=40;
for sf=reshape(sigmas,1,[])
 caseStarted=tic;dl=sf-d.c.sigma_low_S_m;dh=sf-d.c.sigma_high_S_m;
 f=dl*loadL+dh*loadH;U=(dl*UL+dh*UH)/DV^2;
 Z=Z0;AZ=Kf*Z+dl*(Lf*Z)+dh*(Hf*Z);
 history=struct([]);extra=0;backsolves=0;backsolveSeconds=0;reducedSeconds=0;
 stoppingReason='uninitialized';
 while true
  reducedStarted=tic;
  H=Z'*AZ;H=(H+H')/2;q=Z'*f;
  alpha=-(H\q);wfree=Z*alpha;w=zeros(nfull,1);w(free)=wfree;v=u0+w;
  R=(-q'*alpha)/DV^2;raw=U-R;
  backgroundTerm=(w'*(d.K0*w))/DV^2;
  inclusionLowTerm=dl*(v'*(d.DL*v))/DV^2;
  inclusionHighTerm=dh*(v'*(d.DH*v))/DV^2;
  stable=backgroundTerm+inclusionLowTerm+inclusionHighTerm;
  cross=U+(f'*wfree)/DV^2;
  residual=-f-AZ*alpha;
  reducedSeconds=reducedSeconds+toc(reducedStarted);
  backStarted=tic;z=factor\residual;backsolveSeconds=backsolveSeconds+toc(backStarted);backsolves=backsolves+1;
  errorBound=(residual'*z)/DV^2;
  certifiedLower=stable-errorBound; % Do not clip to zero or hide a vacuous bound.
  if certifiedLower>0 && errorBound>=0,relativeBound=errorBound/certifiedLower;
  else,relativeBound=Inf;end
  tolerance=1e-8*max([abs(U) abs(stable) 1e-12]);
  boundViolation=errorBound < -tolerance || stable < -tolerance || stable>U+tolerance;
  record=struct('rank',size(Z,2),'extra_vectors',extra,'upper_deltaG_S',stable,...
   'subtractive_deltaG_S',raw,'cross_deltaG_S',cross,'error_bound_S',errorBound,...
   'lower_deltaG_S',certifiedLower,'certified_relative_bound',relativeBound,...
   'residual_norm_A',norm(residual),'galerkin_projection_residual',norm(Z'*residual,inf),...
   'upper_bound_violation',boundViolation);
  if isempty(history),history=record;else,history(end+1)=record;end %#ok<AGROW>
  if boundViolation
   stoppingReason='numerical_bound_violation_not_clipped';break;
  elseif relativeBound<=relativeTarget
   stoppingReason='background_residual_certificate_met';break;
  elseif extra>=maxExtra
   stoppingReason='additional_basis_budget_reached';break;
  end
  oldRank=size(Z,2);[Z,kept]=append_energy_basis(Z,z,Kf);
  if ~kept,stoppingReason='background_residual_direction_numerically_dependent';break;end
  newZ=Z(:,oldRank+1:end);AZ=[AZ Kf*newZ+dl*(Lf*newZ)+dh*(Hf*newZ)]; %#ok<AGROW>
  extra=extra+size(newZ,2);
 end
 r=struct('case_id',d.c.id,'method','vector4_background_residual_rankguard','segments',4,...
  'candidate_columns',size(P,2),'original_phase_source_groups',nonemptySourceGroups,...
  'initial_rank',initialRank,'initial_columns_discarded',discardedInitial,...
  'initial_gram_eigenvalue_relative_threshold',1e-10,...
  'initial_minimum_retained_gram_eigenvalue',initialMinimumKeptEigenvalue,...
  'initial_maximum_discarded_gram_eigenvalue',initialMaximumDiscardedEigenvalue,...
  'residual_direction_relative_energy_threshold',1e-14,...
  'rank',size(Z,2),'extra_basis_vectors',extra,'additional_basis_limit',maxExtra,...
  'sigma_f_S_m',sf,'G0_S',d.G0,'U_S',U,'predicted_deltaG_S',stable,...
  'subtractive_deltaG_S',raw,'cross_deltaG_S',cross,'R_pred_S',R,...
  'positive_background_term_S',backgroundTerm,'positive_inclusion_low_term_S',inclusionLowTerm,...
  'positive_inclusion_high_term_S',inclusionHighTerm,'residual_error_bound_S',errorBound,...
  'certified_lower_deltaG_S',certifiedLower,'certified_relative_error_bound',relativeBound,...
  'relative_error_bound_target',relativeTarget,'certificate_met',relativeBound<=relativeTarget && ~boundViolation,...
  'upper_bound_violation',boundViolation,'stopping_reason',stoppingReason,...
  'background_endcap_dV_V',diff(d.capV),'background_stationarity_residual_norm_A',backgroundResidualNorm,...
  'background_stationarity_cross_correction_S',2*(wfree'*backgroundResidual)/DV^2,...
  'initial_energy_orthogonality_error',initialOrthogonality,...
  'final_energy_orthogonality_error',norm(Z'*(Kf*Z)-eye(size(Z,2)),inf),...
  'positive_subtractive_relative_disagreement',abs(stable-raw)/max(abs(stable),1e-14),...
  'residual_norm_A',norm(residual),'galerkin_projection_residual',norm(Z'*residual,inf),...
  'background_factorizations_shared_in_function',1,'background_factor_seconds_shared',factorSeconds,...
  'basis_projection_seconds',initialBasisSeconds,'background_residual_backsolves',backsolves,...
  'background_residual_backsolve_seconds',backsolveSeconds,'reduced_evaluation_seconds',reducedSeconds,...
  'target_prediction_seconds',toc(caseStarted),'target_operator_factorizations',0,...
  'units','Conductance S; physical t included; DV=1 V',...
  'development_choice','Toy FEM informed vector response enrichment and residual certificate; fixed 4 groups, normalized Gram cutoff 1e-10, two-pass A0 orthogonalization, residual energy cutoff 1e-14, 40 extra directions, 3 percent certificate. Revision follows background-only rank-dependence diagnostics; no target field used. No target truth used.');
 if isempty(rows),rows=r;else,rows(end+1)=r;end %#ok<AGROW>
 key=['s' strrep(num2str(sf),'.','p')];
 fields.(key)=struct('row',r,'w',w,'alpha',alpha,'history',history,...
  'basis_free',Z,'free',free,'initial_gram_eigenvalues',eigenvalues,'basis_sources','Background inverse only; geometry/phase/coordinate responses plus algebraic target residuals; no target solution.');
end
end

function [Z,kept]=append_energy_basis(Z,v,K0)
% Two-pass modified Gram-Schmidt in physical background energy.
initial=v'*(K0*v);kept=false;
if ~isfinite(initial) || initial<=0,return;end
for pass=1:2
 if ~isempty(Z),v=v-Z*(Z'*(K0*v));end
end
remaining=v'*(K0*v);
if ~isfinite(remaining) || remaining<=max(1e-25,1e-14*initial),return;end
Z(:,end+1)=v/sqrt(remaining);kept=true;
end
