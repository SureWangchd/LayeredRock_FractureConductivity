function info=r6_domination_check(d,sf,layout)
% Read-only finite-matrix sufficient condition for the Round 6 B bound.
% No target field, no COMSOL call, and no modification of r5_engine.
% With B=C(free,nonzero modes), A0=K_M0(free,free), D=t*Yref/N,
% lambda_max(sqrt(D)*B'*A0^-1*B*sqrt(D)) <= 1 implies A0 >= B*D*B'.
% Adding the actual signed modal corrections then preserves domination by
% each changed mode's positive complete target admittance.
started=tic;W=.12;t=1;N=W/2;
free=double(d.free(:));fixed=double(d.fixed(:));
assert(numel(free)+numel(fixed)==numel(d.u0),'Invalid DOF partition');
assert(isempty(intersect(free,fixed)),'Free/fixed DOFs overlap');
assert(norm(d.C(fixed,:),'fro')==0,'Port corrections touch fixed electrodes');
Y=double(layout.Y_reference);
assert(isequal(size(Y),[2 65])&&all(isfinite(Y(:))&Y(:)>0),'Invalid complete reference admittance');
if isnan(sf)
 K=d.K0;
else
 assert(isscalar(sf)&&isfinite(sf)&&sf>0,'Invalid fracture conductivity');
 K=d.K0+(sf-.01)*d.DL+(sf-.1)*d.DH;
end
skew=norm(K-K','fro')/max(norm(K,'fro'),realmin);
assert(skew<1e-10,'Physical M0 matrix is not symmetric');
K=(K+K')/2;A0=K(free,free);
B=d.C(free,[2:65 67:130]);
coeff=t*[Y(1,2:65) Y(2,2:65)]'/N;
clock=tic;F=decomposition(A0,'chol');factorSeconds=toc(clock);
clock=tic;R=full(B'*(F\B));responseSeconds=toc(clock);
Rskew=norm(R-R','fro')/max(norm(R,'fro'),realmin);R=(R+R')/2;
rootD=sqrt(coeff);T=(rootD*rootD').*R;T=(T+T')/2;
ev=eig(T);maximum=max(real(ev));minimum=min(real(ev));margin=1-maximum;
% A strict positive observed margin is stronger than a tiny tolerance pass.
% Floating-point output is a diagnostic, not an interval-arithmetic proof.
info=struct('sigma_f_S_m',sf,'is_background',isnan(sf),...
 'condition','A0 >= B diag(t Yref/Nn) B^T on free DOFs; B=C nonzero-mode columns',...
 'reference_modes',64,'number_of_nonzero_port_columns',128,...
 'width_m',W,'out_of_plane_thickness_m',t,'nonzero_mode_norm_m',N,...
 'coefficient_units','t*Yref/Nn in S/m^2; B in m; B D B^T in S',...
 'maximum_reference_trace_energy_eigenvalue',maximum,...
 'minimum_reference_trace_energy_eigenvalue',minimum,...
 'domination_spectral_margin',margin,...
 'strict_positive_observed_margin',margin>0,...
 'numerical_domination_pass',maximum<=1+1e-10,...
 'physical_matrix_relative_skew',skew,'port_response_relative_skew',Rskew,...
 'ndofs',numel(d.u0),'nfree',numel(free),'factor_seconds',factorSeconds,...
 'response_seconds',responseSeconds,'elapsed_seconds',toc(started),...
 'target_truth_used',false,...
 'scope','Sufficient numerical check for finite-matrix M0/M8 to M32/M64 bounds only; not full-COMSOL/continuum certification');
assert(info.numerical_domination_pass,...
 'M0 finite matrix does not dominate the positive complete reference modal energy');
end
