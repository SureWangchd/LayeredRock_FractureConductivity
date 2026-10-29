function [m,d]=r3_prepare(folder,sourceFolder,reuseSeenTruth)
% Build background-only matrices and responses. No target solve in this entry.
import com.comsol.model.util.*
if ~exist(folder,'dir'),mkdir(folder);end
cache=fullfile(folder,'background_responses.mat');
if isfile(cache)
 d=load(cache);m=mphload(fullfile(folder,'background.mph'),'R3');return;
end
started=tic;
c=jsondecode(fileread(fullfile(sourceFolder,'parameters.json')));
sel=jsondecode(fileread(fullfile(sourceFolder,'selections.json')));
m=mphload(fullfile(sourceFolder,[c.id '.mph']),'R3');
if reuseSeenTruth
 oldU=mphgetu(m);oldX=mphxmeshinfo(m);oldResult=jsondecode(fileread(fullfile(sourceFolder,'SUCCESS.json')));
end
r3_setup_materials(m,sel);
ts=tic;m.sol('sol1').runAll;cost.background_solve_seconds=toc(ts);
u0=mphgetu(m);u0=u0(:);xm=mphxmeshinfo(m);
G0=mphint2(m,'reacf(V)',1,'selection',sel.positive_electrode,'method','summation');
ground=mphint2(m,'reacf(V)',1,'selection',sel.ground_electrode,'method','summation');
power=mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',2,'intorder',8);
assert(abs(G0-power)/G0<1e-6 && abs(G0+ground)/G0<1e-6,'Background balance');
if reuseSeenTruth
 [u10,map]=r3_map_dofs(oldU,oldX,xm);
 save(fullfile(folder,'reused_seen_truth_10.mat'),'u10','oldResult','map','sourceFolder','-v7.3');
end
r3_json(fullfile(folder,'parameters.json'),c);r3_json(fullfile(folder,'selections.json'),sel);
mphsave(m,fullfile(folder,'background.mph'));
bgfields=mpheval(m,{'V','ec.Ex','ec.Ey'},'edim',2,'refine',1,'smooth','none');
save(fullfile(folder,'background_fields.mat'),'bgfields','-v7.3');
capV=zeros(1,2);caps={sel.capminus_domains,sel.capplus_domains};
for j=1:2,capV(j)=mphint2(m,'V',2,'selection',caps{j},'intorder',8)/mphint2(m,'1',2,'selection',caps{j},'intorder',8);end
ts=tic;S=mphmatrix(m,'sol1','out',{'K','N','uscale'},'rowscale','off',...
 'initmethod','sol','initsol','zero','extractafter','v1');
K0=S.K;assert(all(abs(S.uscale(:)-1)<1e-12),'Variable scale must be one');
assert(all(sum(spones(S.N),2)<=1),'Expected only direct Dirichlet constraints');
fixed=find(any(S.N,1));free=setdiff(1:numel(u0),fixed);free=free(:);
assert(norm(K0-K0','fro')/norm(K0,'fro')<1e-10,'Physical K symmetry');
assert(abs(u0'*(K0*u0)-G0)/G0<1e-6,'K energy units/scaling');
coord=xm.dofs.coords';
assert(size(coord,1)==numel(u0),'Scalar V DOFs required');
if max(coord(:))>1,coord=coord*.001;end
m.param.set('r3bgscale','0');m.param.set('r3dlow','1[S/m]');DL=r3_matrix(m);
m.param.set('r3dlow','0[S/m]');m.param.set('r3dhigh','1[S/m]');DH=r3_matrix(m);
m.param.set('r3dhigh','0[S/m]');m.param.set('r3unit','1[S/m]');
cost.operator_assembly_seconds=toc(ts);
% Nested axial blocks with original phase separated; integration at COMSOL quadrature.
nseg=16;nblocks=nseg*2;rhs=zeros(numel(free),3*nblocks);blockinfo=zeros(nblocks,4);
ts=tic;sumD=sparse(size(K0,1),size(K0,2));
for j=1:nseg
 lo=-c.a_mm+2*c.a_mm*(j-1)/nseg;hi=-c.a_mm+2*c.a_mm*j/nseg;
 if j==1,lo=lo-1e-8;end;if j==nseg,hi=hi+1e-8;end
 m.param.set('r3lo',num2str(lo,17));m.param.set('r3hi',num2str(hi,17));
 for phase=0:1
  ib=2*(j-1)+phase+1;m.param.set('r3phase',num2str(phase));B=r3_matrix(m);
  sumD=sumD+B;rhs(:,3*ib-2:3*ib)=-B(free,:)*[u0 coord(:,1) coord(:,2)];
  blockinfo(ib,:)=[j phase lo hi];
 end
 fprintf('BACKGROUND_BLOCKS %s %d/%d\n',c.id,j,nseg);
end
cost.block_assembly_seconds=toc(ts);
assert(norm(sumD-DL-DH,'fro')/norm(DL+DH,'fro')<1e-9,'Segment partition does not close');clear sumD B S;
ts=tic;F=decomposition((K0(free,free)+K0(free,free)')/2,'chol');cost.background_factor_seconds=toc(ts);
ts=tic;responses=F\rhs;cost.auxiliary_backsolve_seconds=toc(ts);
m.param.set('r3bgscale','1');m.param.set('r3unit','0[S/m]');
cost.total_background_preparation_seconds=toc(started);
provenance=struct('target_solve_count_in_prepare',0,'background_solve_count',1,...
 'matrix_extraction','mphmatrix initmethod=sol initsol=zero extractafter=v1 rowscale=off',...
 'source_folder',sourceFolder,'reuse_seen_truth',reuseSeenTruth,...
 'basis_inputs','Background K0, background u0, axial segment/original phase masks; x/z coordinate gradient responses',...
 'units','Matrices include t=1 m; u in V; energies and deltaG numerically W and S at DV=1 V');
save(cache,'K0','DL','DH','u0','free','fixed','coord','xm','rhs','responses','blockinfo','nseg','c','sel','G0','capV','cost','provenance','-v7.3');
r3_json(fullfile(folder,'background_summary.json'),struct('G0_S',G0,'cap_dV_V',diff(capV),'dofs',numel(u0),'free_dofs',numel(free),'cost',cost,'provenance',provenance));
d=load(cache);
end
