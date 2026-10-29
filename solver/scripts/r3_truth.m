function result=r3_truth(m,d,sf,folder,reuseFolder)
% Target validation ONLY. Caller must freeze independent-test prediction first.
if ~exist(folder,'dir'),mkdir(folder);end
stamp=fullfile(folder,'SUCCESS.json');
if isfile(stamp),result=jsondecode(fileread(stamp));return;end
m.param.set('r3bgscale','1');m.param.set('r3unit','0[S/m]');
m.param.set('r3dlow',sprintf('%.17g[S/m]',sf-d.c.sigma_low_S_m));
m.param.set('r3dhigh',sprintf('%.17g[S/m]',sf-d.c.sigma_high_S_m));
ts=tic;reused=false;
if sf==10 && ~isempty(reuseFolder) && isfile(fullfile(reuseFolder,'reused_seen_truth_10.mat'))
 old=load(fullfile(reuseFolder,'reused_seen_truth_10.mat'));
 u=old.u10;G=old.oldResult.sigma_eff_S_m;reused=true;solveSeconds=0;
 source=old.sourceFolder;
 % Postprocess the original already-seen solved MPH; no solver rerun.
 m=mphload(fullfile(source,[d.c.id '.mph']),'Reuse10');
 reuseCleanup=onCleanup(@() com.comsol.model.util.ModelUtil.remove('Reuse10'));
else
 m.sol('sol1').runAll;solveSeconds=toc(ts);ut=mphgetu(m);xt=mphxmeshinfo(m);
 [u,map]=r3_map_dofs(ut,xt,d.xm); %#ok<ASGLU>
 G=mphint2(m,'reacf(V)',1,'selection',d.sel.positive_electrode,'method','summation');source='new_COMSOL_stationary_solve';
end
ts=tic;deltaK=(sf-d.c.sigma_low_S_m)*d.DL+(sf-d.c.sigma_high_S_m)*d.DH;K=d.K0+deltaK;
w=u-d.u0;delta=G-d.G0;U=d.u0'*(deltaK*d.u0);R=w'*(K*w);cross=d.u0'*(deltaK*u);
stable=w'*(d.K0*w)+u'*(deltaK*u);
power=mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',2,'intorder',8);
neg=mphint2(m,'reacf(V)',1,'selection',d.sel.ground_electrode,'method','summation');
result=struct('status','success','case_id',d.c.id,'sigma_f_S_m',sf,'G0_S',d.G0,'G_S',G,...
 'deltaG_S',delta,'relative_gain',delta/d.G0,'U_S',U,'R_S',R,'U_minus_R_S',U-R,...
 'cross_S',cross,'positive_energy_deltaG_S',stable,...
 'identity_relative_error',abs(U-R-delta)/max(abs(delta),1e-12),...
 'cross_relative_error',abs(cross-delta)/max(abs(delta),1e-12),...
 'matrix_energy_relative_error',abs(u'*(K*u)-G)/G,...
 'reaction_flux_relative_error',abs(G+neg)/G,'energy_relative_error',abs(power-G)/G,...
 'solve_seconds',solveSeconds,'reused_seen_sigma10',reused,'source',source);
assert(result.identity_relative_error<1e-4 && result.cross_relative_error<1e-4,'Gain identity accuracy');
assert(result.energy_relative_error<1e-6,'Physical power error');
cc=d.c;cc.sigma_frac_S_m=sf;local=r3_local_metrics(m,cc,d.sel,folder);
result.local=local;result.postprocess_seconds=toc(ts);
fields=mpheval(m,{'V','ec.Ex','ec.Ey','ec.Jx','ec.Jy','ec.Jx*ec.Ex+ec.Jy*ec.Ey'},'edim',2,'refine',1,'smooth','none');
save(fullfile(folder,'fields.mat'),'fields','u','w','-v7.3');
params=d.c;params.sigma_frac_S_m=sf;params.passive_ellipse=false;
r3_json(fullfile(folder,'parameters.json'),params);
if ~reused,mphsave(m,fullfile(folder,'solution.mph'));else,r3_json(fullfile(folder,'reused_model.json'),struct('original_mph',fullfile(source,[d.c.id '.mph'])));end
r3_json(stamp,result);
fprintf('TRUTH %s sf %.3g G %.12g gain %.12g identity %.3g\n',d.c.id,sf,G,delta,result.identity_relative_error);
end
