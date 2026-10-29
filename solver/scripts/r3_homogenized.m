function r3_homogenized(root)
% Common directional-matching baseline for the five FLAT development backgrounds.
% Full homogeneous-anisotropic reference solves, not explicit-target predictions.
out=fullfile(root,'runs','round3_predict','homogenized_reference');if ~exist(out,'dir'),mkdir(out);end
source=fullfile(root,'runs','2026-10-21','homogenized_verified_frac');
sel=jsondecode(fileread(fullfile(source,'selections.json')));
m=mphload(fullfile(source,'homogenized_verified_frac.mph'),'HomogeneousR3');
m.param.set('sigma_xx','0.055[S/m]');m.param.set('sigma_zz',sprintf('%.17g[S/m]',2/(1/.01+1/.1)));
co=m.component('comp1');co.material('frac').selection.set([]);co.material('low').selection.all;
ts=tic;m.sol('sol1').runAll;bgTime=toc(ts);
G0=mphint2(m,'reacf(V)',1,'selection',sel.positive_electrode,'method','summation');
mphsave(m,fullfile(out,'background.mph'));
co.material('low').selection.set(sel.low_domains);co.material('frac').selection.set(sel.fracture_domains);
rows=struct([]);
for sf=[.1 .3 1 3 10 30 100]
 folder=fullfile(out,['sf_' strrep(num2str(sf),'.','p')]);if ~exist(folder,'dir'),mkdir(folder);end
 if isfile(fullfile(folder,'SUCCESS.json'))
  r=jsondecode(fileread(fullfile(folder,'SUCCESS.json')));
 else
  m.param.set('sigma_frac',sprintf('%.17g[S/m]',sf));ts=tic;m.sol('sol1').runAll;st=toc(ts);
  G=mphint2(m,'reacf(V)',1,'selection',sel.positive_electrode,'method','summation');
  power=mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',2,'intorder',8);
  r=struct('sigma_f_S_m',sf,'G0_S',G0,'G_S',G,'predicted_deltaG_S',G-G0,...
   'sigma_xx_S_m',.055,'sigma_zz_S_m',2/(1/.01+1/.1),'solve_seconds',st,...
   'background_solve_seconds_shared',bgTime,'energy_relative_error',abs(G-power)/G,...
   'scope','Common homogenized reference for flat development cases; no target-specific calibration; not a full tensor match');
  mphsave(m,fullfile(folder,'solution.mph'));r3_json(fullfile(folder,'SUCCESS.json'),r);
 end
 if isempty(rows),rows=r;else,rows(end+1)=r;end %#ok<AGROW>
end
writetable(struct2table(rows),fullfile(out,'homogenized_predictions.csv'));
com.comsol.model.util.ModelUtil.remove('HomogeneousR3');
end
