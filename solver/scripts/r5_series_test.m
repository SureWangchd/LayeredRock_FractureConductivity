function result=r5_series_test(root)
% COMSOL and exported-matrix affine series test before target prediction.
out=fullfile(root,'runs','round5_modal_closure','validation');if ~isfolder(out),mkdir(out);end
stamp=fullfile(out,'series_test.json');if isfile(stamp),result=jsondecode(fileread(stamp));return;end
source=fullfile(root,'runs','round4_crossover','periphery','r4_strict_A_swap_far');
sel=jsondecode(fileread(fullfile(source,'selections.json')));c=jsondecode(fileread(fullfile(source,'parameters.json')));
m=mphload(fullfile(source,'background.mph'),'R5Series');co=m.component('comp1');
tags=cell(co.material.tags);for k=1:numel(tags),co.material.remove(tags{k});end
ids=double(sel.all_domains(:)');lay=double(sel.domain_background_layer(:)');edges=[0 c.z_interfaces_mm(:)' 120];
sets={ids(edges(lay+1)<=30),ids(edges(lay)>=30&edges(lay+1)<=90),ids(edges(lay)>=90)};
conductivity=[.01 .02 .04];
for j=1:3
 tag=sprintf('series%d',j);co.material.create(tag,'Common');ma=co.material(tag);ma.selection.set(sets{j});
 s=sprintf('%.17g[S/m]',conductivity(j));ma.propertyGroup('def').set('electricconductivity',{s '0' '0' '0' s '0' '0' '0' s});
 ma.propertyGroup('def').set('relpermittivity',{'1'});
end
m.sol('sol1').runAll;u=mphgetu(m);xm=mphxmeshinfo(m);
pos=mphint2(m,'reacf(V)',1,'selection',sel.positive_electrode,'method','summation');
neg=mphint2(m,'reacf(V)',1,'selection',sel.ground_electrode,'method','summation');
P=mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',2,'intorder',8);mphsave(m,fullfile(out,'uniform_core_series.mph'));
S=mphmatrix(m,'sol1','out',{'K','N','uscale'},'rowscale','off','initmethod','sol','initsol','zero','extractafter','v1');
fixed=find(any(S.N,1));free=setdiff((1:numel(u))',fixed(:));
g=u(fixed);b=-S.K(free,fixed)*g;constant=g'*(S.K(fixed,fixed)*g);v=u;v(free)=S.K(free,free)\b;
expected=.12/(.03/.01+.06/.02+.03/.04);[C,trace]=r5_trace_matrix(xm,64);
matrixG=constant-b'*v(free);modes=reshape(C'*v,65,2)';modes(:,1)=modes(:,1)/.12;modes(:,2:end)=modes(:,2:end)/.06;
assert(abs(pos/expected-1)<1e-6 && abs(P/expected-1)<1e-6,'COMSOL series analytic mismatch');
assert(abs(matrixG/expected-1)<1e-6 && max(abs(u-v))<1e-7,'Affine matrix solve mismatch');
assert(abs(pos+neg)/pos<1e-6 && max(abs(modes(:,2:end)),[],'all')<1e-8,'Series flux or modal mismatch');
result=struct('status','passed','conductivities_lower_core_upper_S_m',conductivity,...
 'COMSOL_G_S',pos,'analytic_G_S',expected,'matrix_G_S',matrixG,'volume_power_W',P,...
 'affine_c_minus_bu_S',matrixG,'reaction_flux_relative_error',abs(pos+neg)/pos,...
 'matrix_COMSOL_max_voltage_difference_V',max(abs(u-v)),'predicted_port_zero_V',modes(:,1),...
 'maximum_nonzero_mode_V',max(abs(modes(:,2:end)),[],'all'),'trace',trace);
r3_json(stamp,result);com.comsol.model.util.ModelUtil.remove('R5Series');disp(result);
end
