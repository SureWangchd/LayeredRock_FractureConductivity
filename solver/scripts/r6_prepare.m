function r6_prepare(root,fine)
% Assemble COMSOL M0 physical operators. No target truth is read.
import com.comsol.model.util.*
out=fullfile(root,'runs','round6_applicability');level='base';if fine,level='fine';end
design=jsondecode(fileread(fullfile(out,'design.json')));
for ig=1:numel(design.source_groups)
 group=design.source_groups(ig);role={group.role};
 folder=fullfile(out,'reference',[group.source_key '_' level]);if ~isfolder(folder),mkdir(folder);end
 cache=fullfile(folder,'reference.mat');if isfile(cache),continue;end
 source=fullfile(out,'reference',[group.source_key '_base'],'source_background');
 idx=find(strcmp({design.cases.id},group.reference_case_id));assert(numel(idx)==1);
 c=design.cases(idx);c.passive_ellipse=true;
 if ~isfolder(source),mkdir(source);end
 if ~isfile(fullfile(source,'SUCCESS.json'))
  r3_json(fullfile(source,'parameters.json'),c);
  timer=tic;r=solve_round2_case(c,source);r.status='success';r.elapsed_seconds=toc(timer);
  r3_json(fullfile(source,'SUCCESS.json'),r);ModelUtil.remove('Rock2D');
 end
 sel=jsondecode(fileread(fullfile(source,'selections.json')));
 m=mphload(fullfile(source,[c.id '.mph']),'R5Reference');
 if fine
  mesh=m.component('comp1').mesh('mesh1');mesh.feature('size').set('hmax','0.8');
  mesh.feature('tri').feature('fs').set('hmax','0.1');mesh.run;
 end
 co=m.component('comp1');tags=cell(co.material.tags);for k=1:numel(tags),co.material.remove(tags{k});end
 ids=double(sel.all_domains(:)');lay=double(sel.domain_background_layer(:)');ph=double(sel.domain_background_material(:)');
 edges=[0 reshape(c.z_interfaces_mm,1,[]) 120];
 lower=edges(lay+1)<=30;upper=edges(lay)>=90;core=~lower&~upper;inside=ismember(ids,sel.inclusion_geometry_domains);
 assert(all(core(inside)),'Ellipse must be within core');
 m.param.set('r5base','1');m.param.set('r5dlow','0[S/m]');m.param.set('r5dhigh','0[S/m]');
 for p=0:1
  sig='sigma_low';delta='r5dlow';if p,sig='sigma_high';delta='r5dhigh';end
  for f=0:1
   selection=ids(core&ph==p&inside==f);if isempty(selection),continue;end
   expr=['r5base*' sig];if f,expr=[expr '+' delta];end
   add_material(co,sprintf('r5p%d_f%d',p,f),selection,expr,expr);
  end
 end
 sx=[.034 .076];sz=[.03/(.022/.01+.008/.1) .03/(.008/.01+.022/.1)];
 add_material(co,'r5lower',ids(lower),sprintf('r5base*%.17g[S/m]',sx(1)),sprintf('r5base*%.17g[S/m]',sz(1)));
 add_material(co,'r5upper',ids(upper),sprintf('r5base*%.17g[S/m]',sx(2)),sprintf('r5base*%.17g[S/m]',sz(2)));
 v=m.sol('sol1').feature('v1');v.set('scalemethod','none');
 ts=tic;m.sol('sol1').runAll;backgroundSeconds=toc(ts);
 u0=mphgetu(m);u0=u0(:);xm=mphxmeshinfo(m);
 G0=mphint2(m,'reacf(V)',1,'selection',sel.positive_electrode,'method','summation');
 power=mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',2,'intorder',8);
 mphsave(m,fullfile(folder,'background_m0.mph'));
 S=mphmatrix(m,'sol1','out',{'K','N','uscale'},'rowscale','off','initmethod','sol','initsol','zero','extractafter','v1');
 K0=S.K;assert(all(abs(S.uscale(:)-1)<1e-12),'Unexpected scaling');
 fixed=find(any(S.N,1));free=setdiff((1:numel(u0))',fixed(:));
 xy=double(xm.dofs.coords);top=find(abs(xy(2,:)-.12)<1e-12);bottom=find(abs(xy(2,:))<1e-12);
 assert(isequal(sort(fixed(:)),sort([top(:);bottom(:)])),'Dirichlet DOFs do not match electrodes');
 assert(abs(G0/(1/(.5/.01+.5/.1))-1)<1e-6,'Strict background analytic check');
 assert(abs(u0'*(K0*u0)-G0)/G0<1e-7 && abs(power-G0)/G0<1e-6,'Reference energy check');
 r3_json(fullfile(folder,'parameters.json'),c);r3_json(fullfile(folder,'selections.json'),sel);
 [ms,md]=mphmeshstats(m,'mesh1');save(fullfile(folder,'mesh.mat'),'ms','md','-v7.3');
 m.param.set('r5base','0');m.param.set('r5dlow','1[S/m]');DL=r3_matrix(m);
 m.param.set('r5dlow','0[S/m]');m.param.set('r5dhigh','1[S/m]');DH=r3_matrix(m);
 m.param.set('r5dhigh','0[S/m]');m.param.set('r5base','1');
 [C,traceinfo]=r5_trace_matrix(xm,64);
 [C40,~]=r6_trace40(xm,64);
 traceinfo.gauss20_vs40_relative_fro=norm(C-C40,'fro')/norm(C40,'fro');
 traceinfo.gauss20_vs40_max_abs_m=full(max(abs(C(:)-C40(:))));
 assert(traceinfo.gauss20_vs40_relative_fro<1e-10,'Modal projection quadrature unresolved');clear C40;
 traceinfo.zero_mode_background_lower_V=(C(:,1)'*u0)/.12;
 traceinfo.zero_mode_background_upper_V=(C(:,66)'*u0)/.12;
 traceinfo.nonzero_background_max_abs_V=max(abs((C(:,setdiff(1:130,[1 66]))'*u0)/.06));
 assert(traceinfo.nonzero_background_max_abs_V<1e-8,'Background has nonuniform trace');
 summary=struct('source_key',group.source_key,'full_length_mm',2*c.a_mm,'role',role{1},'mesh_level',level,'COMSOL_G0_S',G0,'analytic_G0_S',1/(.5/.01+.5/.1),...
  'ndofs',numel(u0),'background_seconds',backgroundSeconds,'trace',traceinfo,'sigma_x_eq_S_m',sx,'sigma_z_eq_S_m',sz,...
  'target_truth_loaded',false,'source_background',source,'method','COMSOL FE reference volume plus analytic finite modal stiffness correction');
 save(cache,'K0','DL','DH','C','u0','xm','fixed','free','top','bottom','G0','c','sel','sx','sz','summary','-v7.3');
 r3_json(fullfile(folder,'reference_summary.json'),summary);
 fprintf('R6_PREPARED %s %s dofs %d G0 %.12g trace_rank %d\n',group.source_key,level,numel(u0),G0,traceinfo.rank);
 ModelUtil.remove('R5Reference');
end
end

function add_material(co,tag,ids,sx,sz)
co.material.create(tag,'Common');ma=co.material(tag);ma.selection.set(ids);
ma.propertyGroup('def').set('electricconductivity',{sx '0' '0' '0' sz '0' '0' '0' sx});
ma.propertyGroup('def').set('relpermittivity',{'1'});
end
