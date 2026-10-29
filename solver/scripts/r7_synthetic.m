function r7_synthetic(root,action)
% Finite 3D COMSOL forward experiments. Previous solvers/data are untouched.
out=fullfile(root,'runs','round7_3d_ct');design=jsondecode(fileread(fullfile(out,'design.json')));
if strcmp(action,'pilot')
 cases=design.cases(1);branch='synthetic3d';sfs=.18;
elseif strcmp(action,'main')
 cases=design.cases;branch='synthetic3d';sfs=[.18 10];
elseif strcmp(action,'m0')
 cases=design.m0_cases;branch='m0';sfs=[.18 10];
elseif strcmp(action,'refine')
 refine_all(out,design);return;
else,error('Unknown R7 action');end
for ic=1:numel(cases)
 c=cases(ic);folder=fullfile(out,branch,c.id);if ~isfolder(folder),mkdir(folder);end
 pending=arrayfun(@(sf)~isfile(fullfile(folder,sf_token(sf),'SUCCESS.json')),sfs);
 if ~any(pending),fprintf('R7_REUSE %s\n',c.id);continue;end
 tag='Rock3D';cleaner=onCleanup(@()remove_model(tag));
 if isfile(fullfile(folder,'background_z','SUCCESS.json'))
  model=mphload(fullfile(folder,'background_z','solution.mph'),tag);
  sel=jsondecode(fileread(fullfile(folder,'selections.json')));
  bg=jsondecode(fileread(fullfile(folder,'background_z','SUCCESS.json')));
 else
  [model,sel]=build_model(c,folder,tag);
  bg=solve_state(model,c,sel,'z',NaN,fullfile(folder,'background_z'),NaN);
  geometry=measure_geometry(model,c,sel,folder);r3_json(fullfile(folder,'actual_geometry.json'),geometry);
 end
 if ~isfile(fullfile(folder,'background_x','SUCCESS.json'))
  solve_state(model,c,sel,'x',NaN,fullfile(folder,'background_x'),NaN);
 end
 if ic==1 && strcmp(branch,'synthetic3d') && strcmp(c.role,'A') && strcmp(c.layout,'P11')
  if ~isfile(fullfile(folder,'background_y','SUCCESS.json'))
   solve_state(model,c,sel,'y',NaN,fullfile(folder,'background_y'),NaN);
  end
 end
 for sf=sfs
  target=fullfile(folder,sf_token(sf));if isfile(fullfile(target,'SUCCESS.json')),continue;end
  solve_state(model,c,sel,'z',sf,target,bg.G_S);
 end
 clear model cleaner;
end
end

function [m,sel]=build_model(c,folder,tag)
import com.comsol.model.util.*
timer=tic;m=ModelUtil.create(tag);m.label([c.id '.mph']);ModelUtil.showProgress(fullfile(folder,'solver_progress.log'));
m.param.set('sigma_low','0.01[S/m]');m.param.set('sigma_high','0.1[S/m]');m.param.set('sigma_frac','0.18[S/m]');m.param.set('frac_on','0');m.param.set('DV','1[V]');
m.component.create('comp1',true);co=m.component('comp1');co.geom.create('geom1',3);g=co.geom('geom1');g.lengthUnit('mm');
isM0=strcmp(c.kind,'synthetic3d_M0');
% Merge only equal-material adjacent blocks; never alter actual material map.
if isM0
 edges=[0 30 reshape(c.physical_z_interfaces_mm,1,[]) 90 120];edges=unique(edges);edges=edges(edges==0|edges==120|(edges>=30&edges<=90));
else
 edges=[0 reshape(c.physical_z_interfaces_mm,1,[]) 120];
end
for j=1:numel(edges)-1
 name=sprintf('layer%d',j);g.create(name,'Block');g.feature(name).set('size',[120 120 edges(j+1)-edges(j)]);g.feature(name).set('pos',[0 0 edges(j)]);
end
g.create('ell','Ellipsoid');el=g.feature('ell');el.set('semiaxes',[c.a_mm c.b_mm c.c_mm]);el.set('pos',[60 60 60]);
el.set('axistype','cartesian');el.set('axis',reshape(c.basis_n,1,[]));el.set('rot',0);el.set('selresult','on');g.run;
ell=reshape(double(co.selection('geom1_ell_dom').entities(3)),1,[]);assert(~isempty(ell));
co.mesh.create('mesh1');me=co.mesh('mesh1');me.create('tet','FreeTet');
me.feature('size').set('custom','on');me.feature('size').set('hmax',5);me.feature('size').set('hmin',.1);me.feature('size').set('hgrad',1.4);me.feature('size').set('hcurve',.25);
me.feature('tet').create('fs','Size');fs=me.feature('tet').feature('fs');fs.selection.geom('geom1',3);fs.selection.set(ell);
fs.set('custom','on');fs.set('hmax',.75);fs.set('hmaxactive',true);fs.set('hmin',.05);fs.set('hminactive',true);fs.set('hgrad',1.3);fs.set('hgradactive',true);
meshTimer=tic;me.run;meshSeconds=toc(meshTimer);
[ms,md]=mphmeshstats(m,'mesh1','qualitymeasure','skewness');save(fullfile(folder,'mesh.mat'),'ms','md','-v7.3');r3_json(fullfile(folder,'mesh_stats.json'),ms);
r7_thinmesh_check(md,ms,c,ell,folder);
it=find(strcmp(ms.types,'tet'));tet=double(md.elem{it})+1;ent=double(md.elementity{it});verts=md.vertex;
scale=1;if max(verts(:))<1,scale=1000;end
centers=(verts(:,tet(1,:))+verts(:,tet(2,:))+verts(:,tet(3,:))+verts(:,tet(4,:)))/4*scale;
ids=unique(ent);phase=zeros(size(ids));layer=phase;zones=phase;
for j=1:numel(ids)
 k=find(ent==ids(j),1);z=centers(3,k);layer(j)=1+sum(z>c.z_interfaces_mm(:));phase(j)=c.layer_materials(layer(j));
 if isM0 && z<30,zones(j)=1;elseif isM0&&z>90,zones(j)=2;end
end
for p=0:1
 sig='sigma_low';if p==1,sig='sigma_high';end
 for f=0:1
  selected=ids(phase==p & ismember(ids,ell)==f & zones==0);if isempty(selected),continue;end
  expression=sig;if f,expression=[sig '+frac_on*(sigma_frac-' sig ')'];end
  addmat(co,sprintf('p%d_f%d',p,f),selected,expression,expression,expression);
 end
end
if isM0
 sx=[.034 .076];sz=[.03/(.022/.01+.008/.1) .03/(.008/.01+.022/.1)];
 for j=1:2,addmat(co,sprintf('outer%d',j),ids(zones==j),sprintf('%.17g[S/m]',sx(j)),sprintf('%.17g[S/m]',sx(j)),sprintf('%.17g[S/m]',sz(j)));end
end
tol=1e-5;faces=cell(3,2);
for ax=1:3
 for side=1:2
  bounds=repmat([-tol 120+tol],3,1);q=(side-1)*120;bounds(ax,:)=[q-tol q+tol];faces{ax,side}=mphselectbox(m,'geom1',bounds,'boundary');assert(~isempty(faces{ax,side}));
 end
end
co.physics.create('ec','ConductiveMedia','geom1');ec=co.physics('ec');
ec.create('pot1','ElectricPotential',2);ec.feature('pot1').selection.set(faces{3,2});ec.feature('pot1').set('V0','DV');
ec.create('gnd1','Ground',2);ec.feature('gnd1').selection.set(faces{3,1});
m.study.create('std1');m.study('std1').create('stat','Stationary');
fracTet=tet(:,ismember(ent,ell));fv=unique(fracTet(:));bbox=[min(verts(:,fv),[],2),max(verts(:,fv),[],2)]*scale;
assert(all(bbox(:,1)>0)&all(bbox(:,2)<120),'Finite ellipsoid reached an external boundary');
sel=struct('all_domains',ids,'inclusion_domains',ell,'background_phase',phase,'background_layer',layer,'zone',zones,...
 'low_domains',ids(phase==0),'high_domains',ids(phase==1),'xminus',faces{1,1},'xplus',faces{1,2},'yminus',faces{2,1},'yplus',faces{2,2},'zminus',faces{3,1},'zplus',faces{3,2});
sel.tetrahedra=size(tet,2);sel.min_mesh_quality=ms.minquality;sel.mean_mesh_quality=ms.meanquality;sel.mesh_seconds=meshSeconds;sel.geometry_and_mesh_seconds=toc(timer);
sel.ellipsoid_mesh_vertex_bbox_mm=bbox;sel.global_hmax_mm=5;sel.ellipse_hmax_mm=.75;sel.geometry_unit='mm';
r3_json(fullfile(folder,'parameters.json'),c);r3_json(fullfile(folder,'selections.json'),sel);
mphsave(m,fullfile(folder,'presolve.mph'));
fprintf('R7_MESH %s tets %d quality %.4g build %.2fs\n',c.id,sel.tetrahedra,sel.min_mesh_quality,toc(timer));
end

function r=solve_state(m,c,sel,load,sf,folder,G0)
if ~isfolder(folder),mkdir(folder);end
ec=m.component('comp1').physics('ec');pos=sel.([load 'plus']);neg=sel.([load 'minus']);
ec.feature('pot1').selection.set(pos);ec.feature('gnd1').selection.set(neg);
if isnan(sf),m.param.set('frac_on','0');else,m.param.set('frac_on','1');m.param.set('sigma_frac',sprintf('%.15g[S/m]',sf));end
m.study('std1').feature('stat').set('usestol',true);m.study('std1').feature('stat').set('stol',1e-8);
timer=tic;m.study('std1').run;solveSeconds=toc(timer);
I=mphint2(m,'reacf(V)',2,'selection',pos,'method','summation');Ig=mphint2(m,'reacf(V)',2,'selection',neg,'method','summation');
q=mphint2(m,'ec.Jx*nx+ec.Jy*ny+ec.Jz*nz',2,'selection',pos,'intorder',6);qg=mphint2(m,'ec.Jx*nx+ec.Jy*ny+ec.Jz*nz',2,'selection',neg,'intorder',6);
power=mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey+ec.Jz*ec.Ez',3,'intorder',6);
A=mphint2(m,'1',2,'selection',pos,'unit','m^2','intorder',6);G=I;assert(G>0&&abs(A/.0144-1)<1e-7);
flux=abs(I+Ig)/G;energy=abs(I-power)/G;assert(max(flux,energy)<1e-5,'3D reaction conservation failed');
if isnan(sf)
 expected=.0066;if strcmp(load,'z'),expected=.12/(.5/.01+.5/.1);end
 assert(abs(G/expected-1)<1e-5,'3D passive background analytic test failed');
end
jynet=mphint2(m,'ec.Jy',3,'intorder',6);jyabs=mphint2(m,'abs(ec.Jy)',3,'intorder',6);jy2=mphint2(m,'ec.Jy^2',3,'intorder',6);j2=mphint2(m,'ec.Jx^2+ec.Jy^2+ec.Jz^2',3,'intorder',6);
u='((x-60[mm])*0.5+(z-60[mm])*sqrt(3)/2)';cap={'minus','plus'};capmetrics=struct();
for k=1:2
 op='<';threshold='-27[mm]';if k==2,op='>';threshold='27[mm]';end
 indicator=['(' u op threshold ')'];vol=mphint2(m,indicator,3,'selection',sel.inclusion_domains,'intorder',10,'unit','mm^3');
 v=mphint2(m,['V*' indicator],3,'selection',sel.inclusion_domains,'intorder',10)/max(vol*1e-9,realmin);
 capmetrics.([cap{k} '_volume_mm3'])=vol;capmetrics.([cap{k} '_mean_V'])=v;
end
capmetrics.definition='Finite 3-mm axial caps, indicator quadrature order10 on curved FE domains; actual quadrature volume recorded, no point tips.';
r=struct('status','success','case_id',c.id,'role',c.role,'layout',c.layout,'kind',c.kind,'load',load,'sigma_f_S_m',sf,'is_background',isnan(sf),...
 'G_S',G,'G0_S',G0,'g_S',G-G0,'relative_gain',(G-G0)/G0,'sigma_app_S_m',G*.12/A,'electrode_area_m2',A,...
 'reaction_positive_A',I,'reaction_ground_A',Ig,'gradient_positive_out_A',q,'gradient_ground_out_A',qg,...
 'reaction_flux_relative_error',flux,'energy_relative_error',energy,'gradient_flux_relative_error',abs(q+qg)/max(abs([q qg])),...
 'volume_dissipation_W',power,'source_power_W',I,'Jy_net_volume_integral_A_m',jynet,'Jy_abs_volume_integral_A_m',jyabs,...
 'Jy_squared_volume_integral_A2_per_m',jy2,'all_J_squared_volume_integral_A2_per_m',j2,'Jy_squared_fraction',jy2/j2,...
 'fracture_dissipation_W',mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey+ec.Jz*ec.Ez',3,'selection',sel.inclusion_domains,'intorder',6),...
 'caps',capmetrics,'tetrahedra',sel.tetrahedra,'min_mesh_quality',sel.min_mesh_quality,'solve_seconds',solveSeconds,...
 'comsol_version','6.3.0.290','dimension',3,'study_stationary_relative_tolerance',1e-8,'units_note','True 3D: G in S; sigma=G*Lz/(Lx*Ly). No out-of-plane thickness.',...
 'completed_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
mphsave(m,fullfile(folder,'solution.mph'));save_slices(m,folder);
c.current_load=load;c.current_fracture_conductivity_S_m=sf;c.passive_ellipsoid=isnan(sf);r3_json(fullfile(folder,'parameters.json'),c);
r3_json(fullfile(folder,'SUCCESS.json'),r);
fprintf('R7_SOLVED %s %s sf %.8g G %.12g gain %.12g flux %.2g energy %.2g Jy2 %.5g time %.2f\n',c.id,load,sf,G,G-G0,flux,energy,jy2/j2,solveSeconds);
end

function geo=measure_geometry(m,c,sel,folder)
ids=sel.all_domains;vol=zeros(size(ids));
for j=1:numel(ids),vol(j)=mphint2(m,'1',3,'selection',ids(j),'intorder',8,'unit','mm^3');end
isell=ismember(ids,sel.inclusion_domains);p=sel.background_phase;
geo=struct('case_id',c.id,'role',c.role,'layout',c.layout,'kind',c.kind,'sample_volume_mm3',sum(vol),...
 'ellipsoid_volume_mm3',sum(vol(isell)),'ellipsoid_analytic_volume_mm3',4*pi*c.a_mm*c.b_mm*c.c_mm/3,...
 'original_low_volume_mm3',sum(vol(p==0)),'original_high_volume_mm3',sum(vol(p==1)),...
 'replaced_low_volume_mm3',sum(vol(isell&p==0)),'replaced_high_volume_mm3',sum(vol(isell&p==1)),...
 'remaining_low_volume_mm3',sum(vol(~isell&p==0)),'remaining_high_volume_mm3',sum(vol(~isell&p==1)),...
 'ellipse_mesh_bbox_mm',sel.ellipsoid_mesh_vertex_bbox_mm,'sample_dimensions_mm',[120 120 120],...
 'volume_definition','Order8 COMSOL physical-volume integration; source-phase bookkeeping in M0 is not actual homogeneous-zone material amount');
geo.ellipsoid_relative_volume_error=geo.ellipsoid_volume_mm3/geo.ellipsoid_analytic_volume_mm3-1;
assert(abs(geo.sample_volume_mm3/120^3-1)<1e-5);assert(abs(geo.ellipsoid_relative_volume_error)<.005,'Curved ellipsoid volume poorly resolved');
T=table(ids(:),p(:),sel.background_layer(:),isell(:),vol(:),'VariableNames',{'domain_id','original_phase_sample','original_layer_sample','is_ellipsoid','volume_mm3'});writetable(T,fullfile(folder,'domain_volumes.csv'));
if contains(lower(c.kind),'m0')
 geo.original_low_volume_mm3=NaN;geo.original_high_volume_mm3=NaN;geo.remaining_low_volume_mm3=NaN;geo.remaining_high_volume_mm3=NaN;
 geo.uniform_lower_volume_mm3=sum(vol(sel.zone==1));geo.uniform_upper_volume_mm3=sum(vol(sel.zone==2));
end
end

function save_slices(m,folder)
[x,z]=meshgrid(linspace(0,120,81));all=struct([]);expressions={'V','ec.Ex','ec.Ey','ec.Ez','ec.Jx','ec.Jy','ec.Jz'};
for k=1:3
 yy=[45 60 75];coords=[x(:)';yy(k)*ones(1,numel(x));z(:)'];vals=cell(1,7);
 [vals{:}]=mphinterp(m,expressions,'coord',coords,'coorderr','off');
 all(k).y_mm=yy(k);all(k).coords_mm=coords;all(k).expressions=expressions;all(k).values=vertcat(vals{:});
 T=array2table([coords' vertcat(vals{:})'],'VariableNames',{'x_mm','y_mm','z_mm','V_V','Ex_V_m','Ey_V_m','Ez_V_m','Jx_A_m2','Jy_A_m2','Jz_A_m2'});
 writetable(T,fullfile(folder,sprintf('slice_y%dmm.csv',yy(k))));
end
save(fullfile(folder,'slices.mat'),'all','-v7.3');
end

function refine_all(out,design)
request=jsondecode(fileread(fullfile(out,'refinement_targets.json')));assert(numel(request.targets)<=4);
for j=1:numel(request.targets)
 t=request.targets(j);c=design.cases(strcmp({design.cases.id},t.case_id));source=fullfile(out,'synthetic3d',c.id);folder=fullfile(out,'precision',[c.id '_' sf_token(t.sigma_f_S_m)]);
 if isfile(fullfile(folder,'target','SUCCESS.json')),continue;end
 if ~isfolder(folder),mkdir(folder);end
 assert(isfile(fullfile(source,sf_token(t.sigma_f_S_m),'SUCCESS.json')));
 com.comsol.model.util.ModelUtil.showProgress(fullfile(folder,'solver_progress.log'));
 tag='Rock3DFine';cleanup=onCleanup(@()remove_model(tag));m=mphload(fullfile(source,'background_z','solution.mph'),tag);
 sel=jsondecode(fileread(fullfile(source,'selections.json')));me=m.component('comp1').mesh('mesh1');
 me.feature('size').set('hmax',3.5);me.feature('tet').feature('fs').set('hmax',.5);timer=tic;me.run;
 [ms,md]=mphmeshstats(m,'mesh1','qualitymeasure','skewness');save(fullfile(folder,'mesh.mat'),'ms','md','-v7.3');r3_json(fullfile(folder,'mesh_stats.json'),ms);
 r7_thinmesh_check(md,ms,c,sel.inclusion_domains,folder);
 it=find(strcmp(ms.types,'tet'));sel.tetrahedra=size(md.elem{it},2);sel.min_mesh_quality=ms.minquality;sel.mean_mesh_quality=ms.meanquality;sel.mesh_seconds=toc(timer);sel.global_hmax_mm=3.5;sel.ellipse_hmax_mm=.5;
 r3_json(fullfile(folder,'selections.json'),sel);c.mesh_level='fine';r3_json(fullfile(folder,'parameters.json'),c);
 bg=solve_state(m,c,sel,'z',NaN,fullfile(folder,'background'),NaN);r3_json(fullfile(folder,'actual_geometry.json'),measure_geometry(m,c,sel,folder));
 result=solve_state(m,c,sel,'z',t.sigma_f_S_m,fullfile(folder,'target'),bg.G_S);base=jsondecode(fileread(fullfile(source,sf_token(t.sigma_f_S_m),'SUCCESS.json')));
 r3_json(fullfile(folder,'comparison.json'),struct('case_id',c.id,'sigma_f_S_m',t.sigma_f_S_m,'base_g_S',base.g_S,'fine_g_S',result.g_S,...
  'gain_change_S',result.g_S-base.g_S,'gain_relative_change',abs(result.g_S-base.g_S)/abs(result.g_S),'base_tets',base.tetrahedra,'fine_tets',result.tetrahedra));
 clear m cleanup;
end
end

function addmat(co,tag,ids,sx,sy,sz)
co.material.create(tag,'Common');ma=co.material(tag);ma.selection.set(ids);ma.propertyGroup('def').set('electricconductivity',{sx '0' '0' '0' sy '0' '0' '0' sz});ma.propertyGroup('def').set('relpermittivity',{'1'});
end
function token=sf_token(sf),token=['sf_' strrep(sprintf('%.10g',sf),'.','p')];end
function remove_model(tag),try,com.comsol.model.util.ModelUtil.remove(tag);catch,end,end
