function r7b_fix_uniform_check(root,importedModelPath,outputFolder,directReferenceFolder)
%R7B_PILOT_ELECTRICAL Uniform-cube interface verification, not CT evidence.
% Existing single COMSOL/LiveLink connection only. No service is started.
% Existing explicit-cut volume mesh is not rebuilt or element-centre recoloured.
% A separate native CAD cube/Ellipsoid provides a direct reference.
if nargin<1||isempty(root),root=fileparts(fileparts(mfilename('fullpath')));end
assert(nargin>=2 && ~isempty(importedModelPath),'Supply the accepted new geometry path');
assert(nargin>=3 && ~isempty(outputFolder),'Supply a separate new uniform-check output folder');
if nargin<4,directReferenceFolder='';end
assert(isfile(importedModelPath),'Actual meshed Union source is missing');
if ~isfolder(outputFolder),mkdir(outputFolder);end
done=fullfile(outputFolder,'PILOT_ELECTRICAL_SUCCESS.json');
if isfile(done),fprintf('R7B_PILOT_ELECTRICAL_REUSE %s\n',done);return;end
cfg=struct('sigma_solid_S_m',1e-5,'sigma_fluid_S_m',.18,'deltaV_V',1,...
 'potential_order',1,'geometry_shape_order','linear','stationary_relative_tolerance',1e-10,...
 'conservation_tolerance',1e-7,'background_analytic_relative_tolerance',1e-6,...
 'gain_comparison_relative_target',.01,'whole_conductance_comparison_relative_target',.005,...
 'formal_CT_target',false,'purpose','Uniform-cube native CAD versus actual explicit-cut interface verification only',...
 'imported_source_model',importedModelPath,'source_geometry_modified',false,...
 'threshold_note','Engineering pilot targets declared before new electrical solves; not universal precision criteria');
r3_json(fullfile(outputFolder,'pilot_configuration.json'),cfg);
types={'imported_localfit','direct_native_CAD'};results=cell(1,2);L=NaN;
for k=1:2
 if k==2 && ~isempty(directReferenceFolder)
  ref=fullfile(directReferenceFolder,'PAIR_SUCCESS.json');results{k}=jsondecode(fileread(ref));
  assert(results{k}.sigma_solid_S_m==cfg.sigma_solid_S_m && results{k}.sigma_fluid_S_m==cfg.sigma_fluid_S_m && abs(results{k}.L_m/L-1)<1e-12);
  assert(isfile(results{k}.target.model_path)&&isfile(results{k}.background.model_path));
  r3_json(fullfile(outputFolder,'reused_direct_reference.json'),struct('reference',ref,'reason','Same already-solved native CAD uniform reference; only imported cube boundary resolution changed'));continue;
 end
 kind=types{k};folder=fullfile(outputFolder,kind);if ~isfolder(folder),mkdir(folder);end
 complete=fullfile(folder,'PAIR_SUCCESS.json');
 if isfile(complete)
  results{k}=jsondecode(fileread(complete));L=results{k}.L_m;continue;
 end
 tag=['R7BPilot_' num2str(k)];clean=onCleanup(@()remove_model(tag)); %#ok<NASGU>
 try
  if k==1
   m=mphload(importedModelPath,tag);co=m.component('comp1');co.sorder('linear');
   me=co.mesh('mesh1');assert(strcmp(char(me.geometricModel()),''),'Imported pilot mesh must define its own topology');
   co.geometricModel('mesh/mesh1');
   vertices=double(me.getVertex());lo=min(vertices,[],2);hi=max(vertices,[],2);L=hi(1)-lo(1);
   assert(L>0&&max(abs(lo))<1e-9*L&&max(abs(hi-L))<1e-9*L,'Union source is not the expected origin-based cube');clear vertices;
   assert(abs(L/(128*8.683e-6)-1)<1e-10,'Uniform pilot has unexpected physical scale');
  else
   [m,me]=direct_cad(tag,L,folder);
  end
  com.comsol.model.util.ModelUtil.showProgress(fullfile(folder,'solver_progress.log'));
  [sel,geom]=existing_domain_check(m,me,L,kind,folder);
  setup_electrical(m,sel,cfg,kind,folder);
  solver=r7_ct_native_solver(m,cfg.stationary_relative_tolerance);
  r3_json(fullfile(folder,'solver_configuration.json'),solver);
  bg=solve_one(m,sel,cfg,L,0,fullfile(folder,'background'),NaN);
  target=solve_one(m,sel,cfg,L,1,fullfile(folder,'target'),bg.G_S);
  results{k}=struct('status','success','kind',kind,'formal_CT_target',false,'L_m',L,...
   'sigma_solid_S_m',cfg.sigma_solid_S_m,'sigma_fluid_S_m',cfg.sigma_fluid_S_m,...
   'G0_S',bg.G_S,'Gf_S',target.G_S,'g_S',target.G_S-bg.G_S,...
   'eta',target.G_S/bg.G_S-1,'G0_analytic_S',cfg.sigma_solid_S_m*L,...
   'background_relative_analytic_error',bg.G_S/(cfg.sigma_solid_S_m*L)-1,...
   'geometry',geom,'background',bg,'target',target,...
   'completed_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
  r3_json(complete,results{k});
 catch ME
  r3_json(fullfile(folder,['FAILURE_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '.json']),...
   struct('status','FAILED','kind',kind,'formal_CT_target',false,'message',ME.message,...
    'identifier',ME.identifier,'report',getReport(ME,'extended','hyperlinks','off')));rethrow(ME);
 end
 clear m clean;
end
imp=results{1};cad=results{2};
gainErr=abs(imp.g_S-cad.g_S)/abs(cad.g_S);wholeErr=abs(imp.Gf_S-cad.Gf_S)/abs(cad.Gf_S);
pass=gainErr<=cfg.gain_comparison_relative_target&&wholeErr<=cfg.whole_conductance_comparison_relative_target;
summary=struct('schema','r7b_uniform_interface_electrical_v1','status','COMPUTED',...
 'formal_CT_target',false,'L_m',L,'imported_localfit',imp,'direct_native_CAD',cad,...
 'imported_minus_CAD_g_S',imp.g_S-cad.g_S,'gain_relative_difference',gainErr,...
 'target_conductance_relative_difference',wholeErr,'comparison_engineering_target_passed',pass,...
 'gain_relative_target',cfg.gain_comparison_relative_target,'whole_conductance_relative_target',cfg.whole_conductance_comparison_relative_target,...
 'interpretation','Same uniform physical parameters, actual partitioned meshes and P1 geometry/potential. Difference includes finite CAD tessellation and FE discretization; this does not certify formal CT geometry or convergence.',...
 'completed_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
r3_json(fullfile(outputFolder,'comparison.json'),summary);
if pass
 summary.status='success';r3_json(done,summary);
else
 summary.status='COMPUTED_COMPARISON_TARGET_NOT_MET';r3_json(fullfile(outputFolder,'PILOT_ELECTRICAL_UNRESOLVED.json'),summary);
end
fprintf('R7B_PILOT_ELECTRICAL_COMPARISON gainRel %.8g wholeRel %.8g passed %d\n',gainErr,wholeErr,pass);
end

function [m,me]=direct_cad(tag,L,folder)
import com.comsol.model.util.*
geometryPath=fullfile(folder,'geometry_meshed.mph');
if isfile(geometryPath)
 m=mphload(geometryPath,tag);me=m.component('comp1').mesh('mesh1');
 fprintf('R7B_PILOT_DIRECT_MESH_REUSE tet %d\n',me.getNumElem('tet'));return;
end
m=ModelUtil.create(tag);m.label('round7b_uniform_direct_native_CAD.mph');
ModelUtil.showProgress(fullfile(folder,'solver_progress.log'));
m.component.create('comp1',true);co=m.component('comp1');co.sorder('linear');
co.geom.create('geom1',3);g=co.geom('geom1');g.lengthUnit('m');
g.create('cube','Block');g.feature('cube').set('size',[L L L]);g.feature('cube').set('pos',[0 0 0]);
g.create('ell_ct','Ellipsoid');e=g.feature('ell_ct');e.set('type','solid');
e.set('semiaxes',L*[.25 .25 .0125]);e.set('pos',L*[.5 .5 .5]);
e.set('axistype','cartesian');e.set('axis',[-sqrt(3)/2 0 .5]);e.set('rot',0);e.set('selresult','on');g.run();
ell=reshape(double(co.selection('geom1_ell_ct_dom').entities(3)),1,[]);assert(numel(ell)==1,'Direct CAD must produce one ellipsoid domain');
co.mesh.create('mesh1','geom1');me=co.mesh('mesh1');me.create('tet','FreeTet');
sz=me.feature('size');sz.set('custom','on');sz.set('hmax',L/20);sz.set('hmin',L/20000);sz.set('hgrad',1.35);sz.set('hcurve',.15);
me.feature('tet').create('ellSize','Size');fs=me.feature('tet').feature('ellSize');
fs.selection.geom('geom1',3);fs.selection.set(ell);fs.set('custom','on');
fs.set('hmax',.0125*L/3);fs.set('hmaxactive',true);
fs.set('hmin',.0125*L/30);fs.set('hminactive',true);fs.set('hgrad',1.25);fs.set('hgradactive',true);
meshTimer=tic;me.run();seconds=toc(meshTimer);
r3_json(fullfile(folder,'mesh_settings.json'),struct('global_hmax_m',L/20,'ellipsoid_hmax_m',.0125*L/3,...
 'ellipsoid_hmin_m',.0125*L/30,'mesh_seconds',seconds,'source','Native CAD Block + solid Ellipsoid'));
mphsave(m,geometryPath);
fprintf('R7B_PILOT_DIRECT_MESH tet %d seconds %.2f\n',me.getNumElem('tet'),seconds);
end

function [sel,rec]=existing_domain_check(m,me,L,kind,folder)
% Classify only the two ALREADY CUT domains. No element entities are set.
v=double(me.getVertex());nt=double(me.getNumElem('tet'));nv=size(v,2);
assert(nt>0,'Pilot mesh contains no tetrahedra');
chunk=100000;count=zeros(0,1);vol=zeros(0,1);fmin=zeros(0,1);fmax=zeros(0,1);insideCount=zeros(0,1);
for first=0:chunk:nt-1
 n=min(chunk,nt-first);t=double(me.getElem('tet',int32(first),int32(n)))+1;
 d=reshape(double(me.getElemEntity('tet',int32(first),int32(n))),[],1);nd=max(d);
 if numel(count)<nd,count(end+1:nd,1)=0;vol(end+1:nd,1)=0;fmin(end+1:nd,1)=Inf;fmax(end+1:nd,1)=-Inf;insideCount(end+1:nd,1)=0;end
 a=v(:,t(1,:));b=v(:,t(2,:));c=v(:,t(3,:));e=v(:,t(4,:));
 volume=abs(sum((b-a).*cross(c-a,e-a,1),1))/6;
 q=(a+b+c+e)/4-L/2;u=[.5 0 sqrt(3)/2]*q;w=[0 1 0]*q;normal=[-sqrt(3)/2 0 .5]*q;
 F=(u/(.25*L)).^2+(w/(.25*L)).^2+(normal/(.0125*L)).^2-1;
 count=count+accumarray(d,ones(n,1),[numel(count) 1],@sum,0);
 vol=vol+accumarray(d,volume(:),[numel(count) 1],@sum,0);
 insideCount=insideCount+accumarray(d,double(F(:)<0),[numel(count) 1],@sum,0);
 for id=reshape(unique(d),1,[]),fmin(id)=min(fmin(id),min(F(d==id)));fmax(id)=max(fmax(id),max(F(d==id)));end
end
ids=find(count>0);
if strcmp(kind,'imported_localfit')
 % The actual CT mesh contains multiple phase/region domains. For THIS
 % homogeneous electrical test ONLY, group the inherited ellipsoid domains
 % as one material. No vertex, element, entity or CT label is changed.
 tags=cell(me.feature('ellvol').componentSelection());inner=[];
 for j=1:numel(tags),inner=union(inner,reshape(double(m.component('comp1').selection(tags{j}).entities(3)),1,[]));end
 inner=inner(:);outer=setdiff(ids,inner);
 assert(~isempty(inner)&&~isempty(outer)&&all(ismember(inner,ids)),'Actual fitted ellipsoid selection missing');
else
 assert(numel(ids)==2,'Direct CAD reference must have two material domains');
 [~,j]=min(vol(ids));inner=ids(j);outer=ids(ids~=inner);
end
ideal=4*pi*(.25*L)^2*(.0125*L)/3;volumeError=sum(vol(inner))/ideal-1;
assert(abs(sum(vol)/L^3-1)<1e-8,'Pilot actual tetrahedron cube volume changed');
assert(abs(volumeError)<.005,'Pilot ellipsoid actual volume differs by more than 0.5 percent');
% A planar inscribed face and the analytic surface enclose a thin shell.
% Save all centre-sign ranges, allowing a declared 0.02 F-band only for
% this domain-label consistency check. Actual thickness/sag acceptance is
% performed independently on the actual boundary by the geometry checks.
assert(max(fmax(inner))<.02&&min(fmin(inner))<-.5,'Small existing domain is inconsistent with ellipsoid interior');
assert(min(fmin(outer))>-.02&&max(fmax(outer))>1,'Outer existing domain intrudes beyond permitted analytic-surface check band');
tol=L*1e-7;top=reshape(mphselectbox(m,'geom1',[-tol L+tol;-tol L+tol;L-tol L+tol],'boundary'),1,[]);
bottom=reshape(mphselectbox(m,'geom1',[-tol L+tol;-tol L+tol;-tol tol],'boundary'),1,[]);
assert(~isempty(top)&&~isempty(bottom)&&isempty(intersect(top,bottom)),'Pilot external electrode selections failed');
stats=mphmeshstats(m,'mesh1');r3_json(fullfile(folder,'mesh_stats.json'),stats);
sel=struct('all_domains',ids','inclusion_domains',inner,'outside_domains',outer,'top',top,'bottom',bottom,...
 'tetrahedra',nt,'vertices',nv,'min_mesh_quality',stats.minquality,'mean_mesh_quality',stats.meanquality);
rec=struct('kind',kind,'actual_domain_count',numel(ids),'domain_ids',ids','domain_volumes_m3',vol(ids)',...
 'domain_tetrahedron_counts',count(ids)','F_min_at_actual_tet_centres',fmin(ids)','F_max_at_actual_tet_centres',fmax(ids)',...
 'negative_F_tetrahedron_count',insideCount(ids)','interior_domain',inner,'exterior_domain',outer,...
 'total_volume_m3',sum(vol),'ellipsoid_actual_volume_m3',sum(vol(inner)),'ellipsoid_ideal_volume_m3',ideal,...
 'ellipsoid_volume_relative_error',volumeError,'centroid_F_consistency_band',.02,...
 'label_method','Inherited explicit-cut ellipsoid selection grouped across original CT domains; no tetrahedron entity reassignment or geometry rebuild',...
 'component_geometric_model',char(m.component('comp1').geometricModel()),'mesh_geometric_model',char(me.geometricModel()),...
 'tetrahedra',nt,'vertices',nv,'geometry_acceptance_still_requires_independent_actual_surface_probes',true);
r3_json(fullfile(folder,'domain_verification.json'),rec);r3_json(fullfile(folder,'selections.json'),sel);
end

function setup_electrical(m,sel,cfg,kind,folder)
co=m.component('comp1');assert(isempty(cell(co.physics.tags))&&isempty(cell(co.material.tags)),'Pilot source must be physics/material-free');
m.param.set('sigma_solid',sprintf('%.17g[S/m]',cfg.sigma_solid_S_m));
m.param.set('sigma_fluid',sprintf('%.17g[S/m]',cfg.sigma_fluid_S_m));m.param.set('frac_on','0');m.param.set('DV','1[V]');
addmat(co,'outer',sel.outside_domains,'sigma_solid');
addmat(co,'inner',sel.inclusion_domains,'sigma_solid+frac_on*(sigma_fluid-sigma_solid)');
co.physics.create('ec','ConductiveMedia','geom1');ec=co.physics('ec');ec.prop('ShapeProperty').set('order_electricpotential',1);
actual=reshape(double(ec.selection.entities(3)),1,[]);
assert(isequal(sort(actual),sort(sel.all_domains)),'Pilot Electric Currents domain selection does not equal actual cube plus ellipsoid');
ec.create('pot1','ElectricPotential',2);ec.feature('pot1').selection.set(sel.top);ec.feature('pot1').set('V0','DV');
ec.create('gnd1','Ground',2);ec.feature('gnd1').selection.set(sel.bottom);
m.study.create('std1');m.study('std1').create('stat','Stationary');m.study('std1').feature('stat').mesh('geom1','mesh1');
r3_json(fullfile(folder,'physics_attachment.json'),struct('kind',kind,'component_geometric_model',char(co.geometricModel()),...
 'physics_geometry_argument','geom1','actual_physics_domains',actual,'physics_covers_all_material_domains',true));
end

function r=solve_one(m,sel,cfg,L,frac,folder,G0)
done=fullfile(folder,'SUCCESS.json');
if isfile(done),r=jsondecode(fileread(done));assert(isfile(r.model_path),'Pilot checkpoint model missing');return;end
if ~isfolder(folder),mkdir(folder);end
m.param.set('frac_on',num2str(frac));timer=tic;m.sol('sol1').runAll();seconds=toc(timer);
I=mphint2(m,'reacf(V)',2,'selection',sel.top,'method','summation');
Ig=mphint2(m,'reacf(V)',2,'selection',sel.bottom,'method','summation');
power=mphint2(m,'ec.Jx*ec.Ex+ec.Jy*ec.Ey+ec.Jz*ec.Ez',3,'intorder',2);
area=mphint2(m,'1',2,'selection',sel.top,'unit','m^2','intorder',2);
G=I;flux=abs(I+Ig)/abs(I);energy=abs(I-power)/abs(I);
modelPath=fullfile(folder,'solution.mph');mphsave(m,modelPath);
r=struct('status','COMPUTED','formal_CT_target',false,'fracture_on',frac,'G_S',G,'G0_S',G0,'g_S',G-G0,...
 'reaction_positive_A',I,'reaction_ground_A',Ig,'dissipation_W',power,'electrode_area_m2',area,...
 'reaction_flux_relative_error',flux,'energy_relative_error',energy,'tetrahedra',sel.tetrahedra,...
 'vertices',sel.vertices,'min_mesh_quality',sel.min_mesh_quality,'potential_order',1,'geometry_order','linear',...
 'solve_seconds',seconds,'model_path',modelPath,'stationary_relative_tolerance',cfg.stationary_relative_tolerance);
r3_json(fullfile(folder,'computed.json'),r);
assert(G>0&&abs(area/L^2-1)<1e-8,'Pilot current/electrode area failed');
assert(max(flux,energy)<=cfg.conservation_tolerance,'Pilot reaction flux or power balance exceeds declared tolerance');
if frac==0,assert(abs(G/(cfg.sigma_solid_S_m*L)-1)<=cfg.background_analytic_relative_tolerance,'Uniform no-fracture analytic solution check failed');end
r.status='success';r.completed_local=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
r3_json(done,r);fprintf('R7B_PILOT_SOLVED frac %d G %.12g g %.12g flux %.3g energy %.3g seconds %.2f\n',frac,G,G-G0,flux,energy,seconds);
end

function addmat(co,tag,ids,expr)
co.material.create(tag,'Common');ma=co.material(tag);ma.selection.set(reshape(ids,1,[]));
ma.propertyGroup('def').set('electricconductivity',{expr '0' '0' '0' expr '0' '0' '0' expr});
end
function remove_model(tag),try,com.comsol.model.util.ModelUtil.remove(tag);catch,end,end
