function r=solve_round2_case(c,folder)
% Round 2 extension of the verified 67-case solver; original solver untouched.
% Geometry coordinates are mm; electrical fields and currents are SI.
% End caps are true ellipse/polygon intersections of axial length 3 mm.
c=defaults(c);zInterfaces=reshape(c.z_interfaces_mm,1,[]);
layerMaterials=reshape(c.layer_materials,1,[]);
assert(numel(layerMaterials)==numel(zInterfaces)+1,'Layer vector lengths disagree');
assert(all(diff([0 zInterfaces 120])>0),'Layer interfaces must increase inside 120 mm');
assert(all(ismember(layerMaterials,[0 1])),'Layer phases must be 0 or 1');
import com.comsol.model.*
import com.comsol.model.util.*
model=ModelUtil.create('Rock2D');model.label([c.id '.mph']);
ModelUtil.showProgress(fullfile(folder,'solver_progress.log'));
model.param.set('L','120[mm]');model.param.set('layer_h',sprintf('%.15g[mm]',c.h_mm));
model.param.set('A',sprintf('%.12g[mm]',c.A_mm));
model.param.set('a',sprintf('%.12g[mm]',c.a_mm));model.param.set('b',sprintf('%.12g[mm]',c.b_mm));
model.param.set('cx',sprintf('%.12g[mm]',c.cx_mm));model.param.set('cz',sprintf('%.12g[mm]',c.cz_mm));
model.param.set('theta',sprintf('%.12g[deg]',c.theta_deg));
model.param.set('d_out',sprintf('%.15g[m]',c.out_of_plane_thickness_m));model.param.set('DV','1[V]');
model.param.set('sigma_low',sprintf('%.12g[S/m]',c.sigma_low_S_m));
model.param.set('sigma_high',sprintf('%.12g[S/m]',c.sigma_high_S_m));
model.param.set('sigma_frac',sprintf('%.12g[S/m]',c.sigma_frac_S_m));
model.component.create('comp1',true);co=model.component('comp1');co.geom.create('geom1',2);
g=co.geom('geom1');g.lengthUnit('mm');g.create('box','Rectangle');g.feature('box').set('size',[120 120]);
amplitudes=c.phase_sign*c.A_mm*ones(size(zInterfaces));
if c.buffer,amplitudes(zInterfaces<=20 | zInterfaces>=100)=0;end
if ~isempty(c.interface_amplitudes_mm),amplitudes=reshape(c.interface_amplitudes_mm,1,[]);end
assert(numel(amplitudes)==numel(zInterfaces),'Interface amplitude vector length mismatch');
assert(all(diff([0 zInterfaces 120])>abs(diff([0 amplitudes 0]))),'Layer interfaces may cross');
if strcmp(c.background,'layered')
 tags=cell(size(zInterfaces));
 for j=1:numel(zInterfaces)
  tag=sprintf('layer%d',j);tags{j}=tag;g.create(tag,'ParametricCurve');
  f=g.feature(tag);f.set('parname','s');f.set('parmin','0');f.set('parmax','120');
  f.set('coord',{'s',sprintf('%.15g+%.15g*sin(2*pi*s/120)',zInterfaces(j),amplitudes(j))});f.set('rtol',1e-8);
 end
 g.create('part','Partition');g.feature('part').selection('input').set({'box'});
 g.feature('part').selection('tool').set(tags);g.feature('part').set('keeptool','off');
end
if c.fracture
 g.create('ell','Ellipse');g.feature('ell').set('semiaxes',[c.a_mm c.b_mm]);
 g.feature('ell').set('pos',[c.cx_mm c.cz_mm]);g.feature('ell').set('rot',c.theta_deg);
 g.feature('ell').set('selresult','on');
end
% Duplicate ellipse intersections only partition the geometry; they do not
% change the original ellipse or prescribe any additional electrical phase.
if c.fracture
 for cap={'capminus','capplus'}
  tag=cap{1};elltag=[tag 'ellipse'];poltag=[tag 'polygon'];
  g.create(elltag,'Ellipse');g.feature(elltag).set('semiaxes',[c.a_mm c.b_mm]);
  g.feature(elltag).set('pos',[c.cx_mm c.cz_mm]);g.feature(elltag).set('rot',c.theta_deg);
  if strcmp(tag,'capminus'),bounds=[-c.a_mm -c.a_mm+c.endcap_length_mm -2*c.b_mm 2*c.b_mm];
  else,bounds=[c.a_mm-c.endcap_length_mm c.a_mm -2*c.b_mm 2*c.b_mm];end
  addpoly(g,poltag,bounds,c);g.create(tag,'Intersection');
  g.feature(tag).selection('input').set({elltag,poltag});g.feature(tag).set('selresult','on');
 end
end
fracTag='ell';isMechanism=~strcmp(c.intervention,'none');
hasRemove=ismember(c.intervention,{'cut','side'});
if hasRemove
 if strcmp(c.intervention,'cut'),bounds=[c.cut_u_mm-c.cut_width_mm/2 c.cut_u_mm+c.cut_width_mm/2 -2*c.b_mm 2*c.b_mm];
 else,bounds=[c.side_u1_mm c.side_u2_mm c.side_width_mm 2*c.b_mm];end
 addpoly(g,'cutter',bounds,c);
 g.create('dfrac','Difference');g.feature('dfrac').selection('input').set({'ell'});
 g.feature('dfrac').selection('input2').set({'cutter'});g.feature('dfrac').set('selresult','on');fracTag='dfrac';
 % A second, identical ellipse/cutter intersection labels the removed region.
 g.create('refell','Ellipse');g.feature('refell').set('semiaxes',[c.a_mm c.b_mm]);
 g.feature('refell').set('pos',[c.cx_mm c.cz_mm]);g.feature('refell').set('rot',c.theta_deg);
 addpoly(g,'rempoly',bounds,c);g.create('removed','Intersection');
 g.feature('removed').selection('input').set({'refell','rempoly'});g.feature('removed').set('selresult','on');
 g.create('compensate','Circle');g.feature('compensate').set('pos',[c.comp_x_mm c.comp_z_mm]);
 g.feature('compensate').set('r',c.comp_radius_mm);g.feature('compensate').set('selresult','on');
end
if isMechanism,addpoly(g,'window',[c.cut_u_mm-3 c.cut_u_mm+3 -3 3],c);g.feature('window').set('selresult','on');end
g.run;
fracIDs=[];compIDs=[];removedIDs=[];windowIDs=[];capminusIDs=[];capplusIDs=[];
if c.fracture,fracIDs=reshape(double(co.selection(['geom1_' fracTag '_dom']).entities(2)),1,[]);end
bodyIDs=fracIDs;
if c.fracture
 capminusIDs=reshape(double(co.selection('geom1_capminus_dom').entities(2)),1,[]);
 capplusIDs=reshape(double(co.selection('geom1_capplus_dom').entities(2)),1,[]);
end
if hasRemove
 compIDs=reshape(double(co.selection('geom1_compensate_dom').entities(2)),1,[]);
 removedIDs=reshape(double(co.selection('geom1_removed_dom').entities(2)),1,[]);
 fracIDs=unique([fracIDs compIDs]);
end
if isMechanism,windowIDs=reshape(double(co.selection('geom1_window_dom').entities(2)),1,[]);end
inclusionIDs=unique([bodyIDs removedIDs]);
% Named output of ellipse follows domain partitioning in the final union.
co.mesh.create('mesh1');me=co.mesh('mesh1');me.create('tri','FreeTri');
globalh=[3 2 1.2];localh=[0.45 0.25 0.15];level=c.mesh_level;
me.feature('size').set('custom','on');me.feature('size').set('hmax',globalh(level));
me.feature('size').set('hmin',0.005);me.feature('size').set('hgrad',1.25);me.feature('size').set('hcurve',0.2);
if ~isempty(fracIDs)
 me.feature('tri').create('fs','Size');fs=me.feature('tri').feature('fs');fs.selection.geom('geom1',2);fs.selection.set(fracIDs);
 fs.set('custom','on');fs.set('hmax',localh(level));fs.set('hmaxactive',true);
 fs.set('hmin',0.002);fs.set('hminactive',true);fs.set('hgrad',1.2);fs.set('hgradactive',true);
end
if isMechanism
 me.feature('tri').create('ws','Size');ws=me.feature('tri').feature('ws');ws.selection.geom('geom1',2);ws.selection.set(windowIDs);
 ws.set('custom','on');ws.set('hmax',localh(level));ws.set('hmaxactive',true);ws.set('hmin',0.002);ws.set('hminactive',true);
end
me.run;
[ms,md]=mphmeshstats(model,'mesh1','qualitymeasure','skewness');
save(fullfile(folder,'mesh.mat'),'ms','md');write_json(fullfile(folder,'mesh_stats.json'),ms);
ti=find(strcmp(ms.types,'tri'));verts=md.vertex;tri=double(md.elem{ti})+1;ent=double(md.elementity{ti});
fprintf('MESH vertices coordinate max %.9g; triangles %d; domains %d\n',max(verts(:)),size(tri,2),numel(unique(ent)));
% Mesh coordinates follow the geometry length unit (mm). Verify at runtime.
if abs(max(verts(:))-120)<1e-5,scale=1;elseif abs(max(verts(:))-0.12)<1e-8,scale=1000;else,error('Unknown mesh coordinate units');end
centers=(verts(:,tri(1,:))+verts(:,tri(2,:))+verts(:,tri(3,:)))/3*scale;
domainIDs=unique(ent);lowIDs=[];highIDs=[];bgLayer=zeros(size(domainIDs));bgMaterial=zeros(size(domainIDs));
if c.passive_ellipse,fracIDs=[];end
for n=1:numel(domainIDs)
 id=domainIDs(n);idx=find(ent==id,1);x=centers(1,idx);z=centers(2,idx);
 lay=1+sum(z>zInterfaces+amplitudes.*sin(2*pi*x/120));bgLayer(n)=lay;
 if ~strcmp(c.background,'layered'),lay=1;bgLayer(n)=lay;end
 bgMaterial(n)=layerMaterials(lay);
 if ~strcmp(c.background,'layered'),bgMaterial(n)=0;end
 if ismember(id,fracIDs),continue;end
 if bgMaterial(n)==0,lowIDs(end+1)=id;else,highIDs(end+1)=id;end
end
addmat(co,'low',lowIDs,'sigma_low');addmat(co,'high',highIDs,'sigma_high');addmat(co,'frac',fracIDs,'sigma_frac');
if strcmp(c.background,'anisotropic')
 model.param.set('sigma_xx',sprintf('%.15g[S/m]',c.sigma_xx_S_m));model.param.set('sigma_zz',sprintf('%.15g[S/m]',c.sigma_zz_S_m));
 co.material('low').propertyGroup('def').set('electricconductivity',{'sigma_xx' '0' '0' '0' 'sigma_zz' '0' '0' '0' 'sigma_xx'});
end
selected=sort([lowIDs highIDs fracIDs]);assert(isequal(selected(:),sort(domainIDs(:))),'Material domains not uniquely covered');
tol=1e-5;
top=mphselectbox(model,'geom1',[-tol 120+tol;120-tol 120+tol],'boundary');
bottom=mphselectbox(model,'geom1',[-tol 120+tol;-tol tol],'boundary');
left=mphselectbox(model,'geom1',[-tol tol;-tol 120+tol],'boundary');
right=mphselectbox(model,'geom1',[120-tol 120+tol;-tol 120+tol],'boundary');
co.physics.create('ec','ConductiveMedia','geom1');ec=co.physics('ec');
ec.prop('d').set('d','d_out');
if strcmp(c.load,'z'),pos=top;neg=bottom;component='ec.Jy';else,pos=right;neg=left;component='ec.Jx';end
ec.create('pot1','ElectricPotential',1);ec.feature('pot1').selection.set(pos);ec.feature('pot1').set('V0','DV');
ec.create('gnd1','Ground',1);ec.feature('gnd1').selection.set(neg);
sel=struct('low_domains',lowIDs,'high_domains',highIDs,'fracture_domains',fracIDs,'all_domains',domainIDs,...
 'domain_background_layer',bgLayer,'domain_background_material',bgMaterial,'top_boundaries',top,'bottom_boundaries',bottom,...
 'left_boundaries',left,'right_boundaries',right,'positive_electrode',pos,'ground_electrode',neg,...
 'interface_amplitudes_mm',amplitudes,'z_interfaces_mm',zInterfaces,'layer_materials',layerMaterials,'mesh_coordinate_to_mm',scale,...
 'compensation_domains',compIDs,'removed_domains',removedIDs,'window_domains',windowIDs,'inclusion_geometry_domains',inclusionIDs,...
 'capminus_domains',capminusIDs,'capplus_domains',capplusIDs);
write_json(fullfile(folder,'selections.json'),sel);
model.study.create('std1');model.study('std1').create('stat','Stationary');
mphsave(model,fullfile(folder,'presolve.mph'));
st=tic;model.study('std1').run;solveSeconds=toc(st);
qpos=mphint2(model,'ec.Jx*nx+ec.Jy*ny',1,'selection',pos,'intorder',8);
qneg=mphint2(model,'ec.Jx*nx+ec.Jy*ny',1,'selection',neg,'intorder',8);
power=mphint2(model,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',2,'intorder',8);
epos=mphint2(model,component,1,'selection',pos,'intorder',8);
len=mphint2(model,'1',1,'selection',pos,'intorder',8,'unit','m');
sigmaGradient=-epos/len*0.12;
reactionPositive=mphint2(model,'reacf(V)',1,'selection',pos,'method','summation');
reactionGround=mphint2(model,'reacf(V)',1,'selection',neg,'method','summation');
sigma=reactionPositive/c.out_of_plane_thickness_m; % L/load-face width=1, DV=1 V.
al=area(model,lowIDs);ah=area(model,highIDs);af=area(model,fracIDs);
repllow=0;replhigh=0;domains_area=zeros(size(domainIDs));diss=zeros(size(domainIDs));
for n=1:numel(domainIDs)
 id=domainIDs(n);domains_area(n)=area(model,id);diss(n)=mphint2(model,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',2,'selection',id,'intorder',8);
 if ismember(id,fracIDs)
  if bgMaterial(n)==0,repllow=repllow+domains_area(n);else,replhigh=replhigh+domains_area(n);end
 end
end
T=table(domainIDs(:),bgLayer(:),bgMaterial(:),ismember(domainIDs(:),fracIDs),domains_area(:),diss(:),...
 'VariableNames',{'domain_id','background_layer','background_material','is_fracture','area_mm2','dissipation_W_per_m'});writetable(T,fullfile(folder,'domain_integrals.csv'));
fields=mpheval(model,{'V','ec.Ex','ec.Ey','ec.Jx','ec.Jy','ec.Jx*ec.Ex+ec.Jy*ec.Ey'},'edim',2,'refine',1,'smooth','none');
save(fullfile(folder,'fields.mat'),'fields','-v7.3');
% Solution-evaluation p follows geometry units (mm), confirmed in prior suite.
writetable(table(fields.p(1,:)'*0.001,fields.p(2,:)'*0.001,fields.d1(:),fields.d2(:),fields.d3(:),fields.d4(:),fields.d5(:),fields.d6(:),...
 'VariableNames',{'x_m','z_m','V_V','Ex_V_m','Ez_V_m','Jx_A_m2','Jz_A_m2','dissipation_W_m3'}),fullfile(folder,'fields.csv'));
r=struct('case_id',c.id,'load',c.load,'sigma_eff_S_m',sigma,'sigma_gradient_S_m',sigmaGradient,'q_positive_out_A_per_m',qpos,...
 'q_ground_out_A_per_m',qneg,'electrode_length_m',len,'power_source_W_per_m',-qpos,...
 'power_volume_W_per_m',power,'flux_relative_error',abs(qpos+qneg)/max(abs([qpos qneg])),...
 'energy_relative_error',abs(power+qpos)/abs(qpos),'area_low_mm2',al,'area_high_mm2',ah,'area_fracture_mm2',af,...
 'replaced_low_mm2',repllow,'replaced_high_mm2',replhigh,'triangles',size(tri,2),...
 'min_mesh_quality',ms.minquality,'mean_mesh_quality',ms.meanquality,'solve_seconds',solveSeconds,...
 'field_coordinate_max_raw',max(fields.p(:)),'mesh_level',level,'comsol_version','6.3.0.290');
r.reaction_positive_A=reactionPositive;r.reaction_ground_A=reactionGround;
r.reaction_flux_relative_error=abs(reactionPositive+reactionGround)/abs(reactionPositive);
r.reaction_power_relative_error=abs(reactionPositive/c.out_of_plane_thickness_m-power)/abs(power);
r.power_source_reaction_W_per_m=reactionPositive/c.out_of_plane_thickness_m;
r.out_of_plane_thickness_m=c.out_of_plane_thickness_m;
r.endcap_length_mm=c.endcap_length_mm;r.h_nominal_mm=c.h_mm;
r.original_ellipse_area_mm2=area(model,inclusionIDs);
bgLowIDs=domainIDs(bgMaterial==0);bgHighIDs=domainIDs(bgMaterial==1);
r.background_low_area_mm2=area(model,bgLowIDs);r.background_high_area_mm2=area(model,bgHighIDs);
r.original_ellipse_replaced_low_mm2=area(model,intersect(inclusionIDs,bgLowIDs));
r.original_ellipse_replaced_high_mm2=area(model,intersect(inclusionIDs,bgHighIDs));
r.capminus_area_mm2=area(model,capminusIDs);r.capplus_area_mm2=area(model,capplusIDs);
r.capminus_high_area_mm2=area(model,intersect(capminusIDs,bgHighIDs));
r.capminus_low_area_mm2=area(model,intersect(capminusIDs,bgLowIDs));
r.capplus_high_area_mm2=area(model,intersect(capplusIDs,bgHighIDs));
r.capplus_low_area_mm2=area(model,intersect(capplusIDs,bgLowIDs));
r.capminus_high_fraction=r.capminus_high_area_mm2/r.capminus_area_mm2;
r.capplus_high_fraction=r.capplus_high_area_mm2/r.capplus_area_mm2;
primitive=@(u)c.b_mm*(u.*sqrt(max(0,1-u.^2/c.a_mm^2))+c.a_mm*asin(u/c.a_mm));
r.original_ellipse_analytic_area_mm2=pi*c.a_mm*c.b_mm;
r.endcap_analytic_area_mm2=primitive(c.a_mm)-primitive(c.a_mm-c.endcap_length_mm);
r.ellipse_geometry_area_error_mm2=r.original_ellipse_area_mm2-r.original_ellipse_analytic_area_mm2;
r.removed_area_mm2=area(model,removedIDs);r.compensation_area_mm2=area(model,compIDs);
r.removed_original_low_mm2=area(model,intersect(removedIDs,bgLowIDs));
r.removed_original_high_mm2=area(model,intersect(removedIDs,bgHighIDs));
r.compensation_original_low_mm2=area(model,intersect(compIDs,bgLowIDs));
r.compensation_original_high_mm2=area(model,intersect(compIDs,bgHighIDs));
r.fracture_dissipation_W_per_m=integ(model,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',fracIDs);
r.capminus_dissipation_W_per_m=integ(model,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',capminusIDs);
r.capplus_dissipation_W_per_m=integ(model,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',capplusIDs);
xy=localcoords([-0.95*c.a_mm 0.95*c.a_mm],[0 0],c);tipV=mphinterp(model,'V','coord',xy);
r.neartip_dV_V=tipV(2)-tipV(1);
if strcmp(c.background,'anisotropic') && (~c.fracture || c.passive_ellipse)
 if strcmp(c.load,'x'),expected=c.sigma_xx_S_m;else,expected=c.sigma_zz_S_m;end
elseif strcmp(c.background,'uniform'),expected=c.sigma_low_S_m;
elseif all(amplitudes==0) && (~c.fracture || c.passive_ellipse)
 thicknesses=diff([0 zInterfaces 120]);fh=sum(thicknesses(layerMaterials==1))/120;
 if strcmp(c.load,'z'),expected=1/((1-fh)/c.sigma_low_S_m+fh/c.sigma_high_S_m);
 else,expected=(1-fh)*c.sigma_low_S_m+fh*c.sigma_high_S_m;end
else,expected=NaN;end
r.analytic_expected_S_m=expected;r.analytic_relative_error=abs(sigma/expected-1);
if c.passive_ellipse
 r.background_ellipse_area_mm2=area(model,inclusionIDs);
 r.prediction_upper_bound_delta_sigma_S_m=integ(model,'(sigma_frac-sigma_low)*(ec.Ex^2+ec.Ey^2)',intersect(inclusionIDs,lowIDs))...
   +integ(model,'(sigma_frac-sigma_high)*(ec.Ex^2+ec.Ey^2)',intersect(inclusionIDs,highIDs));
 xy=localcoords([-0.95*c.a_mm 0.95*c.a_mm],[0 0],c);v0=mphinterp(model,'V','coord',xy);
 r.background_neartip_dV_V=v0(2)-v0(1);r.background_neartip_dV2_V2=(v0(2)-v0(1))^2;
end
if isMechanism
 r.window_area_mm2=area(model,windowIDs);
 r.window_dissipation_W_per_m=integ(model,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',windowIDs);
 r.window_low_dissipation_W_per_m=integ(model,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',intersect(windowIDs,lowIDs));
 r.window_fracture_dissipation_W_per_m=integ(model,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',intersect(windowIDs,fracIDs));
 if hasRemove
  intendedIDs=domainIDs(bgMaterial==c.intervention_phase);
  assert(isempty(setdiff(removedIDs,intendedIDs)),'Removed region not wholly in declared original phase');
  assert(isempty(setdiff(compIDs,intendedIDs)),'Compensation not wholly in declared original phase');
  assert(isempty(setdiff(compIDs,fracIDs)),'Compensation not fracture material');
  assert(abs(r.removed_area_mm2-r.compensation_area_mm2)<max(0.01,1e-4*r.removed_area_mm2),'Removal/compensation area mismatch');
 end
 tt=linspace(-1,1,201);s=[(c.cut_u_mm-3)*ones(size(tt)) (c.cut_u_mm+3)*ones(size(tt))];t=[tt tt];
 xy=localcoords(s,t,c);vp=mphinterp(model,'V','coord',xy);
 r.window_dV_V=trapz(tt,vp(numel(tt)+1:end))/2-trapz(tt,vp(1:numel(tt)))/2;
 % All cases use identical local coordinates for direct current-path comparison.
 [su,tv]=meshgrid(linspace(c.cut_u_mm-6,c.cut_u_mm+6,241),linspace(-6,6,241));xy=localcoords(su(:)',tv(:)',c);
 [v,ex,ez,jx,jz,q]=mphinterp(model,{'V','ec.Ex','ec.Ey','ec.Jx','ec.Jy','ec.Jx*ec.Ex+ec.Jy*ec.Ey'},'coord',xy);
 writetable(table(su(:),tv(:),xy(1,:)'*0.001,xy(2,:)'*0.001,v(:),ex(:),ez(:),jx(:),jz(:),q(:),...
 'VariableNames',{'s_mm','t_mm','x_m','z_m','V_V','Ex_V_m','Ez_V_m','Jx_A_m2','Jz_A_m2','dissipation_W_m3'}),fullfile(folder,'local_window_fields.csv'));
 % Trace through the declared cut, positive long-axis current; distances in m.
 tv=linspace(-6,6,2401);xy=localcoords(c.cut_u_mm*ones(size(tv)),tv,c);
 [jx,jz]=mphinterp(model,{'ec.Jx','ec.Jy'},'coord',xy);js=jx*cosd(c.theta_deg)+jz*sind(c.theta_deg);
 r.section_current_signed_A_per_m=trapz(tv*0.001,js);
 halfWidth=c.b_mm*sqrt(max(0,1-(c.cut_u_mm/c.a_mm)^2));
 inside=abs(tv)<=halfWidth;outsideLeft=tv < -halfWidth;outsideRight=tv > halfWidth;
 r.section_original_ellipse_current_A_per_m=trapz(tv(inside)*0.001,js(inside));
 r.section_background_bypass_current_A_per_m=trapz(tv(outsideLeft)*0.001,js(outsideLeft))+trapz(tv(outsideRight)*0.001,js(outsideRight));
 writetable(table(tv(:),js(:),'VariableNames',{'t_mm','Js_A_m2'}),fullfile(folder,'section_current.csv'));
end
model.result.create('pg1','PlotGroup2D');model.result('pg1').label('Potential V - internal inspection');
model.result('pg1').create('surf1','Surface');model.result('pg1').feature('surf1').set('expr','V');
model.result.create('pg2','PlotGroup2D');model.result('pg2').label('Current density - internal inspection');
model.result('pg2').create('surf1','Surface');model.result('pg2').feature('surf1').set('expr','ec.normJ');
mphsave(model,fullfile(folder,[c.id '.mph']));
try,mphsave(model,fullfile(folder,'model_export.m'));catch ME,fprintf('M export warning %s\n',ME.message);end
fprintf('RESULT sigma %.12g flux %.3g energy %.3g area %.9g\n',sigma,r.flux_relative_error,r.energy_relative_error,al+ah+af);
assert(abs(al+ah+af-14400)<0.02,'Area closure failure');
if c.fracture
 assert(abs(r.ellipse_geometry_area_error_mm2)<0.001,'Original ellipse selection/area failure');
 assert(abs(r.capminus_area_mm2-r.endcap_analytic_area_mm2)<0.001 && abs(r.capplus_area_mm2-r.endcap_analytic_area_mm2)<0.001,'Finite end-cap selection/area failure');
end
assert(r.flux_relative_error<0.01 && r.energy_relative_error<0.01,'Conservation tolerance 1 percent exceeded');
assert(r.reaction_flux_relative_error<1e-6 && r.reaction_power_relative_error<1e-6,'Reaction conservation tolerance 1e-6 exceeded');
if ~isnan(expected),assert(r.analytic_relative_error<1e-5,'Analytic baseline failed');end
end

function c=defaults(c)
d=struct('background','layered','A_mm',0,'phase_sign',1,'buffer',true,'fracture',true,...
 'cx_mm',60,'cz_mm',60,'a_mm',15,'b_mm',1.5,'theta_deg',60,...
 'sigma_low_S_m',0.01,'sigma_high_S_m',0.1,'sigma_frac_S_m',10,'mesh_level',1,...
 'load','z','intervention','none','cut_u_mm',0,'cut_width_mm',1,'side_width_mm',0,...
 'comp_radius_mm',0,'comp_x_mm',20,'comp_z_mm',45,'L_mm',120,'h_mm',10,...
 'out_of_plane_thickness_m',1,'passive_ellipse',false,'sigma_xx_S_m',0,'sigma_zz_S_m',0,...
 'z_interfaces_mm',10:10:110,'layer_materials',mod(0:11,2),'interface_amplitudes_mm',[],...
 'endcap_length_mm',3,'side_u1_mm',-2,'side_u2_mm',2,'intervention_phase',0);
names=fieldnames(d);for k=1:numel(names),if ~isfield(c,names{k}),c.(names{k})=d.(names{k});end,end
assert(c.L_mm==120,'Round 2 is restricted to the 120 mm square');
assert(c.endcap_length_mm>0 && c.endcap_length_mm<c.a_mm,'Invalid finite cap length');
assert(c.out_of_plane_thickness_m==1,'This round explicitly uses 1 m out-of-plane thickness');
end

function addpoly(g,tag,bounds,c)
s=[bounds(1) bounds(2) bounds(2) bounds(1)];t=[bounds(3) bounds(3) bounds(4) bounds(4)];xy=localcoords(s,t,c);
g.create(tag,'Polygon');g.feature(tag).set('x',xy(1,:));g.feature(tag).set('y',xy(2,:));g.feature(tag).set('type','solid');
end
function xy=localcoords(s,t,c)
xy=[c.cx_mm+s*cosd(c.theta_deg)-t*sind(c.theta_deg);c.cz_mm+s*sind(c.theta_deg)+t*cosd(c.theta_deg)];
end
function a=integ(model,expr,ids)
if isempty(ids),a=0;else,a=mphint2(model,expr,2,'selection',ids,'intorder',8);end
end

function addmat(co,tag,ids,sig)
if isempty(ids),return;end
co.material.create(tag,'Common');ma=co.material(tag);ma.selection.set(ids);ma.label([tag ' conductivity']);
ma.propertyGroup('def').set('electricconductivity',{sig '0' '0' '0' sig '0' '0' '0' sig});
ma.propertyGroup('def').set('relpermittivity',{'1'});
end
function a=area(model,ids)
if isempty(ids),a=0;else,a=mphint2(model,'1',2,'selection',ids,'intorder',8,'unit','mm^2');end
end
function write_json(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));fclose(fid);
end
