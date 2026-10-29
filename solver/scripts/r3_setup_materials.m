function r3_setup_materials(m,sel)
% Same geometry/mesh/EC and physical DOFs. Only material expressions change.
co=m.component('comp1');tags=cell(co.material.tags);
for k=1:numel(tags),co.material.remove(tags{k});end
all=reshape(sel.all_domains,1,[]);phase=reshape(sel.domain_background_material,1,[]);
D=reshape(sel.inclusion_geometry_domains,1,[]);
m.param.set('r3bgscale','1');m.param.set('r3dlow','0[S/m]');m.param.set('r3dhigh','0[S/m]');
m.param.set('r3unit','0[S/m]');m.param.set('r3phase','0');
m.param.set('r3lo','-1000');m.param.set('r3hi','1000');
co.variable.create('r3coords');
co.variable('r3coords').set('r3smm','((x-cx)*cos(theta)+(y-cz)*sin(theta))/1[mm]');
for p=0:1
 if p==0,sig='sigma_low';delta='r3dlow';else,sig='sigma_high';delta='r3dhigh';end
 for inside=0:1
  ids=all(phase==p & ismember(all,D)==inside);if isempty(ids),continue;end
  expr=['r3bgscale*' sig];
  if inside,expr=[expr '+' delta '+r3unit*(r3phase==' num2str(p) ')*(r3smm>=r3lo)*(r3smm<r3hi)'];end
  tag=sprintf('r3p%d_d%d',p,inside);co.material.create(tag,'Common');ma=co.material(tag);ma.selection.set(ids);
  ma.propertyGroup('def').set('electricconductivity',{expr '0' '0' '0' expr '0' '0' '0' expr});
  ma.propertyGroup('def').set('relpermittivity',{'1'});
 end
end
v=m.sol('sol1').feature('v1');
v.set('scalemethod','none');
children=cell(v.feature.tags);
for k=1:numel(children)
 try,v.feature(children{k}).set('scalemethod','none');catch,end
end
end
