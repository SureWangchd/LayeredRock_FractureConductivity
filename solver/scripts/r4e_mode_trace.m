function result=r4e_mode_trace(m,c,sel,folder,nmax)
% Read-only posthoc cosine DtN diagnostic of CURRENT target field.
% No geometry/mesh/physics/solver modification; no prediction or target solve.
if nargin<5,nmax=32;end
assert(nmax>=32 && nmax==floor(nmax),'Use at least modes 0..32');
assert(all(abs(c.interface_amplitudes_mm)<1e-14),'Flat peripherals required');
if ~isfolder(folder),mkdir(folder);end
L=c.L_mm*1e-3;t=double(m.param.evaluate('d_out'));DV=double(m.param.evaluate('DV'));
assert(abs(L-.12)<1e-12 && t>0 && DV>0,'Unexpected geometry or loading');
edges=[0 reshape(c.z_interfaces_mm,1,[]) c.L_mm]*1e-3;
phase=reshape(c.layer_materials,1,[]);sigma=c.sigma_low_S_m+phase*(c.sigma_high_S_m-c.sigma_low_S_m);
h=diff(edges);Rtotal=sum(h./sigma);J0=DV/Rtotal; % sigma*dV0/dz, positive; E0z=-J0/sigma
ids=double(sel.all_domains(:)');bgLayer=double(sel.domain_background_layer(:)');
assert(numel(ids)==numel(bgLayer),'Domain layer mapping required');
rows=struct([]);sideResults=struct([]);n=(0:nmax)';intOrder=8;
for sideIndex=1:2
 if sideIndex==1,side='lower';z=.03;layerIndices=find(edges(2:end)<=z+1e-12);
 else,side='upper';z=.09;layerIndices=find(edges(1:end-1)>=z-1e-12);end
 assert(abs(sum(h(layerIndices))-.03)<1e-12,'Core boundary must be a retained partition');
 outerIds=ids(ismember(bgLayer,layerIndices));
 assert(isempty(intersect(outerIds,double(sel.inclusion_geometry_domains(:)'))),'Ellipse intersects peripheral rectangle');
 % Geometry input coordinates and mphselectbox use mm; integrals explicitly use SI.
 tol=1e-6;zmm=z*1e3;
 boundary=double(mphselectbox(m,'geom1',[-tol c.L_mm+tol;zmm-tol zmm+tol],'boundary'));
 assert(~isempty(boundary),'Core interface not found');
 length_m=mphint2(m,'1',1,'selection',boundary,'unit','m','intorder',intOrder,'solnum','end');
 assert(abs(length_m-L)<1e-8,'Selected interface is not one full-width horizontal boundary');
 usedH=max(0,min(edges(2:end),z)-edges(1:end-1));V0=J0*sum(usedH./sigma);
 expressions=cell(1,nmax+1);units=repmat({'V*m'},1,nmax+1);
 for k=0:nmax
  expressions{k+1}=sprintf('(V-(%.17g[V]))*cos(%d*pi*x/(%.17g[m]))',V0,k,L);
 end
 raw=cell(1,nmax+1);
 [raw{:}]=mphint2(m,expressions,1,'selection',boundary,'unit',units,'intorder',intOrder,'solnum','end');
 an=cellfun(@(v)double(v(1)),raw(:))/L;an(2:end)=2*an(2:end);
 assert(all(isfinite(an)),'Nonfinite mode coefficient');
 % Propagate impedance from electrode (delta V=0) towards the core.
 ordered=layerIndices;if sideIndex==2,ordered=fliplr(ordered);end
 Yn=zeros(nmax+1,1);Yn(1)=1/sum(h(ordered)./sigma(ordered));
 for k=1:nmax
  wave=k*pi/L;Z=0;
  for j=ordered
   beta=sigma(j)*wave;q=tanh(wave*h(j));Z=(Z+q/beta)/(1+beta*Z*q);
  end
  Yn(k+1)=1/Z;
 end
 norm=repmat(.5,nmax+1,1);norm(1)=1;modePower=t*L*norm.*Yn.*an.^2;
 directPower=0;
 for ph=0:1
  layerSet=layerIndices(phase(layerIndices)==ph);selected=ids(ismember(bgLayer,layerSet));
  if isempty(selected),continue;end
  sig=c.sigma_low_S_m+ph*(c.sigma_high_S_m-c.sigma_low_S_m);
  expr=sprintf('(%.17g[S/m])*(ec.Ex^2+(ec.Ey+(%.17g[V/m]))^2)',sig,J0/sig);
  directPower=directPower+t*mphint2(m,expr,2,'selection',selected,'unit','W/m','intorder',intOrder,'solnum','end');
 end
 for k=1:nmax+1
  row=struct('side',side,'n',n(k),'k_per_m',n(k)*pi/L,'Y_S_per_m2',Yn(k),'delta_V_cosine_V',an(k),...
   'mode_perturbation_power_W',modePower(k),'cumulative_perturbation_power_W',sum(modePower(1:k)));
  if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
 end
 p4=sum(modePower(1:5));p8=sum(modePower(1:9));p16=sum(modePower(1:17));p32=sum(modePower(1:33));
 summary=struct('side',side,'core_interface_z_mm',zmm,'core_interface_boundary_ids',boundary,...
  'outer_domain_ids',outerIds,'interface_length_m',length_m,'V0_interface_V',V0,'Y0_S_per_m2',Yn(1),...
  'delta_V_mean_V',an(1),'mode0_power_W',modePower(1),'n1to32_power_W',sum(modePower(2:33)),...
  'modes_0to4_power_W',p4,'modes_0to8_power_W',p8,'modes_0to16_power_W',p16,'modes_0to32_power_W',p32,...
  'direct_outer_perturbation_power_W',directPower,'mode32_minus_direct_W',p32-directPower,...
  'mode32_direct_relative_difference',abs(p32-directPower)/max(directPower,1e-30),...
  'mode16_to32_increment_relative_to_direct',(p32-p16)/max(directPower,1e-30),...
  'mode32_fraction_from_nge1',sum(modePower(2:33))/max(p32,1e-30));
 if isempty(sideResults),sideResults=summary;else,sideResults(end+1)=summary;end %#ok<AGROW>
end
writetable(struct2table(rows),fullfile(folder,'mode_coefficients.csv'));
result=struct('status','posthoc_from_saved_COMSOL_target','case_id',c.id,'sigma_f_S_m',c.sigma_frac_S_m,...
 'nmax',nmax,'width_m',L,'thickness_m',t,'DV_V',DV,'background_sigma_dVdz_A_per_m2',J0,...
 'sides',sideResults,'integration_order',intOrder,...
 'energy_definition','Integral over each peripheral rectangle of sigma*|grad(V-V0)|^2, including out-of-plane thickness; this is not the change in raw peripheral Joule dissipation.',...
 'mode_convention','a0=mean(deltaV), an=2/L*integral(deltaV*cos(n*pi*x/L)); P=t*L*(Y0*a0^2+0.5*sum(Yn*an^2)).',...
 'interpretation','Y depends only on the layered periphery; amplitudes use full target field and are posthoc. This is not a new predictor and does not assign unique causal shares of the whole gain difference.',...
 'precision','Uses continuous potential traces, not boundary J. Report n16/n32 truncation and direct volume-energy closure; they do not replace mesh convergence.');
write_json(fullfile(folder,'mode_summary.json'),result);
end

function write_json(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0,'Cannot open %s',path);
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
