function result=r4_local_response(m,c,sel,folder)
% Round 4 post-hoc marginal/8-segment diagnostics of the CURRENT solved field.
% No connection, loading, solving, mesh/geometry changes, MPH saving or prediction.
started=tic;if ~isfolder(folder),mkdir(folder);end
D=row_ids(sel.inclusion_geometry_domains);
allIDs=double(sel.all_domains(:)');bgPhase=double(sel.domain_background_material(:)');
assert(numel(allIDs)==numel(bgPhase)&&numel(unique(allIDs))==numel(allIDs),'Incomplete domain/phase map');
assert(all(ismember(bgPhase,[0,1])),'Original material phases must be 0 or 1');
assert(~isempty(D)&&all(ismember(D,allIDs)),'Missing actual ellipse selection');
if isfield(sel,'removed_domains'),assert(isempty(sel.removed_domains),'Helper expects the complete ellipse');end
if isfield(sel,'compensation_domains'),assert(isempty(sel.compensation_domains),'Compensation is outside this diagnostic');end
t=double(m.param.evaluate('d_out'));DV=double(m.param.evaluate('DV'));
assert(t>0&&abs(DV)>0,'Invalid thickness or electrode voltage');
assert(c.sigma_frac_S_m>0,'The current uniform fracture conductivity must be positive');
intOrder=8;collar=3;nseg=8;
s='((x-cx)*cos(theta)+(y-cz)*sin(theta))';
tr='(-(x-cx)*sin(theta)+(y-cz)*cos(theta))';
edges=linspace(-c.a_mm,c.a_mm,nseg+1);masks={'1'};
for j=1:nseg
 if j==1,mask=sprintf('(%s<%.17g[mm])',s,edges(j+1));
 elseif j==nseg,mask=sprintf('(%s>=%.17g[mm])',s,edges(j));
 else,mask=sprintf('((%s>=%.17g[mm])*(%s<%.17g[mm]))',s,edges(j),s,edges(j+1));end
 masks{end+1}=mask; %#ok<AGROW>
end
e2='(ec.Ex^2+ec.Ey^2)';power='(ec.Jx*ec.Ex+ec.Jy*ec.Ey)';
domainExpr={'1','V','ec.Ex','ec.Ey','ec.Jx','ec.Jy',e2,power};
domainUnits={'mm^2','V*mm^2','V/m*mm^2','V/m*mm^2','A/m^2*mm^2','A/m^2*mm^2','V^2','W/m'};
wholeD=integrals(m,domainExpr,domainUnits,2,D,{'1'},intOrder);
phaseDomain=zeros(nseg+1,8,2);
for ph=0:1
 ids=intersect(D,allIDs(bgPhase==ph));
 phaseDomain(:,:,ph+1)=integrals(m,domainExpr,domainUnits,2,ids,masks,intOrder);
end
Sf=t*wholeD(7)/DV^2;
domainRows=struct([]);
for j=1:nseg
 for ph=0:1
  raw=phaseDomain(j+1,:,ph+1);area=raw(1);means=nan(1,5);
  if area>0,means=raw(2:6)/area;end
  rmsE=NaN;if area>0,rmsE=sqrt(raw(7)/(area*1e-6));end
  row=struct('segment_id',j,'s_lo_mm',edges(j),'s_hi_mm',edges(j+1), ...
   'original_phase',ph,'area_mm2',area,'V_mean_V',means(1), ...
   'Ex_mean_V_per_m',means(2),'Ez_mean_V_per_m',means(3), ...
   'Es_mean_V_per_m',means(2)*cosd(c.theta_deg)+means(3)*sind(c.theta_deg), ...
   'Et_mean_V_per_m',-means(2)*sind(c.theta_deg)+means(3)*cosd(c.theta_deg), ...
   'E_rms_V_per_m',rmsE,'Jx_mean_A_per_m2',means(4),'Jz_mean_A_per_m2',means(5), ...
   'Js_mean_A_per_m2',means(4)*cosd(c.theta_deg)+means(5)*sind(c.theta_deg), ...
   'Jt_mean_A_per_m2',-means(4)*sind(c.theta_deg)+means(5)*cosd(c.theta_deg), ...
   'integral_E2_V2',raw(7),'Sf_contribution_m',t*raw(7)/DV^2,'power_W',t*raw(8));
  domainRows=append_row(domainRows,row);
 end
end
writetable(struct2table(domainRows),fullfile(folder,'r4_local_response.csv'));
% Whole actual finite endcaps, independently selected before any averaging.
capIDs={row_ids(sel.capminus_domains),row_ids(sel.capplus_domains)};
capArea=zeros(1,2);capV=zeros(1,2);
for j=1:2
 q=integrals(m,{'1','V'},{'mm^2','V*mm^2'},2,capIDs{j},{'1'},intOrder);
 assert(q(1)>0,'Empty finite cap');capArea(j)=q(1);capV(j)=q(2)/q(1);
end
caps=struct('minus_area_mm2',capArea(1),'plus_area_mm2',capArea(2), ...
 'minus_mean_V',capV(1),'plus_mean_V',capV(2),'plus_minus_dV_V',diff(capV), ...
 'definition','Area average on actual stored finite cap selections; not tip-point values.');
% A fixed extended-ellipse window. This is NOT a constant-distance collar.
nearMask=sprintf('((%s/(a+%.17g[mm]))^2+(%s/(b+%.17g[mm]))^2<=1)',s,collar,tr,collar);
outerMask=['(1-(' nearMask '))'];
lowerStrip='((y>=20[mm])*(y<30[mm]))';upperStrip='((y>=90[mm])*(y<100[mm]))';
powerMasks={'1',nearMask,outerMask,lowerStrip,upperStrip};powerRegions=struct([]);
powerRegions=append_row(powerRegions,power_row('ellipse','main_partition',-1,wholeD(1),t*wholeD(8)));
matrixPower=zeros(5,2,2);
for ph=0:1
 ids=setdiff(allIDs(bgPhase==ph),D);
 q=integrals(m,{'1',power},{'mm^2','W/m'},2,ids,powerMasks,intOrder);
 matrixPower(:,:,ph+1)=q;
 phaseText={'low','high'};
 powerRegions=append_row(powerRegions,power_row(['near_' phaseText{ph+1}],'main_partition',ph,q(2,1),t*q(2,2)));
 powerRegions=append_row(powerRegions,power_row(['outer_' phaseText{ph+1}],'main_partition',ph,q(3,1),t*q(3,2)));
 powerRegions=append_row(powerRegions,power_row(['operable_lower_' phaseText{ph+1}],'overlapping_strip_diagnostic',ph,q(4,1),t*q(4,2)));
 powerRegions=append_row(powerRegions,power_row(['operable_upper_' phaseText{ph+1}],'overlapping_strip_diagnostic',ph,q(5,1),t*q(5,2)));
end
writetable(struct2table(powerRegions),fullfile(folder,'r4_power_regions.csv'));
totalPower=t*mphint2(m,power,2,'selection',allIDs,'intorder',intOrder,'unit','W/m','solnum','end');
partitionPower=t*(wholeD(8)+sum(reshape(matrixPower(2:3,2,:),[],1)));
% Arc selections use geometric domain adjacency and original MATRIX phase.
boundaries=[];
for k=1:numel(D),boundaries=union(boundaries,row_ids(mphgetadj(m,'geom1','boundary','domain',D(k))));end
boundaries=row_ids(boundaries);byPhase={[],[]};boundaryMap=zeros(0,3);
for k=1:numel(boundaries)
 b=boundaries(k);adj=row_ids(mphgetadj(m,'geom1','domain','boundary',b));outside=setdiff(adj,D);
 if isempty(outside),continue;end
 assert(numel(adj)==2&&numel(outside)==1,'Ellipse boundary must adjoin exactly one matrix domain');
 [known,index]=ismember(outside,allIDs);assert(known,'Unknown adjacent matrix domain');ph=bgPhase(index);
 assert(ismember(ph,[0,1]),'Unknown matrix phase');
 byPhase{ph+1}(end+1)=b;boundaryMap(end+1,:)=[b outside ph]; %#ok<AGROW>
end
assert(~isempty(boundaryMap),'No complete ellipse boundary found');
writetable(array2table(boundaryMap,'VariableNames',{'boundary_id','matrix_domain_id','original_matrix_phase'}), ...
 fullfile(folder,'r4_boundary_domain_map.csv'));
gx=['(' s '*cos(theta)/a^2-' tr '*sin(theta)/b^2)'];
gz=['(' s '*sin(theta)/a^2+' tr '*cos(theta)/b^2)'];
normFactor=['sqrt((' s '/a^2)^2+(' tr '/b^2)^2)'];
nx=['(' gx '/' normFactor ')'];nz=['(' gz '/' normFactor ')'];
ju=['(up(ec.Jx)*' nx '+up(ec.Jy)*' nz ')'];
jd=['(down(ec.Jx)*' nx '+down(ec.Jy)*' nz ')'];
jm=['(0.5*(' ju '+' jd '))'];jump=['(' ju '-' jd ')'];
jout=['max(' jm ',0[A/m^2])'];jin=['max(-' jm ',0[A/m^2])'];
meanV='(0.5*(up(V)+down(V)))';
arcExpr={'1',jm,['abs(' jm ')'],jout,jin,jump,['abs(' jump ')'], ...
 ['(' meanV '*' jout ')'],['(' meanV '*' jin ')'],meanV};
arcUnits={'mm','A/m','A/m','A/m','A/m','A/m','A/m','V*A/m','V*A/m','V*mm'};
arcData=zeros(nseg+1,10,2);arcRows=struct([]);
for ph=0:1
 arcData(:,:,ph+1)=integrals(m,arcExpr,arcUnits,1,row_ids(byPhase{ph+1}),masks,intOrder);
 for j=1:nseg
  q=arcData(j+1,:,ph+1);vout=safe_ratio(q(8),q(4));vin=safe_ratio(q(9),q(5));
  row=struct('segment_id',j,'s_lo_mm',edges(j),'s_hi_mm',edges(j+1), ...
   'original_matrix_phase',ph,'arc_length_mm',q(1),'net_out_A_per_m',q(2), ...
   'absolute_exchange_A_per_m',q(3),'outflow_A_per_m',q(4),'inflow_A_per_m',q(5), ...
   'outflow_A',t*q(4),'inflow_A',t*q(5),'trace_net_difference_A_per_m',q(6), ...
   'trace_absolute_difference_A_per_m',q(7), ...
   'trace_difference_relative_to_absolute_exchange',safe_ratio(q(7),q(3)), ...
   'current_weighted_outflow_V',vout,'current_weighted_inflow_V',vin, ...
   'weighted_in_minus_out_V',vin-vout,'arc_mean_V',safe_ratio(q(10),q(1)));
  arcRows=append_row(arcRows,row);
 end
end
writetable(struct2table(arcRows),fullfile(folder,'r4_boundary_exchange_8segments.csv'));
arcWhole=reshape(sum(arcData(1,:,:),3),1,10);
segmentSf=t*sum(reshape(phaseDomain(2:end,7,:),[],1))/DV^2;
segmentPower=t*sum(reshape(phaseDomain(2:end,8,:),[],1));
arcSegmentSum=reshape(sum(sum(arcData(2:end,:,:),1),3),1,10);
result=struct('status','success','case_id',c.id,'sigma_f_S_m',c.sigma_frac_S_m, ...
 'Sf_m',Sf,'marginal_response_m',Sf,'thickness_m',t,'DV_V',DV, ...
 'ellipse_area_mm2',wholeD(1),'ellipse_integral_E2_V2',wholeD(7), ...
 'ellipse_power_W',t*wholeD(8),'total_power_W',totalPower, ...
 'power_partition_residual_W',partitionPower-totalPower, ...
 'Sf_segment_sum_m',segmentSf,'Sf_segment_closure_residual_m',segmentSf-Sf, ...
 'ellipse_power_segment_closure_residual_W',segmentPower-t*wholeD(8), ...
 'finite_endcaps',caps,'segment_count',nseg,'original_phase_segment_rows',2*nseg, ...
 'whole_inflow_A',t*arcWhole(5),'whole_outflow_A',t*arcWhole(4), ...
 'whole_net_out_A',t*arcWhole(2),'whole_absolute_exchange_A',t*arcWhole(3), ...
 'whole_current_weighted_inflow_V',safe_ratio(arcWhole(9),arcWhole(5)), ...
 'whole_current_weighted_outflow_V',safe_ratio(arcWhole(8),arcWhole(4)), ...
 'arc_trace_absolute_difference_A_per_m',arcWhole(7), ...
 'arc_trace_difference_relative_to_absolute_exchange',safe_ratio(arcWhole(7),arcWhole(3)), ...
 'arc_net_relative_to_absolute_exchange',safe_ratio(abs(arcWhole(2)),arcWhole(3)), ...
 'arc_segment_length_closure_residual_mm',arcSegmentSum(1)-arcWhole(1), ...
 'arc_segment_current_max_abs_closure_residual_A_per_m',max(abs(arcSegmentSum(2:7)-arcWhole(2:7))), ...
 'inflow_outflow_identity_residual_A_per_m',arcWhole(4)-arcWhole(5)-arcWhole(2), ...
 'integration_order',intOrder,'elapsed_seconds',toc(started));
result.Sf_from_uniform_fracture_power_m=t*wholeD(8)/(c.sigma_frac_S_m*DV^2);
result.uniform_fracture_power_identity_relative_error=abs(result.Sf_from_uniform_fracture_power_m-Sf)/max(abs(Sf),1e-30);
result.near_window=struct('expanded_a_mm',c.a_mm+collar,'expanded_b_mm',c.b_mm+collar, ...
 'definition','Expanded ellipse (s/(a+3mm))^2+(t_local/(b+3mm))^2<=1, excluding actual D; NOT an equidistant 3mm collar.');
result.operable_strips=struct('lower_z_mm',[20 30],'upper_z_mm',[90 100], ...
 'definition','Whole-width strips outside the protected z30--90 core and the 20mm electrode buffers; overlap the main near/outer partition, so do not add their powers to that partition.');
result.units=struct('Sf','m = S/(S/m)','G','S','voltage','V','E','V/m','J','A/m^2', ...
 'area','mm^2','power','W including physical out-of-plane thickness','arc_exchange','both A/m per thickness and A including thickness');
result.interpretation='Sf uses the CURRENT solved target field and is a post-hoc derivative, not a background-only predictor. Piecewise Sf sums describe the local field contribution, not unique causal shares. Current magnitude alone is not pathway conductance.';
result.precision='Domain/arc windows use indicator integration at order 8 on the existing FE mesh. Segment/window edges need not follow elements. Closure is not quadrature or mesh convergence. Two-sided current-trace discrepancy is retained separately from global reaction-current conservation.';
write_json(fullfile(folder,'r4_local_response.json'),result);
end

function values=integrals(m,baseExpr,baseUnits,dim,ids,masks,intOrder)
values=zeros(numel(masks),numel(baseExpr));if isempty(ids),return;end
expressions=cell(1,numel(values));units=expressions;counter=0;
for j=1:numel(masks)
 for k=1:numel(baseExpr)
  counter=counter+1;expressions{counter}=['(' baseExpr{k} ')*(' masks{j} ')'];units{counter}=baseUnits{k};
 end
end
raw=cell(1,numel(expressions));
[raw{:}]=mphint2(m,expressions,dim,'selection',ids,'intorder',intOrder,'unit',units,'solnum','end');
raw=cellfun(@(x)double(x(1)),raw);
assert(all(isfinite(raw)),'Nonfinite integral in local response');
assert(max(abs(imag(raw)))<=1e-12*max(1,max(abs(raw))),'Unexpected complex steady DC integral');
values=reshape(real(raw),numel(baseExpr),numel(masks))';
end

function row=power_row(name,partition,ph,area,powerW)
row=struct('region',name,'partition_family',partition,'original_phase',ph,'area_mm2',area,'power_W',powerW);
end

function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end

function v=safe_ratio(numerator,denominator)
v=NaN;if denominator>0,v=numerator/denominator;end
end

function ids=row_ids(values)
ids=unique(double(values(:)'));ids=ids(ids>0);
end

function write_json(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0,'Cannot open %s',path);
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
