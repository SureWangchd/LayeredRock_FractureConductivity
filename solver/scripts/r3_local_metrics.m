function r=r3_local_metrics(model,c,sel,folder)
% Postprocess the passed, already solved model only. No load/connect/solve/save.
% The geometric ellipse selection works for both passive and conducting states.
started=tic;if ~isfolder(folder),mkdir(folder);end
main=row_ids(sel.inclusion_geometry_domains);assert(~isempty(main),'Missing geometric ellipse domains');
assert(isempty(sel.removed_domains)&&isempty(sel.compensation_domains), ...
 'Round 3 local metrics initially support complete ellipses only');
minus=row_ids(sel.capminus_domains);plus=row_ids(sel.capplus_domains);
assert(~isempty(minus)&&~isempty(plus),'Missing finite endcap selections');
[am,vm]=mphint2(model,{'1','V'},2,'selection',minus,'intorder',8, ...
 'unit',{'mm^2','V*mm^2'},'solnum','end');
[ap,vp]=mphint2(model,{'1','V'},2,'selection',plus,'intorder',8, ...
 'unit',{'mm^2','V*mm^2'},'solnum','end');
assert(am>0&&ap>0,'Invalid endcap areas');
endpoint=struct('capminus_area_mm2',am,'capplus_area_mm2',ap, ...
 'capminus_V',vm/am,'capplus_V',vp/ap,'dV_V',vp/ap-vm/am, ...
 'definition','Actual finite endcap area averages; no point-tip substitution.');
write_json(fullfile(folder,'r3_endpoint_potential_averages.json'),endpoint);
% Sample through the major-axis interior. Near-tip exclusion avoids an
% arbitrary boundary-side choice exactly at a geometric vertex.
s=linspace(-c.a_mm*(1-1e-7),c.a_mm*(1-1e-7),801);
xy=[c.cx_mm+s*cosd(c.theta_deg);c.cz_mm+s*sind(c.theta_deg)];
[v,jx,jz,ex,ez]=mphinterp(model,{'V','ec.Jx','ec.Jy','ec.Ex','ec.Ey'}, ...
 'coord',xy,'edim',2,'selection',main,'ext',0,'solnum','end', ...
 'unit',{'V','A/m^2','A/m^2','V/m','V/m'});
js=jx*cosd(c.theta_deg)+jz*sind(c.theta_deg);
jt=-jx*sind(c.theta_deg)+jz*cosd(c.theta_deg);
axisT=table(s(:),xy(1,:)',xy(2,:)',v(:),jx(:),jz(:),js(:),jt(:),ex(:),ez(:), ...
 'VariableNames',{'s_mm','x_mm','z_mm','V_V','Jx_A_m2','Jz_A_m2','Js_A_m2','Jt_A_m2','Ex_V_m','Ez_V_m'});
writetable(axisT,fullfile(folder,'r3_ellipse_axis_fields.csv'));
% Identify all geometric ellipse/exterior-matrix boundaries. The material
% assigned inside the ellipse is irrelevant to this adjacency selection.
allD=row_ids(sel.all_domains);backgroundPhase=double(sel.domain_background_material(:)');
assert(numel(allD)==numel(backgroundPhase),'Background material labels are incomplete');
bnds=[];for k=1:numel(main),bnds=union(bnds,row_ids(mphgetadj(model,'geom1','boundary','domain',main(k))));end
bnds=row_ids(bnds);byphase={[],[]};boundaryRows=zeros(0,3);
for k=1:numel(bnds)
 b=bnds(k);adj=row_ids(mphgetadj(model,'geom1','domain','boundary',b));
 outside=setdiff(adj,main);
 if isempty(outside),continue;end
 assert(numel(adj)==2&&numel(outside)==1,'Unexpected ellipse boundary adjacency');
 [exists,index]=ismember(outside,allD);assert(exists,'Unknown exterior matrix domain');
 ph=backgroundPhase(index);assert(ismember(ph,[0,1]),'Invalid matrix phase');
 byphase{ph+1}(end+1)=b;boundaryRows(end+1,:)=[b outside ph]; %#ok<AGROW>
end
assert(~isempty(boundaryRows),'No ellipse/matrix boundary was identified');
writetable(array2table(boundaryRows,'VariableNames',{'boundary_id','matrix_domain_id','original_matrix_phase'}), ...
 fullfile(folder,'r3_ellipse_arc_domains.csv'));
% Analytic outward normal from the ellipse into the matrix. The same normal
% is used with both current traces; their difference measures discretization.
ss='((x-cx)*cos(theta)+(y-cz)*sin(theta))';
tt='(-(x-cx)*sin(theta)+(y-cz)*cos(theta))';
gx=['(' ss '*cos(theta)/a^2-' tt '*sin(theta)/b^2)'];
gz=['(' ss '*sin(theta)/a^2+' tt '*cos(theta)/b^2)'];
den=['sqrt((' ss '/a^2)^2+(' tt '/b^2)^2)'];
nx=['(' gx '/' den ')'];nz=['(' gz '/' den ')'];
ju=['(up(ec.Jx)*' nx '+up(ec.Jy)*' nz ')'];
jd=['(down(ec.Jx)*' nx '+down(ec.Jy)*' nz ')'];
jm=['(0.5*(' ju '+' jd '))'];jump=['(' ju '-' jd ')'];
baseExpr={'1',jm,['abs(' jm ')'],jump,['abs(' jump ')']};
baseUnits={'mm','A/m','A/m','A/m','A/m'};
edges=linspace(-c.a_mm,c.a_mm,17);expr=baseExpr;units=baseUnits;
for seg=1:16
 % Extend the two outside endpoints: geometric coordinate rounding must not
 % omit a sliver of the closed ellipse. Interior bounds partition exactly.
 if seg==1,mask=sprintf('(%s<%.17g[mm])',ss,edges(seg+1));
 elseif seg==16,mask=sprintf('(%s>=%.17g[mm])',ss,edges(seg));
 else,mask=sprintf('((%s>=%.17g[mm])*(%s<%.17g[mm]))',ss,edges(seg),ss,edges(seg+1));end
 for j=1:5,expr{end+1}=['(' baseExpr{j} ')*' mask];units{end+1}=baseUnits{j};end %#ok<AGROW>
end
whole=zeros(2,5);segments=zeros(16,2,5);
for ph=0:1
 ids=row_ids(byphase{ph+1});if isempty(ids),continue;end
 values=cell(1,numel(expr));
 [values{:}]=mphint2(model,expr,1,'selection',ids,'intorder',8, ...
  'unit',units,'solnum','end');
 values=cellfun(@(v)double(v(1)),values);
 assert(all(isfinite(values)),'Nonfinite ellipse boundary metric');
 whole(ph+1,:)=values(1:5);
 for seg=1:16,segments(seg,ph+1,:)=reshape(values(5+(seg-1)*5+(1:5)),1,1,5);end
end
tableRows=zeros(32,8);counter=0;
for seg=1:16
 for ph=0:1
  counter=counter+1;tableRows(counter,:)=[seg edges(seg) edges(seg+1) ph reshape(segments(seg,ph+1,2:5),1,4)];
 end
end
segmentT=array2table(tableRows,'VariableNames',{'segment_id','s_lo_mm','s_hi_mm','original_matrix_phase', ...
 'net_out_A_per_m','absolute_exchange_A_per_m','trace_net_difference_A_per_m','trace_absolute_difference_A_per_m'});
lengths=zeros(32,1);counter=0;
for seg=1:16,for ph=0:1,counter=counter+1;lengths(counter)=segments(seg,ph+1,1);end,end
segmentT.arc_length_mm=lengths;
writetable(segmentT,fullfile(folder,'r3_ellipse_exchange_16_segments.csv'));
totalMetrics=sum(whole,1);segTotals=reshape(sum(segments,1),2,5);
r=struct('status','success','endcap',endpoint,'axis_samples',numel(s), ...
 'axis_nonfinite_value_count',sum(~isfinite(v(:))), ...
 'segments',16,'phase_segment_rows',32,'selected_arc_count',size(boundaryRows,1), ...
 'low_net_out_A_per_m',whole(1,2),'high_net_out_A_per_m',whole(2,2), ...
 'low_absolute_exchange_A_per_m',whole(1,3),'high_absolute_exchange_A_per_m',whole(2,3), ...
 'total_net_out_A_per_m',totalMetrics(2),'total_absolute_exchange_A_per_m',totalMetrics(3), ...
 'total_trace_absolute_difference_A_per_m',totalMetrics(5), ...
 'trace_difference_relative_to_absolute_exchange',totalMetrics(5)/max(totalMetrics(3),eps), ...
 'net_to_absolute_exchange',abs(totalMetrics(2))/max(totalMetrics(3),eps), ...
 'segment_partition_length_max_abs_residual_mm',max(abs(segTotals(:,1)-whole(:,1))), ...
 'segment_partition_current_max_abs_residual_A_per_m',max(abs(segTotals(:,2:5)-whole(:,2:5)),[],'all'), ...
 'phase_whole_integrals',whole,'integration_order',8,'elapsed_seconds',toc(started));
r.phase_whole_integral_columns={'length_mm','net_A_per_m','absolute_A_per_m','trace_net_difference_A_per_m','trace_absolute_difference_A_per_m'};
r.normal_definition='Outward from the geometric ellipse, common analytic normal for up/down traces.';
r.trace_definition='Mean of two finite-element gradient current traces. Their difference is a local discretization diagnostic, not an electrode-reaction error.';
r.scope='Complete geometric ellipse boundary, including passive backgrounds; a passive partition is not a physical conductivity interface.';
r.units='V in V; axis coordinates and lengths in mm; J in A/m^2; exchange in A/m of out-of-plane thickness.';
r.segmentation_precision='16 indicator-weighted axial intervals on the existing mesh, quadrature order 8. Interval cuts need not coincide with boundary elements; local segment integrals have quadrature/trace error. Partition closure is not a convergence proof.';
r.interpretation='Closed ellipse net flux should approach zero with refinement. Good global reaction-current conservation does not imply equally precise local gradient-trace integrals. Axis samples are local current density, not total through-current.';
if isfield(c,'id'),r.case_id=c.id;end
write_json(fullfile(folder,'r3_local_metrics.json'),r);
end

function ids=row_ids(values)
ids=unique(double(values(:)'));ids=ids(ids>0);
end

function write_json(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0,'Cannot open %s',path);
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
