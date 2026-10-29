function figure_ct_field_matlab(outdir,project_root)
% Final corrected localfit CT geometry and saved target potentials only.
% No COMSOL load, no field solve, no inferred current paths.
if nargin<2||isempty(project_root),project_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));end
if nargin<1||isempty(outdir),outdir=fullfile(project_root,'outputs','JGRSE_20261028_CTField');end
ip=fullfile(outdir,'plot_inputs');sp=fullfile(outdir,'source_notes');tp=fullfile(outdir,'source_tiles');
for p={outdir,ip,sp,tp},if ~isfolder(p{1}),mkdir(p{1});end,end
gi=fullfile(project_root,'outputs','JGRSE_20261028_Fig09_CT','inputs');
fi=fullfile(project_root,'outputs','JGRSE_20261028_CTField_feasibility');
samples={'S1','S5'};fluids=[.18 10];S=cell(1,2);D=cell(1,2);
for j=1:2
 for stem={'geometry_','voxel_display_'}
  name=[stem{1} samples{j} '.mat'];copy_if_absent(fullfile(gi,name),fullfile(ip,name));
 end
 name=[samples{j} '_M1_deltaV_H_minus_T_y050L.csv'];copy_if_absent(fullfile(fi,name),fullfile(ip,name));
 S{j}=load(fullfile(ip,['geometry_' samples{j} '.mat']));
 S{j}.voxel_display=load(fullfile(ip,['voxel_display_' samples{j} '.mat']));
 S{j}.boundary=section_boundary(S{j});
 T=readtable(fullfile(ip,name),'TextType','string');D{j}=grid_fields(T);
end
copy_if_absent(fullfile(fi,'core_deltaV_summary.csv'),fullfile(ip,'core_deltaV_summary.csv'));
copy_if_absent(fullfile(fi,'source_notes.json'),fullfile(ip,'original_slice_sources.json'));
copy_if_absent(fullfile(gi,'ct_main4_pairs.csv'),fullfile(ip,'ct_main4_pairs.csv'));
cmname='ling_L08_F02_stream_colormap.csv';
copy_if_absent(fullfile(project_root,'outputs','JGRSE_20261028_Field3D','source_notes',cmname),fullfile(ip,cmname));
cmap=readmatrix(fullfile(ip,cmname));
W=183;H=199;
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman', ...
 'defaultAxesFontSize',7.3,'defaultTextFontSize',8,'defaultTextColor','black', ...
 'defaultAxesXColor','black','defaultAxesYColor','black','defaultAxesZColor','black', ...
 'defaultAxesLineWidth',.6,'defaultAxesBox','on','defaultAxesTickDir','in', ...
 'defaultAxesTickLength',[.018 .018],'defaultAxesTitleFontWeight','normal', ...
 'defaultAxesLabelFontSizeMultiplier',1.05);
f=figure('Visible','off','Color','white','Units','centimeters','Position',[2 2 W/10 H/10], ...
 'PaperUnits','centimeters','PaperSize',[W H]/10,'PaperPosition',[0 0 W H]/10, ...
 'Renderer','opengl','InvertHardcopy','off');
colormap(f,cmap);av=gobjects(1,8);pos=zeros(8,4);geomcam=struct([]);
% a,b: original unsmoothed CT voxel walls + actual fitted M1 surface.
for j=1:2
 pos(j,:)=[15+(j-1)*90 139 72 51];ax=axesmm(f,pos(j,:),W,H);av(j)=ax;
 draw_ct(ax,S{j});
 % A plane outline is a locator only, not a colored field or current path.
 plot3(ax,[0 1 1 0 0],.5*ones(1,5),[0 0 1 1 0],'--k','LineWidth',.5);
 geomcam(j).sample=samples{j};geomcam(j).position=ax.CameraPosition;
 geomcam(j).target=ax.CameraTarget;geomcam(j).up=ax.CameraUpVector;
 note(f,sprintf('%s, \\sigma_w = %.2g S/m',samples{j},fluids(j)),[16+(j-1)*90 190 72 5],W,H,8.5,'center');
 note(f,sprintf('L = %.6f mm',S{j}.L_um/1000),[21+(j-1)*90 135 62 4],W,H,7.2,'center');
 panel(f,char('a'+j-1),5+(j-1)*90,195,W,H);
end
orientation_triad(f,av(1),[18 141],5.3,W,H);
% c,d: paired target potentials, no target-minus-background interpretation.
for j=1:2
 for r=1:2
  k=2+(j-1)*2+r;xx=17+(j-1)*90+(r-1)*37;
  pos(k,:)=[xx 86 32 32];ax=axesmm(f,pos(k,:),W,H);av(k)=ax;
  V=D{j}.VT;if r==2,V=D{j}.VH;end
  imagesc(ax,D{j}.x,D{j}.z,V);set(ax,'YDir','normal');colormap(ax,cmap);clim(ax,[0 1]);
  contour(ax,D{j}.x,D{j}.z,V,[.25 .5 .75],'LineColor','black','LineWidth',.38);
  draw_boundary(ax,S{j}.boundary,.65);
  if r==2
   for z=[.10 .25 .75 .90],plot(ax,[0 1],[z z],'--k','LineWidth',.35);end
  end
  axis(ax,'equal');xlim(ax,[0 1]);ylim(ax,[0 1]);set(ax,'XTick',[0 .5 1],'YTick',[0 .5 1],'Layer','top','FontSize',6.8);
  xlabel(ax,'x/L','FontSize',7.1);if r==1,ylabel(ax,'z/L','FontSize',7.1);else,set(ax,'YTickLabel',[]);end
  names={'T','H'};title(ax,names{r},'FontSize',8);
 end
 panel(f,char('c'+j-1),5+(j-1)*90,123,W,H);
end
cb=colorbar(av(6),'Location','southoutside');cb.Units='normalized';cb.Position=[59/W 72/H 67/W 1.9/H];
cb.Color='black';cb.FontSize=7;cb.LineWidth=.4;cb.Ticks=[0 .5 1];cb.TickDirection='out';cb.Label.String='';
note(f,'V (V)',[29 68.5 25 5],W,H,7.5,'right');av(6).Position=pos(6,:)./[W H W H];
% e,f: same retained core, full common +/-350 mV scale (no clipping).
for j=1:2
 k=6+j;pos(k,:)=[18+(j-1)*90 24 68 34];ax=axesmm(f,pos(k,:),W,H);av(k)=ax;
 imagesc(ax,D{j}.x,D{j}.z,D{j}.DV);set(ax,'YDir','normal');colormap(ax,cmap);clim(ax,[-350 350]);
 % A zero equipotential-difference contour, not a phase boundary/streamline.
 contour(ax,D{j}.x,D{j}.z,D{j}.DV,[0 0],'LineColor','black','LineStyle',':','LineWidth',.5);
 draw_boundary(ax,S{j}.boundary,.7);
 axis(ax,'equal');xlim(ax,[0 1]);ylim(ax,[.25 .75]);set(ax,'XTick',[0 .5 1],'YTick',[.25 .5 .75],'Layer','top');
 xlabel(ax,'x/L');ylabel(ax,'z/L');title(ax,'Retained core','FontSize',8);
 panel(f,char('e'+j-1),5+(j-1)*90,63,W,H);
end
cb2=colorbar(av(8),'Location','southoutside');cb2.Units='normalized';cb2.Position=[59/W 8/H 67/W 1.9/H];
cb2.Color='black';cb2.FontSize=7;cb2.LineWidth=.4;cb2.Ticks=[-350 0 350];cb2.TickDirection='out';cb2.Label.String='';
note(f,'\deltaV (mV)',[21 4.5 33 5],W,H,7.5,'right');av(8).Position=pos(8,:)./[W H W H];
drawnow;stem=fullfile(outdir,'Fig_CT_physical_fields');savefig(f,[stem '.fig']);
% Preserve high-resolution field/geometry tiles, vector outer axes and text.
% Tile labels a-h are raster indices: they do not redefine figure panels a-f.
tilemap={'a: S1 CT','b: S5 CT','c-left: S1 T','c-right: S1 H', ...
 'd-left: S5 T','d-right: S5 H','e: S1 target difference','f: S5 target difference'};
for k=1:8
 ax=av(k);pp=pos(k,:);
 tmp=figure('Visible','off','Color','white','Units','centimeters','Position',[1 1 pp(3:4)/10], ...
  'PaperUnits','centimeters','PaperSize',pp(3:4)/10,'PaperPosition',[0 0 pp(3:4)/10], ...
  'Renderer','opengl','InvertHardcopy','off');
 ca=copyobj(ax,tmp);colormap(ca,colormap(ax));ca.Position=[0 0 1 1];ca.PositionConstraint='innerposition';
 ca.Title.String='';ca.XLabel.String='';ca.YLabel.String='';ca.ZLabel.String='';axis(ca,'off');
 tile=fullfile(tp,sprintf('panel_%c.png',char('a'+k-1)));
 if k<=2,print(tmp,tile,'-dpng','-r600');else,exportgraphics(tmp,tile,'Resolution',600,'BackgroundColor','white','Padding',0);end
 close(tmp);I=imread(tile);
 if k<=2
  delete(ax);ar=axes(f,'Units','normalized','Position',pp./[W H W H]);image(ar,I);axis(ar,'image');axis(ar,'off');ar.Position=pp./[W H W H];
 else
  xr=ax.XLim;yr=ax.YLim;delete(ax.Children);image(ax,xr,fliplr(yr),I);set(ax,'YDir','normal');
  axis(ax,'equal');xlim(ax,xr);ylim(ax,yr);ax.Position=pp./[W H W H];
 end
end
drawnow;print(f,[stem '.png'],'-dpng','-r300');print(f,[stem '.pdf'],'-dpdf','-painters','-r600');
print(f,[stem '.svg'],'-dsvg','-painters','-r600');print(f,[stem '.tiff'],'-dtiff','-r600');close(f);
meta=struct('figure_mm',[W H],'source','Final corrected R7B localfit M1; existing target slices only', ...
 'samples',{samples},'fluid_sigma_S_m',fluids,'L_mm',[S{1}.L_um S{2}.L_um]/1000, ...
 'potential_limits_V',[0 1],'difference_limits_mV',[-350 350],'difference_definition','V_H,target - V_T,target', ...
 'slice_y_over_L',.5,'retained_core_z_over_L',[.25 .75],'homogenized_remote_z_over_L',[.10 .25;.75 .90], ...
 'actual_surface_face_counts',[size(S{1}.surface_faces1,1),size(S{2}.surface_faces1,1)], ...
 'surface_decimation',false,'original_voxel_labels_smoothed',false,'new_solver_calls',0, ...
 'field_grid_size',[129 129],'field_interpolation','Only graphics interpolation and contour extraction on existing samples; no field resampling or smoothing', ...
 'boundary_source','Free boundary of the saved actual M1 fracture-section polygons, not ideal CAD ellipse', ...
 'top_geometry_display','Original skeleton voxel walls for x_index<28 or y_index<28, original categorical bottom plane; complete actual M1 fracture surface', ...
 'zero_contour_meaning','deltaV=0; not a pore boundary and not a current path','tile_bindings',{tilemap},'cameras',geomcam);
fid=fopen(fullfile(outdir,'render_metadata.json'),'w');fprintf(fid,'%s',jsonencode(meta,PrettyPrint=true));fclose(fid);
helper=fullfile(fileparts(mfilename('fullpath')),'finish_field3d_paraview.py');
if isfile(helper)
 [status,msg]=system(sprintf('"%s" "%s" "%s" --stem Fig_CT_physical_fields --panels abcdefgh',repository_python(),helper,outdir));
 assert(status==0,'High-resolution image embedding failed: %s',msg);fprintf('%s',msg);
end
fprintf('Saved %s\n',stem);
end

function copy_if_absent(src,dst)
if ~isfile(dst),assert(isfile(src),'Missing source %s',src);copyfile(src,dst);end
end
function D=grid_fields(T)
x=unique(T.x_over_L);z=unique(T.z_over_L);assert(numel(x)==129&&numel(z)==129&&height(T)==129^2);
[~,ix]=ismember(T.x_over_L,x);[~,iz]=ismember(T.z_over_L,z);idx=sub2ind([129 129],iz,ix);assert(numel(unique(idx))==129^2);
D=struct('x',x,'z',z,'VT',zeros(129),'VH',zeros(129),'DV',zeros(129));
D.VT(idx)=T.V_V_T;D.VH(idx)=T.V_V_H;D.DV(idx)=T.deltaV_mV;
assert(max(abs(1000*(D.VH-D.VT)-D.DV),[],'all')<1e-8);
core=z>=.25&z<=.75;assert(max(abs(D.DV(core,:)),[],'all')<350);
end
function B=section_boundary(s)
X=s.slice_inside_x_um/s.L_um;Z=s.slice_inside_z_um/s.L_um;v=[];e=[];
for k=1:size(X,1)
 q=[X(k,:)' Z(k,:)'];q=q(all(isfinite(q),2),:);n=size(q,1);offset=size(v,1);
 v=[v;q];e=[e;offset+[(1:n)' [2:n 1]']]; %#ok<AGROW>
end
[u,~,iv]=unique(round(v,12),'rows');e=sort(iv(e),2);[ue,~,ic]=unique(e,'rows');counts=accumarray(ic,1);e=ue(counts==1,:);
degree=accumarray(e(:),1,[size(u,1),1]);assert(all(degree(degree>0)==2),'Actual slice boundary not closed.');
B=struct('vertices',u,'edges',e);
end
function draw_boundary(ax,B,lw)
x=[B.vertices(B.edges(:,1),1),B.vertices(B.edges(:,2),1),nan(size(B.edges,1),1)]';
z=[B.vertices(B.edges(:,1),2),B.vertices(B.edges(:,2),2),nan(size(B.edges,1),1)]';
plot(ax,x(:),z(:),'Color',[0 .45 0],'LineWidth',lw);
end
function draw_ct(ax,s)
hold(ax,'on');axis(ax,'equal');axis(ax,'off');view(ax,115,20);ax.Projection='orthographic';
lab=s.labels_xyz;n=128;k=28;c=k/n;cg=(k:n)/n;
rgb=repmat([.753 .753 .753],size(s.voxel_display.faces1,1),1);rgb(s.voxel_display.face_kind(:)==1,:)=.98;
patch(ax,'Vertices',s.voxel_display.vertices_normalized,'Faces',double(s.voxel_display.faces1), ...
 'FaceVertexCData',rgb,'FaceColor','flat','EdgeColor','none','FaceLighting','flat', ...
 'AmbientStrength',.8,'DiffuseStrength',.3,'SpecularStrength',.03,'BackFaceLighting','reverselit');
[X,Y]=meshgrid(cg,cg);phase=squeeze(lab(k+1:end,k+1:end,1))';g=ones(size(phase));g(phase==1)=.753;g=.96*g;g=g([1:end end],[1:end end]);
surf(ax,X,Y,0*X,repmat(g,1,1,3),'FaceColor','flat','EdgeColor','none','FaceLighting','none');
patch(ax,'Vertices',s.surface_vertices_um/s.L_um,'Faces',double(s.surface_faces1), ...
 'FaceColor',[0 .5 0],'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',.18, ...
 'DiffuseStrength',.72,'SpecularStrength',.75,'SpecularExponent',24,'SpecularColorReflectance',.65,'BackFaceLighting','reverselit');
light(ax,'Position',[1.083 -.958 -.833],'Style','local','Color',[1 1 1]);light(ax,'Position',[1.5 1.25 1],'Style','local','Color',[.28 .28 .28]);
P=[0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];E=[1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;1 5;2 6;4 8];
for j=1:size(E,1),t=P(E(j,:),:);plot3(ax,t(:,1),t(:,2),t(:,3),'Color',[.2 .2 .2],'LineWidth',.4);end
plot3(ax,[c c],[c c],[0 1],'k','LineWidth',.35);plot3(ax,[c c],[1 1],[0 1],'k','LineWidth',.35);plot3(ax,[1 1],[c c],[0 1],'k','LineWidth',.35);
plot3(ax,[0 c c 1],[1 1 c c],[1 1 1 1],'k','LineWidth',.35);
xlim(ax,[-.06 1.13]);ylim(ax,[-.06 1.13]);zlim(ax,[-.07 1.10]);daspect(ax,[1 1 1]);
end
function orientation_triad(f,ax,origin,len,W,H)
forward=ax.CameraTarget-ax.CameraPosition;forward=forward/norm(forward);right=cross(forward,ax.CameraUpVector);right=right/norm(right);up=cross(right,forward);up=up/norm(up);names={'x','y','z'};
for k=1:3,d=[right(k) up(k)];tip=origin+len*d;
 annotation(f,'arrow',[origin(1) tip(1)]/W,[origin(2) tip(2)]/H,'Color','black','LineWidth',.6,'HeadLength',3,'HeadWidth',2.5);
 p=tip+1.6*d/norm(d);note(f,['\it' names{k}],[p(1)-1.5 p(2)-1.5 3 3],W,H,7,'center');
end
end
function ax=axesmm(f,p,W,H)
ax=axes(f,'Units','normalized','Position',p./[W H W H],'PositionConstraint','innerposition');hold(ax,'on');box(ax,'on');ax.Layer='top';
end
function panel(f,s,x,y,W,H)
annotation(f,'textbox',[x/W (y-3.5)/H 12/W 5/H],'String',['(' s ')'],'EdgeColor','none','Margin',0,'FontSize',9,'FontAngle','italic','FontName','Times New Roman','Color','black','VerticalAlignment','middle');
end
function note(f,s,p,W,H,fs,align)
annotation(f,'textbox',p./[W H W H],'String',s,'EdgeColor','none','Margin',0,'FontName','Times New Roman','FontSize',fs,'Color','black','HorizontalAlignment',align,'VerticalAlignment','middle');
end

function executable=repository_python()
executable=getenv('JGRSE_PYTHON');
if isempty(executable),executable='python';end
end
