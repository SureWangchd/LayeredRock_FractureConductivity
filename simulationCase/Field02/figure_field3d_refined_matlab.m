function figure_field3d_refined_matlab(outdir,project_root)
% Real finite-ellipsoid fields from four saved COMSOL 3D solutions.
if nargin<2||isempty(project_root),project_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));end
if nargin<1||isempty(outdir),outdir=fullfile(project_root,'outputs','JGRSE_20261028_Field3D','first_version_refined');end
ip=fullfile(project_root,'outputs','JGRSE_20261028_Field3D','data');sp=fullfile(outdir,'source_notes');if ~isfolder(sp),mkdir(sp);end
source_tiles=fullfile(outdir,'source_tiles');if ~isfolder(source_tiles),mkdir(source_tiles);end
if isfile(fullfile(outdir,'plot_inputs','geometry.json')),ip=fullfile(outdir,'plot_inputs');end
cmapfile=fullfile(sp,'ling_L08_F02_stream_colormap.csv');
if ~isfile(cmapfile),copyfile(fullfile(project_root,'outputs','JGRSE_20261028_FieldSequence','source_notes','ling_L08_F02_stream_colormap.csv'),cmapfile);end
cmap=readmatrix(cmapfile);G=jsondecode(fileread(fullfile(ip,'geometry.json')));
M=readtable(fullfile(ip,'selected_four_states.csv'),'TextType','string');
E=readtable(fullfile(ip,'same_role_peripheral_effects.csv'),'TextType','string');
roles={'A','B'};layouts={'P77','P11'};
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman', ...
 'defaultAxesFontSize',7.5,'defaultTextFontSize',8,'defaultTextColor','black', ...
 'defaultAxesXColor','black','defaultAxesYColor','black','defaultAxesZColor','black', ...
 'defaultAxesLineWidth',.6,'defaultAxesBox','on','defaultAxesTickDir','in', ...
 'defaultAxesTickLength',[.018 .018],'defaultAxesTitleFontWeight','normal');
W=183;H=213;
f=figure('Visible','off','Color','white','Units','centimeters','Position',[2 2 W/10 H/10], ...
 'PaperUnits','centimeters','PaperSize',[W H]/10,'PaperPosition',[0 0 W H]/10, ...
 'Renderer','opengl','InvertHardcopy','off');
colormap(f,cmap);
data=cell(2,2);maxJ=0;minJ=Inf;diffs=cell(2,1);maxD=0;
for r=1:2
 for c=1:2
  name=sprintf('dense_%s_%s_sf10.mat',roles{r},layouts{c});
  data{r,c}=load(fullfile(ip,name));d=data{r,c};
  mag=sqrt(double(d.slice_JX).^2+double(d.slice_JY).^2+double(d.slice_JZ).^2);
  maxJ=max(maxJ,max(mag(:,:,2),[],'all'));minJ=min(minJ,min(mag(:,:,2),[],'all'));
 end
 diffs{r}=1000*(data{r,2}.slice_V(:,:,1)-data{r,1}.slice_V(:,:,1));
 q=data{r,1}.X_mm>=25&data{r,1}.X_mm<=95&data{r,1}.Z_mm>=30&data{r,1}.Z_mm<=90;
 maxD=max(maxD,max(abs(diffs{r}(q))));
end
% Shared limits cover all displayed samples; no per-panel rescaling.
Jlim=[.003 20];
if minJ<Jlim(1),Jlim(1)=10^floor(log10(minJ));end
if maxJ>Jlim(2),Jlim(2)=10^ceil(log10(maxJ));end
Dlim=ceil(maxD/5)*5;
fprintf('Displayed |J| %.8g..%.8g, shared log range %.8g..%.8g; |delta V|max %.8g mV.\n',minJ,maxJ,Jlim,maxD);
av=gobjects(1,6);pos=zeros(6,4);
xs=[7 98];ys=[140 74];
stream_stats=struct([]);
for r=1:2
 for c=1:2
  k=(r-1)*2+c;name=[roles{r} '_' layouts{c}];p=G.cases.(name);d=data{r,c};
  pos(k,:)=[xs(c) ys(r) 79 60];ax=axesmm(f,pos(k,:),W,H);av(k)=ax;
  draw_background(ax,p);
  J=sqrt(double(d.slice_JX(:,:,2)).^2+double(d.slice_JY(:,:,2)).^2+double(d.slice_JZ(:,:,2)).^2);
  surf(ax,d.X_mm,60+zeros(size(d.X_mm)),d.Z_mm,J, ...
   'FaceColor','interp','EdgeColor','none','FaceLighting','none');
  colormap(ax,cmap);set(ax,'ColorScale','log');clim(ax,Jlim);
  cache=fullfile(ip,['paths_' name '_sf10.mat']);
  if isfile(cache),q=load(cache);paths=q.paths;seed=q.seed;
  else
   [X,Y,Z]=meshgrid(d.x_mm,d.y_mm,d.z_mm);
   mag=sqrt(double(d.JX).^2+double(d.JY).^2+double(d.JZ).^2);
   [sx,sy]=meshgrid([35 45 55 65 75 85],[30 45 60 75 90]);sz=104.5+zeros(size(sx));
   seed=[sx(:),sy(:),sz(:)];
   paths=stream3(X,Y,Z,double(d.JX)./mag,double(d.JY)./mag,double(d.JZ)./mag, ...
    sx,sy,sz,[.15 10000]);
   save(cache,'paths','seed','-v7');
  end
  for n=1:numel(paths)
   q=paths{n};if size(q,1)<3,continue;end
   plot3(ax,q(:,1),q(:,2),q(:,3),'Color',[0 .447 .741],'LineWidth',.42);
   if mod(n,4)==1
    [~,ii]=min(abs(q(:,3)-98));jj=min(ii+12,size(q,1));v=q(jj,:)-q(ii,:);
    if norm(v)>1e-10,v=3*v/norm(v);quiver3(ax,q(ii,1),q(ii,2),q(ii,3),v(1),v(2),v(3),0, ...
     'Color',[0 .447 .741],'LineWidth',.42,'MaxHeadSize',.65);end
   end
  end
  draw_ellipsoid(ax,p);
  view(ax,115,20);ax.Projection='orthographic';axis(ax,'equal');axis(ax,'off');
  xlim(ax,[-10 165]);ylim(ax,[-20 145]);zlim(ax,[-20 138]);daspect(ax,[1 1 1]);
  % Coordinate directions are geometric axes; stream arrows indicate J.
  triad(ax,[132 -8 0],17);
  text(ax,65,127,-8,'120 mm','FontSize',7,'HorizontalAlignment','center','Clipping','off');
  camzoom(ax,1.04);
  t=M(M.role==roles{r}&M.layout==layouts{c},:);
  distances={'Far','Near'};
  note(f,sprintf('%s, %s',roles{r},layouts{c}),[xs(c)+15 ys(r)+58 50 5],W,H,8.7);
  note(f,sprintf('g = %.2f \\muS',t.g_uS),[xs(c)+20 ys(r)+1 50 4],W,H,7.3);
  panel(f,char('a'+k-1),xs(c)-2,ys(r)+64,W,H);
  ddy=cellfun(@(q)max(abs(q(:,2)-q(1,2))),paths);
  stream_stats(k).case=name;stream_stats(k).count=numel(paths);stream_stats(k).max_path_y_excursion_mm=max(ddy);
  stream_stats(k).seed_xyz_mm=seed;
 end
end
note(f,'\sigma_f = 10 S/m',[75 205 50 4],W,H,8);
cb=colorbar(av(4),'Location','southoutside');cb.Units='normalized';cb.Position=[58/W 68/H 69/W 1.8/H];
cb.Color='black';cb.FontSize=7;cb.LineWidth=.4;cb.TickDirection='out';
cb.Ticks=unique([Jlim(1) .03 .3 3 Jlim(2)]);cb.Label.String='';
note(f,'|J| (A m^{-2})',[23 66 32 4],W,H,7.5);
av(4).Position=pos(4,:)./[W H W H];
for r=1:2
 k=4+r;d=data{r,1};p=G.cases.([roles{r} '_P77']);
 pos(k,:)=[23+(r-1)*93 12 56 48];ax=axesmm(f,pos(k,:),W,H);av(k)=ax;
 imagesc(ax,d.X_mm(1,:),d.Z_mm(:,1),diffs{r});set(ax,'YDir','normal');
 colormap(ax,cmap);clim(ax,[-Dlim Dlim]);
 % Same unchanged material core; no displaced peripheral band enters here.
 co=p.basis_u(1);si=p.basis_u(3);ar=p.a_mm*sqrt(1-((45-p.cy_mm)/p.b_mm)^2);br=p.c_mm*sqrt(1-((45-p.cy_mm)/p.b_mm)^2);
 t=linspace(0,2*pi,601);ex=p.cx_mm+co*ar*cos(t)-si*br*sin(t);ez=p.cz_mm+si*ar*cos(t)+co*br*sin(t);
 plot(ax,ex,ez,'Color',[0 .5 0],'LineWidth',.8);
 for zz=p.physical_z_interfaces_mm(:).'
  if zz>30&&zz<90
   A=co^2/ar^2+si^2/br^2;B=2*co*si*(1/ar^2-1/br^2)*(zz-p.cz_mm);
   C=(si^2/ar^2+co^2/br^2)*(zz-p.cz_mm)^2-1;disc=B^2-4*A*C;
   if disc>0
    rr=p.cx_mm+sort([(-B-sqrt(disc))/(2*A),(-B+sqrt(disc))/(2*A)]);
    plot(ax,[25 rr(1)],[zz zz],'Color',[.4 .4 .4],'LineWidth',.3);
    plot(ax,[rr(2) 95],[zz zz],'Color',[.4 .4 .4],'LineWidth',.3);
   else,plot(ax,[25 95],[zz zz],'Color',[.4 .4 .4],'LineWidth',.3);end
  end
 end
 axis(ax,'equal');xlim(ax,[25 95]);ylim(ax,[30 90]);set(ax,'XTick',[30 60 90],'YTick',[30 60 90],'Layer','top');
 xlabel(ax,'x (mm)');ylabel(ax,'z (mm)');
 z=E(E.role==roles{r},:);
 title(ax,sprintf('%s: \\delta g = %+.3f \\muS',roles{r},z.delta_g_uS),'FontSize',8.4);
 panel(f,char('e'+r-1),5+(r-1)*93,62,W,H);
end
cb2=colorbar(av(6),'Location','westoutside');cb2.Units='normalized';cb2.Position=[95/W 12/H 2/W 48/H];
cb2.Color='black';cb2.FontSize=7;cb2.LineWidth=.4;cb2.Ticks=[-Dlim 0 Dlim];
cb2.Label.String='\deltaV (mV)';cb2.Label.FontSize=8;
av(6).Position=pos(6,:)./[W H W H];
drawnow;stem=fullfile(outdir,'Fig_3D_current_fields');savefig(f,[stem '.fig']);
% Rasterize only physical field/geometry viewports at 600 dpi. Main labels,
% chart axes and colourbar labels remain vector in the publication PDF.
for k=1:6
 ax=av(k);pp=pos(k,:);
 tmp=figure('Visible','off','Color','white','Units','centimeters','Position',[1 1 pp(3:4)/10], ...
  'PaperUnits','centimeters','PaperSize',pp(3:4)/10,'PaperPosition',[0 0 pp(3:4)/10],'Renderer','opengl');
 ca=copyobj(ax,tmp);colormap(ca,colormap(ax));ca.Position=[0 0 1 1];
 ca.Title.String='';ca.XLabel.String='';ca.YLabel.String='';ca.ZLabel.String='';axis(ca,'off');
 tile=fullfile(source_tiles,sprintf('panel_%c.png',char('a'+k-1)));
 if k<=4
  print(tmp,tile,'-dpng','-r600');
 else
  exportgraphics(tmp,tile,'Resolution',600,'BackgroundColor','white','Padding',0);
 end
 close(tmp);
 I=imread(tile);
 if k<=4
  delete(ax.Children);
  ar=axes(f,'Units','normalized','Position',pp./[W H W H]);
  image(ar,I);axis(ar,'image');axis(ar,'off');ar.Position=pp./[W H W H];
 else
  delete(ax.Children);image(ax,[25 95],[90 30],I);set(ax,'YDir','normal');
  axis(ax,'equal');xlim(ax,[25 95]);ylim(ax,[30 90]);
 end
 ax.Position=pp./[W H W H];
end
drawnow;
print(f,[stem '.png'],'-dpng','-r300');print(f,[stem '.pdf'],'-dpdf','-painters','-r600');
print(f,[stem '.svg'],'-dsvg','-painters','-r600');print(f,[stem '.tiff'],'-dtiff','-r600');close(f);
fid=fopen(fullfile(outdir,'render_metadata.json'),'w');fprintf(fid,'%s',jsonencode(struct('J_limits',Jlim,'deltaV_limit_mV',Dlim,'streams',stream_stats),PrettyPrint=true));fclose(fid);
helper=fullfile(fileparts(mfilename('fullpath')),'finish_field3d_paraview.py');
if isfile(helper),[status,msg]=system(sprintf('"%s" "%s" "%s"',repository_python(),helper,outdir));assert(status==0,'%s',msg);fprintf('%s',msg);end
fprintf('Saved %s\n',stem);
end

function draw_background(ax,p)
hold(ax,'on');
for k=1:numel(p.physical_material_runs)
 r=p.physical_material_runs(k);col=[1 1 1];if r.phase==1,col=[.82 .82 .82];end
 patch(ax,[12 12 12 12],[12 120 120 12],[r.z_min_mm r.z_min_mm r.z_max_mm r.z_max_mm],col, ...
  'EdgeColor',[.35 .35 .35],'LineWidth',.25,'FaceLighting','none');
 patch(ax,[12 120 120 12],[12 12 12 12],[r.z_min_mm r.z_min_mm r.z_max_mm r.z_max_mm],col*.97, ...
  'EdgeColor',[.35 .35 .35],'LineWidth',.25,'FaceLighting','none');
end
P=[0 0 0;120 0 0;120 120 0;0 120 0;0 0 120;120 0 120;120 120 120;0 120 120];
EE=[1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;1 5;2 6;4 8];
for k=1:size(EE,1),q=P(EE(k,:),:);plot3(ax,q(:,1),q(:,2),q(:,3),'Color',[.4 .4 .4],'LineWidth',.3);end
end
function draw_ellipsoid(ax,p)
[xx,yy,zz]=ellipsoid(0,0,0,p.a_mm,p.b_mm,p.c_mm,180);R=[p.basis_u(:) p.basis_v(:) p.basis_n(:)];
X=p.cx_mm+R(1,1)*xx+R(1,2)*yy+R(1,3)*zz;Y=p.cy_mm+R(2,1)*xx+R(2,2)*yy+R(2,3)*zz;Z=p.cz_mm+R(3,1)*xx+R(3,2)*yy+R(3,3)*zz;
surf(ax,X,Y,Z,'FaceColor',[0 .5 0],'EdgeColor','none','FaceAlpha',.24,'FaceLighting','gouraud', ...
 'AmbientStrength',.5,'DiffuseStrength',.6,'SpecularStrength',.25);
view(ax,115,20);light(ax,'Position',[180 150 120],'Style','local');
dv=(ax.CameraPosition-ax.CameraTarget)';dv=dv/norm(dv);abc=[p.a_mm;p.b_mm;p.c_mm];Q=null(((R'*dv)./abc)');
t=linspace(0,2*pi,401);q=[p.cx_mm;p.cy_mm;p.cz_mm]+R*(abc.*(Q*[cos(t);sin(t)]));
plot3(ax,q(1,:),q(2,:),q(3,:),'Color',[0 .35 0],'LineWidth',.65);
end
function triad(ax,o,L)
labs='xyz';
for k=1:3,d=zeros(1,3);d(k)=L;quiver3(ax,o(1),o(2),o(3),d(1),d(2),d(3),0,'k','LineWidth',.55,'MaxHeadSize',.35,'Clipping','off');
 q=o+1.2*d;text(ax,q(1),q(2),q(3),labs(k),'FontSize',7,'FontAngle','italic','Clipping','off');end
end
function ax=axesmm(f,p,W,H)
ax=axes(f,'Units','normalized','Position',p./[W H W H],'PositionConstraint','innerposition');hold(ax,'on');
end
function panel(f,s,x,y,W,H)
annotation(f,'textbox',[x/W (y-4)/H 9/W 5/H],'String',['(' s ')'],'EdgeColor','none','Margin',0, ...
 'FontName','Times New Roman','FontAngle','italic','FontSize',9,'Color','black');
end
function note(f,s,p,W,H,fs)
annotation(f,'textbox',p./[W H W H],'String',s,'HorizontalAlignment','center','VerticalAlignment','middle', ...
 'EdgeColor','none','Margin',0,'FontName','Times New Roman','FontSize',fs,'Color','black');
end

function executable=repository_python()
executable=getenv('JGRSE_PYTHON');
if isempty(executable),executable='python';end
end
