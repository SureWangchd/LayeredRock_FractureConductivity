function figure_field_sequence_matlab(outdir, project_root)
% Existing strict A/B COMSOL fields, displayed without a new solve.
% Original triangle traces are retained for colour and interpolation.
if nargin<2||isempty(project_root),project_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));end
if nargin<1||isempty(outdir),outdir=fullfile(project_root,'outputs','JGRSE_20261028_FieldSequence');end
if ~isfolder(outdir),mkdir(outdir);end
ip=fullfile(outdir,'data');
geom=jsondecode(fileread(fullfile(ip,'geometry.json')));
sp=fullfile(outdir,'source_notes');if ~isfolder(sp),mkdir(sp);end
cmapfile=fullfile(sp,'ling_L08_F02_stream_colormap.csv');
if ~isfile(cmapfile)
 copyfile(fullfile(project_root,'outputs','JGRSE_20261028_Fig02_RevisionACD','source_notes','ling_L08_F02_stream_colormap.csv'),cmapfile);
end
cmap=readmatrix(cmapfile);
W=183;H=188;
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman', ...
 'defaultAxesFontSize',7.5,'defaultTextFontSize',8,'defaultTextColor',[0 0 0], ...
 'defaultAxesXColor',[0 0 0],'defaultAxesYColor',[0 0 0], ...
 'defaultAxesLineWidth',.6,'defaultAxesBox','on','defaultAxesTickDir','in', ...
 'defaultAxesTickLength',[.018 .018],'defaultAxesTitleFontWeight','normal', ...
 'defaultAxesLabelFontSizeMultiplier',1.04,'defaultLegendFontName','Times New Roman', ...
 'defaultLegendFontSize',7,'defaultLegendTextColor',[0 0 0],'defaultLegendBox','off');
f=figure('Visible','off','Color','white','Units','centimeters', ...
 'Position',[2 2 W/10 H/10],'PaperUnits','centimeters','PaperSize',[W H]/10, ...
 'PaperPosition',[0 0 W H]/10,'Renderer','opengl','InvertHardcopy','off');
colormap(f,cmap);
roles={'A','B'};tokens={'0p1','0p159375','0p3'};sf=[.1 .159375 .3];
fieldaxes=gobjects(2,3);fieldpos=zeros(6,4);fieldmeta=struct([]);
xs=[17 74 131];ys=[124 64];pw=45;roi=[20 100];
% A shared physical scale is chosen once for all six views.
maxJ=0;
for r=1:2
 for c=1:3
  d=load(fullfile(ip,sprintf('state_%s_sf_%s.mat',roles{r},tokens{c})));
  q=all(d.p_mm>=roi(1)&d.p_mm<=roi(2),2);
  maxJ=max(maxJ,max(d.Jmag_A_m2(q)));
 end
end
upperJ=ceil(maxJ);lowerJ=.03;
fprintf('Shared logarithmic colour limits: %.9g to %.9g A/m2; data maximum %.9g.\n',lowerJ,upperJ,maxJ);
for r=1:2
 for c=1:3
  k=(r-1)*3+c;name=sprintf('state_%s_sf_%s',roles{r},tokens{c});
  d=load(fullfile(ip,[name '.mat']));
  roleGeom=geom.roles.(roles{r});par=roleGeom.ellipse;
  pos=[xs(c) ys(r) pw pw];fieldpos(k,:)=pos;
  ax=axesmm(f,pos,W,H);fieldaxes(r,c)=ax;
  p=d.p_mm;tt=double(d.triangles1);
  patch(ax,'Faces',tt,'Vertices',[p -ones(size(p,1),1)],'FaceVertexCData',d.Jmag_A_m2(:), ...
   'FaceColor','interp','EdgeColor','none');
  clim(ax,[lowerJ upperJ]);set(ax,'ColorScale','log');
  % Only boundaries separating different materials are outlined.
  interfacez=roleGeom.physical_conductivity_interfaces_z_mm;
  for z=interfacez(:).'
   if z>roi(1)&&z<roi(2)
    dz=z-par.cz_mm;co=cosd(par.theta_deg);si=sind(par.theta_deg);
    aa=co^2/par.a_mm^2+si^2/par.b_mm^2;
    bb=2*co*si*(1/par.a_mm^2-1/par.b_mm^2)*dz;
    cc=(si^2/par.a_mm^2+co^2/par.b_mm^2)*dz^2-1;
    disc=bb^2-4*aa*cc;
    if disc>0
     rr=par.cx_mm+sort([(-bb-sqrt(disc))/(2*aa),(-bb+sqrt(disc))/(2*aa)]);
     plot(ax,[roi(1) rr(1)],[z z],'-','Color',[.35 .35 .35],'LineWidth',.36);
     plot(ax,[rr(2) roi(2)],[z z],'-','Color',[.35 .35 .35],'LineWidth',.36);
    else
     plot(ax,roi,[z z],'-','Color',[.35 .35 .35],'LineWidth',.36);
    end
   end
  end
  cache=fullfile(ip,[name '_display_grid.mat']);
  if isfile(cache)
   g=load(cache);
  else
   % Merge coordinates only for the triangle locator. Original per-triangle
   % vertex values are used below: no averaging across material boundaries.
   tr=triangulation(double(d.locator_triangles1),d.locator_p_mm);
   [X,Z]=meshgrid(linspace(19,101,821));
   [ti,bc]=pointLocation(tr,[X(:) Z(:)]);valid=~isnan(ti);
   jx=nan(numel(X),1);jz=jx;
   orig=tt(ti(valid),:);
   jx(valid)=sum(bc(valid,:).*reshape(d.Jx_A_m2(orig),size(orig)),2);
   jz(valid)=sum(bc(valid,:).*reshape(d.Jz_A_m2(orig),size(orig)),2);
   JX=reshape(jx,size(X));JZ=reshape(jz,size(X));
   assert(all(isfinite(JX(:)))&&all(isfinite(JZ(:))),'Uncovered field sample; inspect saved triangle topology.');
   save(cache,'X','Z','JX','JZ','-v7');g=struct('X',X,'Z',Z,'JX',JX,'JZ',JZ);
  end
  speed=hypot(g.JX,g.JZ);
  % Identical seeds in each state. Lines show direction, not flux weights.
  seedx=linspace(21,99,21);seedz=99.8+zeros(size(seedx));
  lines=stream2(g.X,g.Z,g.JX./speed,g.JZ./speed,seedx,seedz,[.4 10000]);
  for n=1:numel(lines)
   xy=lines{n};if size(xy,1)<3,continue;end
   plot(ax,xy(:,1),xy(:,2),'-','Color',[.08 .08 .08],'LineWidth',.30);
   if mod(n,4)==1
    [~,ii]=min(abs(xy(:,2)-91));jj=min(ii+9,size(xy,1));
    v=xy(jj,:)-xy(ii,:);
    if norm(v)>1e-10
     v=2.0*v/norm(v);quiver(ax,xy(ii,1),xy(ii,2),v(1),v(2),0, ...
      'Color','black','LineWidth',.45,'MaxHeadSize',.9);
    end
   end
  end
  t=linspace(0,2*pi,721);
  ex=par.cx_mm+par.a_mm*cos(t)*cosd(par.theta_deg)-par.b_mm*sin(t)*sind(par.theta_deg);
  ez=par.cz_mm+par.a_mm*cos(t)*sind(par.theta_deg)+par.b_mm*sin(t)*cosd(par.theta_deg);
  plot(ax,ex,ez,'Color',[0 .5 0],'LineWidth',.85);
  if r==1 && c==1
   % The local s coordinate increases along the saved ellipse rotation.
   u=[cosd(par.theta_deg) sind(par.theta_deg)];normal=[-u(2) u(1)];
   arrow0=[par.cx_mm par.cz_mm]-18*u-12*normal;arrowv=14*u;
   quiver(ax,arrow0(1),arrow0(2),arrowv(1),arrowv(2),0, ...
    'Color','black','LineWidth',.75,'MaxHeadSize',.55);
   labelpos=arrow0+arrowv+[4 -3];
   text(ax,labelpos(1),labelpos(2),'+\its','FontSize',7.3,'Color','black', ...
    'BackgroundColor','white','Margin',.5,'HorizontalAlignment','left');
  end
  view(ax,2);axis(ax,'equal');xlim(ax,roi);ylim(ax,roi);set(ax,'XTick',[20 60 100],'YTick',[20 60 100],'Layer','top','SortMethod','childorder');
  if r==1,set(ax,'XTickLabel',[]);else,xlabel(ax,'x (mm)');end
  if c==1,ylabel(ax,'z (mm)');else,set(ax,'YTickLabel',[]);end
  if r==1,title(ax,sprintf('\\sigma_f = %.6g S/m',sf(c)),'FontSize',8.3);end
  text(ax,24,94,roles{r},'FontSize',9,'Color','black','BackgroundColor','white','Margin',1);
  panel(f,char('a'+k-1),xs(c)-10,ys(r)+pw+5,W,H);
  fieldmeta(k).role=roles{r};fieldmeta(k).sigma_f_S_m=sf(c);
  fieldmeta(k).data_max_J_A_m2=max(d.Jmag_A_m2);
  fieldmeta(k).colour_limits_A_m2=[lowerJ upperJ];fieldmeta(k).colour_scale='logarithmic';fieldmeta(k).streamlines=numel(lines);
  fieldmeta(k).positive_s_unit_vector=[cosd(par.theta_deg) sind(par.theta_deg)];
  fieldmeta(k).positive_s_arrow_shown=(r==1 && c==1);
 end
end
cb=colorbar(fieldaxes(2,3),'Location','southoutside');
cb.Units='normalized';cb.Position=[56/W 51/H 78/W 2.2/H];
cb.FontSize=7;cb.Color='black';cb.TickDirection='out';cb.LineWidth=.4;
cb.Ticks=[.03 .1 .3 1 3 5];cb.TickLabels={'0.03','0.1','0.3','1','3','5'};
cb.Label.String='';
annotation(f,'textbox',[24/W 49.5/H 29/W 5/H],'String','|J| (A m^{-2})', ...
 'HorizontalAlignment','right','EdgeColor','none','Margin',0,'FontName','Times New Roman','FontSize',8,'Color','black');
fieldaxes(2,3).Position=[xs(3) ys(2) pw pw]./[W H W H];

% The bottom row ties spatial fields to the measured ordering and local
% conductivity sensitivity. Only existing COMSOL data are used.
plot_response(f,ip,W,H);
drawnow;
stem=fullfile(outdir,'Fig_field_sequence_strict_AB');
savefig(f,[stem '.fig']);
% Keep the editable, actual triangular fields in FIG. The print copy uses
% high-resolution image tiles for just the six dense field backgrounds;
% labels, frames and quantitative charts remain vector in the PDF/SVG.
for r=1:2
 for c=1:3
  k=(r-1)*3+c;ax=fieldaxes(r,c);
  tmp=figure('Visible','off','Color','white','Units','centimeters', ...
   'Position',[1 1 pw/10 pw/10],'Renderer','opengl');
  ca=copyobj(ax,tmp);colormap(tmp,cmap);ca.Position=[0 0 1 1];
  ca.XLabel.String='';ca.YLabel.String='';ca.Title.String='';
  delete(findall(ca,'Type','text'));axis(ca,'off');
  tile=fullfile(sp,sprintf('field_%c_600dpi.png',char('a'+k-1)));
  exportgraphics(tmp,tile,'Resolution',600,'BackgroundColor','white','Padding',0);
  close(tmp);
  I=imread(tile);
  kids=ax.Children;delete(kids(~arrayfun(@(h)strcmp(h.Type,'text'),kids)));
  im=image(ax,roi,fliplr(roi),I);uistack(im,'bottom');
  set(ax,'XColor','black','YColor','black','YDir','normal');
  axis(ax,'equal');xlim(ax,roi);ylim(ax,roi);
  ax.Position=fieldpos(k,:)./[W H W H];
 end
end
drawnow;
print(f,[stem '.png'],'-dpng','-r300');
print(f,[stem '.pdf'],'-dpdf','-painters','-r600');
print(f,[stem '.svg'],'-dsvg','-painters','-r600');
print(f,[stem '.tiff'],'-dtiff','-r600');
close(f);
helper=fullfile(fileparts(mfilename('fullpath')),'finish_field_sequence.py');
if isfile(helper)
 [status,msg]=system(sprintf('"%s" "%s" "%s"',repository_python(),helper,outdir));
 assert(status==0,'PDF/SVG field embedding failed: %s',msg);fprintf('%s',msg);
end
fid=fopen(fullfile(outdir,'render_metadata.json'),'w');fprintf(fid,'%s',jsonencode(fieldmeta,PrettyPrint=true));fclose(fid);
fprintf('Saved %s\n',stem);
end

function plot_response(f,ip,W,H)
T=readtable(fullfile(ip,'original_D_scan.csv'));T.D_S=T.D_B_minus_A_S;
T=T(T.sigma_f_S_m>=.1 & T.sigma_f_S_m<=.3,:);
ax=axesmm(f,[20 12 63 26],W,H);
plot(ax,T.sigma_f_S_m,1e6*T.D_S,'-k','LineWidth',1.15);
yline(ax,0,':k','LineWidth',.6);
sf=[.1 .159375 .3];cols={[0 0 1],[0 0 0],[1 0 0]};marks={'o','s','^'};
for k=1:3
 [err,i]=min(abs(T.sigma_f_S_m-sf(k)));assert(err<1e-8);
 plot(ax,sf(k),1e6*T.D_S(i),'Marker',marks{k},'Color',cols{k}, ...
  'MarkerFaceColor','white','MarkerSize',4.2,'LineWidth',1.05);
end
xlim(ax,[.09 .31]);ylim(ax,[-50 140]);set(ax,'XTick',[.1 .2 .3],'YTick',[-50 0 50 100]);
xlabel(ax,'\sigma_f (S/m)');ylabel(ax,'\itg\rm_B - \itg\rm_A (\muS)');
title(ax,'Gain ordering','FontSize',8.4);panel(f,'g',5,43,W,H);
text(ax,.104,-42,'A > B','FontSize',7,'VerticalAlignment','bottom');
text(ax,.25,12,'B > A','FontSize',7,'VerticalAlignment','bottom');
P=readtable(fullfile(ip,'segment_marginal_B_minus_A.csv'));
ax=axesmm(f,[112 12 64 26],W,H);hs=gobjects(1,3);
for k=1:3
 Q=P(abs(P.sigma_f_S_m-sf(k))<1e-8,:);
 hs(k)=plot(ax,Q.s_center_mm,1000*Q.Sf_B_minus_A_m,['-' marks{k}], ...
  'Color',cols{k},'MarkerFaceColor','white','MarkerSize',3.3,'LineWidth',1.15);
end
yline(ax,0,':k','LineWidth',.5);xlim(ax,[-30 30]);set(ax,'XTick',[-30 0 30]);
xlabel(ax,'\its\rm (mm)');ylabel(ax,'\Delta\itS\rm_{f,i} (mm)');
title(ax,'Segment-wise sensitivity difference','FontSize',8.4);panel(f,'h',97,43,W,H);
lg=legend(ax,hs,{'0.1','0.159375','0.3'},'Location','southwest','NumColumns',3,'FontSize',6.6);
lg.ItemTokenSize=[9 6];
end

function a=axesmm(f,p,W,H)
a=axes(f,'Units','normalized','Position',p./[W H W H],'PositionConstraint','innerposition');hold(a,'on');
end
function panel(f,s,x,y,W,H)
annotation(f,'textbox',[x/W (y-4)/H 10/W 5/H],'String',['(' s ')'],'EdgeColor','none', ...
 'Margin',0,'FontName','Times New Roman','FontAngle','italic','FontSize',9,'Color','black');
end

function executable=repository_python()
executable=getenv('JGRSE_PYTHON');
if isempty(executable),executable='python';end
end
