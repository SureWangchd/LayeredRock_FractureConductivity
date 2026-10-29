function figure09_ct_matlab(outdir,project_root)
% Existing corrected localfit CT results only. No COMSOL call or new solve.
if nargin<2||isempty(project_root),project_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));end
if nargin<1||isempty(outdir),outdir=fullfile(project_root,'outputs','JGRSE_20261028_Fig09_FinalLabels');end
ip=fullfile(outdir,'inputs');
T=readtable(fullfile(ip,'ct_main8_targets.csv'),'TextType','string','Delimiter',',');
D=readtable(fullfile(ip,'ct_preselected6_mesh_pairs.csv'),'TextType','string','Delimiter',',');
assert(height(T)==8&&height(D)==6);assert(max(abs(T.eta_percent-100*T.g_S./T.G0_S))<1e-8);
S1=load(fullfile(ip,'geometry_S1.mat'));S5=load(fullfile(ip,'geometry_S5.mat'));
S1.voxel_display=load(fullfile(ip,'voxel_display_S1.mat'));
S5.voxel_display=load(fullfile(ip,'voxel_display_S5.mat'));
section_phase=squeeze(S1.labels_xyz(:,65,:))';
section_g=ones(size(section_phase));section_g(section_phase==1)=.753;
imwrite(flipud(repmat(section_g,1,1,3)),fullfile(outdir,'source_notes','ct_S1_section_native128.png'));
W=183;H=195;blue=[0 0 1];red=[1 0 0];
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman', ...
 'defaultAxesFontSize',8,'defaultTextFontSize',8,'defaultTextColor',[0 0 0], ...
 'defaultAxesXColor',[0 0 0],'defaultAxesYColor',[0 0 0],'defaultAxesZColor',[0 0 0], ...
 'defaultAxesLineWidth',.6,'defaultAxesBox','on','defaultAxesTickDir','in', ...
 'defaultAxesTickLength',[.018 .018],'defaultAxesTitleFontWeight','normal', ...
 'defaultAxesLabelFontSizeMultiplier',1.05,'defaultLegendFontName','Times New Roman', ...
 'defaultLegendTextColor',[0 0 0],'defaultLegendBox','off');
f=figure('Visible','off','Color','white','Units','centimeters', ...
 'Position',[2 2 W/10 H/10],'PaperUnits','centimeters','PaperSize',[W H]/10, ...
 'PaperPosition',[0 0 W H]/10,'Renderer','opengl');

% Native categorical CT on actual voxel planes, actual exported M1 surface.
ga=axesmm(f,[16 135 72 51],W,H);draw_ct(ga,S1);
gb=axesmm(f,[106 135 72 51],W,H);draw_ct(gb,S5);
orientation_triad(f,ga,[20 137],5.5,W,H);
note(f,'S1',[20 185 65 5],W,H,8.5,'center');
note(f,'S5',[110 185 65 5],W,H,8.5,'center');
note(f,'\itL\rm = 1.111 mm',[24 132 58 5],W,H,7.2,'center');
note(f,'\itL\rm = 0.512 mm',[114 132 58 5],W,H,7.2,'center');
panel(f,'a',5,190,W,H);panel(f,'b',95,190,W,H);

% Actual CT trace at y=L/2, positive-coordinate-side convention.
note(f,'Peripheral simplification',[17 120 75 5],W,H,8.5,'center');
ca=axesmm(f,[20 84 29 29],W,H);draw_section(ca,S1,false);
cb=axesmm(f,[58 84 29 29],W,H);draw_section(cb,S1,true);
title(ca,'T','Color','black','FontSize',8);title(cb,'H','Color','black','FontSize',8);
set(ca,'XTick',[0 .5 1],'YTick',[0 .10 .25 .75 .90 1], ...
 'YTickLabel',{'0','0.10','0.25','0.75','0.90','1'},'FontSize',6.6);
set(cb,'XTick',[0 .5 1],'YTick',[],'FontSize',6.6);
xlabel(ca,'\itx\rm / \itL','FontSize',7);ylabel(ca,'\itz\rm / \itL','FontSize',7);
xlabel(cb,'\itx\rm / \itL','FontSize',7);
note(f,'Background mismatch < 3.3 \times 10^{-9}',[15 73 77 4],W,H,6.6,'center');
panel(f,'c',5,126,W,H);

% All four M1 pairs. Each gain uses that variant's own calibrated G0.
ad=axesmm(f,[110 83 67 35],W,H);hh=gobjects(2,1);
for k=1:4
 a=T(T.plot_group_order==k&T.role=="T",:);b=T(T.plot_group_order==k&T.role=="H",:);
 assert(height(a)==1&&height(b)==1);
 plot(ad,[k-.09 k+.09],[a.eta_percent b.eta_percent],'-k','LineWidth',.55);
 hh(1)=plot(ad,k-.09,a.eta_percent,'o','Color',blue,'MarkerFaceColor','white','MarkerSize',4,'LineWidth',1.15);
 hh(2)=plot(ad,k+.09,b.eta_percent,'^','Color',red,'MarkerFaceColor','white','MarkerSize',4.2,'LineWidth',1.15);
end
xlim(ad,[.6 4.4]);ylim(ad,[0 16]);set(ad,'XTick',1:4, ...
 'XTickLabel',{'S1 / 0.18','S1 / 10','S5 / 0.18','S5 / 10'},'YTick',[0 4 8 12 16],'FontSize',7);
xlabel(ad,'Sample / \sigma_w (S/m)','FontSize',7.5);ylabel(ad,'Conductance gain, \itg\rm / \itG\rm_0 (%)','FontSize',8);
title(ad,'Fracture-induced gain','FontSize',8.5);
lg=legend(ad,hh,{'T','H'},'NumColumns',2,'Location','northeast','FontSize',7,'AutoUpdate','off');
if isprop(lg,'ItemTokenSize'),lg.ItemTokenSize=[9 7];end
panel(f,'d',95,126,W,H);

% Fixed actual M1 interface; M2/M3 refine the FE mesh only.
Q=table();
for j=1:2
 sam="S1";sf=.18;if j==2,sam="S5";sf=10;end
 z=sortrows(D(D.sample==sam&D.water_sigma_S_m==sf,:),'mesh_order');assert(height(z)==3);
 ax=axesmm(f,[20+(j-1)*90 19 66 39],W,H);
 yyT=100*z.gT_S/z.gT_S(1);yyH=100*z.gH_S/z.gH_S(1);
 h1=plot(ax,1:3,yyT,'-o','Color',blue,'MarkerFaceColor','white','MarkerSize',4,'LineWidth',1.15);
 h2=plot(ax,1:3,yyH,'--^','Color',red,'MarkerFaceColor','white','MarkerSize',4.2,'LineWidth',1.15);
 xlim(ax,[.8 3.2]);ylim(ax,[99.65 100.04]);set(ax,'XTick',1:3,'XTickLabel',{'M1','M2','M3'},'YTick',[99.7 99.8 99.9 100]);
 ytickformat(ax,'%.1f');xlabel(ax,'Mesh');ylabel(ax,'\itg\rm / \itg\rm_{M1} (%)');
 title(ax,sprintf('%s, \\sigma_w = %.2g S/m',sam,sf),'FontSize',8.5);
 lg=legend(ax,[h1 h2],{'T','H'},'Location','northeast','NumColumns',2,'FontSize',6.8,'AutoUpdate','off');
 if isprop(lg,'ItemTokenSize'),lg.ItemTokenSize=[11 7];end
 ai=axesmm(f,[31+(j-1)*90 23 27 13],W,H);
 plot(ai,1:3,z.epsilon_percent,'-dk','MarkerFaceColor','white','MarkerSize',2.4,'LineWidth',.75);
 xlim(ai,[.7 3.3]);set(ai,'XTick',1:3,'XTickLabel',{'1','2','3'},'FontSize',6.1);
 if j==1,ylim(ai,[-25.09 -25.04]);set(ai,'YTick',[-25.08 -25.04]);else,ylim(ai,[3.185 3.205]);set(ai,'YTick',[3.19 3.20]);end
 ytickformat(ai,'%.2f');
 title(ai,'\epsilon (%)','FontSize',6.5);
 % M3 pair difference: use the actual M3 row, not the M1 main-batch value.
 text(ax,.975,.055,sprintf('\\epsilon_{M3} = %+.2f%%',z.epsilon_percent(end)), ...
  'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom', ...
  'FontSize',6.7,'Color','black','Interpreter','tex');
 z.gT_over_M1_percent=yyT;z.gH_over_M1_percent=yyH;Q=[Q;z]; %#ok<AGROW>
end
panel(f,'e',5,65,W,H);panel(f,'f',95,65,W,H);
writetable(T,fullfile(outdir,'plotted_main8.csv'));writetable(Q,fullfile(outdir,'plotted_mesh6.csv'));

drawnow;stem=fullfile(outdir,'Fig09_CT_digital_rock');savefig(f,[stem '.fig']);
% Keep native editable 3D in FIG, and rasterize only the two spatial panels
% at 600 dpi to avoid the vector painter's 3D visibility errors.
anns=[findall(f,'Type','textboxshape');findall(f,'Tag','fig09_triad_arrow')];set(anns,'Visible','off');
canvas=fullfile(outdir,'source_notes','geometry_render_canvas.png');print(f,canvas,'-dpng','-r600');set(anns,'Visible','on');
im=imread(canvas);[nh,nw,~]=size(im);
for j=1:2
 if j==1,ax=ga;else,ax=gb;end
 pos=ax.Position;
 x1=max(1,floor(pos(1)*nw)+1);x2=min(nw,ceil((pos(1)+pos(3))*nw));
 y1=max(1,floor((1-pos(2)-pos(4))*nh)+1);y2=min(nh,ceil((1-pos(2))*nh));
 crop=im(y1:y2,x1:x2,:);imwrite(crop,fullfile(outdir,'source_notes',sprintf('geometry_%d_600dpi.png',j)));
 delete(ax);ar=axes(f,'Units','normalized','Position',pos);image(ar,crop);axis(ar,'off');set(ar,'Position',pos,'PositionConstraint','innerposition');
end
exportgraphics(f,[stem '.png'],'Resolution',300,'BackgroundColor','white');
print(f,[stem '.pdf'],'-dpdf','-painters','-r600');print(f,[stem '.tiff'],'-dtiff','-r600');close(f);
helper=fullfile(fileparts(mfilename('fullpath')),'finish_fig09_pdf.py');
if isfile(helper)
 python=getenv('FIG09_PYTHON');
 if isempty(python)
  python=getenv('JGRSE_PYTHON');if isempty(python),python='python';end
 end
 [status,msg]=system(sprintf('"%s" "%s" "%s"',python,helper,outdir));
 assert(status==0,'PDF high-resolution geometry embedding failed: %s',msg);fprintf('%s',msg);
end
fprintf('Saved %s\n',stem);
end

function draw_ct(ax,s)
hold(ax,'on');axis(ax,'equal');axis(ax,'off');view(ax,115,20);ax.Projection='orthographic';
lab=s.labels_xyz; n=128; k=28;c=k/n;cg=(k:n)/n;
% Exact native-voxel faces of retained skeleton slabs; pores are real holes.
% Flat face lighting changes only display shading, not the original labels.
face_rgb=repmat([.753 .753 .753],size(s.voxel_display.faces1,1),1);
face_rgb(s.voxel_display.face_kind(:)==1,:)=.98;
patch(ax,'Vertices',s.voxel_display.vertices_normalized, ...
 'Faces',double(s.voxel_display.faces1),'FaceVertexCData',face_rgb,'FaceColor','flat', ...
 'EdgeColor','none','FaceLighting','flat','AmbientStrength',.8, ...
 'DiffuseStrength',.3,'SpecularStrength',.03,'BackFaceLighting','reverselit');
% The exposed floor is a categorical material section at the original scale.
[X,Y]=meshgrid(cg,cg);tex(ax,X,Y,0*X,squeeze(lab(k+1:end,k+1:end,1))',.96);
v=s.surface_vertices_um/s.L_um;
patch(ax,'Vertices',v,'Faces',double(s.surface_faces1),'FaceColor',[0 .5 0], ...
 'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',.18,'DiffuseStrength',.72, ...
 'SpecularStrength',.75,'SpecularExponent',24,'SpecularColorReflectance',.65,'BackFaceLighting','reverselit');
light(ax,'Position',[1.083 -.958 -.833],'Style','local','Color',[1 1 1]);
light(ax,'Position',[1.5 1.25 1],'Style','local','Color',[.28 .28 .28]);
P=[0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
E=[1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;1 5;2 6;4 8];
for j=1:size(E,1),t=P(E(j,:),:);plot3(ax,t(:,1),t(:,2),t(:,3),'Color',[.2 .2 .2],'LineWidth',.4);end
plot3(ax,[c c],[c c],[0 1],'k','LineWidth',.35);
plot3(ax,[c c],[1 1],[0 1],'k','LineWidth',.35);plot3(ax,[1 1],[c c],[0 1],'k','LineWidth',.35);
plot3(ax,[0 c c 1],[1 1 c c],[1 1 1 1],'k','LineWidth',.35);
xlim(ax,[-.06 1.13]);ylim(ax,[-.06 1.13]);zlim(ax,[-.07 1.10]);daspect(ax,[1 1 1]);
end
function orientation_triad(f,ax,origin,len,W,H)
forward=ax.CameraTarget-ax.CameraPosition;forward=forward/norm(forward);
right=cross(forward,ax.CameraUpVector);right=right/norm(right);
up=cross(right,forward);up=up/norm(up);
names={'x','y','z'};
for k=1:3
 d=[right(k) up(k)];tip=origin+len*d;
 annotation(f,'arrow',[origin(1) tip(1)]/W,[origin(2) tip(2)]/H, ...
  'Color','black','LineWidth',.65,'HeadLength',3,'HeadWidth',2.5,'Tag','fig09_triad_arrow');
 p=tip+1.6*d/norm(d);
 note(f,['\it' names{k}],[p(1)-1.5 p(2)-1.5 3 3],W,H,7.5,'center');
end
end
function tex(ax,X,Y,Z,phase,shade)
g=ones(size(phase));g(phase==1)=.753;g=min(g*shade,1);
g=g([1:end end],[1:end end]);rgb=repmat(g,1,1,3);
surf(ax,X,Y,Z,rgb,'FaceColor','flat','EdgeColor','none','FaceLighting','none');
end
function draw_section(ax,s,isH)
n=128;
phase=squeeze(s.labels_xyz(:,65,:))';g=ones(size(phase));g(phase==1)=.753;
% Native 128-by-128 categorical image avoids PDF stitching between voxel
% quads. Pixel centers and outer edges match the original voxel grid exactly.
image(ax,'XData',[.5/n 1-.5/n],'YData',[.5/n 1-.5/n],'CData',repmat(g,1,1,3));
view(ax,2);set(ax,'YDir','normal');axis(ax,'equal');xlim(ax,[0 1]);ylim(ax,[0 1]);
if isH
 for zz=[.10 .75]
  patch(ax,[0 1 1 0],[zz zz zz+.15 zz+.15],ones(1,4)*.01,[.9 .9 .9],'EdgeColor','none');
  text(ax,.5,zz+.075,.04,'\kappa','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8,'Color','black');
 end
end
xx=s.slice_inside_x_um'/s.L_um;zz=s.slice_inside_z_um'/s.L_um;
% Exported sections currently contain triangles; NaN-padded fourth entries
% must not be passed as vertices, or MATLAB may omit these filled faces.
nv=sum(isfinite(xx),1);
for q=3:4
 sel=nv==q;
 if any(sel)
  patch(ax,xx(1:q,sel),zz(1:q,sel),ones(q,sum(sel))*.02,[0 .5 0],'EdgeColor','none');
 end
end
for v=[.1 .25 .75 .9],plot3(ax,[0 1],[v v],[.03 .03],'--k','LineWidth',.45);end
plot3(ax,[0 1],[0 0],[.03 .03],'-k','LineWidth',1.1);plot3(ax,[0 1],[1 1],[.03 .03],'-k','LineWidth',1.1);
ax.Layer='top';
end
function ax=axesmm(f,mm,W,H)
ax=axes(f,'Units','normalized','Position',mm./[W H W H]);hold(ax,'on');box(ax,'on');ax.PositionConstraint='innerposition';ax.Layer='top';
end
function panel(f,s,x,y,W,H)
annotation(f,'textbox',[x/W (y-3.5)/H 12/W 5/H],'String',['(' s ')'], ...
 'EdgeColor','none','Margin',0,'FontSize',9,'FontAngle','italic','FontName','Times New Roman','Color','black','VerticalAlignment','middle');
end
function note(f,s,mm,W,H,fs,align)
annotation(f,'textbox',mm./[W H W H],'String',s,'EdgeColor','none','Margin',0, ...
 'FontName','Times New Roman','FontSize',fs,'Color','black','HorizontalAlignment',align,'VerticalAlignment','middle');
end
