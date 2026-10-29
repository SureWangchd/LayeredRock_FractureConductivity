function figure08_finite3d_matlab(outdir,project_root)
% Existing synthetic 3D COMSOL results only. No simulation is launched.
if nargin<2||isempty(project_root),project_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));end
if nargin<1||isempty(outdir),outdir=fullfile(project_root,'outputs','JGRSE_20261028_Fig08_Notation');end
ip=fullfile(outdir,'inputs');
T=readtable(fullfile(ip,'synthetic3d_response.csv'),'TextType','string','Delimiter',',');
T=T(T.kind=="synthetic3d_explicit",:);assert(height(T)==16);
C=readtable(fullfile(ip,'controlled_effects.csv'),'TextType','string','Delimiter',',');
D=readtable(fullfile(ip,'fine_joint_comparison.csv'),'TextType','string','Delimiter',',');
P=jsondecode(fileread(fullfile(ip,'r7_strict_B_P11_parameters.json')));
assert(max(abs(T.eta-T.g_S./T.G0_S))<1e-12);
W=183;H=145;blue=[0 0 1];red=[1 0 0];roles=["A","B"];
cols={blue,red};marks={'o','^'};layouts=["P11","P17","P71","P77"];
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

% (a) Exact saved geometry with a TGRS Fig1(b)-style matrix display cutaway.
ag=axesmm(f,[15 79 73 54],W,H);draw_geometry(ag,P);
note(f,'Finite ellipsoid',[27 131 52 5],W,H,8.5,'center');
note(f,'B, P11',[20 126 23 4],W,H,7,'left');
note(f,'60 \times 60 \times 3 mm, 60\circ',[25 77 53 5],W,H,7.1,'center');
panel(f,'a',5,136,W,H);

% (b) All sixteen targets: main = 10 S/m; inset = 0.18 S/m.
ab=axesmm(f,[110 88 66 43],W,H);hb=gobjects(1,2);
for j=1:2
 vals=zeros(1,4);
 for k=1:4,z=T(T.role==roles(j)&T.layout==layouts(k)&T.sigma_f_S_m==10,:);assert(height(z)==1);vals(k)=100*z.eta;end
 hb(j)=plot(ab,1:4,vals,['-' marks{j}],'Color',cols{j},'MarkerFaceColor','white', ...
  'MarkerSize',4,'LineWidth',1.15);
end
xlim(ab,[.7 4.3]);ylim(ab,[0 14.5]);set(ab,'XTick',1:4,'XTickLabel',{'(1,1)','(1,7)','(7,1)','(7,7)'},'YTick',[0 5 10]);
xlabel(ab,'(\itd\rm_-, \itd\rm_+) (mm)');ylabel(ab,'Conductance gain, \itg\rm / \itG\rm_0 (%)');
title(ab,'Fracture-induced gain','FontSize',8.5);
text(ab,4.14,11.6,'\sigma_f = 10 S/m','HorizontalAlignment','right','FontSize',7.0);
lg=legend(ab,hb,{'A','B'},'Location','northeast','NumColumns',2,'FontSize',7,'AutoUpdate','off');
lg.Units='normalized';lg.Position=[112/W 128/H 23/W 3/H];
if isprop(lg,'ItemTokenSize'),lg.ItemTokenSize=[10 7];end
bi=axesmm(f,[122 92 42 12.5],W,H);
for j=1:2
 vals=zeros(1,4);for k=1:4,z=T(T.role==roles(j)&T.layout==layouts(k)&T.sigma_f_S_m==.18,:);vals(k)=100*z.eta;end
 plot(bi,1:4,vals,['-' marks{j}],'Color',cols{j},'MarkerFaceColor','white','MarkerSize',2.6,'LineWidth',.9);
end
xlim(bi,[.7 4.3]);ylim(bi,[2.33 2.51]);set(bi,'XTick',1:4,'XTickLabel',{'(1,1)','(1,7)','(7,1)','(7,7)'}, ...
 'YTick',[2.35 2.45],'FontSize',6.1);
title(bi,'\sigma_f = 0.18 S/m','FontSize',6.7);
panel(f,'b',95,136,W,H);writetable(T,fullfile(outdir,'plotted_all16_gains.csv'));

% (c) Single-side and joint changes have a common, explicitly stated baseline.
ac=axesmm(f,[20 17 66 44],W,H);plot(ac,[.65 3.35],[0 0],':k','LineWidth',.55);
names=["lower_near_minus_far_upper_far","upper_near_minus_far_lower_far","joint_near_minus_far"];
hs=gobjects(1,4);labs=cell(1,4);q=0;Pc=table();
for j=1:2
 for sf=[.18 10]
  yy=zeros(1,3);
  for k=1:3
   z=C(C.role==roles(j)&C.sigma_f_S_m==sf&C.effect_id==names(k),:);assert(height(z)==1);
   yy(k)=100*z.relative_to_reference_gain;Pc=[Pc;z]; %#ok<AGROW>
  end
  q=q+1;if sf==.18,ls='--';face='white';else,ls='-';face=cols{j};end
  hs(q)=plot(ac,1:3,yy,'LineStyle',ls,'Marker',marks{j},'Color',cols{j}, ...
   'MarkerFaceColor',face,'MarkerSize',4,'LineWidth',1.15);
  labs{q}=sprintf('%s, %.2g',roles(j),sf);
 end
end
xlim(ac,[.65 3.35]);ylim(ac,[-2.4 4.2]);set(ac,'XTick',1:3,'XTickLabel',{'Lower','Upper','Both'},'YTick',[-2 0 2 4]);
ylabel(ac,'\delta\itg\rm / \itg\rm_{77} (%)');title(ac,'Peripheral relocation','FontSize',8.5);
text(ac,1.08,.27,'*','FontSize',9,'Color','black');
lg=legend(ac,hs,labs,'NumColumns',2,'Location','northwest','FontSize',6.8,'AutoUpdate','off');
title(lg,'Core, \sigma_f (S/m)','FontWeight','normal','FontSize',6.8);
if isprop(lg,'ItemTokenSize'),lg.ItemTokenSize=[13 7];end
panel(f,'c',5,66,W,H);writetable(Pc,fullfile(outdir,'plotted_peripheral_effects.csv'));

% (d) Only two preselected joint contrasts have direct two-mesh comparisons.
note(f,'Mesh comparison',[112 61 65 5],W,H,8.5,'center');
for j=1:2
 sf=.18;if j==2,sf=10;end
 z=D(D.role==roles(j)&D.sigma_f_S_m==sf,:);assert(height(z)==1);
 ad=axesmm(f,[112+(j-1)*39 20 25 37],W,H);
 yy=1e6*[z.controlled_gain_difference_S,z.fine_controlled_gain_difference_S];
 plot(ad,1:2,yy,['-' marks{j}],'Color',cols{j},'MarkerFaceColor','white','MarkerSize',4,'LineWidth',1.15);
 xlim(ad,[.75 2.25]);set(ad,'XTick',[1 2],'XTickLabel',{'Base','Fine'},'FontSize',7);
 if j==1,ylim(ad,[.54 .57]);set(ad,'YTick',[.54 .55 .56 .57]);else,ylim(ad,[-4.11 -4.08]);set(ad,'YTick',[-4.11 -4.10 -4.09 -4.08]);end
 ytickformat(ad,'%.2f');
 title(ad,sprintf('%s, %.2g S/m',roles(j),sf),'FontSize',7.2);
 if j==1,ylabel(ad,'\itg\rm_{11} - \itg\rm_{77} (\muS)','FontSize',8);end
end
panel(f,'d',95,66,W,H);writetable(D,fullfile(outdir,'plotted_joint_mesh_comparison.csv'));

drawnow;stem=fullfile(outdir,'Fig08_finite_three_dimensional');savefig(f,[stem '.fig']);
% MATLAB's vector painter can draw hidden back-wall lines through the
% ellipsoid. Preserve the editable native 3D FIG above, then rasterize only
% the geometry viewport at 600 dpi. Data axes and figure text stay vector.
anns=findall(f,'Type','textboxshape');set(anns,'Visible','off');
canvas=fullfile(outdir,'source_notes','geometry_render_canvas.png');
print(f,canvas,'-dpng','-r600');set(anns,'Visible','on');
im=imread(canvas);pos=ag.Position;[nh,nw,~]=size(im);
x1=max(1,floor(pos(1)*nw)+1);x2=min(nw,ceil((pos(1)+pos(3))*nw));
y1=max(1,floor((1-pos(2)-pos(4))*nh)+1);y2=min(nh,ceil((1-pos(2))*nh));
crop=im(y1:y2,x1:x2,:);imwrite(crop,fullfile(outdir,'source_notes','geometry_view_600dpi.png'));
delete(ag);ar=axes(f,'Units','normalized','Position',pos);
image(ar,crop);axis(ar,'off');set(ar,'Position',pos,'PositionConstraint','innerposition');
exportgraphics(f,[stem '.png'],'Resolution',300,'BackgroundColor','white');
print(f,[stem '.pdf'],'-dpdf','-painters','-r600');
print(f,[stem '.svg'],'-dsvg','-painters','-r600');print(f,[stem '.tiff'],'-dtiff','-r600');
close(f);
% MATLAB's PDF writer downscales images to screen resolution even with
% -r600. Replace only that image object; native vector text/data are kept.
helper=fullfile(fileparts(mfilename('fullpath')),'finish_fig08_pdf.py');
pyexe=getenv('JGRSE_PYTHON');if isempty(pyexe),pyexe='python';end
if isfile(helper)
 [status,msg]=system(sprintf('"%s" "%s" "%s"',pyexe,helper,outdir));
 assert(status==0,'PDF geometry embedding failed: %s',msg);fprintf('%s',msg);
else
 warning('Run python finish_fig08_pdf.py "bundle_folder" to restore the PDF geometry resolution.');
end
fprintf('Saved %s\n',stem);
end

function draw_geometry(ax,p)
fig08_geometry_tgrs(ax,p,1);
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
