function figure05_framework_matlab(outdir,project_root)
% Figure 5. Read saved predictions and analytical coefficients only.
% Does not run COMSOL, assemble a solver, or use target truth fields.
if nargin<2 || isempty(project_root)
    project_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin<1 || isempty(outdir)
    outdir=fullfile(project_root,'outputs','JGRSE_20261028_Fig05_Framework');
end
idir=fullfile(outdir,'inputs');assert(isfolder(idir),'Prepared inputs are required.');
c=jsondecode(fileread(fullfile(idir,'geometry_A_NN.json')));
pred=jsondecode(fileread(fullfile(idir,'prediction_A_NN_sf0p3.json')));
assert(~pred.target_truth_used && strcmp(pred.role,'A') && strcmp(pred.layout,'NN'));
Y=readtable(fullfile(idir,'NN_port_coefficients.csv'),'TextType','string');
M=readtable(fullfile(idir,'A_NN_sf0p3_predicted_modes.csv'),'TextType','string');
Kstar=8; W=183;H=157; blue=[0 0 1];red=[1 0 0];green=[0 .5 0];
sidecolors={blue,red};sides={'upper','lower'};
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman', ...
    'defaultAxesFontSize',8,'defaultTextFontSize',8,'defaultTextColor',[0 0 0], ...
    'defaultAxesXColor',[0 0 0],'defaultAxesYColor',[0 0 0], ...
    'defaultAxesLineWidth',.6,'defaultAxesBox','on','defaultAxesTickDir','in', ...
    'defaultAxesTickLength',[.018 .018],'defaultAxesTitleFontWeight','normal', ...
    'defaultAxesTitleFontSizeMultiplier',1.05,'defaultAxesLabelFontSizeMultiplier',1.05, ...
    'defaultLegendFontName','Times New Roman','defaultLegendFontSize',7.2, ...
    'defaultLegendTextColor',[0 0 0],'defaultLegendBox','off');
f=figure('Visible','off','Color','white','Units','centimeters', ...
    'Position',[2 2 W/10 H/10],'PaperUnits','centimeters','PaperSize',[W H]/10, ...
    'PaperPosition',[0 0 W H]/10,'InvertHardcopy','off','Renderer','painters');

% (a) The physical layered geometry, including its original electrode buffers.
ax=axesmm(f,[20 87 60 60],W,H);
layers(ax,c,0,120);ellipse(ax,c,green);
interfaces(ax,blue,red);
plot(ax,[0 120],[0 0],'-k','LineWidth',1.3);
plot(ax,[0 120],[120 120],'-k','LineWidth',1.3);
set(ax,'XTick',[0 60 120],'YTick',[0 30 60 90 120]);
xlabel(ax,'\itx\rm (mm)');ylabel(ax,'\itz\rm (mm)');
title(ax,'Layered model','FontSize',8.5);
text(ax,60,110,'Electrode buffer','HorizontalAlignment','center','FontSize',7.2);
text(ax,60,10,'Electrode buffer','HorizontalAlignment','center','FontSize',7.2);
text(ax,13,75,'Core','FontSize',7.5);
text(ax,124,90,'\Gamma_+','Clipping','off','FontSize',8);
text(ax,124,30,'\Gamma_-','Clipping','off','FontSize',8);
text(ax,117,123,'1 V','HorizontalAlignment','right','Clipping','off','FontSize',7);
text(ax,117,3,'0 V','HorizontalAlignment','right','FontSize',7);
panel(f,'a',5,153,W,H);

% (b) Mathematical reference: whole 30 mm exteriors are equivalent tensors.
% The physical buffers are included in those tensors, not kept as extra slabs.
ax=axesmm(f,[110 87 60 60],W,H);
patch(ax,[0 120 120 0],[0 0 30 30],[.94 .94 .94],'EdgeColor','none');
patch(ax,[0 120 120 0],[90 90 120 120],[.94 .94 .94],'EdgeColor','none');
layers(ax,c,30,90);ellipse(ax,c,green);interfaces(ax,blue,red);
plot(ax,[0 120],[0 0],'-k','LineWidth',1.3);
plot(ax,[0 120],[120 120],'-k','LineWidth',1.3);
set(ax,'XTick',[],'YTick',[]);
title(ax,'Self-consistent model','FontSize',8.5);
text(ax,60,110,'Equivalent outer','HorizontalAlignment','center','FontSize',7.4);
text(ax,60,101,'\sigma_{x,+}^{eq},  \sigma_{z,+}^{eq}','HorizontalAlignment','center','FontSize',7.2);
text(ax,60,19,'Equivalent outer','HorizontalAlignment','center','FontSize',7.4);
text(ax,60,9,'\sigma_{x,-}^{eq},  \sigma_{z,-}^{eq}','HorizontalAlignment','center','FontSize',7.2);
text(ax,13,76,{'Explicit','core'},'FontSize',7.5);
text(ax,124,90,'\itV\rm_+(\itx\rm)','Clipping','off','FontSize',7.5);
text(ax,124,30,'\itV\rm_-(\itx\rm)','Clipping','off','FontSize',7.5);
% Short arrows denote current exchange, not prescribed interface voltages.
quiver(ax,94,85,0,12,0,'Color','black','LineWidth',.65,'MaxHeadSize',.45);
quiver(ax,94,35,0,-12,0,'Color','black','LineWidth',.65,'MaxHeadSize',.45);
text(ax,99,96,'\itq\rm_+','FontSize',7.2);
text(ax,99,24,'\itq\rm_-','FontSize',7.2);
text(ax,117,123,'1 V','HorizontalAlignment','right','Clipping','off','FontSize',7);
text(ax,117,3,'0 V','HorizontalAlignment','right','FontSize',7);
captionline(f,[109 74 69 7], ...
    '(\itA\rm_0 + \itH\rm_K)\itu\rm_F = \itb\rm  \rightarrow  \itV\rm_\pm(\itx\rm), \itg',W,H,8);
panel(f,'b',95,153,W,H);

% (c) Information retained: nonzero response orders 1..8 are actual.
% All later orders keep M0, shown at ratio 1 rather than omitted/insulated.
ax=axesmm(f,[20 15 66 44],W,H);
set(ax,'YScale','log');
plot(ax,[0 32],[1 1],'--k','LineWidth',.9);
plot(ax,[8.5 8.5],[.18 5.5],':k','LineWidth',.55);
h=gobjects(1,5);ratios=table();
for s=1:2
    z=sortrows(Y(Y.side==sides{s},:),'n');assert(isequal(z.n,(0:64)'));
    ar=z.Y_actual_S_m2./z.Y_reference_S_m2;hy=ones(65,1);hy(2:9)=ar(2:9);
    sel=z.n<=32;
    h(s)=plot(ax,z.n(sel),ar(sel),'-','Color',sidecolors{s},'LineWidth',1.15);
    if s==1,mi=[0:8 10:4:30];mark='o';else,mi=[0:8 12:4:32];mark='^';end
    h(s+2)=plot(ax,mi,hy(mi+1),mark,'Color',sidecolors{s}, ...
        'LineStyle','none','MarkerSize',3.4,'MarkerFaceColor','white','LineWidth',.9);
    z.actual_over_reference=ar;z.M8_over_reference=hy;
    ratios=[ratios;z]; %#ok<AGROW>
end
h(5)=plot(ax,NaN,NaN,'--k','LineWidth',.9);
xlim(ax,[0 32]);ylim(ax,[.18 5.5]);
set(ax,'XTick',[0 8 16 24 32],'YTick',[.2 .5 1 2 5],'YTickLabel',{'0.2','0.5','1','2','5'});
xlabel(ax,'Mode number, \itn');ylabel(ax,'\itY\rm_n / \itY\rm_n^{eq}');
title(ax,'Retained outer response','FontSize',8.5);
text(ax,8.2,4.35,'\itK\rm = 8','HorizontalAlignment','right','FontSize',7);
text(ax,29,4.15,'Lower','FontSize',7.2,'HorizontalAlignment','right');
text(ax,29,.265,'Upper','FontSize',7.2,'HorizontalAlignment','right');
plot(ax,[17 20],[2.8 2.8],'-k','LineWidth',1.0);text(ax,21.5,2.8,'Layer-resolved','FontSize',7.2);
plot(ax,[17 20],[2.1 2.1],'--k','LineWidth',.9);text(ax,21.5,2.1,'M0','FontSize',7.2);
plot(ax,17.7,1.55,'ok','MarkerFaceColor','white','MarkerSize',3.2,'LineWidth',.9);
plot(ax,19.3,1.55,'^k','MarkerFaceColor','white','MarkerSize',3.2,'LineWidth',.9);
text(ax,21.5,1.55,'M8','FontSize',7.2);
panel(f,'c',5,64,W,H);
writetable(ratios,fullfile(outdir,'plotted_admittance_ratios.csv'));

% (d) Self-consistently predicted nonuniform potential, not target truth.
% K selects response corrections; all 64 saved potential modes are retained.
ax=axesmm(f,[110 15 66 44],W,H);
x=linspace(0,120,601)';curves=table(x,'VariableNames',{'x_mm'});
hc=gobjects(1,2);
for s=1:2
    for K=[0 8]
        rows=sortrows(M(M.side==sides{s} & M.K==K,:),'n');
        assert(height(rows)==65 && isequal(rows.n,(0:64)'));
        v=rows.v_predicted_V;
        method=pred.methods([pred.methods.K]==K);sideidx=3-s;
        assert(max(abs(v-method.port_v(sideidx,:)'))<1e-13);
        nonuniform=cos(x*(1:64)*pi/120)*v(2:65)*1000;
        if K==0,sty='--';else,sty='-';end
        hp=plot(ax,x,nonuniform,sty,'Color',sidecolors{s},'LineWidth',1.15);
        if K==8,hc(s)=hp;end
        curves.(sprintf('%s_M%d_nonuniform_mV',sides{s},K))=nonuniform;
    end
end
plot(ax,[0 120],[0 0],':k','LineWidth',.45);
xlim(ax,[0 120]);ylim(ax,[-55 35]);
set(ax,'XTick',[0 30 60 90 120],'YTick',[-50 -25 0 25]);
xlabel(ax,'\itx\rm (mm)');ylabel(ax,'\itV\rm_\pm - \langle\itV\rm_\pm\rangle (mV)');
title(ax,'Predicted interface potential','FontSize',8.5);
lg=legend(ax,hc,{'Upper','Lower'},'Location','southwest','FontSize',7.1,'AutoUpdate','off');
if isprop(lg,'ItemTokenSize'),lg.ItemTokenSize=[13 8];end
% Independent black line-style key keeps side and method encodings distinct.
plot(ax,[76 88],[28 28],'--k','LineWidth',1.0);text(ax,91,28,'M0','FontSize',7.1);
plot(ax,[76 88],[19 19],'-k','LineWidth',1.0);text(ax,91,19,'M8','FontSize',7.1);
text(ax,116,-50,'A, NN, \sigma_f = 0.3 S/m', ...
    'HorizontalAlignment','right','FontSize',6.8);
panel(f,'d',95,64,W,H);
writetable(curves,fullfile(outdir,'predicted_interface_profiles.csv'));

drawnow;stem=fullfile(outdir,'Fig05_selfconsistent_framework');savefig(f,[stem '.fig']);
exportgraphics(f,[stem '.png'],'Resolution',250,'BackgroundColor','white');
exportgraphics(f,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');
print(f,[stem '.svg'],'-dsvg','-painters');print(f,[stem '.tiff'],'-dtiff','-r600');
close(f);fprintf('Saved %s\n',stem);
end

function ax=axesmm(f,mm,W,H)
ax=axes(f,'Units','normalized','Position',mm./[W H W H]);hold(ax,'on');box(ax,'on');
ax.PositionConstraint='innerposition';ax.Layer='top';
end
function layers(ax,c,lo,hi)
edges=[0;c.z_interfaces_mm(:);120];
for k=1:numel(c.layer_materials)
    a=max(edges(k),lo);b=min(edges(k+1),hi);if b<=a,continue;end
    if c.layer_materials(k)==1,col=[.75294 .75294 .75294];else,col=[1 1 1];end
    patch(ax,[0 120 120 0],[a a b b],col,'EdgeColor','none');
end
xlim(ax,[0 120]);ylim(ax,[0 120]);
end
function ellipse(ax,c,col)
t=linspace(0,2*pi,1001);s=c.a_mm*cos(t);v=c.b_mm*sin(t);
x=c.cx_mm+s*cosd(c.theta_deg)-v*sind(c.theta_deg);
z=c.cz_mm+s*sind(c.theta_deg)+v*cosd(c.theta_deg);
patch(ax,x,z,col,'EdgeColor',col,'LineWidth',.45);
end
function interfaces(ax,blue,red)
plot(ax,[0 120],[90 90],'-','Color',blue,'LineWidth',1.15);
plot(ax,[0 120],[30 30],'-','Color',red,'LineWidth',1.15);
end
function panel(f,s,x,y,W,H)
annotation(f,'textbox',[x/W (y-3.5)/H 12/W 5/H],'String',['(' s ')'], ...
    'EdgeColor','none','Margin',0,'FontSize',9,'FontAngle','italic', ...
    'FontName','Times New Roman','Color','black','VerticalAlignment','middle');
end
function captionline(f,mm,str,W,H,fs)
annotation(f,'textbox',mm./[W H W H],'String',str,'Interpreter','tex', ...
    'EdgeColor','none','Margin',0,'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName','Times New Roman','FontSize',fs,'Color','black');
end
