function figure02_intervention_acd_matlab(outdir, project_root)
% FIGURE02_INTERVENTION_ACD_MATLAB  Targeted revision of panels a, c and d.
% Uses existing fields only. No COMSOL solve and no modification of runs/.
% Call with an output folder containing inputs/ to reproduce a review bundle.
if nargin < 2 || isempty(project_root)
    project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin < 1 || isempty(outdir)
    outdir = fullfile(project_root,'outputs','JGRSE_20261028_Fig02_RevisionACD');
end
if ~isfolder(outdir), mkdir(outdir); end
inputdir = fullfile(outdir,'inputs');
states = {'full','cut','side'};
names = {'Full','Cut','Side removal'};
sourcefiles = {'parameters.json','result.csv','local_window_fields.csv', ...
    'section_current.csv','topology_interfaces.json','mesh_stats.json'};
if ~isfolder(inputdir)
    mkdir(inputdir);
    for k=1:3
        dest=fullfile(inputdir,states{k}); mkdir(dest);
        src=fullfile(project_root,'runs','2026-10-21_round2', ...
            ['strict_h10_L60_B_mech_m3_' states{k}]);
        for j=1:numel(sourcefiles)
            copyfile(fullfile(src,sourcefiles{j}),fullfile(dest,sourcefiles{j}));
        end
    end
    copyfile(fullfile(project_root,'runs','2026-10-21_round2', ...
        'strict_h10_L60_B_bgz','result.csv'), ...
        fullfile(inputdir,'background_result.csv'));
end

set(groot,'defaultAxesFontName','Times New Roman', ...
    'defaultTextFontName','Times New Roman', ...
    'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultAxesXColor',[0 0 0],'defaultAxesYColor',[0 0 0], ...
    'defaultTextColor',[0 0 0],'defaultAxesLineWidth',0.6, ...
    'defaultAxesTickDir','in','defaultAxesBox','on', ...
    'defaultAxesTickLength',[0.018 0.018], ...
    'defaultAxesTitleFontWeight','normal', ...
    'defaultAxesLabelFontSizeMultiplier',1.05, ...
    'defaultAxesTitleFontSizeMultiplier',1.05, ...
    'defaultLegendFontName','Times New Roman', ...
    'defaultLegendFontSize',7.5,'defaultLegendTextColor',[0 0 0], ...
    'defaultLegendBox','off');
blue=[0 0 1]; red=[1 0 0]; green=[0 0.5 0];
casecolors={[0 0 0],red,blue};
markers={'o','^','o'};
bg=readtable(fullfile(inputdir,'background_result.csv'));
sigma0=bg.sigma_eff_S_m(1);
sigma=zeros(3,1); gain=zeros(3,1); drop=zeros(3,1);
fields=cell(1,3); params=cell(1,3); results=cell(1,3);
for k=1:3
    params{k}=jsondecode(fileread(fullfile(inputdir,states{k},'parameters.json')));
    results{k}=readtable(fullfile(inputdir,states{k},'result.csv'));
    fields{k}=readtable(fullfile(inputdir,states{k},'local_window_fields.csv'));
    sigma(k)=results{k}.sigma_eff_S_m(1);
    gain(k)=100*(sigma(k)-sigma0)/sigma0;
    drop(k)=1000*results{k}.window_dV_V(1);
    assert(height(fields{k})==241^2 && all(isfinite(fields{k}.V_V)), ...
        'Expected complete saved 241-by-241 field grid.');
end
p=params{2};
loss=100*(sigma(1)-sigma)/sigma(1);
gainloss=100*(sigma(1)-sigma)/(sigma(1)-sigma0);
metrics=table(string(states(:)),sigma,repmat(sigma0,3,1),gain,loss,gainloss,drop, ...
    'VariableNames',{'state','sigma_eff_S_m','sigma0_S_m', ...
    'gain_percent','total_conductivity_loss_percent', ...
    'incremental_gain_loss_percent','window_dV_mV'});
writetable(metrics,fullfile(outdir,'plotted_metrics.csv'));

W=183; H=158;
f=figure('Visible','off','Color','white','Units','centimeters', ...
    'Position',[2 2 W/10 H/10],'PaperUnits','centimeters', ...
    'PaperSize',[W H]/10,'PaperPosition',[0 0 W H]/10, ...
    'InvertHardcopy','off','Renderer','painters');
% A compact two-column composition: one shared specimen and three local
% views, followed by the macroscopic and local electrical responses.
cmap=readmatrix(fullfile(outdir,'source_notes','ling_L08_F02_stream_colormap.csv'));
colormap(f,cmap);
ax=newaxes(f,[20 82 66 66],W,H);
draw_global(ax,p,'full',green);
title(ax,'Intervention layout','FontSize',9,'FontWeight','normal');
set(ax,'XTick',[0 60 120],'YTick',[0 60 120]);
xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
ang=linspace(0,2*pi,361);
patch(ax,p.comp_x_mm+p.comp_radius_mm*cos(ang), ...
    p.comp_z_mm+p.comp_radius_mm*sin(ang),green,'EdgeColor',green);
% Black geometry arrow follows the material-pointer convention in L08-F03e.
% The arrow ends at the compensation circle's rim, leaving the disk visible.
arrowX=(20+[10 p.comp_x_mm]*66/p.L_mm)/W;
arrowY=(82+[76 p.comp_z_mm+p.comp_radius_mm]*66/p.L_mm)/H;
annotation(f,'arrow',arrowX,arrowY,'Color','black','LineWidth',0.8, ...
    'HeadLength',4,'HeadWidth',4);
text(ax,3,81,{'Add-back','(cut / side)'},'FontSize',7.2, ...
    'HorizontalAlignment','left','VerticalAlignment','middle');
% Local coordinates shown once; the arrows are geometry axes, not currents.
quiver(ax,88,67,9*cosd(60),9*sind(60),0,'Color','black', ...
    'LineWidth',0.6,'MaxHeadSize',0.65);
quiver(ax,88,67,-9*sind(60),9*cosd(60),0,'Color','black', ...
    'LineWidth',0.6,'MaxHeadSize',0.65);
text(ax,94,77,'s','FontAngle','italic','FontSize',8);
text(ax,77,74,'r','FontAngle','italic','FontSize',8);
panel(f,'a',5,154,W,H);
% Matching areas are given in the caption, rather than repeated on the plot.

% (b) Same physical local window, same seeds, actual saved current field.
% Streamlines indicate direction/path, not current magnitude or flux density.
stream_count=zeros(3,1);
potential_profiles=[];
field_y=[126 98 70];
field_axes=gobjects(1,3);
for k=1:3
    ax=newaxes(f,[110 field_y(k) 66 22],W,H);
    field_axes(k)=ax;
    T=fields{k};
    sv=unique(T.s_mm)-p.cut_u_mm; tv=unique(T.t_mm);
    [SS,TT]=meshgrid(sv,tv);
    Jx=reshape(T.Jx_A_m2,241,241);
    Jz=reshape(T.Jz_A_m2,241,241);
    Js=Jx*cosd(p.theta_deg)+Jz*sind(p.theta_deg);
    Jt=-Jx*sind(p.theta_deg)+Jz*cosd(p.theta_deg);
    mag=hypot(Js,Jt);
    draw_local(ax,p,states{k},green);
    % Equal start coordinates in all cases. Normalization changes only the
    % streamline integration speed, never the measured field direction.
    seedS=[linspace(-5.8,5.8,9) repmat(5.98,1,13)];
    seedT=[repmat(1.98,1,9) linspace(-1.9,1.9,13)];
    paths=stream2(SS,TT,Js./mag,Jt./mag,seedS,seedT,[0.5 10000]);
    stream_count(k)=numel(paths);
    for q=1:numel(paths)
        xy=paths{q};
        if size(xy,1)<4, continue; end
        colorJ=interp2(SS,TT,mag,xy(:,1),xy(:,2),'linear');
        surface(ax,[xy(:,1) xy(:,1)],[xy(:,2) xy(:,2)], ...
            zeros(size(xy,1),2),[colorJ colorJ], ...
            'FaceColor','none','EdgeColor','interp','LineWidth',0.6);
        % One small tangent arrow on selected streamlines, away from axes.
        if mod(q,7)==1 && size(xy,1)>30
            ii=round(size(xy,1)*0.28);
            jj=min(ii+4,size(xy,1));
            d=xy(jj,:)-xy(ii,:);
            if norm(d)>1e-8 && abs(xy(ii,1))<5.2 && abs(xy(ii,2))<1.7
                d=0.35*d/norm(d);
                quiver(ax,xy(ii,1),xy(ii,2),d(1),d(2),0, ...
                    'Color','black','LineWidth',0.35,'MaxHeadSize',1.4, ...
                    'ShowArrowHead','on');
            end
        end
    end
    % Re-outline the true material boundary so every gap stays visible.
    boundary_local(ax,p,states{k},green,0.7);
    xlim(ax,[-6 6]); ylim(ax,[-2 2]);
    set(ax,'XTick',[],'YTick',[],'Layer','top','ColorScale','linear');
    clim(ax,[0 15]);
    title(ax,names{k},'FontSize',8.5,'FontWeight','normal');
    V=reshape(T.V_V,241,241);
    [~,it]=min(abs(tv));
    potential_profiles=[potential_profiles V(it,:).']; %#ok<AGROW>
end
panel(f,'b',95,154,W,H);
cb=colorbar(field_axes(3),'Location','southoutside');
cb.Units='normalized'; cb.Position=[110/W 59/H 66/W 2.2/H];
cb.Ticks=0:3:15;
cb.FontName='Times New Roman'; cb.FontSize=7; cb.Color='black';
cb.Label.String='';
annotation(f,'textbox',[110/W 62/H 66/W 4/H], ...
    'String','|J| (A m^{-2})','HorizontalAlignment','center', ...
    'EdgeColor','none','Margin',0,'FontName','Times New Roman', ...
    'FontSize',8,'Color','black');
% colorbar creation can move its axes: restore the precise physical aspect.
field_axes(3).Position=[110 field_y(3) 66 22]./[W H W H];
prof=table(sv,potential_profiles(:,1),potential_profiles(:,2),potential_profiles(:,3), ...
    'VariableNames',{'s_minus_sc_mm','V_full_V','V_cut_V','V_side_V'});
writetable(prof,fullfile(outdir,'plotted_potential_profiles.csv'));

% (c) Actual gain relative to the common saved no-fracture background.
ax=newaxes(f,[20 11 66 34],W,H);
for k=1:3
    % Independent discrete states: no artificial continuous connecting line.
    % L07-F09 uses category bars; narrow unfilled bars retain our lighter style.
    bar(ax,k,gain(k),0.32,'FaceColor','white','EdgeColor',casecolors{k}, ...
        'LineWidth',1.15,'ShowBaseLine','off');
    text(ax,k,gain(k)+0.75,sprintf('%.2f',gain(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',7.5);
end
set(ax,'XTick',1:3,'XTickLabel',names,'YTick',[0 10 20 30]);
xlim(ax,[0.5 3.5]); ylim(ax,[0 30]);
ylabel(ax,'Conductivity gain (%)');
title(ax,'Macroscopic response','FontSize',8.5);
panel(f,'c',5,50,W,H);

% (d) A true line sample through t=0, not a fitted or reconstructed field.
ax=newaxes(f,[110 11 66 34],W,H);
hside=plot(ax,sv,potential_profiles(:,3),'-o','Color',blue,'LineWidth',1.15, ...
    'MarkerIndices',1:60:241,'MarkerSize',3.4,'MarkerFaceColor','white');
hfull=plot(ax,sv,potential_profiles(:,1),'--','Color','black','LineWidth',1.0);
[~,icutleft]=min(abs(sv+p.cut_width_mm/2));
[~,icutright]=min(abs(sv-p.cut_width_mm/2));
cutmarkers=unique([1:40:241 icutleft icutright]);
hcut=plot(ax,sv,potential_profiles(:,2),'-^','Color',red,'LineWidth',1.15, ...
    'MarkerIndices',cutmarkers,'MarkerSize',3.5,'MarkerFaceColor','white');
xlim(ax,[-6 6]); ylim(ax,[0.44 0.605]);
set(ax,'XTick',[-6 -3 0 3 6],'YTick',[0.45 0.50 0.55 0.60]);
xlabel(ax,'\its\rm - \its_c\rm (mm)'); ylabel(ax,'\itV\rm(\its\rm, 0) (V)');
title(ax,'Local potential','FontSize',8.5);
lg=legend(ax,[hfull hcut hside],names,'Location','southeast','FontSize',7);
if isprop(lg,'ItemTokenSize'), lg.ItemTokenSize=[14 8]; end
panel(f,'d',95,50,W,H);
% A small second view reveals the close full/side curves without changing
% their coordinates in the main axes. Borrow only the inset layout from
% L04-F12; that source used rescaling, whereas this is a true local zoom.
axi=newaxes(f,[116 29 23 12],W,H);
sel=find(sv>=-3-1e-9 & sv<=3+1e-9);
plot(axi,sv(sel),potential_profiles(sel,3),'-o','Color',blue, ...
    'LineWidth',0.9,'MarkerIndices',1:30:numel(sel), ...
    'MarkerSize',2.6,'MarkerFaceColor','white');
plot(axi,sv(sel),potential_profiles(sel,1),'--','Color','black','LineWidth',0.85);
xlim(axi,[-3 3]); ylim(axi,[0.5045 0.509]);
set(axi,'XTick',[-3 0 3],'YTick',[0.505 0.507 0.509], ...
    'FontSize',6.2,'LineWidth',0.5,'TickLength',[0.035 0.035]);
axi.YAxis.Exponent=0;
title(axi,'Full / Side','FontSize',6.5,'FontWeight','normal');

drawnow;
stem=fullfile(outdir,'Fig02_key_path_intervention');
savefig(f,[stem '.fig']);
exportgraphics(f,[stem '.png'],'Resolution',250,'BackgroundColor','white');
exportgraphics(f,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');
print(f,[stem '.svg'],'-dsvg','-painters');
print(f,[stem '.tiff'],'-dtiff','-r600');
close(f);
fprintf('Saved %s\n',stem);
disp(metrics);
end

function ax=newaxes(f,mm,W,H)
ax=axes(f,'Units','normalized','Position',mm./[W H W H]);
hold(ax,'on'); box(ax,'on');
ax.PositionConstraint='innerposition';
end

function panel(f,letter,x,y,W,H)
annotation(f,'textbox',[x/W (y-3.5)/H 12/W 5/H], ...
    'String',['(' letter ')'],'EdgeColor','none','Margin',0, ...
    'FontName','Times New Roman','FontSize',9,'FontAngle','italic', ...
    'Color','black','VerticalAlignment','middle');
end

function shape=ellipse_shape(p,state,local)
v=linspace(0,2*pi,2401);
s=p.a_mm*cos(v); t=p.b_mm*sin(v);
shape=polyshape(s,t,'Simplify',false);
if strcmp(state,'cut')
    xx=p.cut_u_mm+[-.5 .5 .5 -.5]*p.cut_width_mm;
    yy=[-3 -3 3 3];
    shape=subtract(shape,polyshape(xx,yy));
elseif strcmp(state,'side')
    xx=[p.side_u1_mm p.side_u2_mm p.side_u2_mm p.side_u1_mm];
    yy=[p.side_width_mm p.side_width_mm 3 3];
    shape=subtract(shape,polyshape(xx,yy));
end
if local
    shape=translate(shape,[-p.cut_u_mm 0]);
else
    shape=rotate(shape,p.theta_deg);
    shape=translate(shape,[p.cx_mm p.cz_mm]);
end
end

function draw_global(ax,p,state,green)
z=[0 p.z_interfaces_mm(:).' p.L_mm];
for j=1:numel(p.layer_materials)
    if p.layer_materials(j)==1, col=[.75 .75 .75]; else, col=[1 1 1]; end
    patch(ax,[0 p.L_mm p.L_mm 0],[z(j) z(j) z(j+1) z(j+1)],col,'EdgeColor','none');
end
shape=ellipse_shape(p,state,false);
plot(ax,shape,'FaceColor',green,'FaceAlpha',1,'EdgeColor',green,'LineWidth',0.45);
if ~strcmp(state,'full')
    ang=linspace(0,2*pi,361);
    patch(ax,p.comp_x_mm+p.comp_radius_mm*cos(ang), ...
        p.comp_z_mm+p.comp_radius_mm*sin(ang),green,'EdgeColor',green,'LineWidth',0.4);
    plot(ax,[p.comp_x_mm 11],[p.comp_z_mm+1.3 78],'-k','LineWidth',0.45);
    text(ax,4,84,'Add-back','FontSize',7,'VerticalAlignment','middle');
end
% Expanded locator only: the unchanged panel (b) remains 12-by-4 mm.
% Do not describe this 18-by-8 mm outline as the exact field crop.
ss=p.cut_u_mm+[-9 9 9 -9 -9]; tt=[-4 -4 4 4 -4];
xx=p.cx_mm+ss*cosd(p.theta_deg)-tt*sind(p.theta_deg);
zz=p.cz_mm+ss*sind(p.theta_deg)+tt*cosd(p.theta_deg);
plot(ax,xx,zz,'--k','LineWidth',1.05);
if strcmp(state,'full')
    text(ax,7,107,'High','FontSize',7.5);
    text(ax,98,88,'Low','FontSize',7.5,'HorizontalAlignment','center');
end
xlim(ax,[0 p.L_mm]); ylim(ax,[0 p.L_mm]); daspect(ax,[1 1 1]);
ax.Layer='top';
end

function draw_local(ax,p,state,green)
% Exact polygon clipping of the shared local window by horizontal CT-free
% material layers. No center-based approximation is used for the ellipse.
window=polyshape([-6 6 6 -6],[-6 -6 6 6]);
z=[0 p.z_interfaces_mm(:).' p.L_mm];
for j=1:numel(p.layer_materials)
    if p.layer_materials(j)~=1, continue; end
    xx=[0 p.L_mm p.L_mm 0]-p.cx_mm;
    zz=[z(j) z(j) z(j+1) z(j+1)]-p.cz_mm;
    ss=xx*cosd(p.theta_deg)+zz*sind(p.theta_deg)-p.cut_u_mm;
    tt=-xx*sind(p.theta_deg)+zz*cosd(p.theta_deg);
    sh=intersect(polyshape(ss,tt),window);
    if sh.NumRegions>0
        plot(ax,sh,'FaceColor',[.75 .75 .75],'FaceAlpha',1,'EdgeColor','none');
    end
end
shape=intersect(ellipse_shape(p,state,true),window);
plot(ax,shape,'FaceColor','white','FaceAlpha',1,'EdgeColor',green,'LineWidth',0.7);
xlim(ax,[-6 6]); ylim(ax,[-6 6]); daspect(ax,[1 1 1]);
end

function boundary_local(ax,p,state,green,width)
shape=ellipse_shape(p,state,true);
[x,y]=boundary(shape);
plot(ax,x,y,'Color',green,'LineWidth',width);
end
