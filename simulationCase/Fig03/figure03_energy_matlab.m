function figure03_energy_matlab(outdir,project_root)
% Saved-COMSOL-field postprocessing only. Does not call COMSOL or solve PDEs.
if nargin<2 || isempty(project_root)
    project_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin<1 || isempty(outdir)
    outdir=fullfile(project_root,'outputs','JGRSE_20261028_Fig03_Energy');
end
if ~isfolder(outdir),mkdir(outdir);end
idir=fullfile(outdir,'inputs');
ids={'strict_h10_L60_A_m3_full','strict_h10_L60_B_m3_full', ...
     'core_h20_L60_B_m3_full','core_h20_L60_C_m3_full'};
if ~isfolder(idir)
    mkdir(idir);
    rr=fullfile(project_root,'runs','round3_predict');
    copyfile(fullfile(rr,'energy_decomposition.csv'),fullfile(idir,'energy_decomposition.csv'));
    for k=1:4
        dest=fullfile(idir,ids{k});mkdir(dest);
        src=fullfile(rr,'development',ids{k});
        for name={'background_fields.mat','background_summary.json','parameters.json'}
            copyfile(fullfile(src,name{1}),fullfile(dest,name{1}));
        end
        copyfile(fullfile(src,'sf_10','r3_ellipse_axis_fields.csv'),fullfile(dest,'target_axis_fields.csv'));
        copyfile(fullfile(src,'sf_10','r3_endpoint_potential_averages.json'),fullfile(dest,'target_endcaps.json'));
    end
end
allEnergy=readtable(fullfile(idir,'energy_decomposition.csv'),'TextType','string');
E=allEnergy([],:);profiles=cell(1,4);caps=zeros(4,1);
for k=1:4
    row=allEnergy(allEnergy.case_id==ids{k} & allEnergy.sigma_f_S_m==10 & allEnergy.scope=="development",:);
    assert(height(row)==1,'Need exactly one saved development state per case.');
    E=[E;row]; %#ok<AGROW>
    src=fullfile(idir,ids{k});
    T=readtable(fullfile(src,'target_axis_fields.csv'));
    bg=load(fullfile(src,'background_fields.mat'),'bgfields');
    b=bg.bgfields;
    P=double(b.p.');tr=double(b.t.')+1; v0=double(b.d1(:));
    if max(P(:))<1,P=P*1000;end
    assert(min(tr(:))==1 && max(tr(:))<=size(P,1));
    % Preserve the saved triangle connectivity. Linear barycentric sampling
    % of the real, piecewise-linear flat-layer background, without a new mesh.
    tri=triangulation(tr,P);
    [ti,bc]=pointLocation(tri,[T.x_mm T.z_mm]);
    assert(all(isfinite(ti)),'Background query outside saved mesh.');
    nod=tr(ti,:);
    V0=sum(bc.*reshape(v0(nod),size(nod)),2);
    assert(all(isfinite(V0)));
    T.background_V_V=V0;
    T.perturbation_V_V=T.V_V-V0;
    profiles{k}=T;
    writetable(T,fullfile(outdir,[ids{k} '_axis_comparison.csv']));
    bgsummary=jsondecode(fileread(fullfile(src,'background_summary.json')));
    caps(k)=bgsummary.cap_dV_V;
end
writetable(E,fullfile(outdir,'selected_energy_states.csv'));
% Pair orientation remains second-minus-first for every signed quantity.
dU=[E.U_S(2)-E.U_S(1);E.U_S(4)-E.U_S(3)];
dR=[E.R_S(2)-E.R_S(1);E.R_S(4)-E.R_S(3)];
dg=[E.deltaG_S(2)-E.deltaG_S(1);E.deltaG_S(4)-E.deltaG_S(3)];
pairs=table(["Strict B-A";"Core C-B"],dU,dR,dg,dU-dR, ...
    'VariableNames',{'pair','delta_U_S','delta_R_S','delta_g_S','delta_U_minus_delta_R_S'});
writetable(pairs,fullfile(outdir,'plotted_pair_changes.csv'));

set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman', ...
    'defaultAxesFontSize',8,'defaultTextFontSize',8,'defaultTextColor',[0 0 0], ...
    'defaultAxesXColor',[0 0 0],'defaultAxesYColor',[0 0 0], ...
    'defaultAxesLineWidth',0.6,'defaultAxesBox','on','defaultAxesTickDir','in', ...
    'defaultAxesTickLength',[.018 .018],'defaultAxesTitleFontWeight','normal', ...
    'defaultAxesTitleFontSizeMultiplier',1.05,'defaultAxesLabelFontSizeMultiplier',1.05, ...
    'defaultLegendFontName','Times New Roman','defaultLegendFontSize',7.5, ...
    'defaultLegendTextColor',[0 0 0],'defaultLegendBox','off');
W=183;H=132;
f=figure('Visible','off','Color','white','Units','centimeters', ...
    'Position',[2 2 W/10 H/10],'PaperUnits','centimeters', ...
    'PaperSize',[W H]/10,'PaperPosition',[0 0 W H]/10, ...
    'InvertHardcopy','off','Renderer','painters');
cols={[0 0 1],[1 0 0],[1 0 0],[0 .5 0]};
marks={'o','^','^','s'}; labels={'A','B','B','C'};
pairnames={'Strict A/B','Core B/C'};
% Top row: actual background profiles, not the already fractured fields.
for j=1:2
    ax=newaxes(f,[20+90*(j-1) 76 66 42],W,H);
    patch(ax,[-30 -27 -27 -30],[0 0 1 1],[.94 .94 .94],'EdgeColor','none');
    patch(ax,[27 30 30 27],[0 0 1 1],[.94 .94 .94],'EdgeColor','none');
    h=gobjects(1,2);
    for n=1:2
        k=2*(j-1)+n;T=profiles{k};
        if j==2 && n==2,sty='--';else,sty='-';end
        h(n)=plot(ax,T.s_mm,T.background_V_V,sty,'Color',cols{k}, ...
            'LineWidth',1.15,'Marker',marks{k},'MarkerSize',3.4, ...
            'MarkerFaceColor','white','MarkerIndices',(1+40*(n-1)):100:height(T));
    end
    xlim(ax,[-30 30]);ylim(ax,[.15 .85]);
    set(ax,'XTick',[-30 -15 0 15 30],'YTick',[.2 .4 .6 .8]);
    xlabel(ax,'\its\rm (mm)');ylabel(ax,'Background potential, \itV\rm_0 (V)');
    title(ax,pairnames{j},'FontSize',8.5);
    legend(ax,h,labels(2*j-1:2*j),'Location','northwest','FontSize',7.5);
    text(ax,28,.275,sprintf('\\Delta\\itV\\rm_{0,cap} = %.4f V',caps(2*j-1)), ...
        'HorizontalAlignment','right','FontSize',7.2);
    text(ax,28,.205,sprintf('\\itU\\rm \\approx %.3f mS',1000*E.U_S(2*j-1)), ...
        'HorizontalAlignment','right','FontSize',7.2);
    panel(f,char('a'+j-1),5+90*(j-1),124,W,H);
end

% Small actual gains get their own scale instead of being hidden by U and R.
ax=newaxes(f,[20 13 66 42],W,H);
x=[1 2 4 5];
for k=1:4
    bar(ax,x(k),E.deltaG_S(k)*1000,.42,'FaceColor','white', ...
        'EdgeColor',cols{k},'LineWidth',1.15,'ShowBaseLine','off');
    text(ax,x(k),E.deltaG_S(k)*1000+.11,sprintf('%.3f',E.deltaG_S(k)*1000), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',7.1);
end
xlim(ax,[.3 5.7]);ylim(ax,[0 6.5]);
set(ax,'XTick',x,'XTickLabel',labels,'YTick',[0 2 4 6]);
text(ax,1.5,6.1,'Strict','HorizontalAlignment','center','FontSize',7.5);
text(ax,4.5,6.1,'Core','HorizontalAlignment','center','FontSize',7.5);
ylabel(ax,'Added conductance, \itg\rm (mS)');title(ax,'Actual gain','FontSize',8.5);
panel(f,'c',5,61,W,H);

% Hero panel: paired contributions, explicitly separating near-zero dU.
ax=newaxes(f,[110 13 66 42],W,H);
bb=[dU -dR dg]*1000;
hb=gobjects(1,2);paircolors={[0 0 1],[1 0 0]};
for j=1:2
    pos=(1:3)+(-.16+.32*(j-1));
    hb(j)=bar(ax,pos,bb(j,:),.28,'FaceColor','white','EdgeColor',paircolors{j}, ...
        'LineWidth',1.15,'ShowBaseLine','off');
    plot(ax,pos(1),bb(j,1),'o','Color',paircolors{j},'MarkerFaceColor','white', ...
        'MarkerSize',3.2,'LineWidth',.9);
end
plot(ax,[.5 3.5],[0 0],':k','LineWidth',.55);
xlim(ax,[.5 3.5]);ylim(ax,[-.10 2.05]);
set(ax,'XTick',1:3,'XTickLabel',{'\delta\itU','-\delta\itR','\delta\itg'},'YTick',[0 .5 1 1.5 2]);
ylabel(ax,'Paired difference (mS)');title(ax,'\delta\itg\rm = \delta\itU\rm - \delta\itR','FontSize',8.5);
lg=legend(ax,hb,{'Strict B - A','Core C - B'},'Location','northwest','FontSize',7.1);
if isprop(lg,'ItemTokenSize'),lg.ItemTokenSize=[12 8];end
panel(f,'d',95,61,W,H);

drawnow;stem=fullfile(outdir,'Fig03_energy_redistribution');
savefig(f,[stem '.fig']);
exportgraphics(f,[stem '.png'],'Resolution',250,'BackgroundColor','white');
exportgraphics(f,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');
print(f,[stem '.svg'],'-dsvg','-painters');print(f,[stem '.tiff'],'-dtiff','-r600');
close(f);disp(E(:,{'case_id','U_S','R_S','deltaG_S'}));disp(pairs);
fprintf('Saved %s\n',stem);
end

function ax=newaxes(f,mm,W,H)
ax=axes(f,'Units','normalized','Position',mm./[W H W H]);
hold(ax,'on');box(ax,'on');ax.PositionConstraint='innerposition';ax.Layer='top';
end

function panel(f,letter,x,y,W,H)
annotation(f,'textbox',[x/W (y-3.5)/H 12/W 5/H],'String',['(' letter ')'], ...
    'EdgeColor','none','Margin',0,'FontName','Times New Roman','FontSize',9, ...
    'FontAngle','italic','Color','black','VerticalAlignment','middle');
end
