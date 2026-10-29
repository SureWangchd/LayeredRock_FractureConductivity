function figure06_independent_matlab(outdir,project_root)
% Plot saved, prospective Round-5 predictions and complete COMSOL targets.
% No COMSOL calls, no new solves, no adjustment of the frozen predictor.
if nargin<2 || isempty(project_root)
    project_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin<1 || isempty(outdir)
    outdir=fullfile(project_root,'outputs','JGRSE_20261028_Fig06_MainText');
end
idir=fullfile(outdir,'inputs');assert(isfolder(idir),'Prepared inputs are required.');
T=readtable(fullfile(idir,'prospective_comparison.csv'),'TextType','string');
R=readtable(fullfile(idir,'prospective_ranking.csv'),'TextType','string');
S=readtable(fullfile(idir,'single_side_64_states.csv'),'TextType','string');
AllS=readtable(fullfile(idir,'prospective_side_effects.csv'),'TextType','string');
Mesh=readtable(fullfile(idir,'selected_mesh_effects.csv'),'TextType','string');
assert(height(T)==144 && all(T.phase=="prospective"));
assert(max(abs(T.g_relative_error_percent-100*abs(T.g_prediction_S-T.g_COMSOL_S)./abs(T.g_COMSOL_S)))<1e-10);
for K=[0 8 32]
    assert(sum(T.K==K)==48 && sum(S.K==K)==64);
end
W=183;H=145;blue=[0 0 1];red=[1 0 0];
colors={[0 0 0],blue,red};Ks=[0 8 32];marks={'s','o','^'};
layouts=["P11","P77","P17","P71","P41","P14"];
sfs=[.1 .18 .3 10];
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman', ...
    'defaultAxesFontSize',8,'defaultTextFontSize',8,'defaultTextColor',[0 0 0], ...
    'defaultAxesXColor',[0 0 0],'defaultAxesYColor',[0 0 0], ...
    'defaultAxesLineWidth',.6,'defaultAxesBox','on','defaultAxesTickDir','in', ...
    'defaultAxesTickLength',[.018 .018],'defaultAxesTitleFontWeight','normal', ...
    'defaultAxesTitleFontSizeMultiplier',1.05,'defaultAxesLabelFontSizeMultiplier',1.05, ...
    'defaultLegendFontName','Times New Roman','defaultLegendFontSize',7, ...
    'defaultLegendTextColor',[0 0 0],'defaultLegendBox','off');
f=figure('Visible','off','Color','white','Units','centimeters', ...
    'Position',[2 2 W/10 H/10],'PaperUnits','centimeters','PaperSize',[W H]/10, ...
    'PaperPosition',[0 0 W H]/10,'InvertHardcopy','off','Renderer','painters');

% (a) New structural effects first: six categorical layouts, all retained.
ax=axesmm(f,[20 88 66 43],W,H);rank=table();x=1:6;yt=zeros(1,6);
for j=1:6
    z=R(R.layout==layouts(j) & R.sigma_f_S_m==.18,:);assert(height(z)==3);
    assert(max(z.D_COMSOL_S)-min(z.D_COMSOL_S)<1e-14);yt(j)=1e6*z.D_COMSOL_S(1);
    rank=[rank;z]; %#ok<AGROW>
end
plot(ax,[.5 6.5],[0 0],':k','LineWidth',.55);
h=gobjects(1,4);
h(1)=plot(ax,x,yt,'-ok','LineWidth',1.15,'MarkerSize',3.0,'MarkerFaceColor','black');
for m=1:3
    vals=zeros(1,6);
    for j=1:6
        z=rank(rank.layout==layouts(j)&rank.K==Ks(m),:);assert(height(z)==1);vals(j)=1e6*z.D_prediction_S;
    end
    if m==1
        h(m+1)=plot(ax,x,vals,'--k','LineWidth',1.05);
    else
        % Categorical offsets expose markers without changing the geometry.
        dx=.045*(2*m-5);
        h(m+1)=plot(ax,x+dx,vals,marks{m},'Color',colors{m},'LineStyle','none', ...
            'LineWidth',.95,'MarkerSize',3.8,'MarkerFaceColor','white');
    end
end
xlim(ax,[.5 6.5]);ylim(ax,[-45 15]);
set(ax,'XTick',1:6,'XTickLabel',{'(1,1)','(7,7)','(1,7)','(7,1)','(4,1)','(1,4)'},'YTick',[-40 -20 0]);
xlabel(ax,'(\itd\rm_-, \itd\rm_+) (mm)');ylabel(ax,'\itg\rm_B - \itg\rm_A (\muS)');
title(ax,'New-layout ranking','FontSize',8.5);
text(ax,6.2,10.5,'\sigma_f = 0.18 S/m','HorizontalAlignment','right','FontSize',7.1);
lg=legend(ax,h,{'COMSOL','M0','M8','M32'},'Location','southwest','NumColumns',2,'FontSize',7.1,'AutoUpdate','off');
if isprop(lg,'ItemTokenSize'),lg.ItemTokenSize=[13 8];end
panel(f,'a',5,136,W,H);writetable(rank,fullfile(outdir,'plotted_rank_sigma018.csv'));

% (b) All 48 target errors for each method, not aggregate-only bars.
ax=axesmm(f,[110 88 66 43],W,H);set(ax,'YScale','log');
hb=gobjects(1,3);plotted=table();
for m=1:3
    for j=1:4
        z=T(T.K==Ks(m)&T.sigma_f_S_m==sfs(j),:);z=sortrows(z,{'role','layout'});assert(height(z)==12);
        xp=j+(m-2)*.22+linspace(-.065,.065,12)';
        hb(m)=plot(ax,xp,z.g_relative_error_percent,marks{m},'Color',colors{m}, ...
            'LineStyle','none','LineWidth',.75,'MarkerSize',2.8,'MarkerFaceColor','white');
        z.display_x=xp;plotted=[plotted;z]; %#ok<AGROW>
    end
end
xlim(ax,[.5 4.5]);ylim(ax,[1e-5 25]);
set(ax,'XTick',1:4,'XTickLabel',{'0.1','0.18','0.3','10'},'YTick',10.^[-5 -3 -1 1]);
xlabel(ax,'\sigma_f (S m^{-1})');ylabel(ax,'Gain error, \itE\rm_g (%)');
title(ax,'All 48 new states','FontSize',8.5);
lg=legend(ax,hb,{'M0','M8','M32'},'Location','northwest','NumColumns',3,'FontSize',7.1,'AutoUpdate','off');
if isprop(lg,'ItemTokenSize'),lg.ItemTokenSize=[7 7];end
panel(f,'b',95,136,W,H);writetable(plotted,fullfile(outdir,'plotted_gain_errors_all48.csv'));

% (c) Magnitudes on logarithmic axes; signed effects remain in source CSV.
ax=axesmm(f,[20 17 66 44],W,H);set(ax,'XScale','log','YScale','log');
plot(ax,[.025 200],[.025 200],'--k','LineWidth',.65);
hc=gobjects(1,4);slot=0;
% Larger red open symbols first, with smaller blue symbols visible inside.
for K=[32 8]
    if K==8,col=blue;ms=3.0;else,col=red;ms=4.0;end
    for role=["gA","gB"]
        if role=="gA",mk='o';else,mk='^';end
        z=S(S.K==K & S.quantity==role,:);assert(height(z)==32);
        xx=1e6*abs(z.COMSOL_S);yy=1e6*abs(z.prediction_S);assert(all(xx>0&yy>0));
        slot=slot+1;hc(slot)=plot(ax,xx,yy,mk,'Color',col,'LineStyle','none', ...
            'LineWidth',.8,'MarkerSize',ms,'MarkerFaceColor','none');
    end
end
xlim(ax,[.025 200]);ylim(ax,[.025 200]);
set(ax,'XTick',[.1 1 10 100],'YTick',[.1 1 10 100]);
xlabel(ax,'COMSOL |\delta\itg\rm| (\muS)');ylabel(ax,'Predicted |\delta\itg\rm| (\muS)');
title(ax,'Single-side effects','FontSize',8.5);
lg=legend(ax,hc([3 4 1 2]),{'A, M8','B, M8','A, M32','B, M32'}, ...
    'Location','southeast','FontSize',6.8,'NumColumns',1,'AutoUpdate','off');
if isprop(lg,'ItemTokenSize'),lg.ItemTokenSize=[7 7];end
panel(f,'c',5,66,W,H);
writetable(S,fullfile(outdir,'plotted_single_side_effects.csv'));

% Inset shows the documented single-side failure on linear microS axes.
zi=S(S.quantity=="gB" & S.sigma_f_S_m==.1 & S.effect=="upper_4_to_1_at_lower_1",:);
assert(height(zi)==3);
ai=axesmm(f,[27 39 22 17],W,H);plot(ai,[1 1.27],[1 1.27],'--k','LineWidth',.55);
for K=[32 8]
    z=zi(zi.K==K,:);if K==8,col=blue;ms=3.7;else,col=red;ms=3.7;end
    plot(ai,z.COMSOL_S*1e6,z.prediction_S*1e6,'^','Color',col,'MarkerSize',ms, ...
        'MarkerFaceColor','white','LineWidth',.9);
end
xlim(ai,[1 1.27]);ylim(ai,[1 1.27]);set(ai,'XTick',[1 1.2],'YTick',[1 1.2],'FontSize',6.2);
title(ai,'B, upper','FontSize',6.7);
bad=zi(zi.K==8,:);text(ai,1.025,1.235,sprintf('%.2f%%',bad.relative_error_percent),'FontSize',6.5);
writetable(zi,fullfile(outdir,'inset_single_side_failure.csv'));

% (d) Net effect and an observed mesh scale, without a sign claim.
ax=axesmm(f,[110 17 66 44],W,H);
net=AllS(AllS.quantity=="gB" & AllS.sigma_f_S_m==.18 & AllS.effect=="joint_P11_minus_P77",:);
mr=Mesh(Mesh.role=="B" & Mesh.sigma_f_S_m==.18 & Mesh.effect=="both_near_minus_both_far",:);
assert(height(net)==3 && height(mr)==1 && lower(string(mr.resolved_on_3x_observed_mesh_scale))=="false");
assert(max(abs(net.COMSOL_S*1e6-mr.coarse_effect_uS))<1e-8);
plot(ax,[.5 4.5],[0 0],':k','LineWidth',.55);
coarse=plot(ax,.94,mr.coarse_effect_uS,'ok','MarkerSize',4.0,'MarkerFaceColor','white','LineWidth',.9);
errorbar(ax,1.06,mr.fine_effect_uS,mr.sum_abs_underlying_gain_mesh_changes_uS, ...
    'Color','black','LineStyle','none','LineWidth',.7,'CapSize',4);
fine=plot(ax,1.06,mr.fine_effect_uS,'ok','MarkerSize',3.3,'MarkerFaceColor','black','LineWidth',.9);
for m=1:3
    z=net(net.K==Ks(m),:);plot(ax,m+1,z.prediction_S*1e6,marks{m}, ...
        'Color',colors{m},'MarkerFaceColor','white','MarkerSize',4.0,'LineWidth',1.0);
end
xlim(ax,[.5 4.5]);ylim(ax,[-.32 .16]);
set(ax,'XTick',1:4,'XTickLabel',{'COMSOL','M0','M8','M32'},'YTick',[-.3 -.2 -.1 0 .1]);
ylabel(ax,'Net change, \delta\itg\rm (\muS)');title(ax,'Near cancellation','FontSize',8.5);
text(ax,1.55,.125,'Unresolved','HorizontalAlignment','left','FontSize',7.0);
% Point only to the COMSOL near-zero net change and its observed mesh scale.
arrow_x=110+([1.60 1.13]-.5)/4*66;
arrow_y=17+([.107 .004]+.32)/.48*44;
annotation(f,'arrow',arrow_x/W,arrow_y/H,'Color','black', ...
    'LineWidth',.55,'HeadWidth',3,'HeadLength',3);
text(ax,.7,-.292,'B, P77 \rightarrow P11, \sigma_f = 0.18 S/m','FontSize',6.7);
lg=legend(ax,[coarse fine],{'Base','Fine'},'Location','northeast','FontSize',7,'AutoUpdate','off');
if isprop(lg,'ItemTokenSize'),lg.ItemTokenSize=[7 7];end
panel(f,'d',95,66,W,H);
writetable(net,fullfile(outdir,'plotted_net_effect_predictions.csv'));
writetable(mr,fullfile(outdir,'plotted_net_mesh_scale.csv'));

drawnow;stem=fullfile(outdir,'Fig06_independent_prediction');savefig(f,[stem '.fig']);
exportgraphics(f,[stem '.png'],'Resolution',250,'BackgroundColor','white');
exportgraphics(f,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');
print(f,[stem '.svg'],'-dsvg','-painters');print(f,[stem '.tiff'],'-dtiff','-r600');
close(f);fprintf('Saved %s\n',stem);
end

function ax=axesmm(f,mm,W,H)
ax=axes(f,'Units','normalized','Position',mm./[W H W H]);hold(ax,'on');box(ax,'on');
ax.PositionConstraint='innerposition';ax.Layer='top';
end
function panel(f,s,x,y,W,H)
annotation(f,'textbox',[x/W (y-3.5)/H 12/W 5/H],'String',['(' s ')'], ...
    'EdgeColor','none','Margin',0,'FontSize',9,'FontAngle','italic', ...
    'FontName','Times New Roman','Color','black','VerticalAlignment','middle');
end
