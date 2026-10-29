function fig08_geometry_tgrs(ax,p,variant)
% Exact saved specimen / finite ellipsoid geometry. Only matrix display is
% cut away; no domain is changed in the COMSOL model or result tables.
if nargin<3,variant=1;end
hold(ax,'on');axis(ax,'equal');axis(ax,'off');
views=[115 20;125 25;105 18];view(ax,views(variant,1),views(variant,2));ax.Projection='orthographic';
% Retained matrix display consists of two rear slabs. Their exposed top
% faces establish depth, as in the TGRS Fig1(b) solid/cutaway presentation.
xc=22;yc=22;
for k=1:numel(p.physical_material_runs)
 r=p.physical_material_runs(k);lo=r.z_min_mm;hi=r.z_max_mm;
 col=[1 1 1];if r.phase==1,col=[.753 .753 .753];end
 % Inner cut faces, outer faces, and the two narrow exposed ends.
 face(ax,[xc xc xc xc],[yc 120 120 yc],[lo lo hi hi],col,.34);
 face(ax,[xc 120 120 xc],[yc yc yc yc],[lo lo hi hi],col*.92,.34);
 face(ax,[0 0 0 0],[0 120 120 0],[lo lo hi hi],col*.88,.22);
 face(ax,[0 120 120 0],[0 0 0 0],[lo lo hi hi],col*.86,.22);
 face(ax,[0 xc xc 0],[120 120 120 120],[lo lo hi hi],col*.96,.22);
 face(ax,[120 120 120 120],[0 yc yc 0],[lo lo hi hi],col*.83,.22);
end
% Real uppermost phase is high-conductivity matrix. The display cut does
% not truncate the fracture: ellipse min x>44.9 mm, min y=30 mm >22 mm.
face(ax,[0 xc xc 0],[0 0 120 120],[120 120 120 120],[.82 .82 .82],.4);
face(ax,[xc 120 120 xc],[0 0 yc yc],[120 120 120 120],[.82 .82 .82],.4);
face(ax,[xc 120 120 xc],[yc yc 120 120],[0 0 0 0],[1 1 1],.35);
[xx,yy,zz]=ellipsoid(0,0,0,p.a_mm,p.b_mm,p.c_mm,180);
R=[p.basis_u(:),p.basis_v(:),p.basis_n(:)];
X=p.cx_mm+R(1,1)*xx+R(1,2)*yy+R(1,3)*zz;
Y=p.cy_mm+R(2,1)*xx+R(2,2)*yy+R(2,3)*zz;
Z=p.cz_mm+R(3,1)*xx+R(3,2)*yy+R(3,3)*zz;
surf(ax,X,Y,Z,'FaceColor',[0 .5 0],'EdgeColor','none', ...
 'FaceLighting','gouraud','AmbientStrength',.18,'DiffuseStrength',.72, ...
 'SpecularStrength',.75,'SpecularExponent',24,'SpecularColorReflectance',.65, ...
 'BackFaceLighting','reverselit');
light(ax,'Position',[130 -115 -100],'Style','local','Color',[1 1 1]);
light(ax,'Position',[180 150 120],'Style','local','Color',[.28 .28 .28]);
% Cube outline; hidden remote edges lighter and dashed. Do not draw an
% artificial edge through the inclusion. True equal physical axis scales.
P=[0 0 0;120 0 0;120 120 0;0 120 0;0 0 120;120 0 120;120 120 120;0 120 120];
E=[1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;1 5;2 6;4 8];
for k=1:size(E,1),t=P(E(k,:),:);plot3(ax,t(:,1),t(:,2),t(:,3),'Color',[.25 .25 .25],'LineWidth',.45);end
% Core limits stay on visible cut faces, with no line across the ellipsoid.
for z=[30 90]
 plot3(ax,[xc xc 120],[120 yc yc],[z z z],'--k','LineWidth',.65);
end
% Keep both annotations in the outside white margin. Their positions do
% not alter the camera/limits or the actual specimen dimensions.
text(ax,64,129,-12,'120 mm','HorizontalAlignment','center','FontSize',7, ...
 'Color','black','Clipping','off');
o=[133 -35 -7];L=21;
for j=1:3
 d=zeros(1,3);d(j)=L;quiver3(ax,o(1),o(2),o(3),d(1),d(2),d(3),0,'k', ...
  'LineWidth',.6,'MaxHeadSize',.3,'Clipping','off');
 label='xyz';pos=o+1.23*d;
 ha={'right','left','center'};va={'top','middle','bottom'};
 text(ax,pos(1),pos(2),pos(3),label(j),'FontAngle','italic','FontSize',7, ...
  'Color','black','HorizontalAlignment',ha{j},'VerticalAlignment',va{j},'Clipping','off');
end
xlim(ax,[-5 165]);ylim(ax,[-5 145]);zlim(ax,[-15 138]);daspect(ax,[1 1 1]);
% Exact apparent contour of the unchanged ellipsoid for an orthographic
% camera. This adds a thin rim, not a wireframe or an inflated thickness.
dv=(ax.CameraPosition-ax.CameraTarget)';dv=dv/norm(dv);
abc=[p.a_mm;p.b_mm;p.c_mm];m=(R'*dv)./abc;Q=null(m');
t=linspace(0,2*pi,400);q=Q*[cos(t);sin(t)];
xyz=[p.cx_mm;p.cy_mm;p.cz_mm]+R*(abc.*q);
plot3(ax,xyz(1,:),xyz(2,:),xyz(3,:),'Color',[0 .20 0],'LineWidth',.4);
end
function face(ax,x,y,z,c,lw)
patch(ax,x,y,z,c,'EdgeColor',[.22 .22 .22],'LineWidth',lw,'FaceLighting','none');
end
