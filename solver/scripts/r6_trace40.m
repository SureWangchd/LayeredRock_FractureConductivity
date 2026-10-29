function [C,info]=r6_trace40(xm,nmax)
% Integrate actual quadratic COMSOL edge shape functions in SI, not samples.
xy=double(xm.dofs.coords);assert(max(xy(:))<1,'Expected SI xmesh coordinates');
ed=double(xm.elements.edg.dofs)+1;assert(size(ed,1)==3,'Quadratic scalar V edge expected');
nd=size(xy,2);C=sparse(nd,2*(nmax+1));W=.12;
ng=40;j=(1:ng-1)';beta=j./sqrt(4*j.^2-1);[Q,D]=eig(diag(beta,1)+diag(beta,-1));
[q,ix]=sort(diag(D));weights=2*Q(1,ix)'.^2;s=(q+1)/2;weights=weights/2;
shape=[2*(s-.5).*(s-1),4*s.*(1-s),2*s.*(s-.5)];
lengths=zeros(1,2);counts=zeros(1,2);midError=0;
for side=1:2
 z=[.03 .09];mask=all(abs(reshape(xy(2,ed),size(ed))-z(side))<1e-12,1);
 edges=ed(:,mask);counts(side)=size(edges,2);
 for e=1:size(edges,2)
  ids=edges(:,e);[xx,order]=sort(xy(1,ids));ids=ids(order);h=xx(3)-xx(1);
  assert(h>0,'Degenerate interface edge');midError=max(midError,abs(xx(2)-mean(xx([1 3]))));
  x=xx(1)+s*h;psi=cos(x*((0:nmax)*pi/W));local=h*(shape'*(weights.*psi));
  cols=(side-1)*(nmax+1)+(1:nmax+1);C(ids,cols)=C(ids,cols)+local;
  lengths(side)=lengths(side)+h;
 end
end
assert(max(abs(lengths-W))<1e-11 && midError<1e-12,'Interface partition is incomplete');
expected0=zeros(2*(nmax+1),1);expected0([1 nmax+2])=W;
constantError=max(abs(C'*ones(nd,1)-expected0));
expected2=zeros(2*(nmax+1),1);
for side=1:2
 for n=0:nmax
  if n==0,v=W^3/3;else,v=2*W*(-1)^n/(n*pi/W)^2;end
  expected2((side-1)*(nmax+1)+n+1)=v;
 end
end
quadraticError=max(abs(C'*(xy(1,:)'.^2)-expected2));
assert(constantError<1e-12 && quadraticError<1e-14,'Trace polynomial integration failed');
info=struct('method','40-point Gauss integration of actual quadratic COMSOL edge shapes',...
 'units','C in m; xmesh in m','nmax',nmax,'lengths_m',lengths,'edges',counts,...
 'constant_trace_absolute_error_m',constantError,'quadratic_trace_absolute_error_m3',quadraticError,...
 'maximum_midnode_error_m',midError,'rank',rank(full(C)));
end
