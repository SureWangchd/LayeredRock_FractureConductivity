function [rows,chosen]=r3_predict(d,sigmas,methods)
% All inputs are background-only. Never read a full target solution here.
if nargin<3,methods={'field','vector'};end
rows=struct([]);chosen=struct;Kf=d.K0(d.free,d.free);
for im=1:numel(methods)
 method=methods{im};
 for ns=[2 4 8 16]
  tp=tic;P=[];
  for j=1:ns
   for ph=0:1
    blocks=find(ceil(d.blockinfo(:,1)/(d.nseg/ns))==j & d.blockinfo(:,2)==ph);
    if strcmp(method,'field'),dirs=1;else,dirs=1:3;end
    for dir=dirs,P(:,end+1)=sum(d.responses(:,3*blocks-3+dir),2);end %#ok<AGROW>
   end
  end
  diagE=sum(P.*(Kf*P),1);valid=diagE>max(diagE)*1e-16;P=P(:,valid);
  P=P./sqrt(diagE(valid));C=P'*(Kf*P);C=(C+C')/2;
  [Q,E]=eig(C,'vector');keep=E>max(E)*1e-10;
  Z=P*(Q(:,keep)./sqrt(E(keep)'));rank=size(Z,2);Zall=zeros(numel(d.u0),rank);Zall(d.free,:)=Z;
  AL=Zall'*(d.DL*Zall);AH=Zall'*(d.DH*Zall);
  qL=Zall'*(d.DL*d.u0);qH=Zall'*(d.DH*d.u0);
  UL=d.u0'*(d.DL*d.u0);UH=d.u0'*(d.DH*d.u0);
  H0=Z'*(Kf*Z);basisSeconds=toc(tp);
  for sf=sigmas
   tp=tic;dl=sf-d.c.sigma_low_S_m;dh=sf-d.c.sigma_high_S_m;
   H=H0+dl*AL+dh*AH;q=dl*qL+dh*qH;alpha=-(H\q);w=Zall*alpha;v=d.u0+w;
   U=dl*UL+dh*UH;R=-q'*alpha;raw=U-R;
   % Equivalent nonnegative form avoids subtracting large U and R.
   stable=w'*(d.K0*w)+dl*(v'*(d.DL*v))+dh*(v'*(d.DH*v));
   r=struct('case_id',d.c.id,'method',method,'segments',ns,'candidate_columns',size(P,2),...
    'rank',rank,'sigma_f_S_m',sf,'G0_S',d.G0,'U_S',U,'predicted_deltaG_S',stable,...
    'subtractive_deltaG_S',raw,'R_pred_S',R,'background_endcap_dV_V',diff(d.capV),...
    'basis_projection_seconds',basisSeconds,'reduced_evaluation_seconds',toc(tp));
   if isempty(rows),rows=r;else,rows(end+1)=r;end %#ok<AGROW>
   if strcmp(method,'vector')&&ns==16
    key=['s' strrep(num2str(sf),'.','p')];chosen.(key)=struct('row',r,'w',w,'alpha',alpha);
   end
  end
 end
end
end
