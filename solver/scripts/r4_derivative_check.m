function r4_derivative_check(root)
% A small same-mesh central-difference check; no predictor changes.
out=fullfile(root,'runs','round4_crossover');
design=jsondecode(fileread(fullfile(out,'derivative_check_design.json')));entries=design.entries;
allRows=struct([]);
for j=1:numel(entries)
 e=entries(j);sf=e.sigma_f_S_m;role=char(e.role);folder=fullfile(out,'original',role);
 d=load(fullfile(folder,'background_responses.mat'));source=fullfile(folder,sf_token(sf));
 assert(isfile(fullfile(source,'solution.mph')),'Choose an already solved NEW original state for derivative check');
 m=mphload(fullfile(source,'solution.mph'),'R4FD');
 t=double(m.param.evaluate('d_out'));DV=double(m.param.evaluate('DV'));
 Sf=t/DV^2*mphint2(m,'ec.Ex^2+ec.Ey^2',2,'selection',d.sel.inclusion_geometry_domains,'intorder',8,'unit','V^2');
 checkRoot=fullfile(out,'diagnostics','finite_difference',role);
 for relativeStep=[.01 .005]
  h=sf*relativeStep;points=[sf-h sf+h];G=zeros(1,2);
  for k=1:2
   target=fullfile(checkRoot,sf_token(points(k)));
   r=r3_truth(m,d,points(k),target,'');G(k)=r.G_S;
   r3_json(fullfile(target,'diagnostic_purpose.json'),struct('purpose','Same-mesh central finite difference of actual COMSOL G, posthoc validation only',...
    'center_sigma_f_S_m',sf,'relative_step',relativeStep,'role',role,'predictor_changed',false));
  end
  fd=(G(2)-G(1))/(2*h);
  row=struct('role',role,'center_sigma_f_S_m',sf,'relative_step',relativeStep,'step_S_m',h,...
   'sigma_minus_S_m',points(1),'sigma_plus_S_m',points(2),'G_minus_S',G(1),'G_plus_S',G(2),...
   'Sf_integral_m',Sf,'central_difference_m',fd,'relative_difference',abs(fd-Sf)/abs(Sf),...
   'same_mesh',true,'source_center',source);
  if isempty(allRows),allRows=row;else,allRows(end+1)=row;end %#ok<AGROW>
  writetable(struct2table(allRows),fullfile(out,'derivative_check.csv'));
 end
 com.comsol.model.util.ModelUtil.clear;clear m d;
end
end

function s=sf_token(sf)
s=['sf_' strrep(sprintf('%.10g',sf),'.','p')];
end
