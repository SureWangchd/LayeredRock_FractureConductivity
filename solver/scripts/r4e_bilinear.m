function summary=r4e_bilinear(root)
% Offline finite-change identities on one common finite-element space.
% Reads existing COMSOL matrices/solutions only. No model loading or solve.
% NN: lower near, upper near; FF: lower far, upper far.
% NF/FN change exactly one of these strip positions relative to FF/NN.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
out=fullfile(root,'runs','round4_extension_20261024','diagnostics');
if ~isfolder(out),mkdir(out);end
stamp=fullfile(out,'bilinear_complete.json');
if isfile(stamp),summary=jsondecode(fileread(stamp));return;end
sfValues=[0.1 0.159375 0.3];layouts={'FF','NF','FN','NN'};
mainRows=struct([]);edgeRows=struct([]);mappingRecords=struct([]);operatorRecords=struct([]);
fieldRecords=struct([]);
for roleCell={'A','B'}
 role=roleCell{1};d=struct;
 for il=1:numel(layouts)
  layout=layouts{il};folder=case_folder(root,role,layout);
  cache=fullfile(folder,'background_responses.mat');assert(isfile(cache),'Missing background cache: %s',cache);
  z=load(cache,'K0','DL','DH','u0','xm','G0','c','free','fixed');
  assert(z.c.out_of_plane_thickness_m==1,'Expected explicitly recorded t=1 m');
  assert(strcmp(z.c.load,'z'),'Expected cross-layer loading');
  z.folder=folder;
  if strcmp(layout,'FF'),targetXmesh=z.xm;end
  [z.u0,map,perm]=r3_map_dofs(z.u0,z.xm,targetXmesh);
  z.K0=z.K0(perm,perm);z.DL=z.DL(perm,perm);z.DH=z.DH(perm,perm);
  z.perm=perm;inverse=zeros(size(perm));inverse(perm)=1:numel(perm);
  z.free=sort(inverse(z.free));z.fixed=sort(inverse(z.fixed));
  record=struct('role',role,'layout',layout,'background_cache',cache,'mapping',map);
  mappingRecords=append_record(mappingRecords,record);
  d.(layout)=rmfield(z,'xm');
 end
 L=d.NF.K0-d.FF.K0;U=d.FN.K0-d.FF.K0;
 closure=d.NN.K0-d.FF.K0-L-U;
 closureRelative=norm(closure,'fro')/max(norm(L,'fro')+norm(U,'fro'),realmin);
 assert(closureRelative<1e-9,'Four-layout operator change does not separate into lower + upper');
 for il=1:numel(layouts)
  layout=layouts{il};z=d.(layout);
  lowRel=norm(z.DL-d.FF.DL,'fro')/max(norm(d.FF.DL,'fro'),realmin);
  highRel=norm(z.DH-d.FF.DH,'fro')/max(norm(d.FF.DH,'fro'),realmin);
  assert(lowRel<1e-10 && highRel<1e-10,'Ellipse phase operators changed across peripheral layouts');
  assert(isequal(z.free,d.FF.free) && isequal(z.fixed,d.FF.fixed),'Electrode constraints changed');
  record=struct('role',role,'layout',layout,'operator_sum_relative_error',closureRelative,...
   'ellipse_low_operator_relative_change',lowRel,'ellipse_high_operator_relative_change',highRel);
  operatorRecords=append_record(operatorRecords,record);
 end
 bgLower=d.FF.u0'*(L*d.NN.u0);bgUpper=d.FF.u0'*(U*d.NN.u0);
 bgActual=d.NN.G0-d.FF.G0;
 for sf=sfValues
  v=struct;g=struct;dg=struct;
  for il=1:numel(layouts)
   layout=layouts{il};z=d.(layout);target=target_folder(z.folder,sf);
   f=load(fullfile(target,'fields.mat'),'u');assert(isfield(f,'u'),'No saved cached-order truth DOF vector');
   u=f.u(:);assert(numel(u)==numel(z.perm),'Target/background DOF count differs');
   u=u(z.perm);truth=jsondecode(fileread(fullfile(target,'SUCCESS.json')));
   assert(abs(truth.sigma_f_S_m-sf)<=1e-13,'Target conductivity mismatch');
   assert(abs(truth.G0_S-z.G0)<1e-11,'Truth uses a different background');
   K=z.K0+(sf-z.c.sigma_low_S_m)*z.DL+(sf-z.c.sigma_high_S_m)*z.DH;
   energy=u'*(K*u);energyError=abs(energy-truth.G_S)/max(abs(truth.G_S),realmin);
   assert(energyError<1e-7,'Saved target vector and assembled operator disagree');
   r=K*u;residual=max(abs(r(z.free)));
   boundaryDifference=max(abs(u(z.fixed)-z.u0(z.fixed)));
   assert(boundaryDifference<1e-10,'Target and background imposed potentials differ');
   record=struct('role',role,'layout',layout,'sigma_f_S_m',sf,...
    'target_folder',target,'matrix_energy_relative_error',energyError,...
    'free_residual_inf_A',residual,'fixed_voltage_difference_V',boundaryDifference);
   fieldRecords=append_record(fieldRecords,record);
   v.(layout)=u;g.(layout)=truth.G_S;dg.(layout)=truth.deltaG_S;
  end
  targetLower=v.FF'*(L*v.NN);targetUpper=v.FF'*(U*v.NN);
  lowerGain=targetLower-bgLower;upperGain=targetUpper-bgUpper;
  actualTotal=g.NN-g.FF;actualGain=dg.NN-dg.FF;
  reconstructed=lowerGain+upperGain;absError=abs(reconstructed-actualGain);
  % Separate identity closure from mesh accuracy. Absolute floor protects a
  % nearly zero effect; all residuals are still retained in the output.
  assert(absError<max(1e-10,1e-4*abs(actualGain)),'Bilinear gain identity failed');
  record=struct('role',role,'sigma_f_S_m',sf,...
   'from_layout','FF','to_layout','NN',...
   'lower_target_bilinear_S',targetLower,'upper_target_bilinear_S',targetUpper,...
   'lower_background_bilinear_S',bgLower,'upper_background_bilinear_S',bgUpper,...
   'lower_increment_decomposition_S',lowerGain,'upper_increment_decomposition_S',upperGain,...
   'sum_increment_decomposition_S',reconstructed,'actual_total_G_difference_S',actualTotal,...
   'actual_background_G_difference_S',bgActual,'actual_increment_difference_S',actualGain,...
   'identity_absolute_error_S',absError,...
   'identity_relative_to_effect',absError/max(abs(actualGain),realmin),...
   'identity_relative_to_absolute_components',absError/max(abs(lowerGain)+abs(upperGain),realmin),...
   'lower_intervention_at_upper_far_S',dg.NF-dg.FF,...
   'lower_intervention_at_upper_near_S',dg.NN-dg.FN,...
   'upper_intervention_at_lower_far_S',dg.FN-dg.FF,...
   'upper_intervention_at_lower_near_S',dg.NN-dg.NF,...
   'finite_interaction_S',dg.NN-dg.NF-dg.FN+dg.FF);
  mainRows=append_record(mainRows,record);
  edges={'FF','NF','lower';'FN','NN','lower';'FF','FN','upper';'NF','NN','upper'};
  for ie=1:size(edges,1)
   p=edges{ie,1};q=edges{ie,2};side=edges{ie,3};dk=d.(q).K0-d.(p).K0;
   targetCross=v.(p)'*(dk*v.(q));backgroundCross=d.(p).u0'*(dk*d.(q).u0);
   observedTotal=g.(q)-g.(p);observedBackground=d.(q).G0-d.(p).G0;
   observedGain=dg.(q)-dg.(p);bilinearGain=targetCross-backgroundCross;
   edgeError=abs(bilinearGain-observedGain);
   assert(edgeError<max(1e-10,1e-4*abs(observedGain)),'Single-side edge identity failed');
   record=struct('role',role,'sigma_f_S_m',sf,'from_layout',p,'to_layout',q,...
    'moved_side',side,'target_bilinear_S',targetCross,'background_bilinear_S',backgroundCross,...
    'increment_bilinear_S',bilinearGain,'actual_total_G_difference_S',observedTotal,...
    'actual_background_G_difference_S',observedBackground,'actual_increment_difference_S',observedGain,...
    'total_identity_absolute_error_S',abs(targetCross-observedTotal),...
    'background_identity_absolute_error_S',abs(backgroundCross-observedBackground),...
    'increment_identity_absolute_error_S',edgeError);
   edgeRows=append_record(edgeRows,record);
  end
 end
 clear d L U closure targetXmesh;
end
writetable(struct2table(mainRows),fullfile(out,'bilinear_decomposition.csv'));
writetable(struct2table(edgeRows),fullfile(out,'bilinear_edge_checks.csv'));
writetable(struct2table(fieldRecords),fullfile(out,'bilinear_field_checks.csv'));
writetable(struct2table(operatorRecords),fullfile(out,'bilinear_operator_checks.csv'));
r3_json(fullfile(out,'bilinear_mapping.json'),mappingRecords);
summary=struct('status','complete','method','posthoc_two_solution_finite_change_identity',...
 'target_solve_count',0,'background_solve_count',0,'interpolation_used',false,...
 'sigma_f_S_m',sfValues,'decomposition_rows',numel(mainRows),'single_side_edge_rows',numel(edgeRows),...
 'maximum_increment_identity_absolute_error_S',max([mainRows.identity_absolute_error_S]),...
 'maximum_single_side_increment_identity_absolute_error_S',max([edgeRows.increment_identity_absolute_error_S]),...
 'units','t=1 m and DV=1 V included in K; bilinear products have conductance units S',...
 'interpretation','The lower/upper terms partition one jointly changed state pair. They are not independent intervention effects, percentage contributions, or predictions. Single-side counterfactual effects and their interaction are reported separately.');
r3_json(stamp,summary);disp(summary);
end

function rows=append_record(rows,row)
% Assign the first scalar record directly: a fieldless empty struct cannot
% reliably accept indexed assignment of a record with a new field schema.
assert(isscalar(row),'Expected one diagnostic record');
if isempty(rows),rows=row;else,rows(end+1)=orderfields(row,rows(1));end
end

function folder=case_folder(root,role,layout)
switch layout
 case 'NN',folder=fullfile(root,'runs','round4_crossover','periphery',['r4_strict_' role '_swap_near']);
 case 'FF',folder=fullfile(root,'runs','round4_crossover','periphery',['r4_strict_' role '_swap_far']);
 otherwise,folder=fullfile(root,'runs','round4_extension_20261024','cases',['r4e_strict_' role '_' layout]);
end
end

function folder=target_folder(parent,sf)
folder=fullfile(parent,['sf_' strrep(sprintf('%.10g',sf),'.','p')]);
assert(isfile(fullfile(folder,'SUCCESS.json')) && isfile(fullfile(folder,'fields.mat')),...
 'Missing complete target: %s',folder);
end
