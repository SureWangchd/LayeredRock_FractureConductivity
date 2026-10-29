function result=r3_refine_case(sourceTargetFolder,outputFolder)
% Targeted finer-mesh check of ONE already validated Round 3 target state.
% Called explicitly by the main task, serially. Connects to an existing server;
% never launches one. No predictor matrices/bases or old-case sweeps are built.
% Both fine target and fine passive background share the newly refined mesh.
sourceTargetFolder=char(sourceTargetFolder);outputFolder=char(outputFolder);
srcCanonical=char(java.io.File(sourceTargetFolder).getCanonicalPath);
outCanonical=char(java.io.File(outputFolder).getCanonicalPath);
assert(~strcmpi(srcCanonical,outCanonical),'Refinement output must differ from source');
if ~isfolder(outputFolder),mkdir(outputFolder);end
sourceModel=fullfile(sourceTargetFolder,'solution.mph');
assert(isfile(sourceModel),'A newly solved target solution.mph is required');
coarse=jsondecode(fileread(fullfile(sourceTargetFolder,'SUCCESS.json')));
parameters=jsondecode(fileread(fullfile(sourceTargetFolder,'parameters.json')));
assert(strcmp(coarse.status,'success'),'Source target is not successfully solved');
stamp=fullfile(outputFolder,'SUCCESS.json');
if isfile(stamp)
 result=jsondecode(fileread(stamp));
 assert(strcmpi(result.source_target_folder,srcCanonical),'Existing output belongs to another source');
 fprintf('R3_REFINE_SKIP %s sf %.6g\n',result.case_id,result.sigma_f_S_m);return;
end
% Load only small provenance/selection variables, never the background matrices.
backgroundFolder=fileparts(sourceTargetFolder);
for level=1:5
 cache=fullfile(backgroundFolder,'background_responses.mat');
 if isfile(cache),break;end
 parent=fileparts(backgroundFolder);assert(~strcmp(parent,backgroundFolder),'Background cache not found');
 backgroundFolder=parent;
end
assert(isfile(cache),'No background_responses.mat ancestor was found');
small=load(cache,'c','sel','G0');c=small.c;sel=small.sel;
assert(strcmp(c.id,coarse.case_id),'Background cache and source target case IDs differ');
sf=coarse.sigma_f_S_m;
assert(abs(parameters.sigma_frac_S_m-sf)<1e-12*max(1,sf),'Source target conductivity mismatch');
base='D:/Program Files/COMSOL/COMSOL63/Multiphysics';addpath(fullfile(base,'mli'));
import com.comsol.model.util.*
owned=false;
try,mphstart('127.0.0.1',2037,base);owned=true;
catch err,try,ModelUtil.tags;catch,rethrow(err);end,end
connectionCleanup=onCleanup(@()disconnect_if_owned(owned)); %#ok<NASGU>
tag=['R3Refine_' strrep(char(java.util.UUID.randomUUID),'-','')];
modelCleanup=onCleanup(@()remove_owned_model(tag)); %#ok<NASGU>
started=tic;
targetCheckpoint=fullfile(outputFolder,'fine_target_metrics.json');
targetMph=fullfile(outputFolder,'fine_target.mph');
backgroundCheckpoint=fullfile(outputFolder,'fine_background_metrics.json');
backgroundMph=fullfile(outputFolder,'fine_background.mph');
if isfile(targetCheckpoint)&&isfile(targetMph)
 target=jsondecode(fileread(targetCheckpoint));
 assert(strcmpi(target.source_target_folder,srcCanonical),'Target checkpoint belongs to another source');
 model=mphload(targetMph,tag);
 meshInfo=jsondecode(fileread(fullfile(outputFolder,'mesh_refinement.json')));
 fprintf('R3_REFINE_RESUME_BACKGROUND %s sf %.6g\n',c.id,sf);
else
 model=mphload(sourceModel,tag);
 assert(model.component('comp1').geom('geom1').getSDim==2,'Only 2D refinement is allowed');
 mesh=model.component('comp1').mesh('mesh1');
 globalSize=mesh.feature('size');localSize=mesh.feature('tri').feature('fs');
 oldGlobal=char(globalSize.getString('hmax'));oldLocal=char(localSize.getString('hmax'));
 [oldStats,~]=mphmeshstats(model,'mesh1','qualitymeasure','skewness');
 oldTriangles=triangle_count(oldStats);
 assert(strcmp(char(model.component('comp1').geom('geom1').lengthUnit),'mm'), ...
  'Mesh sizes here are specified in geometry mm');
 globalSize.set('custom','on');globalSize.set('hmax','0.8');
 localSize.set('custom','on');localSize.set('hmax','0.10');localSize.set('hmaxactive',true);
 meshTimer=tic;mesh.run;meshSeconds=toc(meshTimer);
 [ms,md]=mphmeshstats(model,'mesh1','qualitymeasure','skewness');
 triangles=triangle_count(ms);
 assert(triangles>oldTriangles,'Requested mesh refinement did not increase triangle count');
 save(fullfile(outputFolder,'fine_mesh.mat'),'ms','md','-v7.3');
 r3_json(fullfile(outputFolder,'fine_mesh_stats.json'),ms);
 meshInfo=struct('source_target_folder',srcCanonical,'geometry_unit','mm', ...
  'old_global_hmax_expression',oldGlobal,'old_ellipse_hmax_expression',oldLocal, ...
  'new_global_hmax_mm',0.8,'new_ellipse_hmax_mm',0.10, ...
  'source_triangles',oldTriangles,'fine_triangles',triangles, ...
  'fine_min_mesh_quality',ms.minquality,'fine_mean_mesh_quality',ms.meanquality, ...
  'mesh_seconds',meshSeconds,'geometry_rebuilt',false,'predictor_basis_rebuilt',false);
 r3_json(fullfile(outputFolder,'mesh_refinement.json'),meshInfo);
 set_target(model,c,sf);
 fprintf('R3_REFINE_TARGET_START %s sf %.6g triangles %d -> %d\n',c.id,sf,oldTriangles,triangles);
 solveTimer=tic;model.sol('sol1').runAll;solveSeconds=toc(solveTimer);
 target=measure_state(model,sel,'target',solveSeconds,triangles);
 target.source_target_folder=srcCanonical;target.sigma_f_S_m=sf;
 save_fields(model,fullfile(outputFolder,'fine_target_fields.mat'));
 mphsave(model,targetMph);r3_json(targetCheckpoint,target);
end
if isfile(backgroundCheckpoint)&&isfile(backgroundMph)
 background=jsondecode(fileread(backgroundCheckpoint));
 assert(strcmpi(background.source_target_folder,srcCanonical),'Background checkpoint belongs to another source');
else
 model.param.set('r3bgscale','1');model.param.set('r3unit','0[S/m]');
 model.param.set('r3dlow','0[S/m]');model.param.set('r3dhigh','0[S/m]');
 % No second mesh.run: the background uses the very same target triangulation.
 fprintf('R3_REFINE_BACKGROUND_START %s\n',c.id);
 solveTimer=tic;model.sol('sol1').runAll;solveSeconds=toc(solveTimer);
 background=measure_state(model,sel,'passive_background',solveSeconds,meshInfo.fine_triangles);
 background.source_target_folder=srcCanonical;
 save_fields(model,fullfile(outputFolder,'fine_background_fields.mat'));
 mphsave(model,backgroundMph);r3_json(backgroundCheckpoint,background);
end
fineDelta=target.G_S-background.G_S;
coarseDelta=coarse.deltaG_S;gainDen=max(abs(coarseDelta),1e-12);
Lz=c.L_mm*1e-3;Lx=Lz;
result=struct('status','complete','case_id',c.id,'sigma_f_S_m',sf, ...
 'source_target_folder',srcCanonical,'source_model',sourceModel, ...
 'background_cache',cache,'loaded_cache_variables',{{'c','sel','G0'}}, ...
 'coarse_G_S',coarse.G_S,'coarse_G0_S',coarse.G0_S,'coarse_deltaG_S',coarseDelta, ...
 'fine_G_S',target.G_S,'fine_G0_S',background.G_S,'fine_deltaG_S',fineDelta, ...
 'fine_relative_gain',fineDelta/background.G_S, ...
 'fine_minus_coarse_G_S',target.G_S-coarse.G_S, ...
 'G_relative_change',abs(target.G_S-coarse.G_S)/max(abs(coarse.G_S),1e-30), ...
 'fine_minus_coarse_G0_S',background.G_S-coarse.G0_S, ...
 'G0_relative_change',abs(background.G_S-coarse.G0_S)/max(abs(coarse.G0_S),1e-30), ...
 'fine_minus_coarse_deltaG_S',fineDelta-coarseDelta, ...
 'deltaG_relative_change',abs(fineDelta-coarseDelta)/gainDen, ...
 'deltaG_difference_normalized_by_background',abs(fineDelta-coarseDelta)/max(abs(background.G_S),1e-30), ...
 'near_zero_gain_denominator_floor_S',1e-12, ...
 'fine_sigma_eff_S_m',target.G_S*Lz/(Lx*target.out_of_plane_thickness_m), ...
 'fine_sigma0_eff_S_m',background.G_S*Lz/(Lx*background.out_of_plane_thickness_m), ...
 'source_triangles',meshInfo.source_triangles,'fine_triangles',meshInfo.fine_triangles, ...
 'target',target,'background',background,'mesh',meshInfo, ...
 'predictor_recomputed',false,'elapsed_this_invocation_seconds',toc(started));
result.conservation_pass=target.conservation_pass&&background.conservation_pass;
result.power_balance_pass=target.power_balance_pass&&background.power_balance_pass;
if result.conservation_pass&&result.power_balance_pass,result.status='success';
else,result.status='completed_needs_balance_review';end
result.interpretation='One finer mesh trend for an already chosen target. Fine deltaG uses its own same-mesh background; no coarse/fine DOF subtraction, prediction refit, or new basis construction.';
parameters.refinement=meshInfo;
r3_json(fullfile(outputFolder,'parameters.json'),parameters);
r3_json(fullfile(outputFolder,'selections.json'),sel);
r3_json(fullfile(outputFolder,'refinement_result.json'),result);
row=struct('case_id',c.id,'sigma_f_S_m',sf,'coarse_G_S',coarse.G_S,'fine_G_S',target.G_S, ...
 'coarse_G0_S',coarse.G0_S,'fine_G0_S',background.G_S, ...
 'coarse_deltaG_S',coarseDelta,'fine_deltaG_S',fineDelta, ...
 'G_relative_change',result.G_relative_change,'deltaG_relative_change',result.deltaG_relative_change, ...
 'deltaG_difference_normalized_by_background',result.deltaG_difference_normalized_by_background, ...
 'coarse_triangles',meshInfo.source_triangles,'fine_triangles',meshInfo.fine_triangles, ...
 'target_flux_relative_error',target.reaction_flux_relative_error,'target_energy_relative_error',target.energy_relative_error, ...
 'background_flux_relative_error',background.reaction_flux_relative_error,'background_energy_relative_error',background.energy_relative_error);
writetable(struct2table(row),fullfile(outputFolder,'refinement_comparison.csv'));
if strcmp(result.status,'success'),r3_json(stamp,result);end
fprintf('R3_REFINE_DONE %s sf %.6g deltaG_change %.6g G_change %.6g status %s\n', ...
 c.id,sf,result.deltaG_relative_change,result.G_relative_change,result.status);
end

function set_target(model,c,sf)
model.param.set('r3bgscale','1');model.param.set('r3unit','0[S/m]');
model.param.set('r3dlow',sprintf('%.17g[S/m]',sf-c.sigma_low_S_m));
model.param.set('r3dhigh',sprintf('%.17g[S/m]',sf-c.sigma_high_S_m));
end

function state=measure_state(model,sel,kind,solveSeconds,triangles)
DV=double(model.param.evaluate('DV'));thickness=double(model.param.evaluate('d_out'));
pos=mphint2(model,'reacf(V)',1,'selection',sel.positive_electrode,'method','summation');
neg=mphint2(model,'reacf(V)',1,'selection',sel.ground_electrode,'method','summation');
powerPerM=mphint2(model,'ec.Jx*ec.Ex+ec.Jy*ec.Ey',2,'intorder',8,'unit','W/m');
power=thickness*powerPerM;sourcePower=abs(pos*DV);G=abs(pos/DV);
state=struct('state',kind,'G_S',G,'reaction_positive_A',pos,'reaction_ground_A',neg, ...
 'DV_V',DV,'out_of_plane_thickness_m',thickness,'power_source_W',sourcePower, ...
 'power_volume_W',power,'power_volume_W_per_m',powerPerM, ...
 'reaction_flux_relative_error',abs(pos+neg)/max(abs(pos),1e-30), ...
 'energy_relative_error',abs(power-sourcePower)/max(abs(sourcePower),1e-30), ...
 'triangles',triangles,'solve_seconds',solveSeconds);
state.conservation_pass=state.reaction_flux_relative_error<=1e-6;
state.power_balance_pass=state.energy_relative_error<=1e-6;
end

function save_fields(model,path)
fields=mpheval(model,{'V','ec.Ex','ec.Ey','ec.Jx','ec.Jy','ec.Jx*ec.Ex+ec.Jy*ec.Ey'}, ...
 'edim',2,'refine',1,'smooth','none');
u=mphgetu(model);xm=mphxmeshinfo(model);
field_units=struct('p_coordinates','geometry length unit mm','V','V','Ex_Ez','V/m', ...
 'Jx_Jz','A/m^2','dissipation_density','W/m^3','coordinate_mapping','COMSOL x,y = physical x,z');
save(path,'fields','u','xm','field_units','-v7.3');
end

function n=triangle_count(stats)
index=find(strcmp(stats.types,'tri'),1);assert(~isempty(index),'No triangular elements');
n=double(stats.numelem(index));
end

function remove_owned_model(tag)
try,com.comsol.model.util.ModelUtil.remove(tag);catch,end
end

function disconnect_if_owned(owned)
if owned,try,com.comsol.model.util.ModelUtil.disconnect;catch,end,end
end
