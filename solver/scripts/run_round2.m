function run_round2(group)
% Extend the established EC/MUMPS workflow; preserve the first 67 solutions.
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'runs','2026-10-21_round2');
manifest=fullfile(out,'manifests',[group '.json']);
cases=jsondecode(fileread(manifest));
if ~iscell(cases),cases=num2cell(cases);end
diary(fullfile(out,['suite_' group '_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '.log']));
base='D:/Program Files/COMSOL/COMSOL63/Multiphysics';addpath(fullfile(base,'mli'));
mphstart('127.0.0.1',2037,base);
import com.comsol.model.util.*
cleanup=onCleanup(@() ModelUtil.disconnect);
fprintf('COMSOL %s; MATLAB %s; manifest %s\n',char(ModelUtil.getComsolVersion),version,manifest);
defaults=struct('background','layered','A_mm',0,'phase_sign',1,'buffer',true,...
 'fracture',true,'passive_ellipse',false,'cx_mm',60,'cz_mm',60,'a_mm',30,'b_mm',1.5,'theta_deg',60,...
 'sigma_low_S_m',0.01,'sigma_high_S_m',0.1,'sigma_frac_S_m',10,...
 'mesh_level',1,'load','z','intervention','none','cut_u_mm',0,'cut_width_mm',1,...
 'side_u1_mm',-2,'side_u2_mm',2,'side_width_mm',0,...
 'comp_radius_mm',0,'comp_x_mm',20,'comp_z_mm',45,'intervention_phase',0,...
 'L_mm',120,'h_mm',10,'out_of_plane_thickness_m',1,'sigma_xx_S_m',0,'sigma_zz_S_m',0);
failed=0;
for k=1:numel(cases)
 c=defaults;item=cases{k};names=fieldnames(item);
 for j=1:numel(names),c.(names{j})=item.(names{j});end
 c.z_interfaces_mm=reshape(c.z_interfaces_mm,1,[]);
 c.layer_materials=reshape(c.layer_materials,1,[]);
 if ~isfield(c,'interface_amplitudes_mm')
  c.interface_amplitudes_mm=c.A_mm*ones(size(c.z_interfaces_mm));
  c.interface_amplitudes_mm(c.z_interfaces_mm<=10 | c.z_interfaces_mm>=110)=0;
 else,c.interface_amplitudes_mm=reshape(c.interface_amplitudes_mm,1,[]);end
 assert(numel(c.layer_materials)==numel(c.z_interfaces_mm)+1);
 assert(all(diff([0 c.z_interfaces_mm 120])>0));
 folder=fullfile(out,c.id);
 if exist(fullfile(folder,'SUCCESS.json'),'file')
  prior=jsondecode(fileread(fullfile(folder,'parameters.json')));
  assert(isequaln(orderfields(prior),orderfields(jsondecode(jsonencode(c)))),'Saved case parameters differ; use a new ID.');
  fprintf('SKIP %s\n',c.id);continue;
 end
 if ~exist(folder,'dir'),mkdir(folder);end
 write_json(fullfile(folder,'parameters.json'),c);
 copyfile(fullfile(root,'scripts','solve_round2_case.m'),fullfile(folder,'solve_round2_case_snapshot.m'));
 copyfile(fullfile(root,'scripts','run_round2.m'),fullfile(folder,'run_round2_snapshot.m'));
 fprintf('\nCASE_START %s %s\n',c.id,char(datetime('now')));tt=tic;
 try
  result=solve_round2_case(c,folder);
  result.elapsed_seconds=toc(tt);result.status='success';
  write_json(fullfile(folder,'SUCCESS.json'),result);
  writetable(struct2table(result),fullfile(folder,'result.csv'));
  if exist(fullfile(folder,'FAILURE.txt'),'file'),movefile(fullfile(folder,'FAILURE.txt'),fullfile(folder,'RECOVERED_FAILURE.txt'),'f');end
  fprintf('CASE_SUCCESS %s sigma=%.12g time=%.1fs\n',c.id,result.sigma_eff_S_m,result.elapsed_seconds);
 catch ME
  failed=failed+1;msg=getReport(ME,'extended','hyperlinks','off');fprintf('%s\n',msg);
  fid=fopen(fullfile(folder,'FAILURE.txt'),'w');fprintf(fid,'%s\n',msg);fclose(fid);
 end
 ModelUtil.clear;
end
diary off;
assert(failed==0,'%d cases failed; inspect retained FAILURE.txt.',failed);
end
function write_json(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));fclose(fid);
end
