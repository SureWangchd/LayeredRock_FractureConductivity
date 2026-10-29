function r4e_modes_run()
% Read saved targets only. Connects to existing server 2037; never starts it.
root=fileparts(fileparts(mfilename('fullpath')));out=fullfile(root,'runs','round4_extension_20261024');
base='D:/Program Files/COMSOL/COMSOL63/Multiphysics';addpath(fullfile(base,'mli'));maxNumCompThreads(4);
old=jsondecode(fileread(fullfile(root,'runs','round4_crossover','periphery_manifest.json')));
new=jsondecode(fileread(fullfile(out,'manifest.json')));
oldcases=old.cases;if ~iscell(oldcases),oldcases=num2cell(oldcases);end
newcases=new.cases;if ~iscell(newcases),newcases=num2cell(newcases);end
jobs={};sf=.159375;
for role={'A','B'}
 for layout={'NN','NF','FN','FF'}
  if ismember(layout{1},{'NN','FF'})
   oldLayout='swap_near';if strcmp(layout{1},'FF'),oldLayout='swap_far';end
   ix=find(cellfun(@(c)strcmp(c.configuration_role,role{1})&&strcmp(c.layout,oldLayout),oldcases));
   assert(numel(ix)==1);c=oldcases{ix};folder=fullfile(root,'runs','round4_crossover','periphery',c.id);
  else
   ix=find(cellfun(@(c)strcmp(c.configuration_role,role{1})&&strcmp(c.layout,layout{1}),newcases));
   assert(numel(ix)==1);c=newcases{ix};folder=fullfile(out,'cases',c.id);
  end
  c.sigma_frac_S_m=sf;c.passive_ellipse=false;
  modelPath=fullfile(folder,'sf_0p159375','solution.mph');selectionPath=fullfile(folder,'source_background','selections.json');
  assert(isfile(modelPath)&&isfile(selectionPath),'Finish all target solves before modal postprocessing');
  jobs{end+1}=struct('c',c,'model',modelPath,'selections',selectionPath,'dest',fullfile(out,'modes',[role{1} '_' layout{1}])); %#ok<AGROW>
 end
end
mphstart('127.0.0.1',2037,base);import com.comsol.model.util.*
cleanup=onCleanup(@() ModelUtil.disconnect); %#ok<NASGU>
diary(fullfile(out,['run_modes_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '.log']));
for j=1:numel(jobs)
 job=jobs{j};if isfile(fullfile(job.dest,'mode_summary.json')),continue;end
 m=mphload(job.model,'R4EModes');sel=jsondecode(fileread(job.selections));
 r4e_mode_trace(m,job.c,sel,job.dest,32);ModelUtil.clear;clear m;
 fprintf('R4E_MODES_DONE %s\n',job.c.id);
end
diary off;
end
