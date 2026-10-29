function r7b_review_run(root)
% Recheck the originally incomplete S1/M3/0.18 pair, without a wall-time kill.
out=fullfile(root,'runs','round7b_review_20261029');
addpath(fullfile(root,'scripts'));addpath(fullfile(out,'scripts'));
maxNumCompThreads(26);
base='D:/Program Files/COMSOL/COMSOL63/Multiphysics';addpath(fullfile(base,'mli'));
mphstart('127.0.0.1',2038,base);import com.comsol.model.util.*
cleanup=onCleanup(@disconnect); %#ok<NASGU>
stamp=char(datetime('now','Format','yyyyMMdd_HHmmss'));
diary(fullfile(out,['run_' stamp '.log']));
try
 result=r7b_review_pair_tolerance_check(root,'S1','M3',.18);
 r3_json(fullfile(out,'S1_M3_pair_completed.json'),result);
catch ME
 r3_json(fullfile(out,['FAILURE_' stamp '.json']),struct('status','FAILED','message',ME.message,'report',getReport(ME,'extended','hyperlinks','off')));
 rethrow(ME);
end
diary off;
end
function disconnect(),try,com.comsol.model.util.ModelUtil.disconnect;catch,end;end
