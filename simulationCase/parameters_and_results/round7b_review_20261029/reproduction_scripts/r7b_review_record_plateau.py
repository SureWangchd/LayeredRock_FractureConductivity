from pathlib import Path
import re,json,csv
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs'/'round7b_review_20261029'
state=OUT/'tolerance_checks'/'S1_M3_water_0p18_T_background_tol_1em11'
text=(state/'solver_progress.log').read_text(encoding='utf-8-sig',errors='replace')
rows=[]
for line in text.splitlines():
    m=re.fullmatch(r'\s{10,}(\d+)\s+([0-9.eE+\-]+)\s*',line)
    if m: rows.append({'linear_iteration':int(m[1]),'reported_LinErr':float(m[2])})
assert len(rows)>=300 and rows[-1]['linear_iteration']>=300
recent=rows[-100:]
assert min(r['reported_LinErr'] for r in recent)>=1e-11
with (state/'iteration_history.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
record={'status':'ONE_EM11_NOT_ACCEPTED_AS_CONVERGED_IN_OBSERVED_RUN','last_iteration':rows[-1]['linear_iteration'],
 'last_100_iteration_min_LinErr':min(r['reported_LinErr'] for r in recent),
 'last_100_iteration_max_LinErr':max(r['reported_LinErr'] for r in recent),
 'requested_stationary_tolerance':1e-11,
 'action':'Stop this diagnostic run; perform a separate fixed-mesh four-state check at 2e-11, five times tighter than the saved 1e-10 tolerance. Keep the 1e-11 attempt uncompleted.',
 'source_progress_log':str(state/'solver_progress.log'),'native_solver_nonconvergence_proven':False,
 'interpretation':'The observed reported error repeatedly approaches 1e-11 and rebounds. The rounded log includes a value printed as 1e-11 without solver termination. This trace does not identify the exact cause or establish a rigorous arithmetic error floor.'}
(state/'OBSERVED_PLATFORM_AND_ACTIVE_STOP.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'two_em11_check_plan_before_results.json').write_text(json.dumps({'sample':'S1','mesh':'M3','water_sigma_S_m':.18,'roles':['T_background','H_background','T_target','H_target'],'source_stationary_tolerance':1e-10,'new_stationary_tolerance':2e-11,'tightening_factor':5,'original_1em11_attempt_status':'uncompleted','same_mesh_materials_kappa':True,'background_recalibrated':False,'main_batch_files_changed':False},indent=2),encoding='utf-8')
print(json.dumps(record,ensure_ascii=False))
