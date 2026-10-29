from pathlib import Path
import re,json,csv
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs'/'round7b_review_20261029'
state=OUT/'tolerance_checks'/'S1_M3_water_0p18_T_background_tol_2em11'
assert not (state/'comparison.json').exists()
rows=[]
for line in (state/'solver_progress.log').read_text(encoding='utf-8-sig',errors='replace').splitlines():
 m=re.fullmatch(r'\s{10,}(\d+)\s+([0-9.eE+\-]+)\s*',line)
 if m: rows.append({'linear_iteration':int(m[1]),'reported_LinErr':float(m[2])})
assert rows[-1]['linear_iteration']>=300
recent=rows[-100:]
with (state/'iteration_history.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
report={'status':'GMRES_2EM11_NOT_ACCEPTED_AS_CONVERGED_IN_OBSERVED_RUN',
 'last_iteration':rows[-1]['linear_iteration'],'last_100_min_LinErr':min(r['reported_LinErr'] for r in recent),
 'last_100_max_LinErr':max(r['reported_LinErr'] for r in recent),
 'action':'Active stop after repeated error plateaus. A single native FGMRES/right-preconditioned crosscheck at 2e-11 follows, with the same original saved state and AMG settings.',
 'source_log':str(state/'solver_progress.log'),'native_solver_nonconvergence_proven':False}
(state/'OBSERVED_PLATFORM_AND_ACTIVE_STOP.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
plan={'type':'fixed-mesh solver crosscheck, not a pure tolerance-only comparison',
 'sample':'S1','mesh':'M3','water_sigma_S_m':.18,
 'roles':['T_background','H_background','T_target','H_target'],
 'old_linear_solver':'GMRES/left preconditioning','old_stationary_tolerance':1e-10,
 'check_linear_solver':'COMSOL native FGMRES/right preconditioning','check_stationary_tolerance':2e-11,
 'unchanged':['source mesh','materials','electrodes','frozen kappa','AMG settings','restart=50','maxlinit=2000','rhob=400'],
 'initial_guess':'original accepted saved field; full DOF equality checked',
 'further_parameter_chasing_allowed':False,'old_1em11_and_2em11_GMRES_attempts_status':'uncompleted and retained',
 'official_local_API_reference':'COMSOL6.3 comsol_api_solver.51.31.html; linsolver=fgmres; prefuntype=right'}
(OUT/'FGMRES_crosscheck_plan_before_results.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False))
