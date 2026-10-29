"""Close the bounded native-solver diagnostic without inventing a result."""
from pathlib import Path
import csv,json,re,datetime
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs'/'round7b_review_20261029'
state=OUT/'tolerance_checks'/'S1_M3_water_0p18_T_background_fgmres_tol_2em11'
assert not (state/'comparison.json').exists(), 'A successful state exists; do not stop it.'
rows=[]
for line in (state/'solver_progress.log').read_text(encoding='utf-8-sig',errors='replace').splitlines():
 m=re.fullmatch(r'\s{10,}(\d+)\s+([0-9.eE+\-]+)\s*',line)
 if m: rows.append({'linear_iteration':int(m[1]),'reported_LinErr':float(m[2])})
assert rows[-1]['linear_iteration']>=300
recent=rows[-100:]
with (state/'iteration_history.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
record={'status':'STOPPED_AFTER_OBSERVED_PLATFORM_NO_ACCEPTED_STRICT_SOLUTION',
 'time_local':datetime.datetime.now().isoformat(timespec='seconds'),
 'sample':'S1','mesh':'M3','fluid_sigma_S_m':.18,'state':'T_background',
 'check_type':'fixed-mesh solver crosscheck; original GMRES/left 1e-10 versus FGMRES/right 2e-11',
 'last_iteration':rows[-1]['linear_iteration'],
 'last_100_min_reported_LinErr':min(r['reported_LinErr'] for r in recent),
 'last_100_max_reported_LinErr':max(r['reported_LinErr'] for r in recent),
 'requested_stationary_tolerance':2e-11,
 'native_solver_accepted_strict_solution':False,
 'strict_G_S':None,'four_state_solver_change_scale_S':None,'full_S1_combined_scale_S':None,
 'remaining_three_states_started':False,
 'termination_type':'active stop after the pre-agreed bounded diagnostic; not a native convergence failure return',
 'reason':'The logged error repeatedly approaches the requested tolerance and rebounds; no accepted new solution was returned. No further tolerance or solver changes are made.',
 'exact_cause_of_platform_established':False,
 'saved_main_batch_unchanged':True,'mesh_materials_electrodes_kappa_unchanged':True,
 'sensitivity_scale_is_zero':False,'source_log':str(state/'solver_progress.log')}
(state/'OBSERVED_PLATFORM_AND_ACTIVE_STOP.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'S1_bounded_fixed_mesh_check_final_status.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
statuses=[dict(sample='S1',mesh='M3',fluid_sigma_S_m=.18,status='four-state check incomplete',accepted_states=0,
 original_solver='GMRES/left',original_tolerance=1e-10,attempted_solver='FGMRES/right',attempted_tolerance=2e-11,
 four_state_solver_change_scale_S='',combined_observed_scale_S='',source=str(OUT/'S1_bounded_fixed_mesh_check_final_status.json')),
 dict(sample='S5',mesh='M3',fluid_sigma_S_m=10.,status='four-state tolerance check complete',accepted_states=4,
 original_solver='GMRES/left',original_tolerance=1e-8,attempted_solver='GMRES/left',attempted_tolerance=1e-10,
 four_state_solver_change_scale_S=1.830460425161079e-12,combined_observed_scale_S=4.5272845192786265e-9,
 source=str(ROOT/'runs'/'round7b_localfit'/'tolerance_checks'/'S5_M3_water_10_pair_four_states'/'comparison.json'))]
with (OUT/'CT_fixed_mesh_check_status.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=statuses[0]);w.writeheader();w.writerows(statuses)
print(json.dumps(record,ensure_ascii=False))
