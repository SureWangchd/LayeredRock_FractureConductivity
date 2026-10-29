from pathlib import Path
import json,csv
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'runs'/'round7b_review_20261029'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
records=[]
def add(scheme,tolerance,initial,kind,path,iterations='',lo='',hi='',count=0):
 records.append(dict(sample='S1',mesh='M3',fluid_sigma_S_m=.18,scheme=scheme,stationary_tolerance=tolerance,initial_field=initial,
 classification=kind,observed_linear_iterations=iterations,last_100_min_printed_LinErr=lo,last_100_max_printed_LinErr=hi,
 accepted_new_states=count,required_states=4,remaining_states_started=False,evidence=str(path)))
old=ROOT/'runs'/'round7b_localfit'/'tolerance_checks'/'S1_M3_water_0p18_T_background_tol_1em11'/'OPTIONAL_TOLERANCE_TIMEOUT_20261027_183930.json'
add('original GMRES/left',1e-11,'original initial field','2026-10-27 active stop at 1200-second budget; native nonconvergence not established',old)
add('review GMRES/left',1e-11,'original initial field','active stop after 25-minute diagnostic checkpoint; no per-iteration residual trace',OUT/'attempts'/'zero_start_20261029_221636'/'ACTIVE_TERMINATION.json')
add('saved-field initialization adapter',1e-11,'saved field copy','API property mismatch before solve; useinitsol unavailable on actual Variables node',OUT/'attempts'/'warm_initial_API_20261029_224242'/'ADAPTER_FAILURE.json')
for suffix,scheme in [('T_background_tol_1em11','GMRES/left'),('T_background_tol_2em11','GMRES/left'),('T_background_fgmres_tol_2em11','native FGMRES/right')]:
 path=OUT/'tolerance_checks'/f'S1_M3_water_0p18_{suffix}'/'OBSERVED_PLATFORM_AND_ACTIVE_STOP.json'
 r=read(path)
 low=r.get('last_100_iteration_min_LinErr',r.get('last_100_min_LinErr',r.get('last_100_min_reported_LinErr')))
 high=r.get('last_100_iteration_max_LinErr',r.get('last_100_max_LinErr',r.get('last_100_max_reported_LinErr')))
 add(scheme,1e-11 if suffix.endswith('1em11') else 2e-11,'verified original saved field',
  'active stop after observed error plateau; no accepted stricter solution',path,r['last_iteration'],low,high)
with (OUT/'CT_attempts_summary.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
print(f'{len(records)} classified records')
