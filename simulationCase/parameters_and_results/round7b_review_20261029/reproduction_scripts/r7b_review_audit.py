"""Read-only reconstruction of CT tolerance evidence and numerical scales."""
from pathlib import Path
import json, csv

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / 'runs' / 'round7b_localfit'
OUT = ROOT / 'runs' / 'round7b_review_20261029'
OUT.mkdir(exist_ok=True)
load = lambda p: json.loads(p.read_text(encoding='utf-8-sig'))
def pair(s,m,w):
    return load(OLD/'states'/f'{s}_{m}'/f'water_{w}'/'PAIR_SUCCESS.json')
a,b=pair('S5','M2','10'),pair('S5','M3','10')
tol_path=OLD/'tolerance_checks'/'S5_M3_water_10_pair_four_states'/'comparison.json'
t=load(tol_path)
terms={
 'abs_gT_M3_minus_M2_S':abs(b['gT_S']-a['gT_S']),
 'abs_gH_M3_minus_M2_S':abs(b['gH_S']-a['gH_S']),
 'abs_background_mismatch_M2_S':abs(a['G0_H_S']-a['G0_T_S']),
 'abs_background_mismatch_M3_S':abs(b['G0_H_S']-b['G0_T_S']),
 'sum_abs_four_state_tolerance_changes_S':sum(abs(x) for x in t['conductance_changes_by_role_S']),
}
scale=3*sum(terms.values())
record={'S5_scale_terms':terms,'empirical_multiplier':3,'combined_observed_scale_S':scale,
 'combined_observed_scale_microS':scale*1e6,'delta_M3_S':b['delta_gH_minus_gT_S'],
 'delta_over_scale':abs(b['delta_gH_minus_gT_S'])/scale,
 'formula':'3*(abs(gT3-gT2)+abs(gH3-gH2)+abs(G0H2-G0T2)+abs(G0H3-G0T3)+sum_four_states(abs(Gstrict-Gsaved)))',
 'scale_type':'empirical observed numerical change scale, not a guaranteed error bound or confidence interval',
 'source_pair_M2':str(OLD/'states'/'S5_M2'/'water_10'/'PAIR_SUCCESS.json'),
 'source_pair_M3':str(OLD/'states'/'S5_M3'/'water_10'/'PAIR_SUCCESS.json'),
 'source_tolerance':str(tol_path)}
old_timeout=OLD/'tolerance_checks'/'S1_M3_water_0p18_T_background_tol_1em11'/'OPTIONAL_TOLERANCE_TIMEOUT_20261027_183930.json'
record['old_S1_status']=load(old_timeout)
new_path=OUT/'tolerance_checks'/'S1_M3_water_0p18_pair_four_states_fgmres_tol_2em11'/'comparison.json'
record['new_S1_complete']=new_path.exists()
if new_path.exists():
    n=load(new_path)
    record['new_S1_result']=n
    p2,p3=pair('S1','M2','0p18'),pair('S1','M3','0p18')
    partial=abs(p3['gT_S']-p2['gT_S'])+abs(p3['gH_S']-p2['gH_S'])+abs(p2['G0_H_S']-p2['G0_T_S'])+abs(p3['G0_H_S']-p3['G0_T_S'])
    record['new_S1_combined_scale_S']=3*(partial+n['observed_gain_or_pair_change_scale_S'])
    record['new_S1_delta_over_scale']=abs(p3['delta_gH_minus_gT_S'])/record['new_S1_combined_scale_S']
rows=[]
for s in ('S1','S5'):
    p18,p10=pair(s,'M1','0p18'),pair(s,'M1','10')
    for key in ('G0_T_S','G0_H_S','Gf_T_S','Gf_H_S','gT_S','gH_S'):
        rows.append({'sample':s,'quantity':key,'water_0p18_value_S':p18[key],'water_10_value_S':p10[key],
                     'deviation_from_uniform_all_phase_scaling_percent':100*(p10[key]/(p18[key]*10/.18)-1)})
record['flow_band_interpretation']='Only fluid conductivity changes between 0.18 and 10 S/m; skeleton stays at 1e-5 S/m and H kappa is recalibrated on its background for each fluid. These two states do not satisfy all-phase uniform conductivity scaling.'
with (OUT/'S1_S5_fluid_band_comparison.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
(OUT/'CT_tolerance_source_audit.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
state_rows=[];summary_rows=[];scale_rows=[]
available=[('S5',10.0,t,tol_path)]
if new_path.exists(): available.insert(0,('S1',.18,record['new_S1_result'],new_path))
for sample,water,tr,trpath in available:
    wt='0p18' if water==.18 else '10'
    check_type='fixed-mesh solver crosscheck; GMRES/left 1e-10 to FGMRES/right 2e-11' if sample=='S1' else 'fixed-mesh tolerance check; GMRES/left 1e-8 to GMRES/left 1e-10'
    m2,m3=pair(sample,'M2',wt),pair(sample,'M3',wt)
    for k,role in enumerate(tr['state_roles']):
        state_rows.append(dict(sample=sample,fluid_sigma_S_m=water,mesh='M3',state_role=role,check_type=check_type,
            saved_tolerance=tr['old_stationary_tolerances'][k],checked_tolerance=tr['strict_stationary_tolerances'][k],
            saved_G_S=tr['old_G_by_role_S'][k],checked_G_S=tr['strict_G_by_role_S'][k],
            checked_minus_saved_G_S=tr['conductance_changes_by_role_S'][k],
            source_model=tr['old_model_paths'][k],checked_model=tr['strict_model_paths'][k],
            source_pair_comparison=str(trpath)))
    vals={'T_gain_mesh_change':abs(m3['gT_S']-m2['gT_S']),
          'H_gain_mesh_change':abs(m3['gH_S']-m2['gH_S']),
          'background_mismatch_M2':abs(m2['G0_H_S']-m2['G0_T_S']),
          'background_mismatch_M3':abs(m3['G0_H_S']-m3['G0_T_S']),
          'four_state_solver_check_changes':sum(abs(x) for x in tr['conductance_changes_by_role_S'])}
    combined=3*sum(vals.values())
    for component,value in vals.items():
        scale_rows.append(dict(sample=sample,fluid_sigma_S_m=water,component=component,value_S=value,value_microS=value*1e6,multiplier=3,contribution_to_combined_scale_S=3*value,check_type=check_type))
    summary_rows.append(dict(sample=sample,fluid_sigma_S_m=water,mesh='M3',
        saved_tolerance=tr['old_stationary_tolerances'][0],checked_tolerance=tr['strict_stationary_tolerances'][0],
        kappa_fixed_S_m=tr['kappa_eq_S_m'],background_recalibrated=False,check_type=check_type,
        saved_gT_S=tr['old_gT_S'],checked_gT_S=tr['strict_gT_S'],
        saved_gH_S=tr['old_gH_S'],checked_gH_S=tr['strict_gH_S'],
        saved_pair_difference_S=tr['old_delta_gH_minus_gT_S'],checked_pair_difference_S=tr['strict_delta_gH_minus_gT_S'],
        saved_epsilon_percent=100*tr['old_epsilon_delta_over_gT'],checked_epsilon_percent=100*tr['strict_epsilon_delta_over_gT'],
        four_state_solver_change_scale_S=vals['four_state_solver_check_changes'],
        combined_empirical_scale_S=combined,pair_difference_to_scale_ratio=abs(m3['delta_gH_minus_gT_S'])/combined,
        checked_relative_background_mismatch=tr['strict_relative_background_mismatch'],
        original_background_match_tolerance=tr['required_original_background_match_tolerance'],
        checked_background_match_within_original_tolerance=tr['strict_background_still_within_original_matching_tolerance'],
        source_pair_comparison=str(trpath)))
for name,rows in [('CT_four_state_tolerance_comparison.csv',state_rows),('CT_pair_tolerance_summary.csv',summary_rows),('CT_combined_scale_components.csv',scale_rows)]:
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print(json.dumps({'S5_scale_microS':record['combined_observed_scale_microS'],'S5_delta_over_scale':record['delta_over_scale'],'S1_complete':record['new_S1_complete']},ensure_ascii=False))
