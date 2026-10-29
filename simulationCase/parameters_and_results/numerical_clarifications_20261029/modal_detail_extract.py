from pathlib import Path
import csv,json
import argparse
parser = argparse.ArgumentParser(description="Read saved numerical outputs; do not launch a solver.")
parser.add_argument('--project-root', type=Path, required=True,
    help='Full original computation tree containing runs/. The compact public archive does not include every auxiliary matrix.')
parser.add_argument('--output-dir', type=Path, default=Path.cwd() / 'numerical_check_outputs')
args = parser.parse_args()
if not (args.project_root / 'runs').is_dir():
    parser.error('--project-root must contain the original runs/ directory')
args.output_dir.mkdir(parents=True, exist_ok=True)

R=args.project_root.resolve(); O=args.output_dir.resolve()
def rc(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def rj(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def wj(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def wc(n,x):
    with (O/n).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in x for k in r)));w.writeheader();w.writerows(x)
r5=R/'runs/round5_modal_closure';r6=R/'runs/round6_applicability';r7=R/'runs/round7_3d_ct'
net=[r for r in rc(r6/'side_and_rank_intervals.csv') if r['category']=='joint_net' and r['reference_K']=='8' and r['target_N']=='64']
out=[]
for r in net:
    d={k:r[k] for k in ['full_length_mm','role','sigma_f_S_m','contrast','reference_K','target_N']}
    for k in ['COMSOL_contrast_S','reference_contrast_S','contrast_lower_S','contrast_upper_S','observed_mesh_resolution_scale_S']:d[k.replace('_S','_uS')]=float(r[k])*1e6
    d['resolved_against_observed_mesh_scale']=r['resolved_against_observed_mesh_scale'];out.append(d)
wc('R6_eight_joint_contrasts.csv',out)
near=[r for r in rc(r5/'selected_mesh_effects.csv') if r['role']=='B' and float(r['sigma_f_S_m'])==.18 and r['comparison_kind']=='net_both_sides']
four=[r for r in rc(R/'runs/round4_extension_20261024/fine_factorial_comparison.csv') if r['quantity'] in ('gA','gB') and r['effect'] in ('average_upper','average_lower','joint')]
# Preserve all the small table if naming conventions differ.
if not four:four=[r for r in rc(R/'runs/round4_extension_20261024/fine_factorial_comparison.csv') if r['quantity'] in ('gA','gB')]
wc('R4_factorial_mesh_scales.csv',four)
wj('R5_B018_net_mesh_scale.json',near)
three=[r for r in rc(r7/'controlled_effects.csv') if r['role']=='A' and r['sigma_f_S_m']=='0.18' and r['effect_id']=='lower_near_minus_far_upper_far']
wc('R7_Fig09c_star_scale.csv',three)
times={'R5_freeze':{k:v for k,v in rj(r5/'method_frozen.json').items() if k in ('version','frozen_local','K_star','K_reference')},'R5_prediction_seal':{k:v for k,v in rj(r5/'predictions_ready.json').items() if k!='entries'},'R6_prediction_seal':{k:v for k,v in rj(r6/'predictions_ready.json').items() if k!='entries'}}
times['R5_prediction_count']=len(list((r5/'predictions').glob('*/sf_*/prediction.json')))
times['R6_prediction_count']=len(list((r6/'predictions').glob('*/sf_*/prediction.json')))
wj('prospective_record_evidence.json',times)

mesh=[]
for batch,run in [('R5',r5),('R6',r6)]:
    for p in sorted((run/'precision').glob('*/mesh_refinement.json')):
        d=rj(p);mesh.append({'batch':batch,'case':p.parent.name,'base_global_hmax_mm':float(d['old_global_hmax_expression']),'base_ellipse_hmax_mm':float(d['old_ellipse_hmax_expression']),'fine_global_hmax_mm':d['new_global_hmax_mm'],'fine_ellipse_hmax_mm':d['new_ellipse_hmax_mm'],'base_elements':d['source_triangles'],'fine_elements':d['fine_triangles'],'element':'triangle','potential_order':2,'source':p.as_posix()})
for r in rc(r7/'mesh_comparison.csv'):
    if r.get('evidence_family')!='synthetic_3D':continue
    mesh.append({'batch':'R7 synthetic','case':r['case_id']+'_sf'+r['sigma_f_S_m'],'base_global_hmax_mm':5,'base_ellipse_hmax_mm':.75,'fine_global_hmax_mm':3.5,'fine_ellipse_hmax_mm':.5,'base_elements':int(float(r['base_tetrahedra'])),'fine_elements':int(float(r['fine_tetrahedra'])),'element':'tetrahedron','potential_order':2,'source':(r7/'mesh_comparison.csv').as_posix()})
wc('nonCT_mesh_specifications.csv',mesh)
print(json.dumps({'R6_joint_count':len(out),'R5_near':near,'Fig09c_star':three,'mesh_rows':len(mesh),'times':times},ensure_ascii=False,indent=2))
