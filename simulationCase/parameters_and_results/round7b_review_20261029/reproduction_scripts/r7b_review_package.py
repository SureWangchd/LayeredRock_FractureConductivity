from pathlib import Path
import json,csv,zipfile,shutil
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs'/'round7b_review_20261029'
OLD=ROOT/'runs'/'round7b_localfit'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
provenance=OUT/'source_evidence';provenance.mkdir(exist_ok=True)
sources=[OLD/'mesh_convergence.csv',OLD/'ct_response_localfit.csv',
 OLD/'tolerance_checks'/'S5_M3_water_10_pair_four_states'/'comparison.json',
 OLD/'tolerance_checks'/'S5_M3_water_10_pair_four_states'/'four_state_tolerance_comparison.csv',
 OLD/'tolerance_checks'/'S1_M3_water_0p18_T_background_tol_1em11'/'OPTIONAL_TOLERANCE_TIMEOUT_20261027_183930.json']
index=[]
for sample,water in [('S1','0p18'),('S5','10')]:
 for mesh in ('M2','M3'):
  base=OLD/'states'/f'{sample}_{mesh}'/f'water_{water}'
  sources.extend([base/'PAIR_SUCCESS.json',base/'kappa_frozen.json'])
  if mesh=='M3':
   freeze=read(base/'kappa_frozen.json')
   for role,path in [('T_background',base/'T_background'/'SUCCESS.json'),('H_background',Path(freeze['H_background_success_json'])),('T_target',base/'T_target'/'SUCCESS.json'),('H_target',base/'H_target'/'SUCCESS.json')]:
    data=read(path);sources.append(path)
    index.append(dict(sample=sample,mesh=mesh,fluid_sigma_S_m=data['water_sigma_S_m'],role=role,status='original accepted state',model_path=data['model_path'],result_json=str(path),model_exists=Path(data['model_path']).exists()))
for role in ('T_background','H_background','T_target','H_target'):
 path=OLD/'tolerance_checks'/f'S5_M3_water_10_{role}_tol_1em10'/'comparison.json'
 data=read(path);sources.append(path)
 index.append(dict(sample='S5',mesh='M3',fluid_sigma_S_m=10,role=role,status='accepted old stricter tolerance state',model_path=data['strict_model_path'],result_json=str(path),model_exists=Path(data['strict_model_path']).exists()))
for source in sources:
 dest=provenance/source.relative_to(OLD);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
with (OUT/'actual_model_index.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=index[0]);w.writeheader();w.writerows(index)
scripts=OUT/'reproduction_scripts';scripts.mkdir(exist_ok=True)
for name in ['r7b_review_prepare.py','r7b_review_run.m','r7b_review_launch.ps1','r7b_review_watch.ps1','r7b_review_diagnostic_restart.ps1','r7b_review_audit.py','r7b_review_record_plateau.py','r7b_review_prepare_crosscheck.py','r7b_review_close_checks.py','r7b_review_attempts_summary.py','r7b_review_package.py']:
 shutil.copy2(ROOT/'scripts'/name,scripts/name)
bundle=ROOT/'outputs'/'CT_fixed_mesh_review_20261029.zip'
manifest=[]
for p in OUT.rglob('*'):
 if p.is_file() and p.suffix not in ('.mph','.pid') and p.name!='public_files_manifest.csv':
  manifest.append(dict(relative_path=p.relative_to(OUT).as_posix(),bytes=p.stat().st_size))
with (OUT/'public_files_manifest.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=['relative_path','bytes']);w.writeheader();w.writerows(sorted(manifest,key=lambda x:x['relative_path']))
with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as z:
 for p in OUT.rglob('*'):
  if p.is_file() and p.suffix not in ('.mph','.pid'):
   z.write(p,p.relative_to(OUT))
print(json.dumps({'bundle':str(bundle),'files':len(zipfile.ZipFile(bundle).namelist()),'bytes':bundle.stat().st_size,'models_indexed':len(index),'all_indexed_models_exist':all(x['model_exists'] for x in index),'new_S1_accepted_strict_models':0},ensure_ascii=False))
