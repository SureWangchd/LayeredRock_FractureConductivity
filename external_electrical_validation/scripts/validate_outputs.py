"""Check real data, calculated tables, conservation, geometry and delivery files."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]

def main():
    checks=[]
    def check(name,ok,detail=None):
        checks.append(dict(check=name,passed=bool(ok),detail=detail))
    raw=pd.read_csv(ROOT/'data/Han2022_table1_electrical_porosity.csv')
    data=pd.read_csv(ROOT/'EXPERIMENT_DATA.csv')
    keys=['sample','effective_pressure_MPa']
    a=raw.set_index(keys).sort_index();b=data.set_index(keys).sort_index()
    check('All 60 measured conductivity and porosity entries preserved',len(data)==60 and np.array_equal(a[['conductivity_S_per_m','porosity_fraction']].values,b[['conductivity_S_per_m','porosity_fraction']].values))
    roles=data.use_in_this_test.value_counts().to_dict()
    check('Two calibration, one intact holdout, three fractured targets',roles.get('background calibration')==2 and roles.get('intact direction holdout')==1 and roles.get('fractured test target')==3,roles)
    p=pd.read_csv(ROOT/'MODEL_PREDICTIONS.csv');c=pd.read_csv(ROOT/'COMPARISON.csv')
    check('Finite predictions and positive conductivity',len(p)==len(c)==39 and np.isfinite(p.predicted_S_per_m).all() and (p.predicted_S_per_m>0).all(),len(p))
    check('Error formula',np.allclose(c.signed_error_percent,100*(c.predicted_S_per_m/c.measured_S_per_m-1),rtol=0,atol=1e-10))
    check('Increment formula',np.allclose(c.predicted_fracture_increment_S_per_m,c.predicted_S_per_m-c.background_S_per_m,rtol=0,atol=1e-12))
    target=data[data.effective_pressure_MPa==5].set_index('sample').conductivity_S_per_m
    check('Every compared target is a source measurement',np.array_equal(c.measured_S_per_m.values,c['sample'].map(target).values))
    check('No negative gain from adding high-conductivity crack at fixed background',(c.predicted_fracture_increment_S_per_m>=-1e-10).all())
    diagnostics=[]
    for path in (ROOT/'runs/cell_results').glob('*.json'):
        item=json.loads(path.read_text())
        for direction in ['normal_solve','parallel_solve']:
            if direction in item:diagnostics.append(item[direction])
    specimens=[]
    for path in (ROOT/'runs/specimen_results').glob('*.json'):
        item=json.loads(path.read_text());diagnostics.append(item);specimens.append(item)
    maxres=max(d['relative_residual'] for d in diagnostics)
    maxbal=max(d['relative_electrode_imbalance'] for d in diagnostics)
    maxpow=max(d['relative_power_current_mismatch'] for d in diagnostics)
    check('All actual PDE residuals below 1.01e-10',maxres<1.01e-10,maxres)
    check('Electrode current conservation',maxbal<1e-6,maxbal)
    check('Power equals current times imposed voltage',maxpow<1e-7,maxpow)
    check('Finite cylinder responses lie within energy bounds',all(d['reuss_direction_bound']*(1-1e-6)<=d['conductivity_S_per_m']<=d['voigt_direction_bound']*(1+1e-6) for d in specimens))
    conv=pd.read_csv(ROOT/'MESH_CONVERGENCE.csv')
    check('Four actual refinement levels for each fractured target',len(conv)==12 and all(g.level.tolist()==[0,1,2,3] for _,g in conv.groupby('sample')))
    req=['REPORT_CN.md','DATA_SOURCES.md','EXPERIMENT_DATA.csv','MODEL_PREDICTIONS.csv','COMPARISON.csv',
         'VALIDATION_FIGURE.pdf','VALIDATION_FIGURE.svg','VALIDATION_FIGURE.png','VALIDATION_ASSESSMENT.md',
         'README.md','FIGURE_CAPTION.md','UNCERTAINTY_BUDGET.csv',
         'requirements.txt','ENVIRONMENT.json','INPUT_PROTOCOL.json']
    check('Required delivery files present',all((ROOT/n).is_file() and (ROOT/n).stat().st_size>0 for n in req))
    check('PDF artifact is a PDF',(ROOT/'VALIDATION_FIGURE.pdf').read_bytes().startswith(b'%PDF'))
    freeze=json.loads((ROOT/'runs/input_freeze.json').read_text())
    check('Frozen input protocol unchanged',freeze['protocol_sha256']==hashlib.sha256((ROOT/'INPUT_PROTOCOL.json').read_bytes()).hexdigest())
    out=dict(passed=all(d['passed'] for d in checks),checks=checks,actual_recorded_PDE_solutions=len(diagnostics))
    (ROOT/'checks').mkdir(exist_ok=True)
    (ROOT/'checks/DELIVERY_VALIDATION.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(out,indent=2,ensure_ascii=False))
    if not out['passed']:raise SystemExit(1)

if __name__=='__main__':main()
