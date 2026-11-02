"""Execute a fixed-input, retrospective comparison with Han2022 laboratory data.

No fractured conductivity is used below except in the final comparison table.
Run: python scripts/run_validation.py --stage nominal|sensitivity|all
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(key,'2')
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
import pandas as pd
from mesh_models import cell_mesh, cylinder_mesh
from tet_fem import solve, element_gradients

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs'
RUN.mkdir(exist_ok=True)
DATA=pd.read_csv(ROOT/'data/Han2022_table1_electrical_porosity.csv')
D5=DATA[DATA.effective_pressure_MPa==5].set_index('sample')

def dump(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')

def calibration(sw):
    return {name:float(np.log(D5.loc[name,'conductivity_S_per_m']/sw)/np.log(D5.loc[name,'porosity_fraction'])) for name in ['I-0','I-90']}

def solve_axis(mesh, tensor, axis, is_cell=True):
    xyz=mesh['nodes']*1e-3
    lo,hi=xyz[:,axis].min(),xyz[:,axis].max()
    i0=np.flatnonzero(np.isclose(xyz[:,axis],lo,atol=1e-12,rtol=0))
    i1=np.flatnonzero(np.isclose(xyz[:,axis],hi,atol=1e-12,rtol=0))
    sol=solve(xyz,mesh['cells'],tensor,np.r_[i0,i1],np.r_[np.zeros(len(i0)),np.ones(len(i1))],method='cg',rtol=1e-10)
    current=float(sol['dirichlet_reaction'][sol['fixed_values']==1].sum())
    if is_cell:
        area=np.prod(np.delete(mesh['lengths']*1e-3,axis))
    else:
        # Polygonal mesh area from volume/length: numerical geometry calibration
        # eliminates a trivial circular-boundary discretization factor.
        area=sol['mesh_info']['volume']/(hi-lo)
    apparent=current*(hi-lo)/area
    out=dict(conductivity_S_per_m=apparent,current_A=current,power_W=sol['energy'],
             relative_residual=sol['relative_residual'],
             relative_power_current_mismatch=abs(sol['energy']/current-1),
             relative_electrode_imbalance=sol['current_balance_absolute']/abs(current),
             nodes=sol['mesh_info']['nodes'],tets=sol['mesh_info']['elements'],
             solver=sol['solver_info'],volume_m3=sol['mesh_info']['volume'])
    return out

def effective_tensor(sample, sw=4.69, aperture=1., spacing=3., level=1, closure='bulk'):
    key=f'{sample}_w{sw:g}_a{aperture:g}_s{spacing:g}_l{level}_{closure}'
    path=RUN/'cell_results'/f'{key}.json'
    if path.exists(): return json.loads(path.read_text())
    phi_c=4*np.pi/3*.062*.021*aperture
    phi_t=float(D5.loc[sample,'porosity_fraction'])
    phi_b=phi_t-phi_c
    if closure=='local': phi_b/=1-phi_c
    m=calibration(sw)
    kb0=sw*phi_b**m['I-0']; kb90=sw*phi_b**m['I-90']
    result=dict(sample=sample,brine_S_per_m=sw,aperture_fraction=aperture,spacing_mm=spacing,level=level,
                porosity_closure=closure,phi_total=phi_t,phi_crack_analytical=phi_c,phi_background=phi_b,
                background_normal_S_per_m=kb0,background_parallel_S_per_m=kb90,calibration_exponents=m)
    if aperture==0:
        result.update(normal_S_per_m=kb0,parallel_S_per_m=kb90,numerical_crack_volume_fraction=0.)
    else:
        mesh=cell_mesh(RUN/'meshes'/f'cell_a{.021*aperture:g}_s{spacing:g}_l{level}.npz',alpha=.021*aperture,spacing=spacing,level=level)
        bg=np.diag([kb90,kb90,kb0])
        sig=np.broadcast_to(bg,(len(mesh['cells']),3,3)).copy()
        sig[mesh['labels']==1]=np.eye(3)*sw
        normal=solve_axis(mesh,sig,2); parallel=solve_axis(mesh,sig,0)
        _,vol,_=element_gradients(mesh['nodes'],mesh['cells'])
        result.update(normal_S_per_m=normal['conductivity_S_per_m'],parallel_S_per_m=parallel['conductivity_S_per_m'],
                      numerical_crack_volume_fraction=float(vol[mesh['labels']==1].sum()/vol.sum()),
                      normal_solve=normal,parallel_solve=parallel)
    path.parent.mkdir(exist_ok=True)
    dump(path,result)
    print(json.dumps(dict(stage='cell',key=key,normal=result['normal_S_per_m'],parallel=result['parallel_S_per_m'])),flush=True)
    return result

def specimen(normal,parallel,angle,length=40.,level=1):
    key=hashlib.sha256(json.dumps([normal,parallel,angle,length,level]).encode()).hexdigest()[:20]
    path=RUN/'specimen_results'/f'{key}.json'
    if path.exists():return json.loads(path.read_text())
    mesh=cylinder_mesh(RUN/'meshes'/f'cylinder_L{length:g}_l{level}.npz',length=length,level=level)
    th=np.deg2rad(angle); n=np.array([np.sin(th),0,np.cos(th)])
    tensor=np.eye(3)*parallel+(normal-parallel)*np.outer(n,n)
    out=solve_axis(mesh,tensor,2,False)
    out.update(angle_deg=angle,length_mm=length,mesh_level=level,tensor=tensor.tolist(),
               reuss_direction_bound=1/np.linalg.inv(tensor)[2,2],voigt_direction_bound=tensor[2,2],
               circular_area_factor=out['volume_m3']/(length*1e-3)/(np.pi*(.0254/2)**2))
    path.parent.mkdir(exist_ok=True);dump(path,out)
    return out

def job(spec):
    sample=spec.get('sample','F-45');level=spec.get('level',1)
    eff=effective_tensor(sample,sw=spec.get('sw',4.69),aperture=spec.get('aperture',1.),spacing=spec.get('spacing',3.),level=level,closure=spec.get('closure','bulk'))
    angle=float(D5.loc[sample,'angle_to_bedding_normal_deg'])
    length=spec.get('length',40.)
    out=specimen(eff['normal_S_per_m'],eff['parallel_S_per_m'],angle,length,level)
    bg=specimen(eff['background_normal_S_per_m'],eff['background_parallel_S_per_m'],angle,length,level)
    row=dict(scenario=spec['scenario'],sample=sample,pressure_MPa=5,angle_deg=angle,
             sw_S_per_m=eff['brine_S_per_m'],aperture_fraction=eff['aperture_fraction'],spacing_mm=eff['spacing_mm'],
             porosity_closure=eff['porosity_closure'],length_mm=length,mesh_level=level,
             predicted_S_per_m=out['conductivity_S_per_m'],background_S_per_m=bg['conductivity_S_per_m'],
             predicted_fracture_increment_S_per_m=out['conductivity_S_per_m']-bg['conductivity_S_per_m'],
             sigma_normal=eff['normal_S_per_m'],sigma_parallel=eff['parallel_S_per_m'],
             phi_crack_analytic=eff['phi_crack_analytical'],phi_crack_mesh=eff['numerical_crack_volume_fraction'],
             specimen_nodes=out['nodes'],specimen_tets=out['tets'],
             relative_residual=out['relative_residual'],power_current_mismatch=out['relative_power_current_mismatch'],
             specimen_area_factor=out['circular_area_factor'])
    return row

def run(stage):
    # Save the exact model configuration fingerprint before predictions.
    proto=ROOT/'INPUT_PROTOCOL.json'
    stamp=RUN/'input_freeze.json'
    if not stamp.exists(): dump(stamp,dict(protocol_sha256=hashlib.sha256(proto.read_bytes()).hexdigest(),time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),note='No independent historical blinding claimed.'))
    experimental=DATA.copy()
    experimental['use_in_this_test']='archived; outside fixed 5 MPa test'
    at5=experimental.effective_pressure_MPa==5
    experimental.loc[at5 & experimental['sample'].isin(['I-0','I-90']),'use_in_this_test']='background calibration'
    experimental.loc[at5 & (experimental['sample']=='I-45'),'use_in_this_test']='intact direction holdout'
    experimental.loc[at5 & experimental['sample'].str.startswith('F-'),'use_in_this_test']='fractured test target'
    experimental.to_csv(ROOT/'EXPERIMENT_DATA.csv',index=False)
    specs=[]
    if stage in ('nominal','all'):
        for lev in (0,1,2,3):
            for s in ['F-0','F-45','F-90']:specs.append(dict(scenario='nominal',sample=s,level=lev))
    if stage in ('sensitivity','all'):
        for s in ['F-0','F-45','F-90']:
            for ap in (0.,.5):specs.append(dict(scenario=f'aperture_{ap:g}',sample=s,aperture=ap))
            for sp in (2.5,3.5):specs.append(dict(scenario=f'spacing_{sp:g}',sample=s,spacing=sp))
            specs.append(dict(scenario='brine_5.2',sample=s,sw=5.2))
            specs.append(dict(scenario='local_porosity',sample=s,closure='local'))
            specs.append(dict(scenario='local_porosity',sample=s,closure='local',level=3))
            for leng in (30.,50.):specs.append(dict(scenario=f'length_{leng:g}',sample=s,length=leng))
    rows=[]
    # Sequential execution preserves cache consistency and limits peak memory.
    for spec in specs:
        rows.append(job(spec))
        pd.DataFrame(rows).to_csv(RUN/f'predictions_{stage}_progress.csv',index=False)
    # Genuine unused intact direction: parameters depend only on I0/I90.
    if stage in ('nominal','all'):
        phi=float(D5.loc['I-45','porosity_fraction']);m=calibration(4.69)
        kn=4.69*phi**m['I-0'];kp=4.69*phi**m['I-90']
        hold=[]
        for L in (30.,40.,50.):
            for lev in (0,1,2,3):
                out=specimen(kn,kp,45,L,lev)
                hold.append(dict(sample='I-45',length_mm=L,level=lev,predicted_S_per_m=out['conductivity_S_per_m'],observed_S_per_m=float(D5.loc['I-45','conductivity_S_per_m']),**{k:out[k] for k in ['reuss_direction_bound','voigt_direction_bound','relative_residual']}))
        pd.DataFrame(hold).to_csv(RUN/'intact_holdout.csv',index=False)
        neg=specimen(kn,kp,-45,40,2);pos=specimen(kn,kp,45,40,2)
        dump(RUN/'operator_symmetry_check.json',dict(positive45=pos['conductivity_S_per_m'],negative45=neg['conductivity_S_per_m'],relative_difference=abs(pos['conductivity_S_per_m']/neg['conductivity_S_per_m']-1)))
    pd.DataFrame(rows).to_csv(RUN/f'predictions_{stage}.csv',index=False)
    files=[RUN/f'predictions_{s}.csv' for s in ('nominal','sensitivity') if (RUN/f'predictions_{s}.csv').exists()] if stage!='all' else [RUN/'predictions_all.csv']
    preds=pd.concat([pd.read_csv(p) for p in files],ignore_index=True)
    preds['uncalibrated_polygon_conductivity_S_per_m']=preds.predicted_S_per_m*preds.specimen_area_factor
    preds.to_csv(ROOT/'MODEL_PREDICTIONS.csv',index=False)
    comp=preds.copy();comp['measured_S_per_m']=comp['sample'].map(D5.conductivity_S_per_m)
    comp['signed_error_percent']=100*(comp.predicted_S_per_m/comp.measured_S_per_m-1)
    comp['absolute_error_percent']=abs(comp.signed_error_percent)
    comp['measured_minus_transferred_background_S_per_m']=comp.measured_S_per_m-comp.background_S_per_m
    comp['measurement_tolerance_S_per_m']=.02*comp.measured_S_per_m
    comp.to_csv(ROOT/'COMPARISON.csv',index=False)
    print(comp[['scenario','sample','mesh_level','predicted_S_per_m','measured_S_per_m','signed_error_percent']].to_string(index=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['nominal','sensitivity','all'],default='nominal')
    run(p.parse_args().stage)
