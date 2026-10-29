from pathlib import Path
import csv,json,sys
import numpy as np
import h5py
from scipy.io import loadmat

import argparse
parser = argparse.ArgumentParser(description="Read saved numerical outputs; do not launch a solver.")
parser.add_argument('--project-root', type=Path, required=True,
    help='Full original computation tree containing runs/. The compact public archive does not include every auxiliary matrix.')
parser.add_argument('--output-dir', type=Path, default=Path.cwd() / 'numerical_check_outputs')
args = parser.parse_args()
if not (args.project_root / 'runs').is_dir():
    parser.error('--project-root must contain the original runs/ directory')
args.output_dir.mkdir(parents=True, exist_ok=True)

ROOT=args.project_root.resolve()
OUT=args.output_dir.resolve()
def readj(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def readc(p): return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def writej(name,d): (OUT/name).write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def writec(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def mesh(p):
    if h5py.is_hdf5(p):
        with h5py.File(p) as f:
            xyz=np.asarray(f['md/vertex'])
            types=[''.join(chr(c) for c in f[r][()].ravel()) for r in f['md/types'][()].ravel()]
            elems={t:np.asarray(f[r]).reshape(-1,{'tri':3,'edg':2,'vtx':1}[t]) for t,r in zip(types,f['md/elem'][()].ravel())}
    else:
        d=loadmat(p,simplify_cells=True)['md'];xyz=d['vertex'].T
        elems={t:np.asarray(e).reshape({'tri':3,'edg':2,'vtx':1}[t],-1).T for t,e in zip(d['types'],d['elem'])}
    if xyz.max()<1: xyz=xyz*1000
    return xyz,elems
def canonical_rows(a):
    return set(map(tuple,np.asarray(a).reshape(len(a),-1)))
def core_signature(p):
    xyz,e=mesh(p)
    tri=e['tri']; inside=np.all((xyz[tri,1]>=30-1e-8)&(xyz[tri,1]<=90+1e-8),axis=1)
    core=tri[inside]
    coords=np.rint(xyz*1e8).astype(np.int64) # 1e-8 mm, much below smallest elements
    vertices=canonical_rows(coords[np.unique(core)])
    triangles=set(tuple(sorted(map(tuple,coords[t]))) for t in core)
    edges={}
    for z in (30,90):
        ee=e['edg']; mask=np.all(np.abs(xyz[ee,1]-z)<1e-8,axis=1)
        edges[z]=set(tuple(sorted(map(tuple,coords[t]))) for t in ee[mask])
    return {'vertices':vertices,'triangles':triangles,'edges':edges,'total_vertices':len(xyz),'total_triangles':len(tri)}

rows=[]
for batch, refs, targets in [
 ('R5',[(role,f'{role}_base') for role in ('A','B')],ROOT/'runs/round5_modal_closure'),
 ('R6',[(f'L{L}_{role}',f'L{L}_{role}_base') for L in (40,64) for role in ('A','B')],ROOT/'runs/round6_applicability')]:
    for ident,refname in refs:
        refpath=targets/'reference'/refname/'mesh.mat'; ref=core_signature(refpath)
        pattern=f'r5_strict_{ident}_P*' if batch=='R5' else f'r6_{ident}_P*'
        for case in sorted((targets/'truth').glob(pattern)):
            targetpath=case/'source_background/mesh.mat'
            if not targetpath.exists(): continue
            trg=core_signature(targetpath)
            rows.append({'batch':batch,'reference':refname,'target':case.name,
                'reference_total_vertices':ref['total_vertices'],'target_total_vertices':trg['total_vertices'],
                'reference_total_triangles':ref['total_triangles'],'target_total_triangles':trg['total_triangles'],
                'reference_core_vertices':len(ref['vertices']),'target_core_vertices':len(trg['vertices']),
                'reference_core_triangles':len(ref['triangles']),'target_core_triangles':len(trg['triangles']),
                'core_vertices_equal_1e_8_mm':ref['vertices']==trg['vertices'],
                'core_triangles_equal_1e_8_mm':ref['triangles']==trg['triangles'],
                'reference_only_core_triangles':len(ref['triangles']-trg['triangles']),
                'target_only_core_triangles':len(trg['triangles']-ref['triangles']),
                'lower_interface_edges_equal':ref['edges'][30]==trg['edges'][30],
                'upper_interface_edges_equal':ref['edges'][90]==trg['edges'][90],
                'reference_path':refpath.as_posix(),'target_path':targetpath.as_posix()})
writec('core_mesh_correspondence.csv',rows)
print('Core mesh comparisons:',len(rows),'all triangle matches:',all(r['core_triangles_equal_1e_8_mm'] for r in rows))
print(json.dumps(rows[0],ensure_ascii=False))

def xmesh_signature(p):
    with h5py.File(p) as f:
        xyz=np.asarray(f['xm/dofs/coords'])*1000
        tri=np.asarray(f['xm/elements/tri/dofs'],dtype=int)
        edg=np.asarray(f['xm/elements/edg/dofs'],dtype=int)
    coords=np.rint(xyz*1e8).astype(np.int64)
    keep=np.all((xyz[tri,1]>=30-1e-8)&(xyz[tri,1]<=90+1e-8),axis=1)
    tri=tri[keep]
    ts=set(tuple(sorted(map(tuple,coords[t]))) for t in tri)
    edges={z:set(tuple(sorted(map(tuple,coords[e]))) for e in edg[np.all(np.abs(xyz[edg,1]-z)<1e-8,axis=1)]) for z in (30,90)}
    return {'total_dofs':len(xyz),'core_dofs':len(np.unique(tri)),'triangles':ts,'edges':edges}
cache={};xr=[]
for r in rows:
    refp=Path(r['reference_path']).with_name('reference.mat')
    trgp=Path(r['target_path']).parent.parent/'background_responses.mat'
    if not trgp.is_file():continue
    if refp not in cache: cache[refp]=xmesh_signature(refp)
    ref=cache[refp];trg=xmesh_signature(trgp)
    xr.append({'batch':r['batch'],'reference':r['reference'],'target':r['target'],
        'reference_total_dofs':ref['total_dofs'],'target_total_dofs':trg['total_dofs'],
        'reference_core_dofs':ref['core_dofs'],'target_core_dofs':trg['core_dofs'],
        'core_quadratic_elements_equal':ref['triangles']==trg['triangles'],
        'lower_quadratic_trace_elements_equal':ref['edges'][30]==trg['edges'][30],
        'upper_quadratic_trace_elements_equal':ref['edges'][90]==trg['edges'][90],
        'reference_path':refp.as_posix(),'target_path':trgp.as_posix()})
writec('quadratic_core_trace_correspondence.csv',xr)
print('Quadratic comparisons',len(xr),'R5 equal',sum(r['core_quadratic_elements_equal'] for r in xr if r['batch']=='R5'),'R6 equal',sum(r['core_quadratic_elements_equal'] for r in xr if r['batch']=='R6'))

r5=ROOT/'runs/round5_modal_closure'
effects=readc(r5/'side_effects.csv')
side=[r for r in effects if r['phase']=='prospective' and int(r['K'])==8 and r['quantity'] in ('gA','gB') and r['effect'].startswith(('upper_','lower_'))]
writec('R5_M8_all64_single_side_checks.csv',side)
failure=[r for r in side if float(r['relative_error_percent'])>10]
writej('R5_side_check_summary.json',{'count':len(side),'point_sign_matches':sum(r['point_sign_correct'].lower()=='true' for r in side),'count_exceeds_10_percent':len(failure),'failures':failure})
print('R5 count/fail:',len(side),len(failure))
near=[r for r in readc(r5/'selected_mesh_effects.csv') if r.get('role')=='B' and float(r.get('sigma_f_S_m',0))==.18]
writej('R5_selected_mesh_evidence.json',near)
r6=ROOT/'runs/round6_applicability'
contrasts=readc(r6/'side_and_rank_intervals.csv')
print('R6 categories',sorted(set(r['category'] for r in contrasts)))
net=[r for r in contrasts if r['category']=='joint_net' and str(r.get('reference_K'))=='8']
if not net: net=[r for r in contrasts if 'joint' in r['category']]
writec('R6_net_contrast_evidence.csv',net)
print('R6 net rows',len(net),'keys',list(contrasts[0]))
