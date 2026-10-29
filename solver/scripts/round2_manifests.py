"""Create/extend bounded round-2 COMSOL manifests without touching old cases."""
from pathlib import Path
import json, copy

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'runs' / '2026-10-21_round2'
MANIFESTS = OUT / 'manifests'

def write_group(group, cases):
    MANIFESTS.mkdir(parents=True, exist_ok=True)
    path = MANIFESTS / (group + '.json')
    if path.exists():
        old = json.loads(path.read_text(encoding='utf-8'))
        if old != cases:
            raise RuntimeError(f'Refusing to replace manifest {path}; use a new group.')
    else:
        path.write_text(json.dumps(cases, ensure_ascii=False, indent=2), encoding='utf-8')
    return path

def case_triplet(base, group, mesh=1, suffix=''):
    """Every full inclusion has its own passive-geometry x/z backgrounds."""
    b = copy.deepcopy(base)
    prefix = b.pop('id') + suffix
    b.update(group=group, mesh_level=mesh,
             background_case_z=prefix+'_bgz', background_case_x=prefix+'_bgx')
    cs=[]
    for tail, passive, direction in [('bgz',True,'z'),('bgx',True,'x'),('full',False,'z')]:
        c=copy.deepcopy(b)
        c.update(id=prefix+'_'+tail, passive_ellipse=passive, fracture=True, load=direction,
                 parent_full_case=prefix+'_full')
        cs.append(c)
    return cs

if __name__ == '__main__':
    base=dict(id='smoke_h10_L20',pair_id='smoke',configuration_role='reference',
              h_mm=10,a_mm=10,b_mm=1.5,theta_deg=60,cx_mm=60,cz_mm=60,
              z_interfaces_mm=list(range(10,120,10)),layer_materials=[0,1]*6,
              geometry_family='equal_tiles_reference',equivalence_audit='single reference')
    print(write_group('smoke',case_triplet(base,'smoke')))
