"""Wave extension with inversion-complement phase symmetry.

The centered ellipse is invariant under a half turn and each background's
phase labels exchange under that turn. Its replaced high/low areas therefore
remain equal analytically. Both background conductivities still require solves.
"""
from round2_manifests import *
tiles=json.loads((OUT/'strict_tile_candidates.json').read_text(encoding='utf-8'))['cases']
cases=[]
for original in tiles:
    if original['h_mm']!=10:continue
    for amplitude in [.5,1.5]:
        base=copy.deepcopy(original)
        token=str(amplitude).replace('.','p')
        base.update(id=base['id']+'_wave'+token,pair_id='strict_h10_L60_wave'+token,
                    family_id='strict_h10_L60_wave'+token,A_mm=amplitude,
                    geometry_family='equal_tiles_wave_inversion_complement',
                    endcap_length_mm=3,
                    interface_amplitudes_mm=[amplitude if 20<z<100 else 0 for z in base['z_interfaces_mm']],
                    interpretation_limit='Half-turn plus phase complement ensures ellipse replacement halves remain equal; actual x/z background conductivities are independently tested against fixed tolerance. Endcap states remain lower-high upper-low when verified.')
        triple=case_triplet(base,'tile_wave',3)
        cases.extend(triple)
        for lv in [1,2]:
            c=copy.deepcopy(triple[-1]);c.update(id=base['id']+f'_m{lv}_full',mesh_level=lv)
            c['parent_full_case']=c['id'];cases.append(c)
print(write_group('tile_wave',cases))
