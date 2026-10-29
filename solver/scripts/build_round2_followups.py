"""Bounded follow-ups; choices and geometry remain explicit and reproducible."""
from round2_manifests import *
import sys, math

screen=json.loads((MANIFESTS/'screen.json').read_text(encoding='utf-8'))
full={c['id']:c for c in screen if not c['passive_ellipse']}
group=sys.argv[1]
cases=[]
if group=='refine':
    # All six selected pairs, rather than only the largest effect.
    for lv in [2,3]:
        for source in full.values():
            c=copy.deepcopy(source)
            c.update(id=c['id'].replace('_full',f'_m{lv}_full'),mesh_level=lv,group=group)
            c['parent_full_case']=c['id']
            cases.append(c)
elif group=='mechanism':
    # Short crack, large L/h and thicker critical low layer; both A and B in last pair.
    choices=[('core_h10_L20_B',0),('core_h10_L60_B',0),
             ('core_h20_L60_B',0),('core_h20_L60_A',1)]
    for prefix,phase in choices:
        base=copy.deepcopy(full[prefix+'_full'])
        a=base['a_mm'];b=base['b_mm']
        primitive=lambda u:b*(u*math.sqrt(1-u*u/a/a)+a*math.asin(u/a))
        removed=primitive(.5)-primitive(-.5)
        threshold=((primitive(2)-primitive(-2))/2-removed)/4
        for lv in [1,2,3]:
            for kind in ['full','cut','side']:
                c=copy.deepcopy(base)
                stem=prefix+f'_mech_m{lv}'
                c.update(id=stem+'_'+kind,group=group,mesh_level=lv,intervention=kind,
                         intervention_phase=phase,cut_u_mm=0,cut_width_mm=1,
                         side_u1_mm=-2,side_u2_mm=2,side_width_mm=threshold,
                         comp_radius_mm=math.sqrt(removed/math.pi),comp_x_mm=20,comp_z_mm=60,
                         parent_full_case=stem+'_full',mechanism_family=prefix,
                         geometric_removal_target_mm2=removed)
                cases.append(c)
elif group=='wave':
    for pair in ['core_h20_L40','core_h10_L60']:
        for amplitude in [.5,1.5]:
            for role in ['A','B']:
                c=copy.deepcopy(full[pair+'_'+role+'_full'])
                c['id']=pair+'_'+role+'_wave'+str(amplitude).replace('.','p')
                c.update(pair_id=pair+'_wave'+str(amplitude),family_id=pair+'_wave'+str(amplitude),
                         A_mm=amplitude,geometry_family='wave_extension_of_matched_flat',
                         interpretation_limit='Interfaces shifted by the same sine except common 10mm buffers. Flat central thicknesses and ellipse fixed; actual wave replacement and both background conductivities must be rechecked, never assumed matched.')
                for key in list(c):
                    if key.startswith('expected_'):del c[key]
                cases.extend(case_triplet(c,group,2))
elif group=='contrast':
    for ratio in [3,30]:
        for role in ['A','B']:
            c=copy.deepcopy(full['core_h20_L40_'+role+'_full'])
            c['id']='core_h20_L40_'+role+f'_contrast{ratio}'
            c.update(pair_id=f'core_h20_L40_contrast{ratio}',family_id=f'core_h20_L40_contrast{ratio}',
                     sigma_high_S_m=.01*ratio,layer_contrast=ratio)
            cases.extend(case_triplet(c,group,3))
else:
    raise ValueError(group)
print(write_group(group,cases))
