"""Additional structural controls chosen before their COMSOL response is known."""
from round2_manifests import *
import math
screen=json.loads((MANIFESTS/'screen.json').read_text(encoding='utf-8'))
full={c['id']:c for c in screen if not c['passive_ellipse']}
tiles=json.loads((OUT/'strict_tile_candidates.json').read_text(encoding='utf-8'))['cases']
cases=[];mechanisms=[]
for source in tiles:
    base=copy.deepcopy(source)
    base.update(geometry_family='equal_tile_internal_rearrangement',family_id=base['pair_id'],endcap_length_mm=3)
    triplet=case_triplet(base,'structural_controls')
    cases.extend(triplet)
    original_full=triplet[-1]
    for lv in [2,3]:
        c=copy.deepcopy(original_full);c['id']=base['id']+f'_m{lv}_full'
        c.update(mesh_level=lv,parent_full_case=c['id']);cases.append(c)
    if base['configuration_role']=='B':
        cutz=65 if base['h_mm']==10 else 50
        mechanisms.append((original_full,(cutz-60)/math.sin(math.pi/3),cutz))

# C keeps the same low bridge and the entire in-ellipse phase map as B.
# Only high-conductivity reservoir thickness outside the ellipse is increased.
base=copy.deepcopy(full['core_h20_L60_B_full'])
a=base['a_mm'];R=math.sqrt((a*math.sin(math.pi/3))**2+(.75)**2)
D=R+7;h=20;w=30-D+h/2
base.update(id='core_h20_L60_C',configuration_role='C',
            configuration_description='same_bridge_and_endcaps_as_B_wider_high_reservoir',
            geometry_family='outside_ellipse_reservoir_control',
            z_interfaces_mm=[10,10+w,60-D,60-h/2,60+h/2,60+D,110-w,110],
            layer_materials=[0,1,0,1,0,1,0,1,0],
            interpretation_limit='Same ellipse phase map and central low slab as B. High reservoirs extend 6mm farther from ellipse on each side; remote high bands compensate exactly. Whether this eases bypass is a result to test, not a design fact.')
base['actual_partition_thicknesses_mm']=[b-a for a,b in zip([0]+base['z_interfaces_mm'],base['z_interfaces_mm']+[120])]
triplet=case_triplet(base,'structural_controls');cases.extend(triplet)
for lv in [2,3]:
    c=copy.deepcopy(triplet[-1]);c.update(id=base['id']+f'_m{lv}_full',mesh_level=lv)
    c['parent_full_case']=c['id'];cases.append(c)
mechanisms.append((triplet[-1],0,60))

# Move the remote high bands of both A/B inward by the same 2mm.
# This leaves all material inside the ellipse untouched and tests dependence
# of the principal A/B effect on the particular compensation placement.
for role in ['A','B']:
    c=copy.deepcopy(full['core_h20_L60_'+role+'_full']);zs=c['z_interfaces_mm'][:]
    zs[0]+=2;zs[1]+=2;zs[-2]-=2;zs[-1]-=2
    c.update(id='remote_shift2_'+role,pair_id='remote_shift2',family_id='remote_shift2',
             geometry_family='remote_compensation_position_control',z_interfaces_mm=zs,
             buffer_bottom_mm=12,buffer_top_mm=12,
             interpretation_limit='Only remote compensation bands moved inward symmetrically by 2mm; central local phase geometry unchanged relative to the corresponding core_h20_L60 role.')
    c['actual_partition_thicknesses_mm']=[b-a for a,b in zip([0]+zs,zs+[120])]
    cases.extend(case_triplet(c,'structural_controls',3))

print(write_group('structural_controls',cases))
cases=[]
for source,u,z in mechanisms:
    a=source['a_mm'];b=source['b_mm']
    primitive=lambda s:b*(s*math.sqrt(1-s*s/a/a)+a*math.asin(s/a))
    removed=primitive(u+.5)-primitive(u-.5)
    threshold=((primitive(u+2)-primitive(u-2))/2-removed)/4
    prefix=source['id'].removesuffix('_full')
    for lv in [1,2,3]:
        for kind in ['full','cut','side']:
            c=copy.deepcopy(source);stem=prefix+f'_mech_m{lv}'
            c.update(id=stem+'_'+kind,group='structural_mechanism',mesh_level=lv,
                     intervention=kind,intervention_phase=0,cut_u_mm=u,cut_width_mm=1,
                     side_u1_mm=u-2,side_u2_mm=u+2,side_width_mm=threshold,
                     comp_radius_mm=math.sqrt(removed/math.pi),comp_x_mm=20,comp_z_mm=z,
                     parent_full_case=stem+'_full',mechanism_family=prefix,
                     geometric_removal_target_mm2=removed)
            cases.append(c)
print(write_group('structural_mechanism',cases))
