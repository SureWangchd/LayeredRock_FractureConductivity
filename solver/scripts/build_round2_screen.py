"""Freeze six geometric pairs chosen before seeing any inclusion response."""
from round2_manifests import *

source=json.loads((OUT/'core_design_candidates.json').read_text(encoding='utf-8'))
keys=('id','pair_id','configuration_role','configuration_description','h_mm','a_mm','b_mm',
      'theta_deg','cx_mm','cz_mm','L_mm','A_mm','z_interfaces_mm','layer_materials',
      'actual_partition_thicknesses_mm','central_actual_thickness_mm',
      'expected_background_high_area_mm2','expected_replaced_high_area_mm2',
      'expected_replaced_low_area_mm2','interpretation_limit')
cases=[]
for original in source['cases']:
    base={k:original[k] for k in keys if k in original}
    base.update(geometry_family='central_slab_matched_remote_compensation',
                family_id=base['pair_id'],endcap_length_mm=3,
                buffer_bottom_mm=10,buffer_top_mm=10,buffer_material=0,
                equivalence_audit='Non-equivalent: opposite finite end phases; each flat profile centrally symmetric. Reflections x/z rotate fixed 60-degree ellipse; 180-degree rotation preserves each separately; diagonal maps change horizontal strata to vertical.')
    cases.extend(case_triplet(base,'screen'))
print(write_group('screen',cases))
