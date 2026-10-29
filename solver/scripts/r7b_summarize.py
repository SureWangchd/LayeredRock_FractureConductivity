#!/usr/bin/env python3
"""Summarize actual Round-7b artifacts. No COMSOL solve or old-data audit.

Default writes current CSV/JSON only. --report requires all 16 prescribed
targets, or an explicit, evidenced GEOMETRY_BLOCKED.json terminal record.
--report --partial creates a distinctly named progress report if incomplete.
Unmeasured error components stay null; they are never silently replaced by 0.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from datetime import datetime
from pathlib import Path

JOBS = [(s, 'M1', w) for s in ('S1', 'S5') for w in (.18, 10.)] + [
    ('S1', 'M2', .18), ('S1', 'M3', .18), ('S5', 'M2', 10.), ('S5', 'M3', 10.)]


def read_json(path):
    p = Path(path)
    return json.loads(p.read_text(encoding='utf-8-sig')) if p.is_file() else None


def dump(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def csv_write(path, rows, defaults=()):
    fields = list(dict.fromkeys([*defaults, *(k for row in rows for k in row)]))
    with Path(path).open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v
                          for k, v in row.items()} for row in rows)


def num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def ratio(a, b):
    return a/b if num(a) and num(b) and b != 0 else None


def water_token(w):
    return 'water_'+format(w, '.10g').replace('.', 'p')


def resolve_evidence(path, root, out):
    p = Path(path)
    if p.is_absolute():
        return p
    if (out/p).is_file():
        return out/p
    return root/p


def state_valid(r, sample, level, water, frac, H):
    if not isinstance(r, dict) or r.get('status') != 'success':
        return False
    return (r.get('sample') == sample and r.get('input_variant') == level
            and r.get('water_sigma_S_m') == water and r.get('fracture_on') == frac
            and r.get('H_on') == H and num(r.get('G_S')) and r['G_S'] > 0
            and isinstance(r.get('model_path'), str) and Path(r['model_path']).is_file())


def geometry_rows(out):
    rows = []
    for sample, level in dict.fromkeys((s, m) for s, m, _ in JOBS):
        folder = out/'inputs'/f'{sample}_{level}'
        gate = read_json(folder/'GEOMETRY_SUCCESS.json')
        candidates = sorted(folder.rglob('geometry_verification.json')) if folder.is_dir() else []
        geom = read_json(candidates[0]) if candidates else {}
        geom = geom or {}
        phase = read_json(folder/'native_phase_verification.json') or {}
        pars = read_json(folder/'parameters.json') or {}
        inside = geom.get('inside_tetrahedron_checks') or {}
        components = inside.get('tetrahedron_face_connected_components')
        row = dict(sample=sample, input_variant=level, geometry_ready_mph=str(folder/'geometry_ready.mph'),
                   geometry_ready_mph_exists=(folder/'geometry_ready.mph').is_file(),
                   acceptance_record_exists=gate is not None,
                   actual_surface_measurements_available=bool(geom),
                   geometry_representation=pars.get('ellipsoid_representation'),
                   actual_volume_m3=geom.get('actual_volume_m3'),
                   ideal_volume_m3=geom.get('ideal_volume_m3'),
                   actual_volume_relative_error=geom.get('actual_volume_relative_error'),
                   actual_central_thickness_m=geom.get('actual_central_thickness_m'),
                   actual_central_thickness_relative_error=geom.get('actual_central_thickness_relative_error'),
                   maximum_probe_absolute_relative_thickness_error=geom.get('maximum_probe_absolute_relative_thickness_error'),
                   maximum_sampled_normal_error_over_central_thickness=geom.get('maximum_sampled_normal_error_over_central_thickness'),
                   tetrahedron_face_components=components,
                   all_evaluated_geometric_gates_pass=geom.get('all_evaluated_geometric_gates_pass'),
                   original_phase_verified=phase.get('status') == 'ORIGINAL_PHASE_AND_REGION_VERIFIED',
                   original_phase_all_tetrahedra_checked=phase.get('all_tetrahedra_checked_against_original_labels'),
                   protected_plane_crossings=phase.get('protected_plane_crossings'),
                   original_phase_region_volume_error_over_L3=phase.get('maximum_volume_error_over_L3'),
                   surface_verification_json=str(candidates[0]) if candidates else None,
                   geometry_acceptance_json=str(folder/'GEOMETRY_SUCCESS.json') if gate else None,
                   tetrahedra=phase.get('tetrahedra'), vertices=phase.get('vertices'),
                   geometry_gates=geom.get('gates'),
                   fixed_surface_identifier=pars.get('fixed_surface_identifier'),
                   fixed_surface_source=pars.get('fixed_surface_source'),
                   refinement_preserves_frozen_surface=pars.get('refinement_preserves_frozen_surface'),
                   fixed_geometry_refinement_verified=(gate or {}).get('fixed_geometry_refinement_verified'))
        accepted = gate is not None and row['geometry_ready_mph_exists']
        accepted = accepted and row['all_evaluated_geometric_gates_pass'] is True and row['original_phase_verified']
        accepted = accepted and components == 1 and geom.get('tetrahedron_connectivity_independently_checked') is True
        if gate is not None and gate.get('geometry_passed') is False:
            accepted = False
        row['accepted_for_summary'] = bool(accepted)
        rows.append(row)
    return rows


def pair_rows(out, geoms):
    rows, states, histories = [], [], []
    geom_map = {(g['sample'], g['input_variant']): g for g in geoms}
    for sample, level, water in JOBS:
        wf = out/'states'/f'{sample}_{level}'/water_token(water)
        fr = read_json(wf/'kappa_frozen.json') or {}
        tbg = read_json(wf/'T_background'/'SUCCESS.json')
        hpath = fr.get('H_background_success_json')
        hbg = read_json(hpath) if hpath else None
        tr = read_json(wf/'T_target'/'SUCCESS.json')
        hr = read_json(wf/'H_target'/'SUCCESS.json')
        records = {'T_background': tbg, 'H_background': hbg, 'T_target': tr, 'H_target': hr}
        expected = {'T_background': (0, 0), 'H_background': (0, 1), 'T_target': (1, 0), 'H_target': (1, 1)}
        valid = {k: state_valid(r, sample, level, water, *expected[k]) for k, r in records.items()}
        row = dict(sample=sample, input_variant=level, water_sigma_S_m=water,
                   prescribed_pair=True, pair_status='incomplete',
                   geometry_accepted=geom_map[(sample, level)]['accepted_for_summary'],
                   kappa_eq_S_m=fr.get('kappa_eq_S_m'), frozen_background_exists=bool(fr),
                   frozen_status=fr.get('status'), calibration_trials=fr.get('number_background_trials'),
                   same_mesh_background_required=True, full_tensor_matched=False,
                   fixed_surface_identifier=geom_map[(sample, level)]['fixed_surface_identifier'],
                   refinement_preserves_frozen_surface=geom_map[(sample, level)]['refinement_preserves_frozen_surface'],
                   fixed_geometry_refinement_verified=geom_map[(sample, level)]['fixed_geometry_refinement_verified'],
                   **{k+'_valid': v for k, v in valid.items()})
        for name, r in records.items():
            if not r:
                continue
            row[name+'_G_S'] = r.get('G_S')
            row[name+'_model_path'] = r.get('model_path')
            row[name+'_reaction_flux_relative_error'] = r.get('reaction_flux_relative_error')
            row[name+'_energy_relative_error'] = r.get('energy_relative_error')
            if name.endswith('target'):
                states.append(dict(sample=sample, input_variant=level, water_sigma_S_m=water,
                    role=name[0], accepted_state=valid[name], model_path=r.get('model_path'),
                    G_S=r.get('G_S'), G0_S=r.get('G0_S'), g_S=r.get('g_S'),
                    solve_seconds=r.get('solve_seconds'), tetrahedra=r.get('tetrahedra'),
                    vertices=r.get('vertices'), min_mesh_quality=r.get('min_mesh_quality'),
                    stationary_relative_tolerance=r.get('stationary_relative_tolerance'),
                    source_success_json=str(wf/name/'SUCCESS.json')))
        if all(valid.values()) and fr.get('status') == 'FROZEN_BEFORE_TARGETS':
            G0T, G0H, GfT, GfH = [records[k]['G_S'] for k in ('T_background', 'H_background', 'T_target', 'H_target')]
            gT, gH = GfT-G0T, GfH-G0H
            delta, residual = gH-gT, G0H-G0T
            own_bg = all(num(r.get('G0_S')) and abs(r['G0_S']/g0-1) < 1e-12 for r, g0 in ((tr, G0T), (hr, G0H)))
            freeze_bg = (num(fr.get('G0_T_S')) and num(fr.get('G0_H_S'))
                         and abs(fr['G0_T_S']/G0T-1) < 1e-12 and abs(fr['G0_H_S']/G0H-1) < 1e-12)
            conservation = all(num(r.get(k)) and r[k] <= 1e-7 for r in records.values()
                               for k in ('reaction_flux_relative_error', 'energy_relative_error'))
            bgmatch = abs(residual)/G0T <= 1e-7
            nearzero = abs(gT) <= 1e-12*G0T
            row.update(G0_T_S=G0T, G0_H_S=G0H, Gf_T_S=GfT, Gf_H_S=GfH,
                       gT_S=gT, gH_S=gH, delta_gH_minus_gT_S=delta,
                       epsilon_delta_over_gT=None if nearzero else delta/gT,
                       etaT=gT/G0T, etaH=gH/G0H,
                       whole_Gf_relative_difference=(GfH-GfT)/GfT,
                       background_residual_H_minus_T_S=residual,
                       background_residual_relative_G0T=residual/G0T,
                       abs_background_residual_over_gT=ratio(abs(residual), abs(gT)),
                       abs_background_residual_over_delta=ratio(abs(residual), abs(delta)),
                       own_same_mesh_background_confirmed=own_bg,
                       frozen_background_matches_saved_solutions=freeze_bg,
                       conservation_all_four_states_pass=conservation,
                       background_match_pass=bgmatch, near_zero_gT=nearzero)
            row['pair_status'] = 'complete' if own_bg and freeze_bg and conservation and bgmatch and row['geometry_accepted'] else 'computed_checks_incomplete_or_failed'
        hist = wf/'calibration_history.csv'
        if hist.is_file():
            with hist.open(encoding='utf-8-sig', newline='') as f:
                for h in csv.DictReader(f):
                    histories.append({**h,'sample':sample,'input_variant':level,'water_sigma_S_m':water})
        rows.append(row)
    return rows, states, histories


def tolerance_evidence(out, sample, water):
    """Require explicit scope and existing solved-model evidence; unknown stays null."""
    records = []
    for p in sorted((out/'tolerance_checks').rglob('comparison.json')):
        r = read_json(p) or {}
        if r.get('sample') != sample or r.get('water_sigma_S_m') != water:
            continue
        scale = r.get('observed_gain_or_pair_change_scale_S')
        paths = r.get('evidence_model_paths', [])
        if (r.get('status') == 'success' and num(scale) and scale >= 0
                and r.get('applies_to_pair_difference') is True and len(paths) >= 2
                and all(Path(q).is_file() for q in paths)):
            records.append((scale, str(p), r.get('input_variant'), r.get('scope_note')))
    if not records:
        return None, [], '未测或记录不足；不能用残差、守恒误差或0代替容差复算变化'
    return sum(x[0] for x in records), [x[1] for x in records], '仅按明确声明适用的真实配对/增益容差复算记录计入；不是严格误差界'


def convergence_rows(out, pairs):
    mapping = {(r['sample'], r['input_variant'], r['water_sigma_S_m']): r for r in pairs}
    rows = []
    for sample, water in (('S1', .18), ('S5', 10.)):
        solver, sources, note = tolerance_evidence(out, sample, water)
        for low, high in (('M1', 'M2'), ('M2', 'M3')):
            a, b = mapping[(sample, low, water)], mapping[(sample, high, water)]
            frozen = (bool(a.get('fixed_surface_identifier'))
                      and a.get('fixed_surface_identifier') == b.get('fixed_surface_identifier')
                      and b.get('refinement_preserves_frozen_surface') is True)
            frozen = frozen or b.get('fixed_geometry_refinement_verified') is True
            row = dict(sample=sample, water_sigma_S_m=water, from_mesh=low, to_mesh=high,
                       comparison_status='incomplete', solver_tolerance_change_scale_S=solver,
                       fixed_geometry_refinement_confirmed=frozen,
                       solver_scale_evidence=sources, solver_scale_note=note,
                       required_relative_g_change_max=.01, required_absolute_epsilon_change_max=.01,
                       resolution_multiplier=3, resolution_scale_is_rigorous_bound=False)
            if a['pair_status'] == b['pair_status'] == 'complete':
                dT, dH = abs(b['gT_S']-a['gT_S']), abs(b['gH_S']-a['gH_S'])
                bg = abs(a['background_residual_H_minus_T_S'])+abs(b['background_residual_H_minus_T_S'])
                relT, relH = ratio(dT, abs(b['gT_S'])), ratio(dH, abs(b['gH_S']))
                e = abs(b['epsilon_delta_over_gT']-a['epsilon_delta_over_gT']) if num(a['epsilon_delta_over_gT']) and num(b['epsilon_delta_over_gT']) else None
                partial_scale = 3*(dT+dH+bg)
                full_scale = None if solver is None else 3*(dT+dH+bg+solver)
                precision = all(num(x) and x <= .01 for x in (relT, relH, e))
                row.update(comparison_status='computed_fixed_geometry' if frozen else 'computed_fixed_geometry_not_verified', gT_from_S=a['gT_S'], gH_from_S=a['gH_S'],
                           gT_to_S=b['gT_S'], gH_to_S=b['gH_S'],
                           delta_from_S=a['delta_gH_minus_gT_S'], delta_to_S=b['delta_gH_minus_gT_S'],
                           epsilon_from=a['epsilon_delta_over_gT'], epsilon_to=b['epsilon_delta_over_gT'],
                           absolute_gT_change_S=dT, absolute_gH_change_S=dH,
                           relative_gT_change=relT, relative_gH_change=relH,
                           absolute_epsilon_change=e, background_residual_scale_S=bg,
                           observed_partial_resolution_scale_S=partial_scale,
                           complete_observed_resolution_scale_S=full_scale,
                           g_and_epsilon_engineering_targets_pass=precision,
                           delta_exceeds_partial_observed_scale=abs(b['delta_gH_minus_gT_S']) > partial_scale,
                           delta_direction_resolved_with_all_declared_terms=(abs(b['delta_gH_minus_gT_S']) > full_scale) if full_scale is not None and frozen else None,
                           missing_solver_term_blocks_full_scale=solver is None)
            rows.append(row)
    return rows


def old_comparison(root, pairs):
    rows = []
    for r in pairs:
        if r['input_variant'] != 'M1':
            continue
        sample, water = r['sample'], r['water_sigma_S_m']
        wf = root/'runs'/'round7_3d_ct'/'ct'/f'{sample}_base'/water_token(water)
        old_t, old_h = read_json(wf/'T_target'/'SUCCESS.json'), read_json(wf/'H_target'/'SUCCESS.json')
        row = dict(sample=sample, water_sigma_S_m=water, new_pair_status=r['pair_status'],
                   old_label='round7_base_centroid_staircase', new_label='round7b_M1_boundary_fitted_ellipsoid',
                   comparison_scope='几何表示、网格组织及本网格背景重标定的合并敏感性；不是单一误差来源严格隔离',
                   old_targets_rerun=False, old_target_T_json=str(wf/'T_target'/'SUCCESS.json'),
                   old_target_H_json=str(wf/'H_target'/'SUCCESS.json'))
        if old_t and old_h and num(old_t.get('g_S')) and num(old_h.get('g_S')):
            ot, oh = old_t['g_S'], old_h['g_S']
            row.update(old_gT_S=ot, old_gH_S=oh, old_delta_S=oh-ot, old_epsilon=ratio(oh-ot, ot))
        if r['pair_status'] == 'complete':
            row.update(new_gT_S=r['gT_S'], new_gH_S=r['gH_S'], new_delta_S=r['delta_gH_minus_gT_S'], new_epsilon=r['epsilon_delta_over_gT'])
        rows.append(row)
    return rows


def model_index(out):
    rows = []
    for p in sorted(out.rglob('*.mph')):
        rel = p.relative_to(out)
        solved = p.name == 'solution.mph' and (p.parent/'SUCCESS.json').is_file()
        r = read_json(p.parent/'SUCCESS.json') or {}
        rows.append(dict(relative_path=str(rel), model_path=str(p.resolve()), exists=True, bytes=p.stat().st_size,
                         solved_success_record=solved, formal_CT_target=bool(r.get('formal_CT_target', False)),
                         sample=r.get('sample'), input_variant=r.get('input_variant'), state=r.get('state'),
                         category='formal_state' if rel.parts[0] == 'states' else 'geometry_or_technical_pilot'))
    return rows


def technical_evidence(out):
    """Read this round's actual pilot solves once, including failed comparisons.

    A reused direct-CAD reference may occur in several comparison records;
    count an existing solved MPH only once. None are formal CT targets.
    """
    surfaces, solved, comparisons, sources = [], [], [], []
    frozen = read_json(out/'surface_frozen.json') or {}
    selected_branch = frozen.get('selected_pilot_branch')
    for p in sorted(out.glob('pilot*/**/geometry_verification.json')):
        r = read_json(p) or {}; topology = r.get('topology') or {}
        branch = p.relative_to(out).parts[0]
        surfaces.append(dict(path=str(p), branch=branch,
            check_kind='CAD_derived_surface' if p.parent.name == 'surface_check' else 'actual_volume_mesh_boundary',
            selected_source_branch=branch == selected_branch,
            selected_import_meshsize=frozen.get('selected_import_meshsize') if branch == selected_branch else None,
            formal_CT_target=False, boundary_triangles=topology.get('face_count'),
            closed_two_manifold=topology.get('closed_two_manifold'),
            surface_components=topology.get('surface_edge_connected_components'),
            tetrahedron_connectivity_independently_checked=r.get('tetrahedron_connectivity_independently_checked'),
            **{k:r.get(k) for k in ('actual_volume_m3','actual_volume_relative_error',
                'actual_central_thickness_m','actual_central_thickness_relative_error',
                'maximum_probe_absolute_relative_thickness_error',
                'maximum_sampled_normal_error_over_central_thickness','all_evaluated_geometric_gates_pass')}))
    seen_models = set()
    for p in sorted(out.glob('pilot*/electrical/**/SUCCESS.json')):
        r = read_json(p) or {}; model = r.get('model_path')
        if r.get('status') != 'success' or not model or not Path(model).is_file():
            continue
        key = str(Path(model).resolve()).casefold()
        if key in seen_models:
            continue
        seen_models.add(key)
        valid = (r.get('formal_CT_target') is False and r.get('fracture_on') in (0,1)
                 and num(r.get('G_S')) and r['G_S'] > 0)
        solved.append(dict(branch=p.relative_to(out).parts[0], implementation=p.parent.parent.name,
            state=p.parent.name, valid_uniform_state=valid, formal_CT_target=False,
            source_success_json=str(p), model_path=model,
            **{k:r.get(k) for k in ('fracture_on','G_S','G0_S','g_S','tetrahedra','vertices',
                'min_mesh_quality','potential_order','geometry_order','solve_seconds',
                'reaction_flux_relative_error','energy_relative_error','stationary_relative_tolerance')}))
    for p in sorted(out.glob('pilot*/electrical/comparison.json')):
        r = read_json(p) or {}; imp = r.get('imported_union') or {}; cad = r.get('direct_native_CAD') or {}
        models = [q.get(role,{}).get('model_path') for q in (imp,cad) for role in ('background','target')]
        available = all(isinstance(q,str) and Path(q).is_file() for q in models)
        comparisons.append(dict(branch=p.relative_to(out).parts[0], comparison_json=str(p),
            formal_CT_target=False, all_four_solution_model_paths_exist=available,
            direct_reference_reused=(p.parent/'reused_direct_reference.json').is_file(),
            imported_G0_S=imp.get('G0_S'), imported_Gf_S=imp.get('Gf_S'), imported_g_S=imp.get('g_S'),
            direct_CAD_G0_S=cad.get('G0_S'), direct_CAD_Gf_S=cad.get('Gf_S'), direct_CAD_g_S=cad.get('g_S'),
            imported_tetrahedra=imp.get('geometry',{}).get('tetrahedra'),
            direct_CAD_tetrahedra=cad.get('geometry',{}).get('tetrahedra'),
            comparison_engineering_target_passed=r.get('comparison_engineering_target_passed'),
            gain_relative_difference=r.get('gain_relative_difference'),
            target_conductance_relative_difference=r.get('target_conductance_relative_difference'),
            gain_relative_target=r.get('gain_relative_target'),
            whole_conductance_relative_target=r.get('whole_conductance_relative_target'),
            original_interpretation=r.get('interpretation')))
    for p in sorted((out/'source_import').glob('*/*/source_verification.json')):
        r=read_json(p) or {}; pars=read_json(p.parent/'source_parameters.json') or {}
        sources.append(dict(source_verification_json=str(p), sample=pars.get('sample'), mode=pars.get('mode'),
            status=r.get('status'), source_geometry_unchanged=r.get('source_geometry_unchanged'),
            original_CT_resampling=pars.get('original_CT_resampling'),
            **{k:r.get(k) for k in ('tetrahedra','vertices','actual_domain_count','protected_region_plane_crossings',
                'maximum_phase_region_volume_error_over_L3','ellipse_mask_not_used_as_physical_state','seconds')}))
    return surfaces, solved, comparisons, sources


def build_attempt_rows(out):
    """Expose failures and unfinished retries without inventing terminal status."""
    paths=set(out.rglob('FAILURE*.json'))
    for pattern in ('union_attempt_*.json','mesher_attempt_*.json','UNION_SURFACE_COMPLETE.json'):
        paths.update(out.rglob(pattern))
    rows=[]
    for p in sorted(paths):
        r=read_json(p) or {}; action=r.get('action',''); rel=str(p.relative_to(out))
        message=r.get('message',''); compact=message.split('Stack trace:')[0].strip()
        if 'vector too long' in message:
            stage='Union'; diagnosis='COMSOL内部 vector too long；尚未到体网格或电学求解'
        elif 'Degenerate final tet' in message:
            stage='final_volume_mesh_validation'; diagnosis='体网格生成后发现零体积单元；未通过几何验收'
        elif '自由四面体' in message or 'FreeTet' in message:
            stage='FreeTet'; diagnosis='体网格构建失败；具体域/面/坐标见错误记录'
        elif p.name.startswith('union_attempt') or p.name=='UNION_SURFACE_COMPLETE.json':
            stage='Union'; diagnosis=None
        elif p.name.startswith('mesher_attempt'):
            stage='FreeTet'; diagnosis=None
        else:
            stage='configuration_or_other'; diagnosis=None
        if 'technical32' in rel or 'technical32' in action:
            scope='CT32_technical_only'
        elif 'inputs' in p.relative_to(out).parts or action.startswith(('build_S1','build_S5')):
            scope='formal_CT_geometry'
        else:
            scope='uniform_or_API_pilot'
        failed=r.get('status') in ('FAILED','failed','ERROR') or p.name.startswith('FAILURE')
        rows.append(dict(source_json=str(p), scope=scope, action=action or p.name,
            operation_stage=stage, status=r.get('status'), failed=failed,
            method=r.get('method'), selection_mode=r.get('selection_mode'),
            selected_domain_candidates=r.get('selected_domain_candidates'),
            total_domains_before_union=r.get('total_domains_before_union'),
            source_triangles=r.get('source_triangles'), actual_triangles=r.get('actual_triangles'),
            actual_domains=r.get('actual_domains'), all_source_domains_retained=r.get('all_source_domains_retained'),
            hmax_m=r.get('hmax_m'), hmin_m=r.get('hmin_m'),
            failure_summary=diagnosis, message_before_stack=compact or None,
            formal_CT_electrical_result=False,
            status_interpretation='体网格已生成不等于验收通过' if r.get('status')=='VOLUME_MESH_CREATED' else None))
    return rows


def failure_diagnostic_rows(out):
    rows=[]
    for p in sorted(out.rglob('mesh_failure_diagnostics.json')):
        r=read_json(p) or {}
        rows.append(dict(source_json=str(p),
            **{k:r.get(k) for k in ('branch','status','source_model','source_model_not_modified',
                'PDE_solved','vertices','tetrahedra','boundary_triangles','exact_zero_signed_volume_tetrahedra',
                'negative_input_ordering_tetrahedra','negative_signed_volume_interpretation',
                'near_zero_diagnostic_threshold_m3','tetrahedra_at_or_below_near_zero_diagnostic_threshold',
                'near_zero_threshold_is_not_a_solver_acceptance_rule','minimum_positive_tet_volume_m3',
                'minimum_positive_tet_edge_m','maximum_tet_longest_shortest_edge_ratio',
                'exact_zero_area_triangles','minimum_triangle_shape_quality',
                'requested_face_id','reported_failure_face_id_cross_attempt_identity_not_assumed','limitations')}))
    return rows


def terminal_status(root, out, pairs):
    completed = sum(r['pair_status'] == 'complete' for r in pairs)
    status = 'complete_fixed_batch' if completed == 8 else 'in_progress'
    block_path = out/'GEOMETRY_BLOCKED.json'
    block = read_json(block_path)
    if status != 'complete_fixed_batch' and block:
        paths = block.get('evidence_paths', [])
        evidence = [resolve_evidence(q, root, out) for q in paths]
        verified = (block.get('status') in ('geometry_blocked', 'blocked', 'GEOMETRY_BLOCKED')
                    and bool(block.get('reason')) and block.get('official_native_cad_mesh_route_attempted') is True
                    and bool(evidence) and all(p.is_file() for p in evidence))
        if verified:
            status = 'geometry_blocked'
        block['evidence_files_present'] = verified
    return dict(status=status, prescribed_target_count=16, complete_same_mesh_pairs=completed,
                complete_prescribed_targets=2*completed, remaining_pairs=8-completed,
                geometry_block_record=block, geometry_block_json=str(block_path) if block else None,
                scope='Only prescribed CT representation correction; no 2D/synthetic-3D reruns')


def fmt(x, percent=False):
    if not num(x):
        return '未得到'
    return f'{100*x:.4g}%' if percent else f'{x:.8g}'


def report_text(summary, pairs, geoms, conv, old):
    terminal = summary['status']
    lines = ['# 第七轮 CT 定点收尾：真实椭球实施记录与结果', '',
             f"本报告状态：**{terminal}**。正式预定16个含裂缝目标中，已形成具有本网格背景与有效几何记录的完整配对 {summary['complete_same_mesh_pairs']}/8。",
             f"本轮另有 {summary.get('unique_valid_uniform_pilot_solves',0)} 个实际保存的均匀立方体电学解；正式 CT 含裂缝成功记录为 {summary.get('actual_valid_target_records',0)}/16。复用同一 CAD 解只计一次，技术先导不充作 CT 答案。",
             '未完成、不满足几何验收或不能分辨的小差异均保留。先导是均匀立方体技术检查，不计入 CT 科学目标。', '']
    if terminal == 'geometry_blocked':
        block = summary['geometry_block_record']
        lines += ['本轮已按任务书的真实几何阻塞停止条件结束，不能称为完成真实椭球 CT 结果。', '',
                  '**具体阻塞：** '+block['reason'], '', '实际失败/规模证据：']
        lines += ['- '+str(p) for p in block['evidence_paths']]
        lines += ['']
    pilot=summary.get('technical_geometry_checks',[])
    if pilot:
        lines += ['技术先导实际几何检查（两档 CAD 导出曲面及已有体网格检查分别列出，不代替正式128³ CT）：', '',
                  '|分支/检查|面数|体积误差|中央厚度误差|最大探针厚度误差|最大法向偏离/全厚度|已评价门槛|选定曲面来源|',
                  '|---|---:|---:|---:|---:|---:|---|---|']
        for r in pilot:
            lines.append(f"|{r.get('branch')}/{r.get('check_kind')}|{r.get('boundary_triangles')}|{fmt(r.get('actual_volume_relative_error'),True)}|{fmt(r.get('actual_central_thickness_relative_error'),True)}|{fmt(r.get('maximum_probe_absolute_relative_thickness_error'),True)}|{fmt(r.get('maximum_sampled_normal_error_over_central_thickness'),True)}|{r.get('all_evaluated_geometric_gates_pass')}|{r.get('selected_source_branch')}|")
        lines += ['', '完整数值与原始测量路径见 technical_geometry_checks.csv。这里的“通过”仅指该记录实际评价的门槛；只测表面的记录没有独立体单元连通性检查，不能充作正式 CT 几何验收。', '']
    uniform=summary.get('uniform_pilot_states',[])
    if uniform:
        lines += ['均匀立方体实际已解状态（每个 MPH 仅计一次）：', '',
                  '|分支/实现/状态|体单元数|G (S)|g (S)|反力守恒相对误差|能量相对误差|',
                  '|---|---:|---:|---:|---:|---:|']
        for r in uniform:
            lines.append(f"|{r['branch']}/{r['implementation']}/{r['state']}|{r.get('tetrahedra')}|{fmt(r.get('G_S'))}|{fmt(r.get('g_S'))}|{fmt(r.get('reaction_flux_relative_error'))}|{fmt(r.get('energy_relative_error'))}|")
        lines += ['', '无裂缝状态没有新增增益，g栏空缺不表示求解失败。已解 MPH 与每次成功记录的路径在 uniform_pilot_states.csv 和 model_index.csv。', '']
    comparisons=summary.get('uniform_interface_comparisons',[])
    if comparisons:
        lines += ['相同均匀物性下，实际切分导入网格与原生 CAD 的接口对照：', '',
                  '|导入分支|导入g (S)|直接CAD g (S)|新增增益相对差|整个Gf相对差|预定接口门槛|复用CAD参考|',
                  '|---|---:|---:|---:|---:|---|---|']
        for r in comparisons:
            lines.append(f"|{r['branch']}|{fmt(r.get('imported_g_S'))}|{fmt(r.get('direct_CAD_g_S'))}|{fmt(r.get('gain_relative_difference'),True)}|{fmt(r.get('target_conductance_relative_difference'),True)}|{r.get('comparison_engineering_target_passed')}|{r.get('direct_reference_reused')}|")
        lines += ['', '保留最初未通过的对照，不只展示成功分支。增益相对差≤1%、整体Gf相对差≤0.5%是接口先导的预定工程门槛；它不证明原生CAD参考已达到连续介质精确解，也不替代正式CT的固定几何加密。表面形状检查与含裂缝电学检查分别回答不同问题。', '']
    source_checks=summary.get('original_CT_source_checks',[])
    if source_checks:
        lines += ['原始 CT 导入的实际完成部分：', '',
                  '|样本/范围|原体单元数|顶点数|原材料域数|原相体积误差/L³|保护面穿越数|',
                  '|---|---:|---:|---:|---:|---:|']
        for r in source_checks:
            lines.append(f"|{r.get('sample')}/{r.get('mode')}|{r.get('tetrahedra')}|{r.get('vertices')}|{r.get('actual_domain_count')}|{fmt(r.get('maximum_phase_region_volume_error_over_L3'))}|{r.get('protected_region_plane_crossings')}|")
        lines += ['', '导入完成只证明原 CT 网格及标签已读入并完成所列核对；不表示它与真实椭球的相交体网格或电学目标已经成功。', '']
    attempts=[r for r in summary.get('build_attempts',[]) if r['scope'] in ('formal_CT_geometry','CT32_technical_only')]
    if attempts:
        lines += ['CT 构建实际失败、完成步骤与有限重试（不据运行中记录自行判定最终失败）：', '',
                  '|范围/动作|阶段|记录状态|具体事实|', '|---|---|---|---|']
        for r in attempts:
            facts=r.get('failure_summary') or r.get('status_interpretation') or ''
            if r.get('selected_domain_candidates') is not None:
                facts+=f"；候选域 {r['selected_domain_candidates']}；原表面三角形 {r.get('source_triangles')}"
            if r.get('actual_domains') is not None:
                facts+=f"；实际域 {r['actual_domains']}；表面三角形 {r.get('actual_triangles')}"
            if r.get('method'):
                facts+=f"；方法 {r['method']}"
            lines.append(f"|{r['scope']}/{r['action']}|{r['operation_stage']}|{r['status']}|{facts or '见原始记录'}|")
        lines += ['', 'build_attempts.csv保留具体错误文本、域/面/坐标、选择范围和原始JSON路径；原始失败记录保留完整堆栈。同一次调用可能同时留下入口失败和内部步骤记录，表中记录行数不等于独立尝试次数。`vector too long`是已观察到的COMSOL内部错误，不能只凭这句话断定为物理内存耗尽。`VOLUME_MESH_CREATED`也不等于有效网格，零体积单元仍须拒绝。', '']
    for r in summary.get('mesh_failure_diagnostics',[]):
        lines += [f"实际失败网格只读定位（{r.get('branch')}）：{r.get('tetrahedra')}个体单元中，导出顶点坐标直接计算得到 {r.get('exact_zero_signed_volume_tetrahedra')} 个严格零体积单元；最小正体积 {fmt(r.get('minimum_positive_tet_volume_m3'))} m³，最大最长/最短边比 {fmt(r.get('maximum_tet_longest_shortest_edge_ratio'))}。该网格不能用于正式求解。",
                  '负的顶点顺序有向体积仅反映导出顺序，本身不证明几何翻转；诊断近零阈值也不是新增验收规则。跨重试的同一面编号未假定代表同一几何面。详见mesh_failure_diagnostics.csv及原始局部导出。', '']
    lines += ['## 1. 现在求解的裂缝是不是指定有限椭球？', '',
              '目标始终是原生 COMSOL `Ellipsoid`，a=b>c、c/a=0.05、面内轴倾角60°。方案规定S1/S5使用原中心128³、原孔隙标签和真实体素尺度；实际完成程度以下表为准，未构建的模型不算完成。',
              '正式新几何必须是原生椭球导出的受控表面与 CT 原相界面真正相交，再生成相容体网格；不是按元素中心染色。CAD 目标、离散三角表面和线性有限元是三个不同层次。', '',
              '|样本/网格|实际几何验收|相交后原CT相检查|体积误差|中央厚度误差|ρ≤0.9 最大厚度误差|最大法向偏离/中央厚度|',
              '|---|---|---|---:|---:|---:|---:|']
    for g in geoms:
        lines.append(f"|{g['sample']}/{g['input_variant']}|{'通过' if g['accepted_for_summary'] else '未建立/未通过'}|{'通过' if g['original_phase_verified'] else '未建立'}|{fmt(g['actual_volume_relative_error'],True)}|{fmt(g['actual_central_thickness_relative_error'],True)}|{fmt(g['maximum_probe_absolute_relative_thickness_error'],True)}|{fmt(g['maximum_sampled_normal_error_over_central_thickness'],True)}|")
    lines += ['', '几何门槛：体积误差≤0.5%；中央及ρ≤0.9探针厚度误差≤1%；所采样实际边界最大法向偏离≤中央全厚度0.5%；人工椭球自身只有一个共享面连通分量。表面采样不是全局 Hausdorff 严格上界，体积小误差也不能替代厚度检查。', '',
              '## 2. 原先的厚度偏差和碎块消失了吗？', '',
              '上表只填写实际网格的测量记录；缺值表示尚未获得该证据。原阶梯模型的中央厚度约偏离8%，不能把这项旧问题仅凭更细网格宣布消除。实际已有测量列于 technical_geometry_checks.csv 或正式 geometry_verification.csv 对应的原始记录；均匀先导的厚度和连通性通过不代替CT验收。',
              'M1/M2/M3必须共享已冻结的合格离散椭球表面。原 CT 图像分辨率仍为S1 8.683 μm、S5 3.997 μm；有限元加密不会提高原图的孔隙结构分辨率。', '',
              '## 3. 修正后 S1 效应是否保留，S5 小效应是否分辨？', '',
              '定义 gT=GfT−G0T，gH=GfH−G0H，δ=gH−gT，ε=δ/gT。负ε表示简化外围后的裂缝新增收益更低。ε是新增收益的差异，不能写成整个岩石电导率的下降百分比。', '',
              '|样本/网格/流体 S/m|状态|gT (S)|gH (S)|δ (S)|ε|ηT|ηH|整个Gf相对差|',
              '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in pairs:
        lines.append(f"|{r['sample']}/{r['input_variant']}/{r['water_sigma_S_m']:g}|{r['pair_status']}|{fmt(r.get('gT_S'))}|{fmt(r.get('gH_S'))}|{fmt(r.get('delta_gH_minus_gT_S'))}|{fmt(r.get('epsilon_delta_over_gT'),True)}|{fmt(r.get('etaT'),True)}|{fmt(r.get('etaH'),True)}|{fmt(r.get('whole_Gf_relative_difference'),True)}|")
    lines += ['']
    for sample, water in (('S1',.18),('S5',10.)):
        available=[r for r in pairs if r['sample']==sample and r['water_sigma_S_m']==water and r['pair_status']=='complete']
        if not available:
            lines.append(f'{sample}/{water:g} 尚无满足当前汇总条件的完整真实椭球配对；不能回答效应保留或小效应分辨。')
            continue
        r=max(available,key=lambda q:q['input_variant'])
        direction='降低' if r['delta_gH_minus_gT_S']<0 else '增加' if r['delta_gH_minus_gT_S']>0 else '不变'
        lines.append(f"{sample}/{water:g} 已有最高网格 {r['input_variant']}：简化外围后新增收益{direction}，ε={fmt(r['epsilon_delta_over_gT'],True)}。幅度与方向是否可靠，仍按以下最后两档的变化尺度判断，不以是否接近旧25%作验收。")
    lines += ['', '正式配对要求每个流体、每档网格均从真正无裂缝T背景重新标定共同正κ，再冻结后求目标；以上仅对已经完成并核对通过的配对报告该结果。孔隙与裂缝同流体、重叠处只赋一次导率。仅匹配z加载总G0，不意味着核心局部电场或全方向响应相同。', '',
              '## 4. 哪些结果经过固定几何精度检查？', '',
              '方案预选S1/0.18与S5/10做三档配对；S1/10与S5/0.18预定只做主网格。实际是否完成由下表和前述状态表说明，不能把预定网格写成已算结果，也不套用代表状态的检查结论。', '',
              '|样本/流体|两档|固定几何记录|gT相对变化|gH相对变化|ε绝对变化（百分点）|三项数值变化门槛|δ超过完整观测尺度|',
              '|---|---|---|---:|---:|---:|---|---|']
    for r in conv:
        flag = r.get('delta_direction_resolved_with_all_declared_terms')
        direction = '是' if flag is True else ('否' if flag is False else '尚缺项/未完成')
        lines.append(f"|{r['sample']}/{r['water_sigma_S_m']:g}|{r['from_mesh']}→{r['to_mesh']}|{'已确认' if r['fixed_geometry_refinement_confirmed'] else '未确认'}|{fmt(r.get('relative_gT_change'),True)}|{fmt(r.get('relative_gH_change'),True)}|{fmt(100*r['absolute_epsilon_change']) if num(r.get('absolute_epsilon_change')) else '未得到'}|{'通过' if r.get('g_and_epsilon_engineering_targets_pass') else '未通过/未完成'}|{direction}|")
    lines += ['', '最后两档要求gT、gH各变化≤1%，ε绝对变化≤0.01。方向判断另看经验尺度 3×(两项增益网格变化绝对值之和＋两档背景匹配残差＋实际检查过的求解容差变化尺度)。它不是严格误差界；未测的求解项保留为空，不填0。',
              'mesh_convergence.csv同时保存已观测部分尺度和完整尺度；完整尺度缺项时，即使δ超过已观测部分，也只说明在这些已检查项上稳定，不能宣布所有误差已排除。', '',
              '## 5. CT现在能承担多强的论文结论？', '']
    if terminal == 'complete_fixed_batch':
        resolved = [r for r in conv if r['from_mesh']=='M2' and r.get('g_and_epsilon_engineering_targets_pass') and r.get('delta_direction_resolved_with_all_declared_terms') is True]
        lines.append(f"固定目标与几何记录已完成。预选最后两档中，有{len(resolved)}/2组同时满足增益/ε门槛及含已测求解项的方向分辨条件；每组具体范围以上表为准。")
        lines.append('仅对通过相应检查的状态，可以定量说明：保留相同CT核心、同一有限理想椭球、电极区并匹配无裂缝总G0，仍可能因外围结构简化改变裂缝新增导电收益。未分辨的S5或其他细小差值应作为限度报告。')
    elif terminal == 'geometry_blocked':
        lines.append('本轮CT只能交付真实的几何实施尝试、阻塞位置和已得到的局部证据，尚不能承担满足用户指定真实椭球条件的正式定量结论。原第七轮阶梯结果可作为表示敏感性记录，不能改名当成本轮成功结果。')
    else:
        lines.append('本轮仍在进行中，尚不能依据未完成的固定几何/网格验证给CT定量结论。')
    lines += ['无论数值门槛是否通过，这都是CT真实图像结构上的数值补充：人工裂缝是规定理想椭球，骨架/流体导率是合成设定，不是天然裂缝测量或实验实测验证；不声称二维冻结预测器已在CT验证。',
              '整体G0匹配不强制局部背景场相同，所以不能把T/H差异全部归给“加入裂缝之后的外围重分布”。二维与合成三维已有论证范围保持不变。', '',
              '## 原阶梯表示对照与复算', '',
              'old_staircase_vs_ellipsoid.csv仅引用已保存R7基础网格结果；不重跑旧目标。新旧差异同时包含裂缝表示、体网格组织和本网格κ重标定，不能宣称已严格隔离单一误差来源。',
              '入口：scripts/r7b_launch.ps1 调用唯一服务2037，沿用26线程和资源保护；scripts/r7b_solve.m 的 main/refinements 或固定单例名运行。先构建并验收真实几何，才运行电学入口。',
              '汇总：python scripts/r7b_summarize.py --report；打包：python scripts/r7b_package.py。已解大型MPH留在本机，实际路径见model_index.csv。压缩包包含明确白名单的小型原生椭球源、实际曲面/体网格边界及失败单元导出，清单见bundle_manifest.json；不包含巨型MPH/MAT及原CT。',
              '全部失败原始日志和失败JSON保留，新输出不覆盖前七轮。']
    return '\n'.join(lines)+'\n'


def summarize(root, report=False, partial=False):
    root = Path(root).resolve();out=root/'runs'/'round7b_ct_ellipsoid';out.mkdir(parents=True,exist_ok=True)
    geoms=geometry_rows(out);pairs,states,history=pair_rows(out,geoms)
    conv=convergence_rows(out,pairs);old=old_comparison(root,pairs);models=model_index(out)
    summary=terminal_status(root,out,pairs)
    technical,uniform,comparisons,source_checks=technical_evidence(out)
    attempts=build_attempt_rows(out)
    diagnostics=failure_diagnostic_rows(out)
    summary.update(generated_local=datetime.now().astimezone().isoformat(),
                   actual_valid_target_records=sum(s['accepted_state'] for s in states),
                   existing_mph_files=len(models), geometry_accepted_levels=sum(g['accepted_for_summary'] for g in geoms),
                   no_COMSOL_run_performed_by_summary=True, unmeasured_error_terms_kept_null=True,
                   technical_geometry_checks=technical, uniform_pilot_states=uniform,
                   unique_valid_uniform_pilot_solves=sum(r['valid_uniform_state'] for r in uniform),
                   uniform_interface_comparisons=comparisons,
                   interface_comparisons_passed=sum(r['comparison_engineering_target_passed'] is True
                       and r['all_four_solution_model_paths_exist'] for r in comparisons),
                   original_CT_source_checks=source_checks, build_attempts=attempts,
                   mesh_failure_diagnostics=diagnostics,
                   failed_formal_CT_build_records=sum(r['failed'] and r['scope']=='formal_CT_geometry' for r in attempts),
                   failed_formal_CT_entry_calls=sum(r['failed'] and r['scope']=='formal_CT_geometry'
                       and r['action'].startswith(('build_S1','build_S5')) for r in attempts))
    csv_write(out/'ct_response_ellipsoid.csv',pairs);csv_write(out/'target_states.csv',states,('sample','input_variant','water_sigma_S_m','role'))
    csv_write(out/'mesh_convergence.csv',conv);csv_write(out/'geometry_verification.csv',geoms)
    probes=[]
    for g in geoms:
        if not g['surface_verification_json']:continue
        source=Path(g['surface_verification_json']).parent/'geometry_probes_actual.csv'
        if not source.is_file():continue
        with source.open(encoding='utf-8-sig',newline='') as f:
            for row in csv.DictReader(f):
                probes.append({**row,'sample':g['sample'],'input_variant':g['input_variant'],'source_actual_probe_csv':str(source)})
    csv_write(out/'geometry_probes_actual.csv',probes,('sample','input_variant','source_actual_probe_csv'))
    csv_write(out/'calibration_history.csv',history,('sample','input_variant','water_sigma_S_m','step'))
    csv_write(out/'old_staircase_vs_ellipsoid.csv',old);csv_write(out/'model_index.csv',models,('relative_path','model_path','exists'))
    csv_write(out/'technical_geometry_checks.csv',technical,('branch','check_kind','path'))
    csv_write(out/'uniform_pilot_states.csv',uniform,('branch','implementation','state','model_path'))
    csv_write(out/'uniform_interface_comparisons.csv',comparisons,('branch','comparison_json'))
    csv_write(out/'original_CT_source_checks.csv',source_checks,('sample','mode','source_verification_json'))
    csv_write(out/'build_attempts.csv',attempts,('scope','action','operation_stage','status','source_json'))
    csv_write(out/'mesh_failure_diagnostics.csv',diagnostics,('branch','source_json','source_model'))
    dump(out/'summary.json',summary)
    if report:
        terminal=summary['status'] in ('complete_fixed_batch','geometry_blocked')
        if not terminal and not partial:
            raise RuntimeError('Round7b not terminal: use --partial for an explicitly incomplete progress report; final report not written')
        name='研究报告_第七轮CT椭球收尾_CN.md' if terminal else '研究报告_第七轮CT椭球收尾_CN_阶段稿.md'
        (out/name).write_text(report_text(summary,pairs,geoms,conv,old),encoding='utf-8')
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--report',action='store_true');parser.add_argument('--partial',action='store_true')
    args=parser.parse_args();s=summarize(args.root,args.report,args.partial)
    print(json.dumps(s,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
