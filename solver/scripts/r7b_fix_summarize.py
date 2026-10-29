#!/usr/bin/env python3
"""Summarize only this machine's new local-fit artifacts, without solving.

The supplied offline diagnostics, CT32 interface pilot, uniform electrical
check and formal 128^3 CT targets are different evidence categories. Existing
Round 7/7b outputs are read only when an explicit reused-reference path names
them; they never contribute new formal target counts. No hashes are computed.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime
from pathlib import Path

from r7b_summarize import (
    JOBS, convergence_rows as base_convergence_rows, csv_write, dump, fmt, geometry_rows, num,
    pair_rows, read_json, resolve_evidence, water_token,
)


def geometry_evidence(out):
    rows = geometry_rows(out)
    for row in rows:
        folder = out/'inputs'/f"{row['sample']}_{row['input_variant']}"
        p = read_json(folder/'parameters.json') or {}
        phase = read_json(folder/'native_phase_verification.json') or {}
        sel = read_json(folder/'selections.json') or {}
        imp = read_json(folder/'IMPORT_COMPLETE_NOT_GEOMETRY_ACCEPTED.json') or {}
        prep = read_json(folder/'PREPROCESS_COMPLETE.json') or {}
        refinement = read_json(folder/'fixed_geometry_refinement.json') or {}
        anchor = read_json(out/'inputs'/f"{row['sample']}_M1"/'parameters.json') or {}
        refined = row['input_variant'] in ('M2','M3')
        method = p.get('refinement_method') or refinement.get('refinement_method')
        if not method and (p.get('fixed_midpoint_refinement') is True
                           or refinement.get('fixed_midpoint_refinement') is True):
            method = 'globally_ordered_original_edge_midpoint_bisection_then_MeshData'
        proof_identity = bool(anchor.get('fixed_surface_identifier')) and (
            p.get('fixed_surface_identifier') == refinement.get('fixed_surface_identifier')
            == anchor.get('fixed_surface_identifier'))
        proof_ok = bool(refined and refinement.get('actual_COMSOL_import_verified') is True
                        and p.get('refinement_preserves_frozen_surface') is True and proof_identity)
        all_faces = read_json(folder/'all_faces_verification.json') or {}
        geom = read_json(folder/'geometry_verification.json') or {}
        export = geom.get('source_metadata') or {}
        # The package validator historically used a Boolean-Union provenance
        # field. This repair must explicitly record inherited cut-cell labels.
        cut = phase.get('ellipse_label_from_explicit_conforming_cut')
        if cut is None:
            cut = 'cut' in str(sel.get('ellipse_selection_source', '')).lower()
        spacing = {'S1':8.683e-6,'S5':3.997e-6}[row['sample']]
        L=128*spacing
        def close_vector(values,expected):
            return (isinstance(values,list) and len(values)==len(expected)
                    and all(num(x) and abs(x-y)<=1e-12*max(abs(y),L) for x,y in zip(values,expected)))
        scale_ok = (p.get('shape_xyz') == [128,128,128] and p.get('protection_regions_present') is True
                    and num(p.get('spacing_m')) and abs(p['spacing_m']/spacing-1)<=1e-12
                    and num(p.get('L_m')) and abs(p['L_m']/L-1)<=1e-12
                    and p.get('source_start_index_xyz')==[86,86,86]
                    and p.get('source_end_exclusive_xyz')==[214,214,214])
        definition_ok = (close_vector(p.get('ellipsoid_center_m'),[.5*L]*3)
                         and close_vector(p.get('ellipsoid_semiaxes_m'),[.25*L,.25*L,.0125*L])
                         and close_vector(p.get('basis_u'),[.5,0,math.sqrt(3)/2])
                         and close_vector(p.get('basis_v'),[0,1,0])
                         and close_vector(p.get('basis_n'),[-math.sqrt(3)/2,0,.5]))
        phase_ok = (phase.get('all_tetrahedra_checked_against_original_labels') is True
                    and phase.get('protected_plane_crossings') == 0
                    and num(phase.get('maximum_volume_error_over_L3'))
                    and phase['maximum_volume_error_over_L3'] < 1e-8)
        row.update(scope='formal_128_CT_geometry', shape_xyz=p.get('shape_xyz'),
                   L_m=p.get('L_m'), spacing_m=p.get('spacing_m'),
                   source_mesh_file=p.get('source_mesh_file'), source_raw=p.get('source_raw'),
                   original_crop_start=p.get('source_start_index_xyz'),
                   original_crop_end_exclusive=p.get('source_end_exclusive_xyz'),
                   ellipsoid_center_m=p.get('ellipsoid_center_m'),
                   ellipsoid_semiaxes_m=p.get('ellipsoid_semiaxes_m'),
                   node_levelset_snap_over_L=p.get('node_levelset_snap_over_L'),
                   node_coordinates_moved_by_snap=False if p.get('node_levelset_snap_over_L') is not None else None,
                   preprocess_record_exists=bool(prep), import_status=imp.get('status'),
                   explicit_cut_phase_provenance=cut,
                   source_scale_and_protection_recorded=scale_ok,
                   original_ellipsoid_definition_preserved=definition_ok,
                   actual_COMSOL_boundary_export=export.get('COMSOL_export',phase.get('actual_boundary_export_file')),
                   geometry_check_uses_preimport_only=export.get('input_is_preimport_mesh'),
                   exported_inside_tetrahedra=phase.get('ellipse_tetrahedra_exported'),
                   exported_inside_vertices=phase.get('ellipse_local_vertices_exported'),
                   all_volume_face_topology_check_available=bool(all_faces),
                   all_volume_faces_checked=all_faces.get('all_volume_faces_checked'),
                   all_volume_face_topology_pass=all_faces.get('all_face_topology_checks_passed'),
                   all_volume_face_records=all_faces.get('face_records'),
                   internal_unmatched_faces=all_faces.get('internal_unmatched_faces'),
                   more_than_two_owner_faces=all_faces.get('more_than_two_owner_faces'),
                   repeated_owner_faces=all_faces.get('repeated_owner_faces'),
                   all_volume_face_check_source=all_faces.get('source'),
                   original_CT_all_element_and_volume_checks_pass=bool(phase_ok),
                   native_definition_mph=str(folder/'native_ellipsoid_definition.mph'),
                   min_mesh_quality=sel.get('min_mesh_quality'),
                   mean_mesh_quality=sel.get('mean_mesh_quality'),
                   refinement_method=method,
                   refinement_proof_required=refined,
                   refinement_actual_COMSOL_import_verified=refinement.get('actual_COMSOL_import_verified'),
                   refinement_proof_matches_M1_identifier=proof_identity if refined else None,
                   fixed_geometry_proof_accepted=proof_ok if refined else None,
                   refinement_actual_COMSOL_min_mesh_quality=refinement.get('actual_COMSOL_min_mesh_quality'),
                   refinement_actual_maximum_overlap_volume_change_m3=refinement.get('maximum_phase_region_ellipse_volume_change_m3'),
                   refinement_actual_relative_ellipse_volume_change=refinement.get('relative_ellipse_volume_change'),
                   refinement_record=refinement or None,
                   label='new_localfit_M1_piecewise_planar_surface_not_previous_Union_surface')
        export_ok=(export.get('input_is_preimport_mesh') is False
                   and bool(row['actual_COMSOL_boundary_export'])
                   and Path(row['actual_COMSOL_boundary_export']).is_file())
        quality_ok=num(sel.get('min_mesh_quality')) and sel['min_mesh_quality']>0
        faces_ok=(all_faces.get('all_volume_faces_checked') is True
                  and all_faces.get('all_face_topology_checks_passed') is True)
        row['accepted_for_summary'] = bool(row['accepted_for_summary'] and scale_ok and definition_ok and phase_ok and cut and export_ok and quality_ok and faces_ok and (not refined or proof_ok))
    return rows


def convergence_rows(out,pairs,geoms=None):
    """Keep numerical formulas, require actual common-M1 refinement evidence.

    The shared previous-round helper accepts a standalone fixed flag as an
    alternative to matching IDs. This round requires both the common M1
    identity and every participating refined mesh's accepted COMSOL proof.
    """
    rows=base_convergence_rows(out,pairs)
    mapping={(g['sample'],g['input_variant']):g for g in (geoms if geoms is not None else geometry_evidence(out))}
    for row in rows:
        sample=row['sample']; levels=(row['from_mesh'],row['to_mesh'])
        anchor=mapping[(sample,'M1')]
        compared=[mapping[(sample,level)] for level in levels]
        identifier=anchor.get('fixed_surface_identifier')
        identity=bool(identifier) and all(g.get('fixed_surface_identifier')==identifier for g in compared)
        accepted=anchor.get('accepted_for_summary') is True and all(g.get('accepted_for_summary') is True for g in compared)
        proof=all(g.get('fixed_geometry_proof_accepted') is True for g in compared if g['input_variant'] in ('M2','M3'))
        fixed=bool(identity and accepted and proof)
        row.update(fixed_geometry_refinement_confirmed=fixed,
                   fixed_M1_surface_identifier=identifier,
                   matching_M1_surface_identifiers=bool(identity),
                   compared_geometries_actually_accepted=bool(accepted),
                   compared_refined_geometries_have_actual_COMSOL_proof=bool(proof),
                   refinement_methods={g['input_variant']:g.get('refinement_method') for g in compared},
                   fixed_geometry_confirmation_basis='Common accepted M1 surface ID and actual accepted COMSOL import proof for every compared M2/M3 mesh')
        if row['comparison_status'].startswith('computed_'):
            row['comparison_status']='computed_fixed_geometry' if fixed else 'computed_fixed_geometry_not_verified'
            if not fixed:
                row['delta_direction_resolved_with_all_declared_terms']=None
    return rows


def uniform_evidence(out):
    folder = out/'uniform_S1'
    cfg = read_json(folder/'pilot_configuration.json') or {}
    summary = read_json(folder/'comparison.json') or {}
    gate = read_json(folder/'PILOT_ELECTRICAL_SUCCESS.json') or {}
    rows = []
    for kind in ('imported_localfit','direct_native_CAD'):
        pair = summary.get(kind) or read_json(folder/kind/'PAIR_SUCCESS.json') or {}
        for name in ('background','target'):
            r = pair.get(name) or read_json(folder/kind/name/'SUCCESS.json') or {}
            if not r:
                continue
            path = r.get('model_path')
            valid = (r.get('status') == 'success' and r.get('formal_CT_target') is False
                     and num(r.get('G_S')) and r['G_S'] > 0 and bool(path) and Path(path).is_file())
            rows.append(dict(scope='uniform_electrical_interface_only', kind=kind, state=name,
                valid_uniform_state=bool(valid), formal_CT_target=False,
                model_path=path, G_S=r.get('G_S'), g_S=r.get('g_S'),
                reaction_flux_relative_error=r.get('reaction_flux_relative_error'),
                energy_relative_error=r.get('energy_relative_error'),
                tetrahedra=r.get('tetrahedra'), vertices=r.get('vertices'),
                min_mesh_quality=r.get('min_mesh_quality'), solve_seconds=r.get('solve_seconds'),
                stationary_relative_tolerance=r.get('stationary_relative_tolerance'),
                reused_reference=(kind == 'direct_native_CAD' and bool(path)
                                  and not Path(path).is_relative_to(out))))
    all_four = len(rows) == 4 and all(r['valid_uniform_state'] for r in rows)
    gain_error = summary.get('gain_relative_difference')
    total_error = summary.get('target_conductance_relative_difference')
    accepted = (bool(gate) and gate.get('status') == 'success' and all_four
                and summary.get('comparison_engineering_target_passed') is True
                and num(gain_error) and gain_error <= .01
                and num(total_error) and total_error <= .005)
    for r in rows:
        if not all(num(r.get(k)) and r[k] <= 1e-7
                   for k in ('reaction_flux_relative_error','energy_relative_error')):
            accepted = False
    record = dict(scope='uniform_electrical_interface_only', formal_CT_evidence=False,
                  comparison_exists=bool(summary), accepted_for_summary=bool(accepted),
                  source_json=str(folder/'comparison.json') if summary else None,
                  gain_relative_difference=gain_error, target_conductance_relative_difference=total_error,
                  original_thresholds=dict(gain=.01, total=.005, conservation=1e-7),
                  all_four_real_saved_states=all_four, configuration=cfg or None,
                  reused_direct_reference=read_json(folder/'reused_direct_reference.json'))
    return record, rows


def technical_evidence(out):
    folder = out/'technical32'
    p = read_json(folder/'parameters.json') or {}
    g = read_json(folder/'geometry_verification.json') or {}
    phase = read_json(folder/'native_phase_verification.json') or {}
    prep = read_json(folder/'PREPROCESS_COMPLETE.json') or {}
    imp = read_json(folder/'IMPORT_COMPLETE_NOT_GEOMETRY_ACCEPTED.json') or {}
    return dict(scope='existing_32_CT_interface_pilot_only', formal_CT_evidence=False,
        parameters_exist=bool(p), shape_xyz=p.get('shape_xyz'), L_m=p.get('L_m'),
        preprocessing_status=prep.get('status'), imported_model_exists=(folder/'geometry_ready.mph').is_file(),
        import_status=imp.get('status'), original_phase_status=phase.get('status'),
        geometry_checks_available=bool(g),
        geometric_gates_pass=g.get('all_evaluated_geometric_gates_pass'),
        tetrahedron_connectivity_independently_checked=g.get('tetrahedron_connectivity_independently_checked'),
        inside_tetrahedron_checks=g.get('inside_tetrahedron_checks'),
        actual_volume_relative_error=g.get('actual_volume_relative_error'),
        maximum_probe_absolute_relative_thickness_error=g.get('maximum_probe_absolute_relative_thickness_error'),
        maximum_sampled_normal_error_over_central_thickness=g.get('maximum_sampled_normal_error_over_central_thickness'),
        source_geometry_json=str(folder/'geometry_verification.json') if g else None,
        source_model=str(folder/'geometry_ready.mph') if (folder/'geometry_ready.mph').is_file() else None)


def checked_pairs(out, geoms, uniform):
    pairs, targets, histories = pair_rows(out, geoms)
    for row in pairs:
        wf=out/'states'/f"{row['sample']}_{row['input_variant']}"/water_token(row['water_sigma_S_m'])
        pr=read_json(wf/'PAIR_SUCCESS.json') or {}
        row['localfit_uniform_gate_pass'] = uniform['accepted_for_summary']
        row['pair_success_record_exists'] = pr.get('status') == 'success'
        agreement = True
        for k in ('gT_S','gH_S','delta_gH_minus_gT_S'):
            x,y=row.get(k),pr.get(k)
            agreement = agreement and num(x) and num(y) and abs(x-y) <= 1e-12*max(abs(x),abs(y),1e-20)
        row['pair_record_matches_recomputed_gains'] = bool(agreement)
        for role in ('T','H'):
            target=read_json(wf/f'{role}_target'/'SUCCESS.json') or {}
            row[f'{role}_target_formal_flag'] = target.get('formal_CT_target') is True
        if row['pair_status'] == 'complete' and not (
                uniform['accepted_for_summary'] and row['pair_success_record_exists'] and agreement
                and row['T_target_formal_flag'] and row['H_target_formal_flag']):
            row['pair_status']='computed_checks_incomplete_or_failed'
    for target in targets:
        wf=out/'states'/f"{target['sample']}_{target['input_variant']}"/water_token(target['water_sigma_S_m'])
        record=read_json(wf/f"{target['role']}_target"/'SUCCESS.json') or {}
        target['accepted_state']=bool(target['accepted_state'] and record.get('formal_CT_target') is True)
    return pairs,targets,histories


def failure_rows(out):
    rows=[]
    for p in sorted(out.rglob('*.json')):
        name=p.name.upper()
        if not any(token in name for token in ('FAIL','BLOCK','UNRESOLVED','RESOURCE_FLOOR','MEMORY_FLOOR','TIMEOUT','OPTIONAL_DISPOSITION')):
            continue
        r=read_json(p) or {}; rel=p.relative_to(out)
        if ('tolerance_checks' in rel.parts or 'tolpair' in str(r.get('stage',''))
                or 'tolpair' in str(r.get('action',''))):
            scope='optional_fixed_state_tolerance'
        elif 'technical32' in rel.parts or 'technical32' in str(r.get('action','')):
            scope='CT32_technical_only'
        elif 'uniform_S1' in rel.parts or 'uniform' in str(r.get('action','')):
            scope='uniform_electrical_only'
        elif 'states' in rel.parts or str(r.get('action','')).startswith('solve_'):
            scope='formal_CT_electrical'
        elif 'inputs' in rel.parts or str(r.get('action','')).startswith(('import_S','refine_S')):
            scope='formal_CT_geometry'
        else:
            scope='execution_or_preprocess'
        rows.append(dict(source_json=str(p), scope=scope, status=r.get('status'),
                         action=r.get('action'), sample=r.get('sample'),
                         input_variant=r.get('input_variant',r.get('level')),
                         message=(r.get('message') or r.get('reason') or r.get('note')
                                  or (r['failure'].get('message') if isinstance(r.get('failure'),dict)
                                      else r.get('failure'))), identifier=r.get('identifier'),
                         time_local=r.get('time_local',r.get('created_local'))))
    return rows


def indices(out):
    models=[]; data=[]
    for p in sorted(out.rglob('*')):
        if not p.is_file():
            continue
        rel=p.relative_to(out)
        if p.suffix.lower() == '.mph':
            r=read_json(p.parent/'SUCCESS.json') or {}
            scope=('formal_CT_state' if rel.parts[0]=='states' else
                   'CT32_technical_only' if rel.parts[0]=='technical32' else
                   'uniform_electrical_only' if rel.parts[0]=='uniform_S1' else
                   'failed_attempt_reference' if rel.parts[0]=='failed_attempts' else
                   'analytic_ellipsoid_parameter_definition' if p.name=='native_ellipsoid_definition.mph' else 'geometry_or_checkpoint')
            solved=(p.name == 'solution.mph' and r.get('status')=='success'
                    and r.get('model_path') and Path(r['model_path']).resolve()==p.resolve())
            if p.name == 'solution_strict.mph':
                r=read_json(p.parent/'comparison.json') or {}
                scope='optional_fixed_state_tolerance'
                solved=(r.get('status')=='success' and r.get('strict_model_path')
                        and Path(r['strict_model_path']).resolve()==p.resolve())
            models.append(dict(model_path=str(p.resolve()), relative_path=str(rel), bytes=p.stat().st_size,
                scope=scope, solved_success_record=bool(solved), formal_CT_target=bool(solved and r.get('formal_CT_target')),
                sample=r.get('sample'), input_variant=r.get('input_variant'), state=r.get('state'),
                water_sigma_S_m=r.get('water_sigma_S_m'), G_S=r.get('strict_G_S') if p.name=='solution_strict.mph' else r.get('G_S')))
        elif p.suffix.lower() in ('.mat','.bin','.npz','.npy'):
            data.append(dict(path=str(p.resolve()), relative_path=str(rel), bytes=p.stat().st_size,
                             kind='actual_COMSOL_boundary_export' if p.name=='actual_boundary.mat' else
                                  'preimport_cut_mesh' if p.name in ('preimport_actual_boundary.mat','vertices_m.bin','tets0.bin','entity_code.bin') else 'raw_auxiliary_or_field'))
    return models,data


def terminal_status(root,out,pairs):
    complete=sum(p['pair_status']=='complete' for p in pairs)
    status='complete_fixed_batch' if complete==8 else 'in_progress'
    blocked=None
    path=out/'LOCALFIT_BLOCKED.json'
    if path.is_file():
        blocked=read_json(path) or {}
        evidence=[resolve_evidence(p,root,out) for p in blocked.get('evidence_paths',[])]
        valid=(bool(blocked.get('reason')) and blocked.get('localfit_route_actually_attempted') is True
               and blocked.get('status') in ('blocked','localfit_blocked','LOCALFIT_BLOCKED')
               and bool(evidence) and all(p.is_file() for p in evidence))
        blocked['evidence_files_present']=valid
        if status!='complete_fixed_batch' and valid:
            status='localfit_blocked'
    return dict(status=status, prescribed_target_count=16, complete_same_mesh_pairs=complete,
                complete_prescribed_targets=complete*2, remaining_pairs=8-complete,
                explicit_localfit_block_record=blocked)


def report_text(summary,pairs,geoms,conv):
    finished=summary['status']=='complete_fixed_batch'
    lines=['# 第七轮B：局部贴体 CT 椭球修复与实算', '',
        f"状态：**{summary['status']}**。完成具有本档背景标定、几何和电学检查的正式配对 **{summary['complete_same_mesh_pairs']}/8**，对应 **{summary['complete_prescribed_targets']}/16** 个含裂缝目标。",
        f"额外四状态严格容差比较完成 **{summary['complete_optional_tolerance_pairs']}/2**；该计数与正式目标分开。未完成的容差误差项保持未知，不因主批完成而变为通过。",
        '','本报告只汇总本机 `runs/round7b_localfit` 的实际新结果。修复包作者的离线测试、32³接口先导、均匀介质对照均不计为正式 CT 目标。',
        '', '## 先回答当前最重要的问题', '',
        f"1. **新网格能否导入并求解？** 实际正式几何验收 {summary['accepted_formal_geometries']}/6 档；新网格均匀电学关口"+
        ('已通过。' if summary['uniform_check']['accepted_for_summary'] else '尚未通过或证据未齐。')+
        f"正式 CT 成功求解记录 {summary['actual_valid_target_records']}/16；只有完成同网格 T/H 背景与检查的配对才用于效应解释。",
        '2. **椭球厚度和连通性是否可靠？** 下表列实际 COMSOL 导出边界的体积、65条厚度线和连通性结果；未测项目保留空缺。通过离散几何门槛不能单独证明电学收敛。',
        '3. **原 CT 结构是否保留？** 每个正式网格分别核对原128³标签、原尺度、保护分区、逐单元体素一致性及各相体积；其结果与几何验收分开记录。新单元继承原单元标签，不以新质心重染材料。',
        '4. **S1/S5效应是否稳定？** 仅从下列完整新配对及固定M1界面的 M2/M3 趋势判断。缺少精度检查时只能报告该网格上的变化，不能宣布稳定符号。',
        '5. **CT能承担多强的结论？** '+
        ('已完成规定数据集，结论的可信范围仍受下方几何、增益收敛和小配对差分辨结果约束。' if finished else
         '当前只能支持已实际完成且通过检查的项目；未完成的S1/S5目标或精度趋势不能由技术先导、旧阶梯结果或包内数字代替。'),
        '', '## 实际几何与原标签', '',
        '|样本/网格|正式验收|四面体|体积偏差|最大厚度偏差|法向距离/中央厚度|椭球面连通分量|原标签全单元检查|',
        '|---|---|---:|---:|---:|---:|---:|---|']
    for g in geoms:
        lines.append('|'+ '|'.join([g['sample']+'/'+g['input_variant'],str(g['accepted_for_summary']),
            str(g.get('tetrahedra') or '未得到'),fmt(g.get('actual_volume_relative_error'),True),
            fmt(g.get('maximum_probe_absolute_relative_thickness_error'),True),
            fmt(g.get('maximum_sampled_normal_error_over_central_thickness'),True),
            str(g.get('tetrahedron_face_components') if g.get('tetrahedron_face_components') is not None else '未得到'),
            str(g['original_CT_all_element_and_volume_checks_pass'])])+'|')
    lines += ['', 'M1是新的受控分片平面椭球界面，并非上次失败Union的CAD三角面。几何目标仍为 a=b=0.25L、c=0.0125L、中心(0.5L,0.5L,0.5L)、60°倾斜。原生Ellipsoid另存定义模型，实际求解截面取导入体网格。节点有符号距离允许在记录阈值内置零，坐标不随之移动；实际几何门槛没有放宽。',
        '', '曲面近似、P1有限元离散、原CT图像分辨率是三项不同限制。M2/M3固定M1的实际曲面做线性h细化；不会恢复CT中未解析的细孔或获得数学上光滑的椭球。',
        '', '全体积网格公共面核对单列于all_volume_faces.csv：每个内部面必须恰有两个不同单元拥有；单侧面仅允许在试样六个外表面。该离线拓扑检查注明输入二进制来源，COMSOL导入后还需通过实际顶点/单元映射、原标签和椭球导出检查。它不等同于电学求解成功。',
        '', '## 技术先导与均匀电学关口', '',
        f"32³实际导入模型存在：{summary['technical32']['imported_model_exists']}；实际几何门槛：{summary['technical32']['geometric_gates_pass']}。仅用于接口检查。",
        f"新网格均匀电学与原生CAD的新增电导相对差：{fmt(summary['uniform_check']['gain_relative_difference'],True)}；总电导相对差：{fmt(summary['uniform_check']['target_conductance_relative_difference'],True)}。",
        '均匀对照临时统一原孔隙和骨架电性，不能解释真实CT的T/H效应。相同原生CAD参考若复用，单独标明，不重复计算为新解。',
        '', '## 正式T/H结果', '',
        'G单位S。gT=Gf,T−G0,T；gH=Gf,H−G0,H；δ=gH−gT；ε=δ/gT。H的上下外围共用κ，仅以本流体、本网格的无裂缝跨层背景标定；不宣称匹配完整方向电导张量。',
        '', '|样本/网格/流体S·m⁻¹|状态|gT (S)|gH (S)|δ (S)|ε|背景相对残差|',
        '|---|---|---:|---:|---:|---:|---:|']
    for p in pairs:
        lines.append('|'+ '|'.join([f"{p['sample']}/{p['input_variant']}/{p['water_sigma_S_m']:g}",p['pair_status'],
            fmt(p.get('gT_S')),fmt(p.get('gH_S')),fmt(p.get('delta_gH_minus_gT_S')),
            fmt(p.get('epsilon_delta_over_gT'),True),fmt(p.get('background_residual_relative_G0T'),True)])+'|')
    lines += ['', '孔隙和裂缝采用同档流体；椭球内原来的孔隙/骨架标签保留，关闭裂缝可恢复原CT。裂缝启用时替换材料，不重复叠加孔隙与裂缝导电。电极接触缓冲与核心保持原定义。',
        '', '## 固定界面精度与小差值分辨', '',
        '|样本/流体|网格比较|固定界面|gT变化|gH变化|ε变化（百分点）|工程精度门槛|已测经验尺度下δ方向|',
        '|---|---|---|---:|---:|---:|---|---|']
    for c in conv:
        lines.append('|'+ '|'.join([f"{c['sample']}/{c['water_sigma_S_m']:g}",c['from_mesh']+'→'+c['to_mesh'],
            str(c['fixed_geometry_refinement_confirmed']),fmt(c.get('relative_gT_change'),True),
            fmt(c.get('relative_gH_change'),True),fmt(100*c['absolute_epsilon_change']) if num(c.get('absolute_epsilon_change')) else '未得到',
            str(c.get('g_and_epsilon_engineering_targets_pass','未得到')),
            str(c.get('delta_direction_resolved_with_all_declared_terms')) if c.get('delta_direction_resolved_with_all_declared_terms') is not None else '未分辨/误差项未齐'])+'|')
    lines += ['', '工程门槛为M2→M3的gT、gH各变化≤1%，ε变化≤1个百分点。小差值方向另与3×(两侧增益网格变化＋背景残差＋实际求解容差复算变化)比较；这是观察到的经验尺度，不是严格数学误差界。未做容差复算时该项为null，不能用0、线性残差或电流守恒误差代替。',
        '本轮四状态容差项只在各预选配对的M3测量；M1→M2行使用该末档经验尺度作趋势参考，不表示M1、M2各自已完成严格容差复算。最终局部分辨判断以M2→M3及M3四状态证据为准。',
        '', '## 失败与复算', '',
        f"本目录保留 {summary['failure_records']} 条失败/未决记录，见 failure_records.csv。重复记录可能来自同一错误的内部记录和上层入口，不能当作独立失败次数。",
        '新模型路径见model_index.csv；网格和场数组见data_index.csv；标定试验见background_calibration.csv。大量模型和原CT文件保留在项目原路径，轻量评审包中的索引不等于模型内嵌。',
        '断点续算入口为scripts/r7b_fix_continue.py；已通过的接口、几何、均匀关口和正式配对直接复用。必要的手动动作通过scripts/r7b_fix_launch.ps1执行，顺序及附加超时处置见README。',
        '本轮未重跑二维或合成三维，未更新冻结预测器。旧约25%阶梯结果属于历史表示敏感性信息，不作为修正后效应必须保留的目标。']
    block=summary.get('explicit_localfit_block_record')
    if block:
        lines += ['', '实际阻塞记录：'+str(block.get('reason')),
                  '证据：'+', '.join(str(x) for x in block.get('evidence_paths',[]))]
    return '\n'.join(lines)+'\n'


def summarize(root,report=False,partial=False):
    root=Path(root).resolve(); out=root/'runs'/'round7b_localfit'; out.mkdir(parents=True,exist_ok=True)
    geoms=geometry_evidence(out); uniform,ustates=uniform_evidence(out)
    pairs,targets,histories=checked_pairs(out,geoms,uniform)
    conv=convergence_rows(out,pairs,geoms); failures=failure_rows(out); models,data=indices(out)
    summary=terminal_status(root,out,pairs)
    tolerance_pairs=[]
    for sample,water in (('S1',.18),('S5',10.)):
        path=out/'tolerance_checks'/f'{sample}_M3_{water_token(water)}_pair_four_states'/'comparison.json'
        record=read_json(path) or {}
        complete=bool(record.get('status')=='success' and record.get('all_four_state_checks_completed') is True
                      and record.get('applies_to_pair_difference') is True
                      and len(record.get('strict_model_paths',[]))==4
                      and all(Path(p).is_file() for p in record.get('strict_model_paths',[])))
        tolerance_pairs.append(dict(sample=sample,water_sigma_S_m=water,complete=complete,
            comparison_json=str(path) if path.is_file() else None,
            observed_scale_S=record.get('observed_gain_or_pair_change_scale_S') if complete else None))
    summary.update(schema='r7b_localfit_evidence_summary_v1',generated_local=datetime.now().astimezone().isoformat(),
        optional_tolerance_pairs=tolerance_pairs,complete_optional_tolerance_pairs=sum(r['complete'] for r in tolerance_pairs),
        accepted_formal_geometries=sum(g['accepted_for_summary'] for g in geoms),
        actual_valid_target_records=sum(t['accepted_state'] for t in targets),
        actual_valid_target_note='Saved successful target states; complete_prescribed_targets additionally requires all pair checks',
        existing_mph_files=len(models), existing_solved_model_records=sum(m['solved_success_record'] for m in models),
        failure_records=len(failures), technical32=technical_evidence(out), uniform_check=uniform,
        uniform_pilot_states=ustates,
        supplied_package_offline_results_are_new_CT_evidence=False, old_target_results_used=False,
        checksums_generated=False, old_results_modified=False,
        scope='Same 128^3 CT and prescribed 8 main plus 8 refined targets; explicit local cut mesh; no global Union/FreeTet retry')
    csv_write(out/'geometry_verification.csv',geoms,('sample','input_variant','accepted_for_summary'))
    facekeys=('sample','input_variant','all_volume_face_topology_check_available','all_volume_faces_checked',
              'all_volume_face_topology_pass','all_volume_face_records','internal_unmatched_faces',
              'more_than_two_owner_faces','repeated_owner_faces','all_volume_face_check_source')
    csv_write(out/'all_volume_faces.csv',[{k:g.get(k) for k in facekeys} for g in geoms],facekeys)
    csv_write(out/'ct_response_localfit.csv',pairs,('sample','input_variant','water_sigma_S_m','pair_status'))
    csv_write(out/'formal_target_states.csv',targets,('sample','input_variant','water_sigma_S_m','role'))
    csv_write(out/'background_calibration.csv',histories,('sample','input_variant','water_sigma_S_m','step'))
    csv_write(out/'mesh_convergence.csv',conv,('sample','water_sigma_S_m','from_mesh','to_mesh'))
    csv_write(out/'uniform_pilot_states.csv',ustates,('kind','state','valid_uniform_state','model_path'))
    csv_write(out/'technical32_status.csv',[summary['technical32']],('scope','imported_model_exists'))
    csv_write(out/'failure_records.csv',failures,('source_json','scope','status','message'))
    csv_write(out/'model_index.csv',models,('model_path','relative_path','scope','solved_success_record'))
    csv_write(out/'data_index.csv',data,('path','relative_path','kind','bytes'))
    dump(out/'summary.json',summary)
    if report:
        terminal=summary['status'] in ('complete_fixed_batch','localfit_blocked')
        if not terminal and not partial:
            raise RuntimeError('Formal report requires all prescribed checked pairs or an evidenced LOCALFIT_BLOCKED record; use --partial for progress')
        name='研究报告_第七轮B局部贴体修复_CN.md' if terminal else '研究报告_第七轮B局部贴体修复_CN_阶段稿.md'
        (out/name).write_text(report_text(summary,pairs,geoms,conv),encoding='utf-8')
    return summary


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--report',action='store_true'); ap.add_argument('--partial',action='store_true')
    a=ap.parse_args(); r=summarize(a.root,a.report,a.partial)
    print(json.dumps({k:r[k] for k in ('status','complete_same_mesh_pairs','actual_valid_target_records','accepted_formal_geometries')},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
