#!/usr/bin/env python3
"""Build a new, explicit conforming cut mesh from the ORIGINAL CT mesh.
This bypasses COMSOL mesh Union/FreeTet, not COMSOL electrical solving.
Never reads a target solution and never changes CT labels or ellipsoid size.
"""
from __future__ import annotations
import argparse,json,sys,time
from datetime import datetime
from pathlib import Path
import numpy as np
from scipy.io import loadmat,savemat
from r7b_fix_localfit import fit,tetvol

def read_mat(path):
 names=['vertices_m','tets0','tet_entity','labels_xyz','L_m','spacing_m']
 try:return loadmat(path,variable_names=names)
 except (NotImplementedError,ValueError):
  import h5py
  with h5py.File(path,'r') as f:
   # MATLAB v7.3 stores reversed dimensions; restore the original MATLAB arrays.
   return {k:np.array(f[k]).transpose(tuple(reversed(range(f[k].ndim)))) for k in names}

def run(root:Path,sample:str,scope:str,tol:float,max_tets:int):
 root=root.resolve()
 src=root/'runs'/'round7_3d_ct'/('ct_import_pilot' if scope=='technical32' else f'ct_inputs/{sample}_base')
 dst=root/'runs'/'round7b_localfit'/('technical32' if scope=='technical32' else f'inputs/{sample}_M1')
 dst.mkdir(parents=True,exist_ok=True)
 done=dst/'PREPROCESS_COMPLETE.json'
 if done.exists():
  d=json.loads(done.read_text(encoding='utf-8'))
  for name,size in d['file_sizes_bytes'].items():
   if not (dst/name).is_file() or (dst/name).stat().st_size!=size:raise RuntimeError(f'Frozen file missing or resized: {name}')
  print(f'Reuse checked preprocessing: {dst}');return
 if any((dst/n).exists() for n in ['vertices_m.bin','tets0.bin','entity_code.bin']):
  raise RuntimeError(f'Partial checkpoint in {dst}. Preserve/rename this directory before a new geometry attempt, never overwrite.')
 start=time.perf_counter();d=read_mat(src/'mesh_input.mat');p0=json.loads((src/'parameters.json').read_text(encoding='utf-8-sig'))
 labels=np.asarray(d['labels_xyz'],np.uint8);n=32 if scope=='technical32' else 128
 if labels.shape!=(n,n,n) or not np.isin(labels,[0,1]).all():raise ValueError('Original CT label shape/values changed')
 L=float(np.asarray(d['L_m']).item());h=float(np.asarray(d['spacing_m']).item())
 if not np.isclose(L/h,n,rtol=1e-12):raise ValueError('CT physical scale inconsistent')
 if not np.isclose(h,{'S1':8.683e-6,'S5':3.997e-6}[sample],rtol=1e-12):raise ValueError('Native voxel size changed')
 v=np.asarray(d['vertices_m'],float).T.copy()/L;t=np.asarray(d['tets0'],np.int32).T.copy()
 rawcode=np.asarray(d['tet_entity']).ravel().astype(np.int64)
 if not np.isin(rawcode,np.arange(1,13)).all():raise ValueError('Unexpected original entity code; do not silently remap')
 code=((rawcode-1)%6+1).astype(np.uint8)
 del rawcode
 if t.min()<0 or t.max()>=len(v):raise ValueError('Source is not 0-based tetrahedra')
 del d
 orig_vol=np.abs(tetvol(v,t))*L**3;expected=np.bincount(code,weights=orig_vol,minlength=7)[1:7]
 p={k:p0[k] for k in ['sample','source','source_raw','source_header','source_start_index_xyz','source_end_exclusive_xyz','shape_xyz'] if k in p0}
 if not all(k in p for k in ['source','source_raw']):raise ValueError('Missing traceable original CT source metadata')
 p.update(schema='r7b_localfit_new_M1_v1',sample=sample,mode='formal' if n==128 else 'technical32',shape_xyz=[n]*3,
          spacing_m=h,L_m=L,total_volume_m3=L**3,background_pore_volume_m3=float(expected[[0,2,4]].sum()),source_mesh_file=str(src/'mesh_input.mat'),source_parameters_file=str(src/'parameters.json'),
          source_read_only=True,protection_regions_present=n==128,technical_pilot_only=n!=128,formal_CT_target=n==128,
          potential_order=1,skeleton_sigma_S_m=1e-5,fluid_and_fracture_sigmas_S_m=[.18,10.],input_variant='M1',
          mesh_tag='mesh1',geometry_tag='geom1',ellipsoid_center_m=[L*.5]*3,ellipsoid_semiaxes_m=[L*.25,L*.25,L*.0125],
          basis_u=[.5,0.,3**.5/2],basis_v=[0.,1.,0.],basis_n=[-3**.5/2,0.,.5],
          ellipsoid_representation='Native Ellipsoid analytic target; adaptively fitted piecewise-planar zero-distance interface; explicit conforming cut tetrahedra, NOT centroid tagging',
          new_surface_version_not_old_CAD_tessellation=True,geometry_preprocess_distance_tolerance_over_L=tol,
          node_levelset_snap_over_L=tol*.1,refinement_preserves_frozen_surface=False,
          all_original_phase_labels_preserved=True,CT_resampling='none')
 (dst/'source_parameters_read_only.json').write_text(json.dumps(p0,ensure_ascii=False,indent=2),encoding='utf-8')
 savemat(dst/'original_labels_and_parameters.mat',{'labels_xyz':labels,'L_m':L,'spacing_m':h,'p_json':json.dumps(p,ensure_ascii=False),
          'expected_phase_region_volume_m3':expected},do_compression=False)
 del orig_vol
 v,t,code,parent,surface,report=fit(v,t,code,tolerance=tol,max_rounds=11,max_tets=max_tets,verbose=True)
 del parent
 # Actual generated geometry gets the independent, existing thickness tests too.
 from r7b_geometry_check import normal_line_probes,BASIS
 probes=normal_line_probes((v-.5)@BASIS,surface,.25,.25,.0125,1.)
 if not all(x['pass_one_interval'] and x['pass_thickness_one_percent'] for x in probes):
  raise RuntimeError('Original 65-ray thickness requirement NOT met; no formal mesh output accepted')
 report.update(actual_probe_count=65,actual_max_probe_relative_thickness_error=max(abs(x['actual_relative_thickness_error']) for x in probes),
               actual_center_relative_thickness_error=probes[0]['actual_relative_thickness_error'],
               input_is_actual_user_CT_mesh=True,COMSOL_PDE_executed_here=False,
               remaining_validation='Actual COMSOL import, original CT SAT, original geometry checker, and target electric solves remain required')
 v*=L
 np.asarray(v,dtype='<f8',order='C').tofile(dst/'vertices_m.bin')
 np.asarray(t,dtype='<i4',order='C').tofile(dst/'tets0.bin')
 np.asarray(code,dtype='u1').tofile(dst/'entity_code.bin')
 # Compact interface export for an independent check and visual inspection.
 inside=t[code>6];used,ii=np.unique(inside,return_inverse=True)
 savemat(dst/'preimport_actual_boundary.mat',{'vertices_m':v[used], 'tetrahedra_inside':ii.reshape(-1,4).astype(np.int32),
          'original_phase_inside':((code[code>6]-1)%2).astype(np.uint8),'index_base':0},do_compression=False)
 sizes={name:(dst/name).stat().st_size for name in ['vertices_m.bin','tets0.bin','entity_code.bin','preimport_actual_boundary.mat','original_labels_and_parameters.mat']}
 p['fixed_surface_identifier']='localfit_'+sample+'_'+scope+'_M1_'+datetime.now().strftime('%Y%m%d_%H%M%S');p['vertices']=len(v);p['tetrahedra']=len(t)
 (dst/'parameters.json').write_text(json.dumps(p,ensure_ascii=False,indent=2),encoding='utf-8')
 (dst/'preprocess_geometry_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 rec=dict(status='PREPROCESS_ONLY_NOT_COMSOL_ACCEPTANCE',sample=sample,scope=scope,vertices=len(v),tetrahedra=len(t),L_m=L,
          vertices_format='little-endian float64, interleaved x y z, zero-based rows',
          tetrahedra_format='little-endian int32, four 0-based vertex IDs per element',
          entity_format='uint8: 1+original_phase+2*region+6*inside_explicit_cut_ellipsoid',file_sizes_bytes=sizes,
          fixed_surface_identifier=p['fixed_surface_identifier'],checksums_used=False,
          seconds=time.perf_counter()-start,formal_COMSOL_completed=False)
 done.write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps(rec,ensure_ascii=False,indent=2))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--sample',choices=['S1','S5'],default='S1')
 ap.add_argument('--scope',choices=['formal','technical32'],default='formal');ap.add_argument('--tolerance',type=float,default=5e-5)
 ap.add_argument('--max-tets',type=int,default=30_000_000);a=ap.parse_args()
 if a.scope=='technical32' and a.sample!='S1':ap.error('only the recorded S1 technical32 source exists')
 if not 0<a.tolerance<=5e-5:ap.error('do not relax the geometry tolerance above 5e-5 L')
 run(a.root,a.sample,a.scope,a.tolerance,a.max_tets)
