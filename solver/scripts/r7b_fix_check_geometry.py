#!/usr/bin/env python3
"""Run the original independent checker on COMSOL-EXPORTED tets, then accept.
This script never accepts a preimport-only mesh as a completed COMSOL geometry.
"""
import argparse,json,subprocess,sys
from pathlib import Path
from r7b_geometry_check import load_input,verify_geometry

def main():
 ap=argparse.ArgumentParser();ap.add_argument('folder',type=Path);ap.add_argument('--technical',action='store_true');a=ap.parse_args()
 folder=a.folder.resolve();p=json.loads((folder/'parameters.json').read_text(encoding='utf-8-sig'))
 for name in ['geometry_ready.mph','actual_boundary.mat','native_phase_verification.json','selections.json']:
  if not (folder/name).is_file():raise SystemExit(f'Missing COMSOL-generated artifact: {name}')
 selections=json.loads((folder/'selections.json').read_text(encoding='utf-8-sig'))
 if selections.get('min_mesh_quality',0)<=0:raise SystemExit('Nonpositive actual mesh quality: no geometry acceptance or PDE allowed')
 if not a.technical and p.get('shape_xyz')!=[128,128,128]:raise SystemExit('Technical32 cannot be accepted as formal CT')
 r=verify_geometry(**load_input(folder/'actual_boundary.mat'),sample=p['sample'],output_dir=folder,
                  make_sections=True,parameters=p,metadata={'COMSOL_export':str(folder/'actual_boundary.mat'),
                  'input_is_preimport_mesh':False,'formal_CT_solved':False})
 if not r['all_evaluated_geometric_gates_pass']:raise SystemExit('Actual geometry failed original gates; no SUCCESS written')
 if a.technical:
  print('TECHNICAL IMPORT GEOMETRY PASSED ONLY; no formal GEOMETRY_SUCCESS or CT target');return
 subprocess.run([sys.executable,str(Path(__file__).with_name('r7b_accept_geometry.py')),str(folder)],check=True)
if __name__=='__main__':main()
