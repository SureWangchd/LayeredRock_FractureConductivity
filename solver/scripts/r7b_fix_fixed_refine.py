#!/usr/bin/env python3
"""Fixed-interface midpoint h refinement of the accepted local-fit CT mesh.

No ellipsoid distance, projection, material classification, smoothing,
COMSOL mesher, or solution field is used. All new vertices are exact edge
midpoints of the preceding frozen tetrahedral mesh. Chosen original edges
are globally sorted and every incident tetrahedron processes them in that
same order, including unmarked neighbours, to preserve complete shared faces.
"""
from __future__ import annotations
import argparse
from datetime import datetime
import json
from pathlib import Path
import shutil
import time

import numpy as np
from numba import njit

EDGES = np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]], dtype=np.int64)
MAX_TETS = 50_000_000
MAX_VERTICES = 9_000_000
CHUNK = 100_000


@njit(cache=True)
def signed_volume(v, t):
    ax=v[t[1],0]-v[t[0],0];ay=v[t[1],1]-v[t[0],1];az=v[t[1],2]-v[t[0],2]
    bx=v[t[2],0]-v[t[0],0];by=v[t[2],1]-v[t[0],1];bz=v[t[2],2]-v[t[0],2]
    cx=v[t[3],0]-v[t[0],0];cy=v[t[3],1]-v[t[0],1];cz=v[t[3],2]-v[t[0],2]
    return (ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)+az*(bx*cy-by*cx))/6.


@njit(cache=True)
def select_longest_edges(v,t,code,box,with_nearfield):
    nv=len(v);keys=np.empty(len(t),np.int64);number=0;ellipse=0;nearfield=0
    for i in range(len(t)):
        selected=code[i]>6
        if selected:ellipse+=1
        elif with_nearfield:
            selected=True
            for j in range(4):
                for a in range(3):
                    q=v[t[i,j],a]
                    if q<box[a,0] or q>box[a,1]:selected=False
            if selected:nearfield+=1
        if not selected:continue
        best=-1.;keybest=np.int64(9223372036854775807)
        for j in range(6):
            a=np.int64(t[i,EDGES[j,0]]);b=np.int64(t[i,EDGES[j,1]])
            if a>b:a,b=b,a
            dx=v[a,0]-v[b,0];dy=v[a,1]-v[b,1];dz=v[a,2]-v[b,2]
            d=dx*dx+dy*dy+dz*dz;key=a*nv+b
            if d>best or (d==best and key<keybest):best=d;keybest=key
        keys[number]=keybest;number+=1
    return keys[:number],ellipse,nearfield


@njit(cache=True)
def split_one(tet, chosen, nv):
    """Split only chosen ORIGINAL edges, in their common global key order."""
    keys=np.empty(6,np.int64);aa=np.empty(6,np.int64);bb=np.empty(6,np.int64);mm=np.empty(6,np.int64);nk=0
    for e in range(6):
        a=np.int64(tet[EDGES[e,0]]);b=np.int64(tet[EDGES[e,1]])
        if a>b:a,b=b,a
        key=a*nv+b;at=np.searchsorted(chosen,key)
        if at<len(chosen) and chosen[at]==key:
            pos=nk
            while pos>0 and keys[pos-1]>key:
                keys[pos]=keys[pos-1];aa[pos]=aa[pos-1];bb[pos]=bb[pos-1];mm[pos]=mm[pos-1];pos-=1
            keys[pos]=key;aa[pos]=a;bb[pos]=b;mm[pos]=nv+at;nk+=1
    children=np.empty((8,4),np.int32);children[0]=tet;n=1
    for e in range(nk):
        before=n
        for j in range(before):
            ia=-1;ib=-1
            for p in range(4):
                if children[j,p]==aa[e]:ia=p
                if children[j,p]==bb[e]:ib=p
            if ia>=0 and ib>=0:
                if n>=8:raise ValueError('Original-edge bisection exceeded eight children')
                children[n]=children[j]
                children[j,ib]=mm[e];children[n,ia]=mm[e];n+=1
    return children,n


@njit(cache=True)
def count_children(t, chosen, nv):
    counts=np.empty(len(t),np.uint8)
    for i in range(len(t)):
        # Most cells are untouched; avoid allocating a child list for them.
        hit=False
        for e in range(6):
            a=np.int64(t[i,EDGES[e,0]]);b=np.int64(t[i,EDGES[e,1]])
            if a>b:a,b=b,a
            key=a*nv+b;at=np.searchsorted(chosen,key)
            if at<len(chosen) and chosen[at]==key:hit=True;break
        if hit:
            _,n=split_one(t[i],chosen,nv);counts[i]=n
        else:counts[i]=1
    return counts


@njit(cache=True)
def write_children(v,t,code,root,chosen,nv,counts):
    total=np.int64(0)
    for n in counts:total+=n
    out=np.empty((total,4),np.int32);codes=np.empty(total,np.uint8);roots=np.empty(total,np.int32)
    source_vol=np.zeros(13);child_vol=np.zeros(13);minimum=1e300;maxerr=0.;negative_reordered=0;offset=0
    for i in range(len(t)):
        oldvol=abs(signed_volume(v,t[i]))
        if not np.isfinite(oldvol) or oldvol<=0:raise ValueError('Nonpositive or nonfinite source tetrahedron; no deletion permitted')
        if counts[i]==1:
            children=np.empty((8,4),np.int32);children[0]=t[i];n=1
        else:children,n=split_one(t[i],chosen,nv)
        if n!=counts[i]:raise ValueError('Dry count and actual subdivision disagree')
        vsum=0.
        for j in range(n):
            vol=signed_volume(v,children[j])
            if not np.isfinite(vol) or vol==0:raise ValueError('Degenerate midpoint child; not deleted')
            if vol<0:
                temporary=children[j,0];children[j,0]=children[j,1];children[j,1]=temporary
                vol=-vol;negative_reordered+=1
            minimum=min(minimum,vol);vsum+=vol
            out[offset+j]=children[j];codes[offset+j]=code[i];roots[offset+j]=root[i]
        err=abs(vsum/oldvol-1.)
        if err>2e-9:raise ValueError('Midpoint children fail original parent-volume partition')
        maxerr=max(maxerr,err);source_vol[code[i]]+=oldvol;child_vol[code[i]]+=vsum;offset+=n
    return out,codes,roots,minimum,maxerr,negative_reordered,source_vol,child_vol


def dump(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8')


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def freeze_design(out,padding):
    path=out/'fixed_midpoint_refinement_design.json'
    design=dict(schema='r7b_fixed_midpoint_h_refinement_v1',samples=['S1','S5'],levels=['M2','M3'],
                padding_native_voxels=padding,nearfield_enabled=padding>0,
                selected_cell_rule='All inherited ellipsoid cells, plus optional cells wholly in a fixed ellipsoid bounding box padded by declared native voxels',
                selected_edge_rule='One longest ORIGINAL edge per selected cell; exact length ties use smallest canonical global edge key',
                subdivision_rule='Globally sort unique original edges; every incident cell splits surviving children on those edges in that same global order',
                propagation_rule='Unselected neighbours sharing chosen original edges are subdivided to maintain shared complete faces',
                new_vertex_rule='Arithmetic midpoint of the preceding actual mesh edge only; original vertices unchanged',
                geometry_or_material_refit=False,solution_fields_read=False,source_phase_labels_reclassified=False,
                max_tetrahedra=MAX_TETS,max_vertices=MAX_VERTICES,
                parent_volume_relative_tolerance=2e-9,
                supersedes='Failed native COMSOL Refine plan; native refinement changed frozen overlap volumes and produced negative element quality before any M2 target was solved',
                scientific_scope='Local h-refinement around the finite ellipsoid and conforming neighbours; not a claim of global CT continuum-limit convergence')
    if path.exists():
        saved=read(path)
        for key,value in design.items():
            if saved.get(key)!=value:raise ValueError(f'Fixed midpoint rule already frozen; mismatch in {key}')
        return saved
    design['frozen_local']=datetime.now().isoformat(timespec='seconds')
    dump(path,design);return design


def fixed_box(p,padding):
    basis=np.column_stack([p['basis_u'],p['basis_v'],p['basis_n']])
    ext=np.sqrt((basis*basis)@(np.asarray(p['ellipsoid_semiaxes_m'])**2))+padding*p['spacing_m']
    center=np.asarray(p['ellipsoid_center_m']);length=p['L_m']
    return np.column_stack([np.maximum(0,center-ext),np.minimum(length,center+ext)])


def fill_midpoint_vertices(source, destination, chosen, chunk=CHUNK):
    """Copy complete source blocks, then append only exact edge midpoints."""
    nv=len(source)
    if len(destination)!=nv+len(chosen):raise ValueError('Unexpected output vertex count')
    for first in range(0,nv,chunk):
        last=min(nv,first+chunk)
        destination[first:last]=source[first:last]
    for first in range(0,len(chosen),chunk):
        edges=chosen[first:first+chunk]
        destination[nv+first:nv+first+len(edges)]=.5*(source[edges//nv]+source[edges%nv])


def prepare(root,sample,level,padding=0,dry_count=False):
    root=Path(root).resolve();out=root/'runs'/'round7b_localfit';out.mkdir(parents=True,exist_ok=True)
    if sample not in ('S1','S5') or level not in ('M2','M3') or padding not in (0,2):raise ValueError('Only prescribed samples/levels and padding 0 or 2 are supported')
    design=freeze_design(out,padding)
    previous='M1' if level=='M2' else 'M2';src=out/'inputs'/f'{sample}_{previous}';dst=out/'inputs'/f'{sample}_{level}'
    p=read(src/'parameters.json');rec=read(src/'PREPROCESS_COMPLETE.json')
    base=read(out/'inputs'/f'{sample}_M1'/'parameters.json')
    if p['sample']!=sample or p['shape_xyz']!=[128]*3 or not p['protection_regions_present']:raise ValueError('Original formal CT identity/protection mismatch')
    if p['fixed_surface_identifier']!=base['fixed_surface_identifier']:raise ValueError('Frozen M1 surface identity changed')
    if previous=='M2' and not p.get('refinement_preserves_frozen_surface'):raise ValueError('M3 input is not a frozen-interface M2 mesh')
    if not dry_count:
        gate=read(src/'GEOMETRY_SUCCESS.json')
        if not gate.get('geometry_passed'):raise ValueError('Previous mesh must pass actual COMSOL geometry before preparing the next final level')
    if (dst/'PREPROCESS_COMPLETE.json').exists() and not dry_count:
        old=read(dst/'PREPROCESS_COMPLETE.json')
        if old.get('fixed_surface_identifier')!=base['fixed_surface_identifier'] or old.get('fixed_midpoint_refinement') is not True:raise ValueError('Different preprocessing already exists')
        for name,size in old['file_sizes_bytes'].items():
            if (dst/name).stat().st_size!=size:raise ValueError(f'Frozen output missing/resized: {name}')
        print(json.dumps(dict(status='REUSE_COMPLETE_PREPROCESSING',folder=str(dst))),flush=True);return old
    if not dry_count and any((dst/name).exists() for name in ('vertices_m.bin','tets0.bin','entity_code.bin','parent_M1_tet0.bin')):
        raise RuntimeError('Partial new-level output exists; preserve it before retry, never silently overwrite')
    nv=int(rec['vertices']);nt=int(rec['tetrahedra']);L=float(p['L_m'])
    expected={'vertices_m.bin':24*nv,'tets0.bin':16*nt,'entity_code.bin':nt}
    for name,size in expected.items():
        if (src/name).stat().st_size!=size:raise ValueError(f'Source byte count mismatch: {name}')
    if nv>MAX_VERTICES or nt>MAX_TETS:raise MemoryError('Input already exceeds fixed resource limits')
    v=np.memmap(src/'vertices_m.bin',mode='r',dtype='<f8',shape=(nv,3))
    t=np.memmap(src/'tets0.bin',mode='r',dtype='<i4',shape=(nt,4))
    code=np.memmap(src/'entity_code.bin',mode='r',dtype='u1',shape=(nt,))
    rootfile=src/'parent_M1_tet0.bin';parents=None
    if rootfile.exists():
        if rootfile.stat().st_size!=4*nt:raise ValueError('Inherited M1-parent array size mismatch')
        parents=np.memmap(rootfile,mode='r',dtype='<i4',shape=(nt,))
    elif previous!='M1':raise ValueError('M2 is missing the inherited M1-parent map')
    started=time.perf_counter();box=fixed_box(p,padding);parts=[];nellipse=0;nnear=0
    pre=out/'fixed_refinement_preflight'/f'{sample}_{level}'/datetime.now().strftime('%Y%m%d_%H%M%S_%f');pre.mkdir(parents=True)
    try:
        for first in range(0,nt,CHUNK):
            last=min(nt,first+CHUNK)
            if np.min(t[first:last])<0 or np.max(t[first:last])>=nv:raise ValueError('Invalid source vertex IDs')
            if not np.all((code[first:last]>=1)&(code[first:last]<=12)):raise ValueError('Invalid inherited material code')
            keys,ne,nn=select_longest_edges(v,t[first:last],code[first:last],box,padding>0)
            if len(keys):parts.append(keys)
            nellipse+=int(ne);nnear+=int(nn)
        if not parts:raise ValueError('No actual inherited ellipsoid cells selected')
        chosen=np.unique(np.concatenate(parts));del parts
        nvout=nv+len(chosen)
        np.save(pre/'chosen_original_edges.npy',chosen)
        if nvout>MAX_VERTICES:raise MemoryError(f'Midpoint vertex count {nvout} exceeds {MAX_VERTICES}')
        counts=np.empty(nt,np.uint8);nout=0;affected=0;splitellipse=0
        for first in range(0,nt,CHUNK):
            last=min(nt,first+CHUNK);cnt=count_children(t[first:last],chosen,nv);counts[first:last]=cnt
            nout+=int(cnt.sum());affected+=int(np.count_nonzero(cnt>1));splitellipse+=int(np.count_nonzero((code[first:last]>6)&(cnt>1)))
            if first%(CHUNK*10)==0 or last==nt:
                print(json.dumps(dict(stage='dry_count',sample=sample,level=level,tetrahedra_read=last,total_source_tetrahedra=nt,predicted_children_so_far=nout,seconds=time.perf_counter()-started)),flush=True)
        np.save(pre/'child_counts.npy',counts)
        if splitellipse!=nellipse:raise ValueError('Not every actual ellipsoid cell will be split')
        preflight=dict(status='DRY_COUNT_COMPLETE_NOT_COMSOL_ACCEPTANCE',sample=sample,level=level,previous_level=previous,
                       input_folder=str(src),source_vertices=nv,source_tetrahedra=nt,predicted_vertices=nvout,predicted_tetrahedra=nout,
                       selected_ellipsoid_tetrahedra=nellipse,selected_additional_nearfield_tetrahedra=nnear,
                       all_ellipsoid_tetrahedra_split=True,affected_source_tetrahedra=affected,chosen_original_edges=len(chosen),
                       padding_native_voxels=padding,nearfield_box_m=box.tolist(),nearfield_enabled=padding>0,
                       resource_limits_passed=nout<=MAX_TETS and nvout<=MAX_VERTICES,max_tetrahedra=MAX_TETS,max_vertices=MAX_VERTICES,
                       fixed_surface_identifier=p['fixed_surface_identifier'],no_solution_fields_read=True,
                       previous_geometry_accepted=(src/'GEOMETRY_SUCCESS.json').exists(),prospective_only=dry_count,
                       preflight_folder=str(pre),seconds=time.perf_counter()-started)
        dump(pre/'preflight.json',preflight)
        print(json.dumps(preflight),flush=True)
        if nout>MAX_TETS:raise MemoryError(f'Actual dry count {nout} exceeds {MAX_TETS}; no output mesh written')
        if dry_count:return preflight
        dst.mkdir(parents=True,exist_ok=True);dump(dst/'FIXED_REFINEMENT_STARTED.json',preflight)
        # Only these new output files are writable. Source vertices and
        # source tetrahedra remain read-only even when orientation is fixed.
        vv=np.memmap(dst/'vertices_m.bin',mode='w+',dtype='<f8',shape=(nvout,3))
        fill_midpoint_vertices(v,vv,chosen)
        vv.flush();minimum=float('inf');maxerr=0.;reordered=0;source_vol=np.zeros(13);child_vol=np.zeros(13)
        try:
            with (dst/'tets0.bin').open('xb') as ft,(dst/'entity_code.bin').open('xb') as fc,(dst/'parent_M1_tet0.bin').open('xb') as fr:
                written=0
                for first in range(0,nt,CHUNK):
                    last=min(nt,first+CHUNK)
                    rootids=np.arange(first,last,dtype=np.int32) if parents is None else parents[first:last]
                    tt,cc,rr,mn,er,nneg,sv,cv=write_children(vv,t[first:last],code[first:last],rootids,chosen,nv,counts[first:last])
                    tt.tofile(ft);cc.tofile(fc);rr.tofile(fr);written+=len(tt)
                    minimum=min(minimum,float(mn));maxerr=max(maxerr,float(er));reordered+=int(nneg);source_vol+=sv;child_vol+=cv
                    if first%(CHUNK*10)==0 or last==nt:
                        print(json.dumps(dict(stage='write_children',sample=sample,level=level,source_cells_done=last,child_cells_written=written,seconds=time.perf_counter()-started)),flush=True)
            if written!=nout:raise ValueError('Written child count differs from dry prediction')
            if np.max(np.abs(child_vol-source_vol))/L**3>1e-10:raise ValueError('Fixed material/ellipsoid overlap volume changed')
        finally:
            vv.flush();vv._mmap.close()
        shutil.copy2(src/'original_labels_and_parameters.mat',dst/'original_labels_and_parameters.mat')
        pnew=dict(p);pnew.update(input_variant=level,vertices=nvout,tetrahedra=nout,
            refinement_preserves_frozen_surface=True,refinement_geometry_rebuilt=False,original_CT_resampled=False,
            fixed_surface_source=str(out/'inputs'/f'{sample}_M1'/'geometry_ready.mph'),
            refinement_preceding_input=str(src/'geometry_ready.mph'),refinement_preceding_binary_folder=str(src),
            refinement_plan_json=str(out/'fixed_midpoint_refinement_design.json'),refinement_facerep='exact inherited planar facets; no COMSOL Refine',
            ellipsoid_representation='Same frozen localfit M1 piecewise-planar ellipsoid; conforming midpoint subdivision only, no surface projection/refit or material reclassification',
            fixed_midpoint_refinement=True,refinement_root_parent_reference=str(out/'inputs'/f'{sample}_M1'/'tets0.bin'))
        dump(dst/'parameters.json',pnew)
        report=dict(preflight,status='REFINED_PREIMPORT_NOT_COMSOL_ACCEPTANCE',fixed_midpoint_refinement=True,
                    geometry_not_rebuilt=True,source_vertices_unchanged=True,new_vertices_are_only_original_edge_midpoints=True,
                    original_CT_labels_inherited_without_reclassification=True,zero_volume_tetrahedra=0,
                    all_output_signed_volumes_positive=True,minimum_tetrahedron_volume_m3=minimum,
                    maximum_parent_partition_relative_error=maxerr,orientation_reordered_children=reordered,
                    orientation_reordering_changes_connectivity_order_only=True,no_bad_cells_removed=True,
                    actual_phase_region_ellipse_volumes_source_m3=source_vol[1:].tolist(),
                    actual_phase_region_ellipse_volumes_output_m3=child_vol[1:].tolist(),
                    maximum_phase_region_ellipse_volume_change_over_L3=float(np.max(abs(child_vol-source_vol))/L**3),
                    all_volume_faces_check_pending=True,actual_COMSOL_geometry_check_pending=True,PDE_solved=False,
                    total_seconds=time.perf_counter()-started,rule=design)
        dump(dst/'fixed_geometry_refinement.json',report)
        names=['vertices_m.bin','tets0.bin','entity_code.bin','parent_M1_tet0.bin','original_labels_and_parameters.mat']
        complete=dict(status='PREPROCESS_ONLY_NOT_COMSOL_ACCEPTANCE',sample=sample,scope='formal',input_variant=level,
                      vertices=nvout,tetrahedra=nout,L_m=L,file_sizes_bytes={name:(dst/name).stat().st_size for name in names},
                      fixed_surface_identifier=p['fixed_surface_identifier'],fixed_midpoint_refinement=True,
                      refinement_preserves_frozen_surface=True,checksum_used=False,formal_COMSOL_completed=False,
                      vertices_format='little-endian float64 interleaved xyz',tetrahedra_format='little-endian int32 four 0-based IDs',
                      entity_format='uint8 inherited 1+phase+2*region+6*inside; unchanged classification',
                      parent_format='int32 zero-based original frozen M1 tetrahedron ID',seconds=time.perf_counter()-started)
        dump(dst/'PREPROCESS_COMPLETE.json',complete);print(json.dumps(complete),flush=True);return complete
    except Exception as exc:
        dump(pre/'FAILED.json',dict(status='FAILED',sample=sample,level=level,error=repr(exc),no_source_modified=True,no_bad_cells_deleted=True,seconds=time.perf_counter()-started))
        raise
    finally:
        v._mmap.close();t._mmap.close();code._mmap.close()
        if parents is not None:parents._mmap.close()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--sample',choices=['S1','S5'],required=True);ap.add_argument('--level',choices=['M2','M3'],required=True)
    ap.add_argument('--padding-native-voxels',type=int,choices=[0,2],default=0)
    ap.add_argument('--dry-count',action='store_true');args=ap.parse_args()
    prepare(args.root,args.sample,args.level,args.padding_native_voxels,args.dry_count)


if __name__=='__main__':main()
