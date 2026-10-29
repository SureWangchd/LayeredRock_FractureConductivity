#!/usr/bin/env python3
"""Check every face of a prepared local-fit volume mesh, without COMSOL.

This checks the NEW mesh only. Canonical triangles are partitioned by their
smallest vertex ID modulo a bucket count; all copies of a face enter the
same bucket. Each bucket is sorted independently, bounding resident memory.
The source vertices/connectivity are opened read-only and never modified.
No checksum, cell deletion, geometry recolouring, or mesh repair is used.

Usage (owned COMSOL service should be stopped during a large check):
  python scripts/r7b_fix_faces_check.py runs/round7b_localfit/inputs/S1_M1
  ... --workers 8 --buckets 128 --chunk-tets 200000

An interrupted check can resume with --work-dir <printed work directory>.
The identity check uses source paths, sizes and modification times only.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

FACE_LOCAL = np.asarray([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])
RECORD_DTYPE = np.dtype([('a', '<u4'), ('b', '<u4'), ('c', '<u4'), ('owner', '<u4')])
BOUNDARY_TOLERANCE_OVER_L = 1e-12
MAX_FAILURE_EXAMPLES = 100


def write_json(path: Path, value: dict):
    """Commit one small result atomically; never alter input binaries."""
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
    os.replace(temp, path)


def read_json(path: Path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def source_spec(folder: Path) -> dict:
    manifest_path = folder / 'PREPROCESS_COMPLETE.json'
    rec = read_json(manifest_path)
    nv, nt, length = int(rec['vertices']), int(rec['tetrahedra']), float(rec['L_m'])
    if not (0 < nv < 2**31 and 0 < nt < 2**32 and np.isfinite(length) and length > 0):
        raise ValueError('Invalid vertex/tetrahedron counts or physical length')
    expected = {'vertices_m.bin': nv*24, 'tets0.bin': nt*16}
    sources = {}
    for name, size in expected.items():
        path = folder / name
        stat = path.stat()
        if stat.st_size != size:
            raise ValueError(f'Input size differs from manifest: {name}')
        sources[name] = dict(path=str(path.resolve()), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
    return dict(folder=str(folder.resolve()), vertices=nv, tetrahedra=nt, L_m=length,
                sources=sources, checksum_used=False)


def make_buckets(spec: dict, work: Path, bucket_count: int, chunk_tets: int) -> dict:
    """Stream every tet exactly once; each writes four face-owner records."""
    marker = work / 'PARTITION_COMPLETE.json'
    if marker.exists():
        done = read_json(marker)
        if done['source'] != spec or done['buckets'] != bucket_count:
            raise ValueError('Existing partition belongs to a different input or bucket rule')
        for b, count in enumerate(done['bucket_face_records']):
            if (work / f'bucket_{b:04d}.bin').stat().st_size != count*16:
                raise ValueError('Incomplete or resized partition file')
        return done
    if any(work.glob('bucket_*.bin')):
        raise RuntimeError('Interrupted partial partition is preserved. Use a new work directory; no files were overwritten.')
    nt, nv = spec['tetrahedra'], spec['vertices']
    tets = np.memmap(spec['sources']['tets0.bin']['path'], mode='r', dtype='<i4', shape=(nt, 4))
    handles = [(work / f'bucket_{b:04d}.bin').open('xb') for b in range(bucket_count)]
    counts = np.zeros(bucket_count, dtype=np.int64)
    start = time.perf_counter()
    next_progress = start
    try:
        for first in range(0, nt, chunk_tets):
            last = min(nt, first+chunk_tets)
            cells = np.asarray(tets[first:last])
            if cells.min() < 0 or cells.max() >= nv:
                raise ValueError(f'Invalid source vertex ID in tetrahedra [{first}, {last})')
            sorted_cells = np.sort(cells, axis=1)
            repeated = np.any(sorted_cells[:, 1:] == sorted_cells[:, :-1], axis=1)
            if np.any(repeated):
                bad = first+np.flatnonzero(repeated)
                write_json(work / 'REPEATED_TET_VERTEX_FAILED.json', dict(
                    tetrahedron_count=int(np.count_nonzero(repeated)),
                    tetrahedron_indices_0based=bad[:MAX_FAILURE_EXAMPLES].tolist(),
                    source_tetrahedra=tets[bad[:MAX_FAILURE_EXAMPLES]].tolist(),
                    no_cells_removed=True))
                raise ValueError('A source tetrahedron repeats a vertex; no cells removed')
            records = np.empty((last-first)*4, dtype=RECORD_DTYPE)
            faces = np.sort(cells[:, FACE_LOCAL].reshape(-1, 3), axis=1)
            records['a'], records['b'], records['c'] = faces.T
            records['owner'] = np.repeat(np.arange(first, last, dtype=np.uint32), 4)
            bins = np.asarray(records['a'] % bucket_count, dtype=np.intp)
            order = np.argsort(bins, kind='stable')
            records = records[order]
            sizes = np.bincount(bins, minlength=bucket_count)
            offsets = np.r_[0, np.cumsum(sizes)]
            for b in np.flatnonzero(sizes):
                records[offsets[b]:offsets[b+1]].tofile(handles[b])
            counts += sizes
            now = time.perf_counter()
            if now >= next_progress or last == nt:
                print(json.dumps(dict(stage='partition', tetrahedra_done=last, tetrahedra_total=nt,
                                      face_records_done=int(counts.sum()), seconds=now-start)), flush=True)
                next_progress = now+15
    finally:
        for handle in handles:
            handle.close()
        # NumPy views can keep this read-only mapping alive in an exception
        # traceback. Explicitly release the Windows file handle as well.
        tets._mmap.close()
        del tets
    assert int(counts.sum()) == nt*4
    done = dict(source=spec, buckets=bucket_count, bucket_rule='minimum canonical vertex ID modulo bucket count',
                bucket_face_records=counts.tolist(), total_face_records=int(counts.sum()),
                maximum_bucket_face_records=int(counts.max()), seconds=time.perf_counter()-start)
    write_json(marker, done)
    return done


def failure_example(records, start, stop, vertices, reason):
    face = np.array([records['a'][start], records['b'][start], records['c'][start]], dtype=np.int64)
    owners = records['owner'][start:stop]
    return dict(reason=reason, face_vertex_ids_0based=face.tolist(), vertices_m=vertices[face].tolist(),
                centroid_m=np.mean(vertices[face], axis=0).tolist(), owner_count=int(stop-start),
                owner_tetrahedron_ids_0based=owners[:MAX_FAILURE_EXAMPLES].astype(np.int64).tolist(),
                owner_list_truncated=bool(len(owners)>MAX_FAILURE_EXAMPLES))


def check_bucket(b: int, spec: dict, work: Path) -> dict:
    result_path = work / f'bucket_{b:04d}.json'
    if result_path.exists():
        return read_json(result_path)
    started = time.perf_counter()
    records = np.fromfile(work / f'bucket_{b:04d}.bin', dtype=RECORD_DTYPE)
    nv, length = spec['vertices'], spec['L_m']
    vertices = np.memmap(spec['sources']['vertices_m.bin']['path'], mode='r', dtype='<f8', shape=(nv, 3))
    result = dict(bucket=b, face_records=int(len(records)), unique_faces=0, single_owner_faces=0,
                  two_owner_faces=0, more_than_two_owner_faces=0, repeated_owner_faces=0,
                  internal_unmatched_faces=0, maximum_owner_count=0, external_boundary_faces_by_plane=[0]*6,
                  external_boundary_area_m2_by_plane=[0.0]*6, failure_examples=[])
    if len(records):
        # All equal face IDs are adjacent regardless of owner or orientation.
        records.sort(order=['a', 'b', 'c', 'owner'], kind='quicksort')
        same = ((records['a'][1:] == records['a'][:-1]) &
                (records['b'][1:] == records['b'][:-1]) &
                (records['c'][1:] == records['c'][:-1]))
        starts = np.r_[0, np.flatnonzero(~same)+1]
        ends = np.r_[starts[1:], len(records)]
        multiplicities = ends-starts
        repeated_pairs = np.flatnonzero(same & (records['owner'][1:] == records['owner'][:-1]))
        repeated_groups = np.unique(np.searchsorted(starts, repeated_pairs, side='right')-1)
        result.update(unique_faces=int(len(starts)), single_owner_faces=int(np.count_nonzero(multiplicities == 1)),
                      two_owner_faces=int(np.count_nonzero(multiplicities == 2)),
                      more_than_two_owner_faces=int(np.count_nonzero(multiplicities > 2)),
                      repeated_owner_faces=int(len(repeated_groups)), maximum_owner_count=int(multiplicities.max()))
        examples = result['failure_examples']
        for group in np.flatnonzero(multiplicities > 2)[:MAX_FAILURE_EXAMPLES]:
            examples.append(failure_example(records, starts[group], ends[group], vertices, 'more_than_two_owners'))
        for group in repeated_groups[:max(0, MAX_FAILURE_EXAMPLES-len(examples))]:
            examples.append(failure_example(records, starts[group], ends[group], vertices, 'same_owner_repeated'))
        single_indices = starts[multiplicities == 1]
        counts = np.zeros(6, dtype=np.int64)
        areas = np.zeros(6, dtype=np.float64)
        tol = length*BOUNDARY_TOLERANCE_OVER_L
        for first in range(0, len(single_indices), 100000):
            indices = single_indices[first:first+100000]
            ids = np.column_stack([records['a'][indices], records['b'][indices], records['c'][indices]])
            q = vertices[ids]
            lower = np.all(np.abs(q) <= tol, axis=1)
            upper = np.all(np.abs(q-length) <= tol, axis=1)
            planes = np.column_stack([lower[:, 0], upper[:, 0], lower[:, 1], upper[:, 1], lower[:, 2], upper[:, 2]])
            valid = np.any(planes, axis=1) & np.all(np.isfinite(q), axis=(1, 2))
            result['internal_unmatched_faces'] += int(np.count_nonzero(~valid))
            cross = np.cross(q[:, 1]-q[:, 0], q[:, 2]-q[:, 0])
            triangle_area = .5*np.linalg.norm(cross, axis=1)
            counts += np.sum(planes, axis=0)
            areas += np.sum(planes*triangle_area[:, None], axis=0)
            for local in np.flatnonzero(~valid)[:max(0, MAX_FAILURE_EXAMPLES-len(examples))]:
                i = int(indices[local])
                examples.append(failure_example(records, i, i+1, vertices, 'single_owner_face_inside_cube'))
        result['external_boundary_faces_by_plane'] = counts.tolist()
        result['external_boundary_area_m2_by_plane'] = areas.tolist()
    result['seconds'] = time.perf_counter()-started
    write_json(result_path, result)
    vertices._mmap.close()
    del vertices
    return result


def run_check(folder: Path, *, buckets=128, chunk_tets=200000, workers=8, work_dir=None) -> dict:
    folder = folder.resolve()
    if not 1 <= buckets <= 2048 or not 1 <= chunk_tets <= 1000000 or not 1 <= workers <= 26:
        raise ValueError('Invalid resource settings (workers must be 1..26)')
    spec = source_spec(folder)
    if work_dir is None:
        work = folder / 'all_faces_checks' / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    else:
        work = Path(work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    print(json.dumps(dict(stage='start', work_dir=str(work), source=spec, workers=workers, buckets=buckets)), flush=True)
    started = time.perf_counter()
    try:
        partition = make_buckets(spec, work, buckets, chunk_tets)
        results = []
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(check_bucket, b, spec, work): b for b in range(buckets)}
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                if len(results) % max(1, buckets//16) == 0 or len(results) == buckets:
                    print(json.dumps(dict(stage='count_owners', buckets_done=len(results), buckets_total=buckets,
                                          seconds=time.perf_counter()-started)), flush=True)
        # Detect accidental replacement/modification during this check without a checksum.
        if source_spec(folder) != spec:
            raise RuntimeError('Input file metadata changed during face check; do not accept this result')
        sums = ['face_records', 'unique_faces', 'single_owner_faces', 'two_owner_faces',
                'more_than_two_owner_faces', 'repeated_owner_faces', 'internal_unmatched_faces']
        stats = {key: sum(result[key] for result in results) for key in sums}
        if stats['face_records'] != spec['tetrahedra']*4:
            raise RuntimeError('Face enumeration count does not equal four times tetrahedron count')
        passed = not any(stats[key] for key in ('more_than_two_owner_faces', 'repeated_owner_faces', 'internal_unmatched_faces'))
        area = np.sum([r['external_boundary_area_m2_by_plane'] for r in results], axis=0)
        failure_examples = [example for r in sorted(results, key=lambda x: x['bucket']) for example in r['failure_examples']][:MAX_FAILURE_EXAMPLES]
        report = dict(schema='r7b_localfit_all_volume_faces_v1', status='PASS' if passed else 'FAILED',
                      all_volume_faces_checked=True, all_face_topology_checks_passed=passed,
                      source=spec, **stats, maximum_owner_count=max(r['maximum_owner_count'] for r in results),
                      external_boundary_plane_order=['x=0', 'x=L', 'y=0', 'y=L', 'z=0', 'z=L'],
                      external_boundary_faces_by_plane=np.sum([r['external_boundary_faces_by_plane'] for r in results], axis=0).tolist(),
                      external_boundary_area_m2_by_plane=area.tolist(),
                      external_boundary_relative_area_error_by_plane=(area/spec['L_m']**2-1).tolist(),
                      area_is_additional_diagnostic_not_acceptance_substitute=True,
                      boundary_coordinate_tolerance_over_L=BOUNDARY_TOLERANCE_OVER_L,
                      owner_rule='Exactly two distinct tetrahedra per internal canonical triangular face; a single owner is allowed only on one cube exterior plane',
                      physical_interface_faces_are_internal=True, no_bad_elements_deleted=True,
                      source_mesh_modified=False, COMSOL_executed=False,
                      formal_electrical_result=False, geometry_success_written=False,
                      workers=workers, bucket_count=buckets, chunk_tetrahedra=chunk_tets,
                      maximum_bucket_face_records=partition['maximum_bucket_face_records'],
                      temporary_bucket_bytes=spec['tetrahedra']*4*16, work_directory=str(work),
                      failure_examples=failure_examples, failure_examples_limit=MAX_FAILURE_EXAMPLES,
                      total_seconds=time.perf_counter()-started, completed_local=datetime.now().isoformat(timespec='seconds'))
        write_json(work / 'all_faces_verification.json', report)
        write_json(folder / 'all_faces_verification.json', report)
        if not passed:
            write_json(work / 'ALL_FACES_FAILED.json', report)
        return report
    except Exception as exc:
        write_json(work / 'CHECK_EXECUTION_FAILED.json', dict(status='CHECK_NOT_COMPLETED', error=repr(exc),
                   source=spec, no_source_modified=True, work_directory=str(work), seconds=time.perf_counter()-started))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_folder', type=Path)
    parser.add_argument('--buckets', type=int, default=128)
    parser.add_argument('--chunk-tets', type=int, default=200000)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--work-dir', type=Path, help='Exact existing work directory for resume, or a new empty directory')
    args = parser.parse_args()
    report = run_check(args.input_folder, buckets=args.buckets, chunk_tets=args.chunk_tets,
                       workers=args.workers, work_dir=args.work_dir)
    print(json.dumps({key: report[key] for key in ('status', 'face_records', 'unique_faces', 'internal_unmatched_faces',
                                                  'more_than_two_owner_faces', 'total_seconds', 'work_directory')}, ensure_ascii=False), flush=True)
    return 0 if report['all_face_topology_checks_passed'] else 2


if __name__ == '__main__':
    sys.exit(main())
