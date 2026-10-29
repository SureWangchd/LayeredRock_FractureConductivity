"""Publish a prepared new version, reusing all inherited model archives."""
from datetime import datetime, timezone
import http.client
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import time
import urllib.parse
import urllib.request
import zipfile

REPOSITORY = 'SureWangchd/LayeredRock_FractureConductivity'
MODELS = {
    'representative_models_2D.zip': 19018921,
    'representative_models_synthetic3D.zip': 859176971,
    'representative_models_CT_S1_M3.zip': 1325648149,
    'representative_models_CT_S5_M3.zip': 1293473200,
}


def main():
    token = os.environ['ZENODO_TOKEN']
    record_id = os.environ['DRAFT_RECORD']
    tag = os.environ['RELEASE_TAG']
    if not record_id.isdigit() or not re.fullmatch(r'v\d+\.\d+\.\d+', tag):
        raise ValueError('Invalid draft ID or release tag')
    deposit = 'https://zenodo.org/api/deposit/depositions/' + record_id
    expected_names = ['LayeredRock_FractureConductivity-' + tag + '.zip',
                      'model_inventory_manuscript_20261029.csv']

    def request(url, method='GET', payload=None):
        if not (url.startswith('https://zenodo.org/api/deposit/depositions/') or
                url.startswith('https://api.github.com/repos/' + REPOSITORY + '/')):
            raise ValueError('Unexpected API endpoint')
        headers = {'User-Agent': 'LayeredRock versioned archive', 'Accept': 'application/json'}
        if url.startswith('https://zenodo.org/'):
            headers['Authorization'] = 'Bearer ' + token
        data = None if payload is None else json.dumps(payload).encode()
        if payload is not None:
            headers['Content-Type'] = 'application/json'
        with urllib.request.urlopen(urllib.request.Request(url, data=data, method=method, headers=headers), timeout=180) as r:
            body = r.read()
            return json.loads(body) if body else None

    def receipt(record):
        obj = {'published_utc': datetime.now(timezone.utc).isoformat(),
            'deposition_id': record['id'], 'version': record['metadata'].get('version'),
            'doi': record.get('doi') or record['metadata'].get('prereserve_doi', {}).get('doi'),
            'submitted': record.get('submitted', False),
            'release_url': 'https://github.com/' + REPOSITORY + '/releases/tag/' + tag,
            'files': [{'name': f['filename'], 'bytes': f['filesize']} for f in record.get('files', [])],
            'inherited_model_archive_bytes': sum(MODELS.values()),
            'inherited_models_reuploaded': False, 'previous_published_version_modified': False}
        Path('zenodo_newversion_publication_receipt.json').write_text(json.dumps(obj, indent=2) + '\n')
        print(json.dumps(obj, indent=2), flush=True)

    record = request(deposit)
    if record['metadata'].get('version') != tag[1:]:
        raise RuntimeError('Prepared Zenodo version differs from the GitHub release version')
    if record.get('submitted'):
        receipt(record)
        return
    existing = {f['filename']: f for f in record.get('files', [])}
    if any(existing.get(n, {}).get('filesize') != size for n, size in MODELS.items()):
        raise RuntimeError('Inherited model archives do not match the approved unchanged inventory')
    release = request('https://api.github.com/repos/' + REPOSITORY + '/releases/tags/' + tag)
    if release.get('draft'):
        raise RuntimeError('Expected a published GitHub release')
    assets = {a['name']: a for a in release['assets']}
    if any(n not in assets or assets[n]['size'] <= 0 for n in expected_names):
        raise RuntimeError('Updated code/data package or current model inventory is missing')
    bucket = record['links']['bucket']
    if not bucket.startswith('https://zenodo.org/api/files/'):
        raise ValueError('Unexpected upload bucket')
    # Delete obsolete code/data packages only from this distinct unpublished draft.
    # The original published record and all model ZIP files remain untouched.
    for name, entry in existing.items():
        if name.startswith('LayeredRock_FractureConductivity-v') and name.endswith('.zip') and name not in expected_names:
            request(deposit + '/files/' + entry['id'], method='DELETE')
    with tempfile.TemporaryDirectory(prefix='layeredrock-version-') as tmp:
        for name in expected_names:
            size = assets[name]['size']
            if existing.get(name, {}).get('filesize') == size:
                print('Already staged: ' + name, flush=True)
                continue
            if name in existing:
                request(deposit + '/files/' + existing[name]['id'], method='DELETE')
            path = Path(tmp) / name
            with urllib.request.urlopen(assets[name]['browser_download_url'], timeout=600) as src, path.open('wb') as dst:
                shutil.copyfileobj(src, dst, length=4 * 1024 * 1024)
            if path.stat().st_size != size:
                raise RuntimeError('Incomplete release asset download: ' + name)
            if name.endswith('.zip'):
                prefix = 'LayeredRock_FractureConductivity-' + tag[1:] + '/'
                with zipfile.ZipFile(path) as z:
                    summary = json.loads(z.read(prefix + 'package_summary_v1.1.0.json'))
                    if (summary.get('status') != 'publication_ready' or
                            summary.get('version') != tag[1:] or not summary.get('CT_supplement')):
                        raise RuntimeError('Expected the coordinated package including the CT supplement')
                    if (summary.get('S1_new_accepted_strict_states') != 0 or
                            summary.get('S1_numerical_check_complete') is not False or
                            summary.get('S5_accepted_stricter_states') != 4):
                        raise RuntimeError('The package must preserve the incomplete S1 and completed S5 check scope')
                    canonical = prefix + 'simulationCase/parameters_and_results/runs/'
                    count = sum(n.startswith(canonical) and '/predictions/' in n and n.endswith('/prediction.json') for n in z.namelist())
                    if count != 80:
                        raise RuntimeError('The canonical prospective prediction batches are incomplete')
            for attempt in range(3):
                u = urllib.parse.urlsplit(bucket + '/' + urllib.parse.quote(name))
                connection = http.client.HTTPSConnection(u.hostname, timeout=600, blocksize=1024 * 1024)
                try:
                    with path.open('rb') as data:
                        connection.request('PUT', u.path, body=data, headers={
                            'Authorization': 'Bearer ' + token, 'User-Agent': 'LayeredRock versioned archive',
                            'Content-Type': 'application/octet-stream', 'Content-Length': str(size)})
                        res = connection.getresponse(); res.read()
                        if res.status >= 300:
                            raise RuntimeError('Zenodo upload HTTP ' + str(res.status))
                    print('Staged: ' + name + ' (' + str(size) + ' bytes)', flush=True)
                    break
                except (OSError, http.client.HTTPException):
                    if attempt == 2:
                        raise
                    time.sleep(15 * (attempt + 1))
                finally:
                    connection.close()
    record = request(deposit)
    actual = {f['filename']: f['filesize'] for f in record.get('files', [])}
    if any(actual.get(n) != size for n, size in MODELS.items()) or any(actual.get(n) != assets[n]['size'] for n in expected_names):
        raise RuntimeError('Draft files are incomplete; leave unpublished')
    metadata = dict(record['metadata'])
    description = metadata.get('description', '').replace(
        'COMSOL Multiphysics 6.3 with LiveLink for MATLAB.',
        'COMSOL Multiphysics 6.3.0.290 with LiveLink for MATLAB R2025b.')
    ct_note = ('<p>The versioned source/data package also contains the bounded CT solver-tolerance '
        'review tables, actual execution logs and reproduction entry points. S1 returned no '
        'solver-accepted new stricter solution in the bounded attempts (0 of 4 states); its '
        'four-state solver-change and full combined scales remain unknown, not zero. S5 has '
        '4 of 4 previously completed stricter states, with combined observed scale '
        '0.004527284519 microS. Completing this archive does not imply that the S1 stricter check passed. '
        'The original thirteen '
        'representative MPH files remain in the inherited model archives; no additional multi-GB '
        'model archive is required for this revision.</p>')
    if 'bounded CT solver-tolerance review tables' not in description:
        description += ct_note
    metadata['description'] = description
    if 'GPL-3.0-only' not in description or 'CC BY 4.0' not in description:
        raise RuntimeError('The code/data license distinction must be explicit before publication')
    request(deposit, method='PUT', payload={'metadata': metadata})
    record = request(deposit + '/actions/publish', method='POST')
    if not record.get('submitted'):
        raise RuntimeError('Publication was not confirmed')
    receipt(record)


if __name__ == '__main__':
    main()
