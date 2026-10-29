"""Archive the approved v1.0.0 release in its prepared Zenodo deposition."""
import json
import os
from pathlib import Path
import shutil
import tempfile
import time
import urllib.request

REPOSITORY = 'SureWangchd/LayeredRock_FractureConductivity'
DEPOSIT = 'https://zenodo.org/api/deposit/depositions/23755034'
TOKEN = os.environ['ZENODO_TOKEN']
EXPECTED = {
    'LayeredRock_FractureConductivity-v1.0.0.zip': 205163991,
    'representative_models_2D.zip': 19018921,
    'representative_models_synthetic3D.zip': 859176971,
    'representative_models_CT_S1_M3.zip': 1325648149,
    'representative_models_CT_S5_M3.zip': 1293473200,
    'model_inventory.csv': 2196,
}

def json_request(url, method='GET'):
    headers = {'User-Agent': 'LayeredRock research archive'}
    if url.startswith('https://zenodo.org/'):
        headers['Authorization'] = 'Bearer ' + TOKEN
    req = urllib.request.Request(url, method=method, headers=headers)
    with urllib.request.urlopen(req, timeout=180) as res:
        return json.load(res)

def main():
    release = json_request('https://api.github.com/repos/' + REPOSITORY + '/releases/tags/v1.0.0')
    assets = {a['name']: a for a in release['assets']}
    for name, size in EXPECTED.items():
        if name not in assets or assets[name]['size'] != size:
            raise RuntimeError('Release file missing or incomplete: ' + name)
    record = json_request(DEPOSIT)
    if record.get('submitted'):
        print('Already published: https://doi.org/10.5281/zenodo.23755034')
        return
    bucket = record['links']['bucket']
    if not bucket.startswith('https://zenodo.org/api/files/'):
        raise RuntimeError('Unexpected Zenodo upload location')
    existing = {f['filename']: f['filesize'] for f in record['files']}
    with tempfile.TemporaryDirectory(prefix='layeredrock-') as directory:
        for name, size in EXPECTED.items():
            if existing.get(name) == size:
                print('Already archived: ' + name, flush=True)
                continue
            path = Path(directory) / name
            print('Downloading release asset: ' + name, flush=True)
            with urllib.request.urlopen(assets[name]['browser_download_url'], timeout=600) as src, path.open('wb') as dst:
                shutil.copyfileobj(src, dst, length=4*1024*1024)
            if path.stat().st_size != size:
                raise RuntimeError('Incomplete download: ' + name)
            print('Uploading to Zenodo: ' + name, flush=True)
            start = time.time()
            with path.open('rb') as data:
                req = urllib.request.Request(bucket + '/' + name, data=data, method='PUT', headers={'Authorization':'Bearer ' + TOKEN, 'Content-Type':'application/octet-stream', 'Content-Length':str(size)})
                with urllib.request.urlopen(req, timeout=1800) as res:
                    json.load(res)
            print('Archived in ' + str(round(time.time()-start,1)) + ' s: ' + name, flush=True)
            path.unlink()
    record = json_request(DEPOSIT)
    actual = {f['filename']:f['filesize'] for f in record['files']}
    if any(actual.get(name) != size for name,size in EXPECTED.items()):
        raise RuntimeError('Archive is incomplete; deposit remains unpublished')
    record = json_request(DEPOSIT + '/actions/publish', method='POST')
    if not record.get('submitted'):
        raise RuntimeError('Zenodo did not confirm publication')
    print('Published: https://doi.org/10.5281/zenodo.23755034', flush=True)

if __name__ == '__main__':
    main()
