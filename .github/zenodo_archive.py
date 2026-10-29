"""Archive the approved v1.0.0 release in its prepared Zenodo deposition."""
import json
import http.client
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import shutil
import tempfile
import time
import urllib.request
import urllib.parse

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

class ProgressReader:
    def __init__(self, stream, name, size):
        self.stream, self.name, self.size = stream, name, size
        self.sent = 0
        self.last_report = time.time()

    def read(self, count):
        block = self.stream.read(count)
        self.sent += len(block)
        now = time.time()
        if now - self.last_report >= 30 or not block:
            print('Transfer progress: ' + self.name + ' ' +
                  str(round(100*self.sent/self.size, 1)) + '%', flush=True)
            self.last_report = now
        return block

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
        def upload(item):
            name, size = item
            if existing.get(name) == size:
                print('Already archived: ' + name, flush=True)
                return
            path = Path(directory) / name
            print('Downloading release asset: ' + name, flush=True)
            with urllib.request.urlopen(assets[name]['browser_download_url'], timeout=600) as src, path.open('wb') as dst:
                shutil.copyfileobj(src, dst, length=4*1024*1024)
            if path.stat().st_size != size:
                raise RuntimeError('Incomplete download: ' + name)
            print('Uploading to Zenodo: ' + name, flush=True)
            start = time.time()
            for attempt in range(3):
                target=urllib.parse.urlsplit(bucket+'/'+name)
                connection=http.client.HTTPSConnection(target.hostname,timeout=300,blocksize=1024*1024)
                try:
                    with path.open('rb') as data:
                        connection.request('PUT',target.path,body=ProgressReader(data,name,size),headers={'Authorization':'Bearer '+TOKEN,'User-Agent':'LayeredRock research archive (GitHub SureWangchd/LayeredRock_FractureConductivity)','Accept':'application/json','Content-Type':'application/octet-stream','Content-Length':str(size)})
                        res=connection.getresponse()
                        body=res.read()
                        if res.status>=300:
                            raise RuntimeError('Zenodo upload HTTP '+str(res.status)+': '+body.decode(errors='replace').replace(TOKEN,'[redacted]')[:500])
                    break
                except (OSError,http.client.HTTPException) as error:
                    print('Transfer interrupted ('+type(error).__name__+'), attempt '+str(attempt+1)+': '+name,flush=True)
                    if attempt==2:raise
                    time.sleep(15*(attempt+1))
                finally:
                    connection.close()
            print('Archived in ' + str(round(time.time()-start,1)) + ' s: ' + name, flush=True)
            path.unlink()
        with ThreadPoolExecutor(max_workers=3) as executor:
            list(executor.map(upload,EXPECTED.items()))
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
