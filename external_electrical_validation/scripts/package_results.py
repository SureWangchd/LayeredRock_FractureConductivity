"""Create a local, checksum-verifiable delivery archive (no upload)."""
from pathlib import Path
import json, hashlib, zipfile, argparse

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args()
    out=Path(args.output).resolve();out.parent.mkdir(parents=True,exist_ok=True)
    files=sorted(f for f in ROOT.rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.name not in ['SHA256SUMS.json'] and f.resolve()!=out and f.suffix not in ['.pyc','.zip'])
    hashes={str(f.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    manifest=ROOT/'SHA256SUMS.json';manifest.write_text(json.dumps(hashes,indent=2,ensure_ascii=False),encoding='utf-8')
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for f in files+[manifest]:z.write(f,ROOT.name+'/'+str(f.relative_to(ROOT)))
    with zipfile.ZipFile(out) as z:
        bad=z.testzip()
        assert bad is None,bad
        for rel,digest in hashes.items():
            assert hashlib.sha256(z.read(ROOT.name+'/'+rel)).hexdigest()==digest,rel
    print(json.dumps(dict(archive=str(out),files=len(files)+1,bytes=out.stat().st_size,
                         archive_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),integrity='CRC and all per-file SHA256 passed'),ensure_ascii=False,indent=2))

if __name__=='__main__':main()
