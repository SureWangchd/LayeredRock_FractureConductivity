"""Create isolated fixed-mesh tolerance-check adapters from accepted r7b code."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'runs' / 'round7b_review_20261029'
DST = OUT / 'scripts'
DST.mkdir(parents=True, exist_ok=True)
changes = []
for name in ('r7b_fix_tolerance_check', 'r7b_fix_pair_tolerance_check'):
    source = ROOT / 'scripts' / (name + '.m')
    text = source.read_text(encoding='utf-8-sig')
    new_name = name.replace('r7b_fix_', 'r7b_review_')
    text = text.replace(name, new_name)
    text = text.replace('r7b_fix_tolerance_check(', 'r7b_review_tolerance_check(')
    text = text.replace("folder=fullfile(out,'tolerance_checks',", "folder=fullfile(root,'runs','round7b_review_20261029','tolerance_checks',")
    target = DST / (new_name + '.m')
    if target.exists():
        raise RuntimeError(f'Adapter already exists: {target}')
    target.write_text(text, encoding='utf-8')
    changes.append({'source': str(source), 'adapter': str(target),
                    'changes': 'function names and destination directory only; source geometry, solution, frozen kappa, solver, tolerance, and validation retained'})
(OUT / 'adapter_manifest.json').write_text(json.dumps(changes, ensure_ascii=False, indent=2), encoding='utf-8')
print(str(OUT))
