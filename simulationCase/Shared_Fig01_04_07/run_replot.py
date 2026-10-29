"""Replot manuscript Figures 1, 5 and 8 from saved data."""
from pathlib import Path
import importlib.util
import sys
import subprocess
root=Path(__file__).resolve().parent
scripts=root/'scripts/jgrse_figures'
sys.path.insert(0,str(scripts))
for script,folder in [('figure02_phenomenon_v3.py','Figure01'),('figure09_applicability_v8.py','Figure08')]:
    subprocess.run([sys.executable,str(scripts/script),'--outdir',str(root/'replots'/folder)],check=True)
spec=importlib.util.spec_from_file_location('mechanism',scripts/'figure06_mechanism_v4.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
mod.OUT=root/'replots/Figure05';mod.main()
