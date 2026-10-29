"""Select a small, declared set of saved fields for round-4 diagnostics."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs" / "round4_crossover"


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def token(sf):
    return "sf_" + format(sf, ".10g").replace(".", "p")


entries = []
brackets = read(OUT / "periphery_brackets.json")
near = next(x for x in brackets if x["layout"] == "swap_near")
for layout, sigmas in [("original", [.1, .159375, .3]),
                       ("swap_near", [.1, .159375, near["lower"], .3])]:
    for role in ("A", "B"):
        case = (OUT / "original" / role if layout == "original" else
                OUT / "periphery" / f"r4_strict_{role}_swap_near")
        for sf in sigmas:
            state = case / token(sf)
            result = read(state / "SUCCESS.json")
            if result.get("round4_reused"):
                source = Path(result["round4_source_folder"])
            else:
                source = state
            model = source / "solution.mph"
            if not model.is_file():
                model = Path(read(source / "reused_model.json")["original_mph"])
            entries.append({
                "layout": layout, "role": role, "sigma_f_S_m": sf,
                "model": model.as_posix(),
                "parameters": (source / "parameters.json").as_posix(),
                "selections": (case / "selections.json").as_posix(),
                "output_folder": (OUT / "local" / f"{layout}_{role}_{token(sf)}").as_posix(),
                "selection_reason": "Original below/near/above crossing; near-exchange common conductivities plus its own lower bracket. Both peripheral layouts remain in the global tables.",
            })
(OUT / "local_targets.json").write_text(json.dumps(entries, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"Selected {len(entries)} saved target fields; no solve or predictor modification.")
