"""Reconstruct the archived S5 empirical scale; no external project or solver."""
from pathlib import Path
import json
root = Path(__file__).resolve().parent
load = lambda p: json.loads(p.read_text(encoding="utf-8-sig"))
src = root / "source_evidence"
m2 = load(src / "states/S5_M2/water_10/PAIR_SUCCESS.json")
m3 = load(src / "states/S5_M3/water_10/PAIR_SUCCESS.json")
check = load(src / "tolerance_checks/S5_M3_water_10_pair_four_states/comparison.json")
terms = {
    "abs_gT_M3_minus_M2_S": abs(m3["gT_S"] - m2["gT_S"]),
    "abs_gH_M3_minus_M2_S": abs(m3["gH_S"] - m2["gH_S"]),
    "abs_background_mismatch_M2_S": abs(m2["G0_H_S"] - m2["G0_T_S"]),
    "abs_background_mismatch_M3_S": abs(m3["G0_H_S"] - m3["G0_T_S"]),
    "sum_abs_four_state_tolerance_changes_S": sum(abs(v) for v in check["conductance_changes_by_role_S"]),
}
scale = 3 * sum(terms.values())
assert abs(scale - 4.5272845192786265e-9) < 1e-20
print(json.dumps({"terms": terms, "empirical_multiplier": 3, "S5_combined_observed_scale_S": scale,
    "S5_combined_observed_scale_microS": scale * 1e6, "S1_new_accepted_strict_states": 0,
    "S1_four_state_scale_known": False, "S1_scale_is_zero": False}, indent=2))
