"""Reproduce Round 6 package diagnostics from saved Round 5 scalars/ports.

No COMSOL calls, no target-field input, no writes in Round 5.  This is a
retrospective implementation check, not an independent validation.
"""
from pathlib import Path
import json
import subprocess
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PACK = ROOT / 'round6_taskpack'
OUT = ROOT / 'runs/round6_applicability/round5_reanalysis'
sys.path.insert(0, str(PACK))
from modal_applicability import modal_interval, contrast_interval


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    test = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s',
        str(PACK), '-p', 'test_modal_applicability.py', '-v'], capture_output=True,
        text=True, encoding='utf8')
    (OUT/'algebraic_selftest.txt').write_text(test.stdout + test.stderr, encoding='utf8')
    if test.returncode:
        raise RuntimeError('Algebraic tests failed; see saved log')
    replay = subprocess.run([sys.executable, str(PACK/'check_round5.py'),
        '--results', str(ROOT/'runs/round5_modal_closure'), '--out', str(OUT)],
        capture_output=True, text=True, encoding='utf8')
    (OUT/'check_round5_stdout.txt').write_text(replay.stdout + replay.stderr, encoding='utf8')
    if replay.returncode:
        raise RuntimeError('Round 5 reproduction failed; see saved log')
    old = ROOT/'runs/round5_modal_closure'
    c = pd.read_csv(old/'closure_comparison.csv')
    c = c[c.phase.eq('prospective')].set_index(['role','layout','sigma_f_S_m','K'])
    p = pd.read_csv(old/'predicted_port_modes.csv')
    p = p[p.phase.eq('new') & p.mesh.eq('base')].set_index(
        ['role','layout','sigma_f_S_m','K','side','n'])
    y = pd.read_csv(old/'port_coefficients.csv').set_index(['layout','side','n'])
    sides = ['lower','upper']

    def yy(layout, actual=True):
        col = 'Y_actual_S_m2' if actual else 'Y_reference_S_m2'
        return np.array([[y.loc[(layout,side,n),col] for n in range(33)] for side in sides])

    def trace(role, layout, sf, order):
        return np.array([[p.loc[(role,layout,sf,order,side,n),'v_predicted_V']
                          for n in range(33)] for side in sides])

    tail = []
    frozen = []
    for (role, layout, sf, k), target in c.iterrows():
        if k != 32:
            continue
        base = c.loc[(role,layout,sf,8)]
        yt = yy(layout)
        yr = yy(layout, False)
        yr[:, :9] = yt[:, :9]
        r = modal_interval(float(base.g_prediction_S), trace(role,layout,sf,8), yr, yt)
        g8, g32, gt = map(float, [base.g_prediction_S,target.g_prediction_S,target.g_COMSOL_S])
        corrected = g8+r.first_variation_S
        tail.append(dict(role=role,layout=layout,sigma_f_S_m=sf,M8_g_S=g8,
            M32_g_S=g32,truth_g_S=gt,tail_Q_S=r.first_variation_S,
            tail_B_S=r.remainder_bound_S,tail_abs_terms_S=r.absolute_mode_term_sum_S,
            corrected_g_S=corrected,lower_model_S=r.target_gain_lower_S,
            upper_model_S=r.target_gain_upper_S,correction_model_error_S=corrected-g32,
            model_bound_ok=r.target_gain_lower_S-1e-12<=g32<=r.target_gain_upper_S+1e-12,
            COMSOL_bound_ok=r.target_gain_lower_S-1e-12<=gt<=r.target_gain_upper_S+1e-12,
            M8_error_percent=100*abs(g8-gt)/abs(gt),
            corrected_error_percent=100*abs(corrected-gt)/abs(gt),
            M8_model_tail_bound_S=r.absolute_reference_error_bound_S,
            corrected_model_bound_S=r.remainder_bound_S))
        if layout == 'P77':
            continue
        ref = c.loc[(role,'P77',sf,32)]
        rr = modal_interval(float(ref.g_prediction_S),trace(role,'P77',sf,32),yy('P77'),yt)
        dmodel = float(target.g_prediction_S-ref.g_prediction_S)
        dtruth = float(target.g_COMSOL_S-ref.g_COMSOL_S)
        q = rr.first_variation_S
        point = float(ref.g_prediction_S)+q
        frozen.append(dict(role=role,sigma_f_S_m=sf,layout=layout,Q_S=q,
            remainder_bound_S=rr.remainder_bound_S,model_delta_S=dmodel,
            COMSOL_delta_S=dtruth,pred_total_S=point,COMSOL_total_S=gt,
            total_error_pct=100*abs(point-gt)/abs(gt),
            effect_error_pct=100*abs(q-dtruth)/abs(dtruth),
            remainder_actual_S=q-dmodel,bound_lower_S=rr.correction_lower_S,
            bound_upper_S=rr.correction_upper_S,
            model_inside_bound=rr.correction_lower_S-1e-12<=dmodel<=rr.correction_upper_S+1e-12,
            COMSOL_inside_bound=rr.correction_lower_S-1e-12<=dtruth<=rr.correction_upper_S+1e-12,
            # Reproduce the package's retrospective POINT-sign comparison.
            # This column is not an interval certificate; B/P11/.1 crosses 0.
            sign_predicted=(np.sign(q)==np.sign(dtruth)),
            abs_term_sum_S=rr.absolute_mode_term_sum_S))
    df = pd.DataFrame(tail).sort_values(['role','layout','sigma_f_S_m'])
    df.to_csv(OUT/'M8_tail_diagnostic_48.csv', index=False)
    pd.DataFrame(frozen).to_csv(OUT/'rejected_frozen_reference_point_rule.csv', index=False)
    contrasts = {
      'lower_7_to_1_at_upper_1':{'P11':1,'P71':-1},
      'lower_7_to_1_at_upper_7':{'P17':1,'P77':-1},
      'upper_7_to_1_at_lower_1':{'P11':1,'P17':-1},
      'upper_7_to_1_at_lower_7':{'P71':1,'P77':-1},
      'lower_7_to_4_at_upper_1':{'P41':1,'P71':-1},
      'lower_4_to_1_at_upper_1':{'P11':1,'P41':-1},
      'upper_7_to_4_at_lower_1':{'P14':1,'P17':-1},
      'upper_4_to_1_at_lower_1':{'P11':1,'P14':-1},
      'joint_P11_minus_P77':{'P11':1,'P77':-1},
      'interaction_1_7':{'P11':1,'P17':-1,'P71':-1,'P77':1}}
    idx=df.set_index(['role','layout','sigma_f_S_m']); rows=[]
    # The supplied table uses D = gB - gA. State the convention explicitly.
    for quantity, roles in [('gA',{'A':1}),('gB',{'B':1}),('D',{'A':-1,'B':1})]:
        for effect, layouts in contrasts.items():
            for sf in sorted(df.sigma_f_S_m.unique()):
                keys = [(role,layout,sf) for role in roles for layout in layouts]
                a = np.array([roles[role]*layouts[layout] for role,layout,_ in keys])
                selected = idx.loc[keys]
                low, high = contrast_interval(selected.lower_model_S.to_numpy(),
                                               selected.upper_model_S.to_numpy(), a)
                truth = float(a@selected.truth_g_S.to_numpy())
                original = float(a@selected.M8_g_S.to_numpy())
                corrected = float(a@selected.corrected_g_S.to_numpy())
                rows.append(dict(quantity=quantity,effect=effect,sigma_f_S_m=sf,
                    truth_S=truth,original_pred_S=original,
                    original_error_percent=100*abs(original-truth)/abs(truth),
                    corrected_pred_S=corrected,
                    corrected_error_percent=100*abs(corrected-truth)/abs(truth),
                    low_S=low,high_S=high,truth_inside=low-1e-12<=truth<=high+1e-12,
                    sign_certified=(low>0 or high<0),positive_if_certified=low>0))
    pd.DataFrame(rows).to_csv(OUT/'M8_tail_contrasts.csv',index=False)
    comparison = {'summary_exact_match':json.loads((OUT/'summary.json').read_text(encoding='utf8'))
        == json.loads((PACK/'reanalysis/summary.json').read_text(encoding='utf8')),
        'scope':'saved R5 data only; no old solves; no target-field input',
        'tables':{}}
    for file, keys in [
        ('M0_applicability_48.csv',['role','layout','sigma_f_S_m']),
        ('M8_tail_diagnostic_48.csv',['role','layout','sigma_f_S_m']),
        ('M8_tail_contrasts.csv',['quantity','effect','sigma_f_S_m']),
        ('rejected_frozen_reference_point_rule.csv',['role','layout','sigma_f_S_m'])]:
        actual=pd.read_csv(OUT/file).set_index(keys).sort_index()
        expected=pd.read_csv(PACK/'reanalysis'/file).set_index(keys).sort_index()
        if not actual.index.equals(expected.index):
            raise ValueError('Reference rows differ: '+file)
        cols=expected.select_dtypes(include='number').columns
        diff=(actual[cols]-expected[cols]).abs()
        boolcols=expected.select_dtypes(include='bool').columns
        comparison['tables'][file]={'rows':len(actual),
            'numeric_allclose_atol_1e-12_rtol_1e-8':bool(np.allclose(actual[cols],expected[cols],atol=1e-12,rtol=1e-8)),
            'max_absolute_difference_by_column':diff.max().to_dict(),
            'boolean_mismatches_by_column':(actual[boolcols]!=expected[boolcols]).sum().to_dict()}
    comparison['M8_tail_summary']={'finite_model_covered':int(df.model_bound_ok.sum()),
        'full_COMSOL_inside_finite_model_interval':int(df.COMSOL_bound_ok.sum()),
        'states':len(df),'warning':'COMSOL containment is empirical, not a promised bound'}
    (OUT/'independent_reproduction.json').write_text(json.dumps(comparison,indent=2,
        ensure_ascii=False),encoding='utf8')
    print(json.dumps(comparison,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
