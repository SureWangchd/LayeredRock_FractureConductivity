"""Reproduce the assistant's retrospective R5 diagnostic, with no COMSOL calls.
Usage: python check_round5.py --results F:/岩层裂缝/runs/round5_modal_closure --out reanalysis
No target field is used to form Q or B. Target SCALAR answers are read only
for the comparisons below. This is retrospective, NOT an independent test.
"""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
from modal_applicability import modal_interval, safe_for_relative_gain


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    c = pd.read_csv(args.results/'closure_comparison.csv')
    p = pd.read_csv(args.results/'predicted_port_modes.csv')
    y = pd.read_csv(args.results/'port_coefficients.csv')
    c = c[c.phase.eq('prospective')].copy()
    if c.duplicated(['role','layout','sigma_f_S_m','K']).any():
        raise ValueError('Duplicate comparison states')
    c = c.set_index(['role','layout','sigma_f_S_m','K'])
    y = y.set_index(['layout','side','n'])
    rows = []
    sides = ['lower', 'upper']; order = 32
    for (role, layout, sf, k), val in c.iterrows():
        if k != 32: continue
        # All M0 layouts for fixed role/sf are identical. Use P77 explicitly.
        base = c.loc[(role, 'P77', sf, 0)]
        modes = p[(p.phase=='new') & (p.mesh=='base') & (p.role==role) &
                  (p.layout=='P77') & np.isclose(p.sigma_f_S_m, sf, rtol=0, atol=1e-12) &
                  (p.K==0) & (p.n<=order)].set_index(['side','n'])
        if len(modes) != 2*(order+1): raise ValueError('Missing or duplicate reference modes')
        v = np.array([[modes.loc[(side,n)].v_predicted_V for n in range(order+1)] for side in sides])
        yr = np.array([[y.loc[(layout,side,n)].Y_reference_S_m2 for n in range(order+1)] for side in sides])
        yt = np.array([[y.loc[(layout,side,n)].Y_actual_S_m2 for n in range(order+1)] for side in sides])
        r = modal_interval(float(base.g_prediction_S), v, yr, yt)
        diff = float(val.g_prediction_S-base.g_prediction_S)
        actualerr = abs(float(base.g_prediction_S-val.g_COMSOL_S))/abs(float(val.g_COMSOL_S))
        rec = dict(role=role,layout=layout,sigma_f_S_m=float(sf),reference='M0',target='M32',
                   reference_gain_S=float(base.g_prediction_S),M32_gain_S=float(val.g_prediction_S),
                   COMSOL_gain_S=float(val.g_COMSOL_S),model_delta_S=diff,
                   empirical_M0_relative_error=actualerr,
                   bound_holds_for_M32=(r.correction_lower_S-1e-12<=diff<=r.correction_upper_S+1e-12),
                   **r.to_dict())
        for tol in [0.01,0.02,0.05,0.1]:
            tag=str(int(tol*100))+'pct'
            rec['model_accepted_'+tag]=safe_for_relative_gain(r,tol)
            rec['COMSOL_actual_pass_'+tag]=actualerr<=tol
        rows.append(rec)
    df=pd.DataFrame(rows).sort_values(['role','sigma_f_S_m','layout'])
    df.to_csv(args.out/'M0_applicability_48.csv',index=False)
    summaries={}
    for tol in [1,2,5,10]:
        accepted=df['model_accepted_'+str(tol)+'pct']; actual=df['COMSOL_actual_pass_'+str(tol)+'pct']
        summaries[str(tol)+'pct']={
          'model_accepted':int(accepted.sum()),'needs_more_information':int((~accepted).sum()),
          'COMSOL_actual_pass':int(actual.sum()),'retrospective_false_accept':int((accepted&~actual).sum()),
          'conservative_flags_that_COMSOL_passed':int((~accepted&actual).sum())}
    # Cancellation example using the unrounded COMSOL scalar results only.
    b = c.xs(('B',0.18,32),level=['role','sigma_f_S_m','K'])['g_COMSOL_S']
    lower=.5*((b['P17']-b['P77'])+(b['P11']-b['P71']))
    upper=.5*((b['P71']-b['P77'])+(b['P11']-b['P17']))
    net=b['P11']-b['P77']
    summary=dict(
       purpose='retrospective criterion diagnostic; no new COMSOL solve or blind validation',
       reference_modes='M0 self-consistent reference trace, not a target field',
       bound_target='M32 mixed finite-modal model, not continuum/full-COMSOL error certification',
       n_states=len(df),bound_holds_for_M32=int(df.bound_holds_for_M32.sum()),
       tolerances=summaries,
       cancellation_example=dict(role='B',sigma_f_S_m=.18,lower_effect_uS=float(lower*1e6),
          upper_effect_uS=float(upper*1e6),net_uS=float(net*1e6),
          algebraic_amplification=float((abs(lower)+abs(upper))/abs(net)),
          warning='Net sign not mesh-resolved in R5; amplification is a descriptive ratio, not a certified value'))
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf8')
    print(json.dumps(summary,indent=2,ensure_ascii=False))

if __name__=='__main__': main()
