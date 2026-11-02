"""Deterministic propagation of reported measurement tolerances, without fitting.

The interval is an extrema calculation conditional on fixed nominal geometry and
porosity closure, not a probability interval. All three input conductivities and
porosities are monotone in the chosen corner directions.
"""
import json
import numpy as np
import pandas as pd
from run_validation import ROOT,RUN,D5,cell_mesh,solve_axis,specimen,dump

def main():
    mesh=cell_mesh(RUN/'meshes/cell_a0.021_s3_l2.npz',level=2)
    pc=4*np.pi/3*.062*.021
    rows=[]
    for sample in ['F-0','F-45','F-90']:
      for closure in ['bulk','local']:
        for sign in [-1,1]:
            # Larger sigma_I, smaller phi_I, larger phi_F give the maximum
            # inferred background for the published porosity ordering.
            m={}
            for name in ['I-0','I-90']:
                sig=float(D5.loc[name,'conductivity_S_per_m'])*(1+sign*.02)
                por=float(D5.loc[name,'porosity_fraction'])-sign*.0001
                m[name]=np.log(sig/4.69)/np.log(por)
            phi=float(D5.loc[sample,'porosity_fraction'])+sign*.0001-pc
            if closure=='local': phi/=1-pc
            kn=4.69*phi**m['I-0'];kp=4.69*phi**m['I-90']
            sig=np.broadcast_to(np.diag([kp,kp,kn]),(len(mesh['cells']),3,3)).copy()
            sig[mesh['labels']==1]=np.eye(3)*4.69
            en=solve_axis(mesh,sig,2)['conductivity_S_per_m']
            ep=solve_axis(mesh,sig,0)['conductivity_S_per_m']
            angle=float(D5.loc[sample,'angle_to_bedding_normal_deg'])
            out=specimen(en,ep,angle,40,2);bg=specimen(kn,kp,angle,40,2)
            rows.append(dict(sample=sample,porosity_closure=closure,corner='minimum' if sign<0 else 'maximum',
                             predicted_S_per_m=out['conductivity_S_per_m'],
                             background_S_per_m=bg['conductivity_S_per_m'],
                             predicted_increment_S_per_m=out['conductivity_S_per_m']-bg['conductivity_S_per_m'],
                             normal_m=m['I-0'],parallel_m=m['I-90'],phi_background=phi))
    df=pd.DataFrame(rows);df.to_csv(ROOT/'MEASUREMENT_PROPAGATION.csv',index=False)
    results=[]
    for (sample,closure),g in df.groupby(['sample','porosity_closure'],sort=False):
        obs=float(D5.loc[sample,'conductivity_S_per_m'])
        results.append(dict(sample=sample,porosity_closure=closure,measured_S_per_m=obs,
            measured_min=.98*obs,measured_max=1.02*obs,
            predicted_min=g.predicted_S_per_m.min(),predicted_max=g.predicted_S_per_m.max(),
            transferred_background_min=g.background_S_per_m.min(),transferred_background_max=g.background_S_per_m.max(),
            inferred_increment_min=.98*obs-g.background_S_per_m.max(),
            inferred_increment_max=1.02*obs-g.background_S_per_m.min(),
            intervals_overlap=bool(max(.98*obs,g.predicted_S_per_m.min())<=min(1.02*obs,g.predicted_S_per_m.max()))))
    pd.DataFrame(results).to_csv(ROOT/'UNCERTAINTY_BUDGET.csv',index=False)
    # Ordinary dielectric displacement at bulk water epsilon_r=80 is tiny;
    # this does not bound electrode or interfacial polarization.
    dump(RUN/'dielectric_scale_check.json',dict(frequency_Hz=1000,assumed_relative_permittivity=80,
         conductivity_lower_S_per_m=.25,omega_epsilon_over_sigma=2*np.pi*1000*8.8541878128e-12*80/.25,
         interpretation='Bulk dielectric estimate only; not a measured phase spectrum or an electrode-polarization test.'))
    print(pd.DataFrame(results).to_string(index=False))

if __name__=='__main__':main()
