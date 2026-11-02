"""Create publication figure and result summaries from completed PDE runs only."""
from pathlib import Path
import json, platform, importlib.metadata
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
BLUE='#0000FF'; RED='#FF0000'; BLACK='#000000'
mpl.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix',
 'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7.5,'ytick.labelsize':7.5,
 'legend.fontsize':7,'axes.linewidth':.65,'lines.linewidth':1.1,'xtick.direction':'in',
 'ytick.direction':'in','xtick.top':True,'ytick.right':True,'text.color':'black',
 'axes.labelcolor':'black','xtick.color':'black','ytick.color':'black',
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})

def main():
    c=pd.read_csv(ROOT/'COMPARISON.csv')
    nominal=c[(c.scenario=='nominal')&(c.mesh_level==3)].set_index('sample')
    local=c[(c.scenario=='local_porosity')&(c.mesh_level==3)].set_index('sample')
    assert len(nominal)==len(local)==3, 'Figure requires real completed predictions.'
    hold=pd.read_csv(ROOT/'runs/intact_holdout.csv')
    h=hold[(hold.length_mm==40)&(hold.level==3)].iloc[0]
    hold['signed_error_percent']=100*(hold.predicted_S_per_m/hold.observed_S_per_m-1)
    hold.to_csv(ROOT/'HOLDOUT_COMPARISON.csv',index=False)
    u=pd.read_csv(ROOT/'UNCERTAINTY_BUDGET.csv').set_index(['sample','porosity_closure'])
    samples=['F-0','F-45','F-90']
    conv=[]
    for s in samples:
        gr=c[(c.scenario=='nominal')&(c['sample']==s)].sort_values('mesh_level')
        prev=None
        for _,r in gr.iterrows():
            mesh=json.loads((ROOT/f'runs/cell_results/{s}_w4.69_a1_s3_l{int(r.mesh_level)}_bulk.json').read_text())
            row=dict(sample=s,level=int(r.mesh_level),cell_nodes=mesh['normal_solve']['nodes'],
                     cell_tets=mesh['normal_solve']['tets'],specimen_nodes=int(r.specimen_nodes),specimen_tets=int(r.specimen_tets),
                     total_S_per_m=r.predicted_S_per_m,increment_S_per_m=r.predicted_fracture_increment_S_per_m,
                     crack_volume_relative_error_percent=100*(r.phi_crack_mesh/r.phi_crack_analytic-1))
            if prev is not None:
                row['total_change_from_previous_percent']=100*(r.predicted_S_per_m/prev.predicted_S_per_m-1)
                row['increment_change_from_previous_percent']=100*(r.predicted_fracture_increment_S_per_m/prev.predicted_fracture_increment_S_per_m-1)
            conv.append(row);prev=r
    pd.DataFrame(conv).to_csv(ROOT/'MESH_CONVERGENCE.csv',index=False)
    summary=[]
    for s in samples:
        n=nominal.loc[s];l=local.loc[s]
        summary.append(dict(sample=s,measured=n.measured_S_per_m,bulk_prediction=n.predicted_S_per_m,
              bulk_error_percent=n.signed_error_percent,local_prediction=l.predicted_S_per_m,
              local_error_percent=l.signed_error_percent,
              inferred_bulk_increment=n.measured_minus_transferred_background_S_per_m,
              predicted_bulk_increment=n.predicted_fracture_increment_S_per_m,
              inferred_local_increment=l.measured_minus_transferred_background_S_per_m,
              predicted_local_increment=l.predicted_fracture_increment_S_per_m))
    summary=pd.DataFrame(summary)
    summary.to_csv(ROOT/'MAIN_COMPARISON_SUMMARY.csv',index=False)
    # Figure contract: quantitative grid; no target is removed; all three F
    # samples at the preselected minimum pressure plus the held-out I45.
    fig,axs=plt.subplots(2,2,figsize=(7.0866,5.55)) # 180 mm wide
    fig.subplots_adjust(left=.095,right=.985,bottom=.09,top=.94,wspace=.31,hspace=.43)
    ax=axs[0,0];x=np.arange(4)
    obs=np.r_[h.observed_S_per_m,nominal.loc[samples,'measured_S_per_m']]
    pn=np.r_[h.predicted_S_per_m,nominal.loc[samples,'predicted_S_per_m']]
    pl=np.r_[h.predicted_S_per_m,local.loc[samples,'predicted_S_per_m']]
    ax.errorbar(x-.15,obs,yerr=.02*obs,fmt='ko',mfc='white',ms=4,capsize=2,lw=.8,label='Measured (±2%)')
    ax.plot(x,pn,'s',color=BLUE,mfc='white',ms=4,label='Bulk pore fraction')
    ax.plot(x+.15,pl,'^',color=RED,mfc='white',ms=4.7,label='Matrix pore fraction')
    ax.set(xticks=x,xticklabels=['I-45','F-0','F-45','F-90'],ylabel=r'Conductivity (S m$^{-1}$)',ylim=(.265,.382),title='Electrical response at 5 MPa')
    ax.legend(loc='upper left',frameon=False,handletextpad=.5,borderpad=.1,labelspacing=.3)
    ax=axs[0,1];x=np.arange(3)
    for closure,df,color,marker,offset,label in [('bulk',nominal,BLUE,'o',-.13,'Bulk'),('local',local,RED,'^',.13,'Matrix')]:
        rr=df.loc[samples]
        inf=rr.measured_minus_transferred_background_S_per_m.to_numpy()*1000
        mn=np.array([u.loc[(s,closure),'inferred_increment_min'] for s in samples])*1000
        mx=np.array([u.loc[(s,closure),'inferred_increment_max'] for s in samples])*1000
        ax.errorbar(x+offset,inf,yerr=np.vstack([inf-mn,mx-inf]),fmt=marker,color=color,mfc='white',ms=4,capsize=2,lw=.8,label=f'{label}: inferred')
        ax.plot(x+offset,rr.predicted_fracture_increment_S_per_m*1000,marker=marker,ls='none',color=color,ms=4,label=f'{label}: predicted')
    ax.axhline(0,color=BLACK,lw=.65,ls=':')
    ax.set(xticks=x,xticklabels=samples,ylabel=r'Increment (mS m$^{-1}$)',ylim=(-25,50),title='Background-corrected increments')
    ax.legend(frameon=False,loc='upper left',ncol=2,columnspacing=.6,handletextpad=.35,labelspacing=.3,fontsize=6.6)
    ax=axs[1,0]
    for s,color,mk in zip(samples,[BLUE,RED,BLACK],['o','^','s']):
        v=[]
        for scenario in ['aperture_0','aperture_0.5','nominal']:
            r=c[(c.scenario==scenario)&(c['sample']==s)&(c.mesh_level==1)].iloc[0]
            v.append(100*r.predicted_S_per_m/r.measured_S_per_m)
        ax.plot([0,.5,1],v,'-',color=color,marker=mk,mfc='white',ms=4,label=s)
    ax.axhline(100,color=BLACK,ls=':',lw=.65)
    ax.set(xlabel=r'Aperture / nominal aperture',ylabel='Predicted / measured (%)',xticks=[0,.5,1],ylim=(98,107),title='Porosity-partition sensitivity')
    ax.legend(frameon=False,loc='upper right')
    ax=axs[1,1]
    for s,color,mk in zip(samples,[BLUE,RED,BLACK],['o','^','s']):
        gr=c[(c.scenario=='nominal')&(c['sample']==s)].sort_values('mesh_level')
        val=100*(gr.predicted_S_per_m/gr.predicted_S_per_m.iloc[-1]-1)
        ax.plot(gr.mesh_level+1,val,'-',color=color,marker=mk,mfc='white',ms=4,label=s)
    ax.set(xticks=[1,2,3,4],xticklabels=['G1','G2','G3','G4'],ylabel=r'Change from G4 (%)',xlabel='Mesh refinement',ylim=(-.006,.16),title='Total-conductivity convergence')
    ax.legend(frameon=False,loc='upper right')
    for label,ax in zip('abcd',axs.flat):
        ax.text(-.18,1.04,f'({label})',transform=ax.transAxes,fontsize=10,fontstyle='italic')
        ax.tick_params(which='both',length=3,pad=3)
    fig.savefig(ROOT/'VALIDATION_FIGURE.pdf')
    fig.savefig(ROOT/'VALIDATION_FIGURE.svg')
    fig.savefig(ROOT/'VALIDATION_FIGURE.png',dpi=600)
    plt.close(fig)
    # Unrounded numerical tables remain the source of every displayed result.
    metrics=dict(holdout_I45_prediction=float(h.predicted_S_per_m),holdout_I45_error_percent=float(100*(h.predicted_S_per_m/h.observed_S_per_m-1)),
       nominal_RMSE_S_per_m=float(np.sqrt(np.mean((nominal.predicted_S_per_m-nominal.measured_S_per_m)**2))),
       local_RMSE_S_per_m=float(np.sqrt(np.mean((local.predicted_S_per_m-local.measured_S_per_m)**2))),
       summary=summary.to_dict(orient='records'),comparison_rows=len(c),
       largest_final_total_grid_change_percent=float(pd.DataFrame(conv).query('level==3').total_change_from_previous_percent.abs().max()),
       largest_final_increment_grid_change_percent=float(pd.DataFrame(conv).query('level==3').increment_change_from_previous_percent.abs().max()))
    (ROOT/'RESULT_SUMMARY.json').write_text(json.dumps(metrics,indent=2),encoding='utf-8')
    env={'python':platform.python_version(),'platform':platform.platform()}
    for name in ['numpy','scipy','pandas','matplotlib','gmsh','pyamg']:
        env[name]=importlib.metadata.version(name)
    (ROOT/'ENVIRONMENT.json').write_text(json.dumps(env,indent=2),encoding='utf-8')
    (ROOT/'requirements.txt').write_text('\n'.join(f'{k}=={v}' for k,v in env.items() if k not in ['python','platform'])+'\n',encoding='utf-8')
    print(json.dumps(metrics,indent=2))

if __name__=='__main__':main()
