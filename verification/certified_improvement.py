"""Robust positions and a normal-mixture sequential gain certificate."""
from pathlib import Path
import json
import math
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]


def position(m,b,a,c,e,L,U):
    signal=m-a*b
    d=np.sign(signal)*max(0.,abs(signal)-c-e)/a
    return float(b+np.clip(d,L-b,U-b))


def boundary(V,alpha=.05,v0=1.):
    return np.sqrt((V+v0)*np.log((V+v0)/(v0*alpha**2)))


def main():
    rng=np.random.default_rng(5259);max_scalar_error=0.;cases=0
    for _ in range(250):
        b=float(rng.uniform(-.3,.3));m=float(rng.uniform(-.8,.8))
        a=float(rng.uniform(.2,3));c=float(rng.uniform(0,.12));e=float(rng.uniform(0,.15))
        L,U=-.4,.4
        objective=lambda w:(m-a*b)*(w-b)-.5*a*(w-b)**2-(c+e)*abs(w-b)
        w=position(m,b,a,c,e,L,U)
        independent=minimize_scalar(lambda q:-objective(q),bounds=(L,U),method='bounded',options={'xatol':1e-13})
        independent_value=max(objective(b),objective(L),objective(U),-independent.fun)
        max_scalar_error=max(max_scalar_error,abs(objective(w)-independent_value))
        assert abs(objective(w)-independent_value)<1e-10
        endpoints=[(mu-a*b)*(w-b)-.5*a*(w-b)**2-c*abs(w-b) for mu in [m-e,m+e]]
        assert abs(min(endpoints)-objective(w))<1e-12
        assert objective(w)>=.5*a*(w-b)**2-1e-12
        cases+=1
    # Independent constrained minimization of a 3D ellipsoid support problem.
    support_error=0.; solver_status_warnings=0
    for _ in range(30):
        S=rng.normal(size=(3,3));d=rng.normal(size=3);rho=.2
        opt=minimize(lambda v:rho*d@S@v,np.zeros(3),jac=lambda v:rho*S.T@d,
                     constraints={'type':'ineq','fun':lambda v:1-v@v,'jac':lambda v:-2*v},
                     method='SLSQP',options={'ftol':1e-9,'maxiter':500})
        # A line-search status is not an accuracy certificate. Project the returned
        # iterate to the feasible ball and check its objective independently.
        solver_status_warnings += int(not opt.success)
        feasible=opt.x/max(1.,np.linalg.norm(opt.x))
        assert feasible@feasible <= 1+1e-12
        err=abs(rho*d@S@feasible+rho*np.linalg.norm(S.T@d));support_error=max(support_error,err)
        assert err<1e-8
    # Numerical integration validates the mixture identity, not the same formula twice.
    integral_error=0.
    for M,V,v0 in [(0.,0.,1.),(.4,.3,.7),(-.8,1.2,.2),(1.1,2.,.4)]:
        integrand=lambda lam:math.sqrt(v0/(2*math.pi))*math.exp(lam*M-.5*(V+v0)*lam*lam)
        numerical=quad(integrand,-np.inf,np.inf,epsabs=1e-12)[0]
        analytic=math.sqrt(v0/(V+v0))*math.exp(M*M/(2*(V+v0)))
        integral_error=max(integral_error,abs(numerical-analytic))
        assert abs(numerical-analytic)<1e-11
    rows=[]
    for e in [0.,.02,.05,.08,.10,.15]:
        w=position(.12,0.,1.,.02,e,-.2,.2)
        gain=.12*w-.5*w*w-(.02+e)*abs(w)
        rows.append({'mean_radius':e,'nominal_position':.1,'robust_position':w,'guaranteed_gain':gain})
    horizon=800;paths=1000;alpha=.05;v0=.1;sigma=.2
    tt=np.arange(horizon);m=.12*np.sin(tt/35)+.07*np.cos(tt/17)
    e=np.full(horizon,.025);mu=m+.8*e*np.sin(tt/23)
    w=np.array([position(mm,0.,1.,.02,ee,-.2,.2) for mm,ee in zip(m,e)])
    charge=.5*w*w+.02*np.abs(w)
    g=m*w-e*np.abs(w)-charge
    expected=mu*w-charge
    assert np.all(expected>=g-1e-14) and np.all(g>=-1e-14)
    noise=rng.normal(0,sigma,size=(paths,horizon))
    residual=noise*w
    M=np.cumsum(residual,axis=1);V=np.cumsum(sigma*sigma*w*w)
    allowance=boundary(V,alpha,v0)
    crossings=np.any(np.abs(M)>allowance,axis=1)
    realized=np.cumsum((mu+noise)*w-charge,axis=1)
    lower=np.cumsum(g)-allowance
    # If the noise boundary is respected, the reported lower curve must be respected.
    assert np.all(realized[~crossings]>=lower-1e-12)
    report={'status':'all_passed','data':'Synthetic Gaussian returns; penalized objective, not cash P&L.',
            'seed':5259,'scalar_optimizer_checks':cases,'maximum_scalar_error':max_scalar_error,
            'ellipsoid_optimizer_checks':30,'maximum_support_error':support_error,
            'ellipsoid_solver_status_warnings_with_verified_objective':solver_status_warnings,
            'mixture_integral_checks':4,'maximum_mixture_integral_error':integral_error,
            'horizon':horizon,'paths':paths,'alpha':alpha,'mean_coverage_failure_probability':0.,
            'noise_sigma':sigma,'mixture_v0':v0,'observed_boundary_crossings':int(crossings.sum()),
            'boundary_crossing_fraction':float(crossings.mean()),
            'cumulative_guaranteed_conditional_gain':float(g.sum()),
            'example_realized_penalized_gain':float(realized[0,-1]),
            'example_time_uniform_lower_gain':float(lower[-1]),'rows':rows}
    (ROOT/'verification/certified_improvement_results.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'tables').mkdir(exist_ok=True)
    (ROOT/'tables/certified_improvement_rows.tex').write_text(''.join(
        f"{r['mean_radius']:.2f} & {r['nominal_position']:.3f} & {r['robust_position']:.3f} & {r['guaranteed_gain']:.5f} \\\\\n" for r in rows))
    fig,ax=plt.subplots(1,2,figsize=(10.8,3.8),layout='constrained')
    ax[0].plot([r['mean_radius'] for r in rows],[r['robust_position'] for r in rows],'o-',color='#0057a6',label='Robust position')
    ax[0].axhline(.1,color='gray',ls='--',label='Nominal position')
    ax[0].set(xlabel='Conditional-mean uncertainty radius',ylabel='Position');ax[0].legend(frameon=False)
    ax[1].plot(realized[0],color='#0057a6',label='Realized penalized gain')
    ax[1].plot(lower,color='#b05e00',label='Time-uniform lower bound')
    ax[1].plot(np.cumsum(g),color='gray',ls='--',label='Conditional certificate')
    ax[1].set(xlabel='Nonoverlapping evaluation period',ylabel='Cumulative objective difference')
    ax[1].legend(frameon=False,fontsize=8)
    for a in ax:a.spines[['top','right']].set_visible(False)
    fig.savefig(ROOT/'figures/06_certified_improvement.png',dpi=220,facecolor='white');plt.close(fig)
    print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2))


if __name__=='__main__':main()
