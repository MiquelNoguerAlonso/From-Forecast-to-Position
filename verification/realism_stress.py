"""Random-execution moment and optimization checks, using synthetic laws."""
from pathlib import Path
from fractions import Fraction as F
import json
import numpy as np
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[1]

def independent_max(fun,lo,hi):
 candidates=[lo,hi,0.]
 for left,right in [(lo,min(0.,hi)),(max(0.,lo),hi)]:
  if right>left:
   res=minimize_scalar(lambda u:-fun(u),bounds=(left,right),method='bounded',options={'xatol':1e-13})
   assert res.success
   candidates.append(res.x)
 return max(candidates,key=fun)

def soft(x,c):return np.sign(x)*max(abs(x)-c,0.)

def main():
 rng=np.random.default_rng(25092026);identity_error=0.;optimizer_error=0.;gain_error=0.
 for _ in range(240):
  weights=rng.dirichlet(np.ones(8));phi=rng.random(8);returns=rng.normal(0,.4,8)
  b=rng.uniform(-1,1);a=rng.uniform(.05,2);c=rng.uniform(0,.1);k=rng.uniform(0,.3)
  lo=-rng.uniform(.01,2);hi=rng.uniform(.01,2);eps=rng.uniform(0,.1)
  p=weights@phi;s2=weights@(phi**2);mf=weights@(phi*returns);h=a*s2+k
  def direct(u):
   final=b+phi*u
   return weights@((final-b)*returns-.5*a*(final**2-b*b)-c*np.abs(phi*u))-.5*k*u*u
  for u in [lo,hi,.2*lo,.3*hi]:
   formula=u*(mf-a*b*p)-c*p*abs(u)-.5*h*u*u
   identity_error=max(identity_error,abs(direct(u)-formula))
  for uncertainty in [0.,eps]:
   # The robust interval is checked by its two linear endpoint perturbations.
   fun=lambda u:min(direct(u)+uncertainty*u,direct(u)-uncertainty*u)
   closed=float(np.clip(soft(mf-a*b*p,c*p+uncertainty)/h,lo,hi))
   opt=independent_max(fun,lo,hi)
   optimizer_error=max(optimizer_error,abs(closed-opt));gain_error=max(gain_error,abs(fun(closed)-fun(opt)))
   assert fun(closed)>=-1e-12
 assert identity_error<1e-12 and gain_error<1e-11 and optimizer_error<2e-6
 vector_error=0.
 for _ in range(40):
  n=3;weights=rng.dirichlet(np.ones(7));fills=rng.random((7,n));ret=rng.normal(size=(7,n));b=rng.normal(size=n);u=rng.normal(size=n)
  raw=rng.normal(size=(n,n));risk=raw.T@raw+np.eye(n);charge=.2*np.eye(n)
  h=charge.copy();eta=np.zeros(n);direct=0.
  for w,f,r in zip(weights,fills,ret):
   d=np.diag(f);h+=w*d.T@risk@d;eta+=w*d.T@(r-risk@b)
   final=b+d@u
   direct+=w*((final-b)@r-.5*(final@risk@final-b@risk@b))
  direct-=.5*u@charge@u
  vector_error=max(vector_error,abs(direct-(u@eta-.5*u@h@u)))
 assert vector_error<1e-11
 a,c=F(1,10),F(1,100);p=s2=F(1,2);h=a*s2;baseline_attempt=F(3,10)
 records=[]
 for theta in [F(0),F(3,10),F(1)]:
  mf=F(2,100)-F(5,100)*theta;u=max(F(0),(mf-c*p)/h)
  g=baseline_attempt*mf-c*p*baseline_attempt-h*baseline_attempt**2/2
  assert u==max(F(0),F(3,10)-theta)
  records.append((theta,mf,u,g))
 assert records[0][3]==F(225,100000) and records[-1][3]==F(-1275,100000)
 labels=['Independent','Mixture, $\\theta=0.30$','Adverse coupling']
 (ROOT/'tables/realism_rows.tex').write_text(''.join(f'{label} & {float(mf):.3f} & {float(u):.2f} & {float(g):.5f} \\\\\n' for label,(_,mf,u,g) in zip(labels,records)))
 result={'synthetic_only':True,'scalar_joint_laws':240,'optimizer_comparisons':480,'vector_joint_laws':40,'max_scalar_identity_error':identity_error,'max_optimizer_position_error':optimizer_error,'max_optimizer_value_error':gain_error,'max_vector_identity_error':vector_error,'independent_optimal_attempt':.3,'independent_gain':.00225,'adverse_gain_of_original_attempt':-.01275,'adverse_optimal_buy_attempt':0.,'stop_buying_theta':.3}
 (ROOT/'verification/realism_stress.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
