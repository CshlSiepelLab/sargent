import numpy as np, tskit, os, json
cell=os.environ['CELL']
p=json.load(open(cell+'/case.params.json')); ntr=p['n_indiv']; nte=p['n_test']
ts=tskit.load(cell+'/case.full.trees'); n=ts.num_individuals
y=np.loadtxt(cell+'/case.full.phen'); y=(y-y.mean())/(y.std() or 1)
G=ts.genotype_matrix(); dos=(G[:,0::2]+G[:,1::2]).astype(float)
af=dos.mean(1)/2.0; maf=np.minimum(af,1-af)
mu=dos.mean(1,keepdims=True); sd=dos.std(1,keepdims=True); keep=sd[:,0]>0
Z=np.zeros_like(dos); Z[keep]=(dos[keep]-mu[keep])/sd[keep]
mt=np.asarray(ts.tables.mutations.time); ms=np.asarray(ts.tables.mutations.site)
sa=np.full(ts.num_sites,np.nan)
for k in range(ts.num_mutations):
    s=ms[k]
    if np.isnan(sa[s]): sa[s]=mt[k]
amed=np.nanmedian(sa[keep])
old=(sa>amed)&keep; young=(sa<=amed)&keep; common=(maf>=0.05)&keep; rare=(maf<0.05)&keep
def grm(m):
    Zb=Z[m]; K=(Zb.T@Zb)/max(1,m.sum()); return K/(np.mean(np.diag(K)) or 1)
def r2(pred,act):
    a=act-act.mean(); return float(np.corrcoef(pred,a)[0,1]**2) if pred.std()>1e-9 else 0.0
def blup(K,tr,te,h2=0.5):
    Ktr=K[np.ix_(tr,tr)]; yt=y[tr]-y[tr].mean(); lam=(1-h2)/h2
    sol=np.linalg.solve(Ktr+lam*np.eye(len(tr)), yt); return K[np.ix_(te,tr)]@sol
A=np.arange(ntr); B=np.arange(ntr,n)
def cross(K): return r2(blup(K,A,B), y[B])
def within(K):
    rng=np.random.default_rng(0); idx=rng.permutation(A); pr=np.zeros(ntr)
    for f in np.array_split(idx,5):
        tr=np.setdiff1d(A,f); pr[f]=blup(K,tr,f)
    return r2(pr,y[A])
S={'baseline':keep,'old':old,'young':young,'rare':rare,'common':common,
   'old_rare':old&rare,'young_rare':young&rare,'old_common':old&common,'young_common':young&common}
with open(cell+'/portability.tsv','w') as f:
    for nm,m in S.items():
        if m.sum()<10: f.write('%s\t%d\tNA\tNA\n'%(nm,m.sum())); continue
        K=grm(m); f.write('%s\t%d\t%.4f\t%.4f\n'%(nm,int(m.sum()),within(K),cross(K)))
