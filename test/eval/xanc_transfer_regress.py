import numpy as np, tskit, os, json
cell=os.environ['CELL']
p=json.load(open(cell+'/case.params.json')); ntr=p['n_indiv']; nA=p.get('n_train_A',ntr) or ntr
ts=tskit.load(cell+'/case.full.trees'); n=ts.num_individuals
G=ts.genotype_matrix()
fA=G[:,0:2*nA].mean(1); fB=G[:,2*ntr:2*n].mean(1)
mt=np.asarray(ts.tables.mutations.time); ms=np.asarray(ts.tables.mutations.site)
sa=np.full(ts.num_sites,np.nan)
for k in range(ts.num_mutations):
    s=ms[k]
    if np.isnan(sa[s]): sa[s]=mt[k]
seg=(fA>0)&(fA<1)&(~np.isnan(sa))
fA=fA[seg]; fB=fB[seg]; age=sa[seg]
T=(fB*(1-fB))/(fA*(1-fA))                 # transfer prerequisite: variance ratio B/A
lage=np.log10(age+1)
zc=lambda x:(x-x.mean())/(x.std() or 1)
zfA, zage, zT = zc(fA), zc(lage), zc(T)
def pcorr(a,b,ctrl):                      # partial corr(a,b | ctrl)
    ra=a-np.polyval(np.polyfit(ctrl,a,1),ctrl); rb=b-np.polyval(np.polyfit(ctrl,b,1),ctrl)
    return np.corrcoef(ra,rb)[0,1]
pc_age=pcorr(zage,zT,zfA)                 # age vs T, controlling freq
pc_frq=pcorr(zfA,zT,zage)                 # freq vs T, controlling age
X=np.column_stack([np.ones_like(zT),zfA]); bf,*_=np.linalg.lstsq(X,zT,rcond=None)
r2f=1-np.var(zT-X@bf)/np.var(zT)
X2=np.column_stack([np.ones_like(zT),zfA,zage]); b2,*_=np.linalg.lstsq(X2,zT,rcond=None)
r2fa=1-np.var(zT-X2@b2)/np.var(zT)
print('%s pc_age=%.4f pc_freq=%.4f r2_freq=%.4f r2_freq+age=%.4f dR2=%.4f nvar=%d'%(
  os.path.basename(cell),pc_age,pc_frq,r2f,r2fa,r2fa-r2f,seg.sum()))
