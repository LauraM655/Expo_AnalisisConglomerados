import numpy as np, matplotlib
from scipy.spatial.distance import pdist
import os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.datasets import make_moons, make_blobs, make_circles
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ari
V,R,A,Y,AC,VC="#94B43B","#A61C31","#35408F","#F3B73C","#2C6DAF","#40B08F"
COL=[V,R,A,Y,AC,VC]
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False})
def sv(n,w=None): plt.savefig(f"graphics/{n}.png",dpi=250,bbox_inches="tight",transparent=False,facecolor="white"); plt.close()
def hc(X,m,k): return fcluster(linkage(X,m),k,"maxclust")
def sc(ax,X,l,t):
    for i,c in enumerate(np.unique(l)): ax.scatter(*X[l==c].T,s=14,color=COL[i%6],alpha=.8)
    ax.set_title(t,fontsize=10); ax.set_xticks([]); ax.set_yticks([])

# 1. Criterios de enlace (esquema)
rng=np.random.default_rng(3)
Ca=rng.normal([0,0],.55,(14,2)); Cb=rng.normal([4,1.2],.55,(14,2))
fig,axs=plt.subplots(1,4,figsize=(13,3.4))
D=np.linalg.norm(Ca[:,None]-Cb[None],axis=2)
i,j=np.unravel_index(D.argmin(),D.shape); i2,j2=np.unravel_index(D.argmax(),D.shape)
ca,cb=Ca.mean(0),Cb.mean(0)
tit=["Simple: min","Completo: max","Promedio: media de\ntodos los pares","Ward: aumento de SCW\n(centroides)"]
for ax,t in zip(axs,tit):
    ax.scatter(*Ca.T,color=A,s=22);ax.scatter(*Cb.T,color=R,s=22);ax.set_title(t,fontsize=10.5);ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect("equal")
axs[0].plot([Ca[i,0],Cb[j,0]],[Ca[i,1],Cb[j,1]],color=V,lw=3)
axs[1].plot([Ca[i2,0],Cb[j2,0]],[Ca[i2,1],Cb[j2,1]],color=V,lw=3)
for a in Ca[::2]:
    for b in Cb[::2]: axs[2].plot([a[0],b[0]],[a[1],b[1]],color=V,lw=.5,alpha=.6)
axs[3].scatter(*ca,marker="*",s=250,color=Y,edgecolor="k");axs[3].scatter(*cb,marker="*",s=250,color=Y,edgecolor="k")
axs[3].plot([ca[0],cb[0]],[ca[1],cb[1]],color=V,lw=3,ls="--")
axs[3].text((ca[0]+cb[0])/2,(ca[1]+cb[1])/2+.5,r"$\frac{n_hn_k}{n_h+n_k}\|\bar X_h-\bar X_k\|^2$",ha="center",fontsize=9)
plt.tight_layout();sv("A_criterios_enlace")

# 2. Donde falla cada enlace (como la figura de k-means)
Xm,_=make_moons(300,noise=.06,random_state=1)
Xb=np.vstack([make_blobs(150,centers=[[0,0]],cluster_std=.5,random_state=1)[0],make_blobs(150,centers=[[5,0]],cluster_std=.5,random_state=2)[0],np.c_[np.linspace(.8,4.2,12),np.random.default_rng(0).normal(0,.05,12)]])
Xd=np.vstack([make_blobs(150,centers=[[0,0]],cluster_std=1.8,random_state=3)[0],make_blobs(40,centers=[[6.5,0]],cluster_std=.25,random_state=4)[0]])
Xc,_=make_circles(300,factor=.3,noise=.04,random_state=0)
Xe=np.vstack([make_blobs(100,centers=[[0,0]],cluster_std=.5,random_state=5)[0],make_blobs(100,centers=[[4,0]],cluster_std=.5,random_state=6)[0],[[2,5],[-3,4.5],[8,-4]]])
Xg=make_blobs(300,centers=[[0,0],[4,0],[2,3.5]],cluster_std=.6,random_state=7)[0]
cases=[(Xg,"ward",3,"Ward: grupos compactos (bien)"),
 (Xm,"single",2,"Simple: formas no convexas (bien)"),
 (Xm,"ward",2,"Ward: formas no convexas (mal)"),
 (Xb,"single",2,"Simple: efecto cadena por puente (mal)"),
 (Xe,"ward",3,"Ward: outliers distorsionan (mal)"),
 (Xe,"single",3,"Simple: aisla outliers como grupos (mal)"),
 (Xd,"complete",2,"Completo/Ward: densidades distintas (parcial)"),
 (Xc,"single",2,"Simple: anillos concentricos (bien)"),
 (Xc,"ward",2,"Ward: anillos concentricos (mal)")]
fig,axs=plt.subplots(3,3,figsize=(11,10))
for ax,(X,m,k,t) in zip(axs.flat,cases): sc(ax,X,hc(X,m,k),t)
plt.tight_layout();sv("B_cuando_falla")

# 3. Outliers y dendrogramas
fig,axs=plt.subplots(1,3,figsize=(13,3.6))
Xo=Xe[:200]; Xo2=Xe
sc(axs[0],Xo2,hc(Xo2,"ward",2),"Datos con 3 outliers (Ward, k=2)")
for ax,X,t in zip(axs[1:],[Xo,Xo2],["Dendrograma Ward sin outliers","Dendrograma Ward con outliers"]):
    dendrogram(linkage(X,"ward"),ax=ax,no_labels=True,color_threshold=0,above_threshold_color=A);ax.set_title(t,fontsize=10);ax.set_yticks([])
plt.tight_layout();sv("C_outliers")

# 4. Estabilidad: jerarquico determinista vs k-means semillas / remuestreo
X=make_blobs(300,centers=[[0,0],[3,0],[1.5,2.6],[6,3]],cluster_std=1.1,random_state=11)[0]
ref_h=hc(X,"ward",4); ref_k=KMeans(4,n_init=1,random_state=0).fit_predict(X)
a_hs=[ari(ref_h,hc(X,"ward",4)) for _ in range(30)]
a_ks=[ari(ref_k,KMeans(4,n_init=1,random_state=s).fit_predict(X)) for s in range(1,31)]
rng=np.random.default_rng(0);bh=[];bk=[]
for _ in range(100):
    idx=rng.choice(300,300,replace=True)
    bh.append(ari(ref_h[idx],hc(X[idx],"ward",4)));bk.append(ari(ref_k[idx],KMeans(4,n_init=10,random_state=0).fit_predict(X[idx])))
fig,axs=plt.subplots(1,2,figsize=(11,3.8))
axs[0].boxplot([a_hs,a_ks],labels=["Ward\n(repetir codigo)","K-means\n(1 semilla distinta)"],patch_artist=True,boxprops=dict(facecolor=VC))
axs[0].set_ylabel("ARI vs. corrida de referencia");axs[0].set_title("Repetir el codigo",fontsize=10.5);axs[0].set_ylim(-.05,1.05)
axs[1].boxplot([bh,bk],labels=["Ward","K-means\n(n_init=10)"],patch_artist=True,boxprops=dict(facecolor=Y))
axs[1].set_title("Remuestreo bootstrap (otra muestra)",fontsize=10.5);axs[1].set_ylim(-.05,1.05)
plt.tight_layout();sv("D_estabilidad")
print("repetir: ward",np.mean(a_hs),"kmeans",np.mean(a_ks),"| bootstrap ward",np.mean(bh),"kmeans",np.mean(bk))

# 5. Aglomerativo vs divisivo + dendrograma con corte
fig,ax=plt.subplots(1,2,figsize=(11,3.6))
X5=make_blobs(40,centers=[[0,0],[3,0],[1.5,2.5]],cluster_std=.5,random_state=2)[0]
L=linkage(X5,"ward");dendrogram(L,ax=ax[0],no_labels=True,color_threshold=.7*L[-1,2],above_threshold_color="gray")
ax[0].axhline(.7*L[-1,2],color=R,ls="--");ax[0].set_title("Dendrograma y linea de corte (k=3)",fontsize=10.5)
ax[1].axis("off");ax[1].text(0,.9,"Aglomerativo (abajo-arriba):\n n grupos  ->  1 grupo\n Es el que usamos (hclust)\n\nDivisivo (arriba-abajo):\n 1 grupo  ->  n grupos\n Mas costoso (diana, bisecting)",va="top",fontsize=11)
plt.tight_layout();sv("E_dendro_corte")

