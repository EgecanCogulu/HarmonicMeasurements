# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt

CMAP="seismic"

alldata=np.loadtxt(r"C:\Users\Egecan\Desktop\J_tm_Gr_Gi.txt")
Jt=alldata.T[0]
tmt=alldata.T[1]
Gr=alldata.T[2]
Gi=alldata.T[3]

x=np.unique(Jt)
y=np.unique(tmt)

X,Y=np.meshgrid(y,x)
GR=Gr.reshape(len(x),len(y))

fig,ax=plt.subplots()
plt.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
plt.colorbar()
ax.set_ylabel("$J/J_{m}$",size=14)
ax.set_xlabel("$J/t$",size=14)

GI=Gi.reshape(len(x),len(y))


fig,ax=plt.subplots()
plt.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
plt.colorbar()
ax.set_ylabel("$J/J_{m}$",size=14)
ax.set_xlabel("$J/t$",size=14)


fig,[ax1,ax2]=plt.subplots(2,1)
im1=ax1.pcolormesh(Y,X,GI,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
im2=ax2.pcolormesh(Y,X,GR,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
fig.colorbar(im1, cax=cbar_ax)

ax1.set_ylabel("$t_{m}/t$",size=14)
ax2.set_xlabel("$J/t$",size=14)
ax2.set_ylabel("$t_{m}/t$",size=14)
