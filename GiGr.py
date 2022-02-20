# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

Test change.
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LogNorm




#CMAP="Spectral"
#CMAP="jet"
#CMAP="seismic"
CMAP="jet"


alldata=np.loadtxt(r"C:\Users\Egecan\Desktop\GitHub\HarmonicMeasurements\J_short_tm_Gr_Gi.txt")
# alldata=np.loadtxt(r"C:\Users\Egecan\Desktop\GitHub\HarmonicMeasurements\J_tm_Gr_Gi.txt")

Jt=alldata.T[0]
tmt=alldata.T[1]
Gr=alldata.T[2]
Gi=alldata.T[3]

x=np.unique(Jt)
y=np.unique(tmt)

X,Y=np.meshgrid(x,y)
GR=(Gr.reshape(len(x),len(y)).T)*(2/np.pi**2)
GI=(Gi.reshape(len(x),len(y)).T)*(2/np.pi**2)
#fig,ax=plt.subplots()
##plt.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
#plt.pcolormesh(X,Y,GR,cmap=CMAP,shading="nearest")
#plt.colorbar()
#ax.set_ylabel("$t_{m}/t$",size=14)
#ax.set_xlabel("$J/t$",size=14)
#ax.set_title("$G_{R}$", size=15)
#
#
#
#
#fig,ax=plt.subplots()
##plt.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
#plt.pcolormesh(X,Y,GI,cmap=CMAP,shading="nearest")
#plt.colorbar()
#ax.set_ylabel("$t_{m}/t$",size=14)
#ax.set_xlabel("$J/t$",size=14)
#ax.set_title("$G_{i}$", size=15)

fig,[ax1,ax2]=plt.subplots(2,1)
#im1=ax1.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
#im2=ax2.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
#
im1=ax1.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=(np.min(GI)),cmap=CMAP,shading="nearest")
plt.close()
divnorm=colors.TwoSlopeNorm(vmin=np.min(GI), vcenter=0., vmax=(np.max(GR)))

fig,[ax1,ax2]=plt.subplots(2,1)
#im1=ax1.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
#im2=ax2.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
#
# im1=ax1.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=(np.min(GI)),cmap=CMAP,shading="nearest")
# im2=ax2.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=(np.min(GI)),cmap=CMAP,shading="nearest")
# #
#im1=ax1.pcolormesh(X,Y,GI,norm=divnorm,cmap=CMAP,shading="nearest")
#im2=ax2.pcolormesh(X,Y,GR,norm=divnorm,cmap=CMAP,shading="nearest")

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.125, 0.05, 0.75])
fig.colorbar(im1, cax=cbar_ax)
ax1.contour(X,Y,np.abs(GI)/np.abs(GR),antialiased=True,linestyles=["dashed"],linewitdhs=[1],levels=[1],colors=["black"],alpha=0.9)
ax1.set_ylabel("$t_{m}/t$",size=14)
ax2.contour(X,Y,np.abs(GI)/np.abs(GR),antialiased=True,linestyles=["dashed"],linewitdhs=[1],levels=[1],colors=["black"],alpha=0.9)
ax2.set_xlabel("$J/t$",size=14)
ax2.set_ylabel("$t_{m}/t$",size=14)
# plt.savefig(r"C:\Users\Egecan\Desktop\GiGr.pdf")
# plt.savefig(r"C:\Users\Egecan\Desktop\GiGr.png",dpi=600)
# divnorm=colors.TwoSlopeNorm(vmin=np.min(GI), vcenter=0., vmax=(np.max(GR)))
# CMAP="seismic"
# fig,[ax1,ax2]=plt.subplots(2,1)
# #im1=ax1.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
# #im2=ax2.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
# #
# im1=ax1.pcolormesh(X,Y,GI,norm=divnorm,cmap=CMAP,shading="nearest")
# im2=ax2.pcolormesh(X,Y,GR,norm=divnorm,cmap=CMAP,shading="nearest")
# #
# #im1=ax1.pcolormesh(X,Y,GI,norm=divnorm,cmap=CMAP,shading="nearest")
# #im2=ax2.pcolormesh(X,Y,GR,norm=divnorm,cmap=CMAP,shading="nearest")

# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
# fig.colorbar(im1, cax=cbar_ax)

# ax1.set_ylabel("$t_{m}/t$",size=14)
# ax2.set_xlabel("$J/t$",size=14)
# ax2.set_ylabel("$t_{m}/t$",size=14)



CMAP="jet"

fig,ax=plt.subplots(figsize=(7,3.5))

im1=plt.pcolormesh(X,Y,np.abs(GI)/np.abs(GR),cmap=CMAP,norm=colors.LogNorm(0.1, 1000),shading="nearest")

plt.colorbar(im1,fraction=0.046, pad=0.04)
plt.contour(X,Y,np.abs(GI)/np.abs(GR),antialiased=True,linestyles=["dashed"],linewitdhs=[1],levels=[1],colors=["black"],alpha=0.9)
ax.set_ylabel("$t_{m}/t$",size=14)
ax.set_xlabel("$J/t$",size=14)
ax.set_title("$|G_{i}|/|G_{r}|$", size=15)
plt.tight_layout()
# plt.savefig(r"C:\Users\Egecan\Desktop\GiGrRatio.pdf")
# plt.savefig(r"C:\Users\Egecan\Desktop\GiGrRatio.png",dpi=600)
