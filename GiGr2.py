# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 14:01:35 2021

@author: Egecan
"""

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
CMAP="RdGy_r"

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

ticklabelsize=16

def plotPNG(show=True,save=False):
    fig,ax=plt.subplots(figsize=(5,4))
    plt.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=(np.min(GI)),cmap=CMAP,shading="nearest")
    plt.contour(X,Y,np.abs(GI)/np.abs(GR),antialiased=True,linestyles=["dashed"],linewitdhs=[1],levels=[1],colors=["black"],alpha=0.9)
    # plt.colorbar()
    ax.set_xticks([0,0.5,1])
    ax.set_yticks([0,0.5,1])
    ax.set_xticklabels([0,0.5,1],fontsize=ticklabelsize)
    ax.set_yticklabels([0,0.5,1],fontsize=ticklabelsize)
    # ax.set_ylabel("$t_{m}/t$",size=14)
    # ax.set_xlabel("$J/t$",size=14)
    # ax.set_title("$G_{R}$", size=15)
    plt.tight_layout()

    if save==True:
        plt.savefig(r"C:\Users\Egecan\Desktop\Gr.png",dpi=600)

    if show==False:
        plt.close()
        
    fig,ax=plt.subplots(figsize=(5,4))
    plt.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=(np.min(GI)),cmap=CMAP,shading="nearest")
    plt.contour(X,Y,np.abs(GI)/np.abs(GR),antialiased=True,linestyles=["dashed"],linewitdhs=[1],levels=[1],colors=["black"],alpha=0.9)

    # ax.set_ylabel("$t_{m}/t$",size=14)
    # ax.set_xlabel("$J/t$",size=14)
    # ax.set_title("$G_{i}$", size=15)
    ax.set_xticks([0,0.5,1])
    ax.set_yticks([0,0.5,1])
    ax.set_xticklabels([0,0.5,1],fontsize=ticklabelsize)
    ax.set_yticklabels([0,0.5,1],fontsize=ticklabelsize)
    plt.tight_layout()
    if save==True:
        plt.savefig(r"C:\Users\Egecan\Desktop\Gi.png",dpi=600)
    if show==False:
        plt.close()
    return(None)   

def plotPDF(show=True,save=False):
    fig,ax=plt.subplots(figsize=(5,4))
    im=plt.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=(np.min(GI)),cmap=CMAP,shading="nearest",rasterized=True)
    # plt.colorbar(im)
    # ax.set_ylabel("$t_{m}/t$",size=14)
    # ax.set_xlabel("$J/t$",size=14)
    # ax.set_title("$G_{R}$", size=15)
    ax.set_xticks([0,0.25,0.5,0.75,1])
    ax.set_yticks([0,0.25,0.5,0.75,1])
    ax.set_xticklabels([0,"",0.5,"",1],fontsize=ticklabelsize)
    ax.set_yticklabels([0,"",0.5,"",1],fontsize=ticklabelsize)
    plt.tight_layout()

    if save==True:
        plt.savefig(r"C:\Users\Egecan\Desktop\Gr.pdf")

    if show==False:
        plt.close()
        
    plt.colorbar
    fig,ax=plt.subplots(figsize=(5,4))
    im=plt.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=(np.min(GI)),cmap=CMAP,shading="nearest",rasterized=True)
    ax.set_xticks([0,0.5,1])
    ax.set_yticks([0,0.5,1])
    ax.set_xticklabels([0,0.5,1],fontsize=ticklabelsize)
    ax.set_yticklabels([0,0.5,1],fontsize=ticklabelsize)
    plt.tight_layout()
    if save==True:
        plt.savefig(r"C:\Users\Egecan\Desktop\Gi.pdf")
    if show==False:
        plt.close()
    return(None)   




def ratioPDF():
       
    CMAP="jet"
    
    fig,ax=plt.subplots(figsize=(10,4))
    
    im1=plt.pcolormesh(X,Y,np.abs(GI)/np.abs(GR),cmap=CMAP,norm=colors.LogNorm(0.1, 10),shading="nearest",rasterized=False)
    
    plt.contour(X,Y,np.abs(GI)/np.abs(GR),antialiased=True,linestyles=["dashed"],linewitdhs=[1],levels=[1],colors=["black"],alpha=0.9)

    ax.set_xticks([0,0.25,0.5,0.75,1])
    ax.set_yticks([0,0.25,0.5,0.75,1])
    ax.set_xticklabels([0,"",0.5,"",1],fontsize=ticklabelsize)
    ax.set_yticklabels([0,"",0.5,"",1],fontsize=ticklabelsize)
    plt.tight_layout()

    # plt.savefig(r"C:\Users\Egecan\Desktop\GiGrRatio.pdf") 
    plt.savefig(r"C:\Users\Egecan\Desktop\GiGrRatio.png",dpi=900) 

# fig,[ax1,ax2]=plt.subplots(2,1)
# #im1=ax1.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
# #im2=ax2.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
# #
# im1=ax1.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=(np.min(GI)),cmap=CMAP,shading="nearest")
# plt.close()
# divnorm=colors.TwoSlopeNorm(vmin=np.min(GI), vcenter=0., vmax=(np.max(GR)))

# fig,[ax1,ax2]=plt.subplots(2,1)
# #im1=ax1.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
# #im2=ax2.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
# #
# # im1=ax1.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=(np.min(GI)),cmap=CMAP,shading="nearest")
# # im2=ax2.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=(np.min(GI)),cmap=CMAP,shading="nearest")
# # #
# #im1=ax1.pcolormesh(X,Y,GI,norm=divnorm,cmap=CMAP,shading="nearest")
# #im2=ax2.pcolormesh(X,Y,GR,norm=divnorm,cmap=CMAP,shading="nearest")

# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.85, 0.125, 0.05, 0.75])
# fig.colorbar(im1, cax=cbar_ax)
# ax1.contour(X,Y,np.abs(GI)/np.abs(GR),antialiased=True,linestyles=["dashed"],linewitdhs=[1],levels=[1],colors=["black"],alpha=0.9)
# ax1.set_ylabel("$t_{m}/t$",size=14)
# ax2.contour(X,Y,np.abs(GI)/np.abs(GR),antialiased=True,linestyles=["dashed"],linewitdhs=[1],levels=[1],colors=["black"],alpha=0.9)
# ax2.set_xlabel("$J/t$",size=14)
# ax2.set_ylabel("$t_{m}/t$",size=14)
# # plt.savefig(r"C:\Users\Egecan\Desktop\GiGr.pdf")
# # plt.savefig(r"C:\Users\Egecan\Desktop\GiGr.png",dpi=600)
# # divnorm=colors.TwoSlopeNorm(vmin=np.min(GI), vcenter=0., vmax=(np.max(GR)))
# # CMAP="seismic"
# # fig,[ax1,ax2]=plt.subplots(2,1)
# # #im1=ax1.pcolormesh(X,Y,GI,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
# # #im2=ax2.pcolormesh(X,Y,GR,vmax=(np.max(GR)),vmin=-1*(np.max(GR)),cmap=CMAP,shading="nearest")
# # #
# # im1=ax1.pcolormesh(X,Y,GI,norm=divnorm,cmap=CMAP,shading="nearest")
# # im2=ax2.pcolormesh(X,Y,GR,norm=divnorm,cmap=CMAP,shading="nearest")
# # #
# # #im1=ax1.pcolormesh(X,Y,GI,norm=divnorm,cmap=CMAP,shading="nearest")
# # #im2=ax2.pcolormesh(X,Y,GR,norm=divnorm,cmap=CMAP,shading="nearest")

# # fig.subplots_adjust(right=0.8)
# # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
# # fig.colorbar(im1, cax=cbar_ax)

# # ax1.set_ylabel("$t_{m}/t$",size=14)
# # ax2.set_xlabel("$J/t$",size=14)
# # ax2.set_ylabel("$t_{m}/t$",size=14)


