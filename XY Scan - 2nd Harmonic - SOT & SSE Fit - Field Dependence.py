# -*- coding: utf-8 -*-
"""
Created on Sun Dec 30 23:42:13 2018

@author: Egecan Cogulu
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from os.path import dirname, abspath
from scipy.optimize import curve_fit
from matplotlib import cm
from  matplotlib.colors import rgb2hex


step_angle=5
Scan="XY"


cmap = cm.get_cmap('viridis',6)    # PiYG
colorlist=[]
for i in range(cmap.N):
    rgba = cmap(i)
    # rgb2hex accepts rgb or rgba
    colorlist.append(rgb2hex(rgba))

def func2(x, a, b, c, d,e,f,g):
    return a * np.sin(b*x+c) + d*np.sin(e*x+f)+g

def func1(x, a, b, c, g):
    return a * np.sin(b*x+c) + g


def func3(x,a,b,c,d,e,f,g):
    
    return (a*np.cos(2*b*x+c)*np.cos(b*x+c)+d*np.cos(e*x+f)+g)

def linear(x,a):
    return(a*x)

def expo1(x,a):
    
    return (a/x)

def expo2(x,a,b):
    
    return (a/x+b)

def create_datapaths(currentfolder,temps):
    datapaths=[]
    datapath=os.getcwd()+"\\"+currentfolder
    for i in temps:
        datapaths.append(datapath+i)
    return datapaths

def XYScan(phiHs,A,B):
    #Parameters
    SigmaAC=5E-5
    HAC=5E-3

    H=3
    HE=1800
    DMI=2
    Kz=1E-2
    thetaD=2*(np.pi/180)
    phiD=135*(np.pi/180)
    thetaH=90*(np.pi/180)
    # phiH=1*(np.pi/180)
    phiHs=np.asfarray(phiHs)*(np.pi/180)
    
    #Defitions of X, Y and Z unit vectors
    xvector=[1,0,0]
    yvector=[0,1,0]
    zvector=[0,0,1]
    
    
    FirstHarm=[]
    SecondHarm=[]
    for phiH in phiHs:

        #DMI Vector Defition
        
        DMIvector= np.asfarray([np.sin(thetaD)*np.cos(phiD),np.sin(thetaD)*np.sin(phiD),np.cos(thetaD)])
        
        #External Field Defition
        
        Hvector= np.asfarray([np.sin(thetaH)*np.cos(phiH),np.sin(thetaH)*np.sin(phiH),np.cos(thetaH)])
        
        # Definitions of X, Y and Z L vectors
        xLvector=np.cross(Hvector,DMIvector)
        xLvector=xLvector/np.linalg.norm(xLvector)
        zLvector = DMIvector
        yLvector = np.cross(zLvector,xLvector)
        
        #HACLx
        HACLx = HAC * np.dot(yvector,xLvector)
        
        #HexLy and HexLz
        HexLy = H * np.dot(Hvector,yLvector)
        HexLz = H * np.dot(Hvector,zLvector)
        
        #SigmaACLy and Lz
        SigmaACLy = SigmaAC*np.dot(yvector,yLvector)
        SigmaACLz = SigmaAC*np.dot(yvector,zLvector)
        
        
        dthetaoverH=np.asfarray([[0,0,0],[-1/HexLy,0,0]])
        
        
        prefactor=HE/(DMI**2*HexLy+2*HE*HexLy*(Kz+DMI**2/(2*HE))+DMI*(H**2+2*HE*Kz+DMI**2))
        dthetaoverSigma=prefactor*np.asfarray([[0,DMI+HexLy,HexLz],[0,HexLz,(DMI*HexLy+HexLz**2+2*HE*Kz+DMI**2)/HexLy]])
        
        
        DELTAtheta=dthetaoverSigma[0][1]*SigmaACLy+dthetaoverSigma[0][2]*SigmaACLz
        DELTAphi=dthetaoverH[1][0]*HACLx+dthetaoverSigma[1][1]*SigmaACLy+dthetaoverSigma[1][2]*SigmaACLz
        
        firstharm=np.dot(xLvector,xvector)*np.dot(xLvector,yvector)
        secondharm=np.dot(xLvector,xvector)*(DELTAphi*np.dot(yLvector,yvector)-DELTAtheta*np.dot(zLvector,yvector)) + ((np.dot(xLvector,yvector))*(DELTAphi*np.dot(yLvector,xvector)-DELTAtheta*np.dot(zLvector,xvector)))
        
        SecondHarm.append(A*secondharm+B*np.cos(phiH))
    
    return(np.asfarray(SecondHarm))


def read(filename):
    """reads the dat file outputs the list row by row. Returns a list"""
    x=[]
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            x.append(row)
    return (x)

def sort_average(rawdata):
    angles=rawdata[:,0]
    voltages=rawdata[:,1]
    numberofsweeps=len(voltages)/73
    
    v=np.split(voltages,numberofsweeps)
    a=np.split(angles,numberofsweeps)
    ave_voltage=[]
    total=np.zeros(73)
    for i in range(int(numberofsweeps)):
        if i%2==0:
            total=total+v[i]
        else:
            total=total+np.flip(v[i])
    ave_voltage=total/(numberofsweeps)
    
    return (ave_voltage,voltages,angles,a[0])



def plotter(datapath,savepath,temp="Unknown Temp"):
    filename=datapath +"\\V vs Angle - Y.txt"
    
    rawdata=np.asfarray(read(filename))

    (ave_voltage,voltages,angles,angle)=sort_average(rawdata)
    ave_voltage=ave_voltage-np.average(ave_voltage)
    ave_resistance=((ave_voltage)/1.5E-3)
    voltages=voltages-np.average(voltages)
    xdata = np.linspace(0, 360,361)
    # ave_voltage=np.append(ave_voltage[18:],ave_voltage[:18])
    # popt, pcov = curve_fit(func2, angle, ave_voltage,[0.1,0.017,1,0,0.017,1,0.1],bounds=( [-1,-1,-1,-1,-1,-1,-1], [100,1,10,0.5,0.018,3.14,10]))

    popt, pcov = curve_fit(XYScan, angle, ave_resistance, [-0.16,-0.16],bounds=( [-1,-1], [1,1]),maxfev=1000000)
    # initial_guess=[0.0003,0.017,-0.06,0.00001,0.017,1,0.1]
    # np.save(r"C:\Users\Egecan\Desktop\Numerical Fitting\Np Arrays\XY2ndHarm",ave_resistance)
    # popt, pcov = curve_fit(func3, angle, ave_resistance,initial_guess,bounds=( [-1,-1,-1,-1,-1,-1,-1], [100,1,10,0.5,0.018,3.14,10]),maxfev=100000)
    print(temp,np.round(popt,5))
    if save==True:
        fig, ax1 = plt.subplots()#Plotting Sine FIT and Average Data
        # ax1.plot(angle[start:start+step],ave_voltage[start:start+step],"o-",linewidth=0.5,markersize=5)
        ax1.plot(angle,ave_resistance*1E3,"o:",markersize=5,label="Data")
        ax1.plot(xdata, XYScan(xdata, *popt)*1E3, 'r-',label="Fit")#Fit
        # ax1.plot(xdata, func3(xdata, *popt)*1E3, 'r-',label="Fit")#Fit
        ax1.set_xlabel("Angle $\\alpha$ ($\degree$)",size=16)
        ax1.set_ylabel("$\\Delta R_{xy}^{2\omega}$ (m$\Omega$)",size=16)
        ax1.set_title("XY Scan - 2$^{nd}$ Harmonic Resistance",size=14)
        ax1.set_xticks([0,90,180,270,360])
        ax1.tick_params(axis='both', which='major', labelsize=15)
        # ax1.set_ylim(0.05,0.3)
        # ax1.set_ylim(0.13,0.23)
        plt.legend()
        plt.tight_layout()
        # plt.savefig(savepath+"Average - "+temp+".png",dpi=600)
        # plt.savefig(savepath+"2nd Harm - XY Scan - "+temp+".pdf",dpi=600)
        plt.close()
        
        fig, ax2 = plt.subplots()#Plotting All data and Average Data
        ax2.set_xlabel("Angle $\\alpha$ ($\degree$)",size=16)
        ax2.set_ylabel("$V_{xy}^{2\omega}$ ($\mu$V)",size=16)
        ax2.set_title("XY Scan - 2$^{nd}$ Harmonic Voltage vs. Angle",size=14)
        halflength=int(len(voltages)/2)
        ax2.plot(angles[:halflength],voltages[:halflength],"o-",markersize=5,label="Forwards")
        ax2.plot(angles[halflength:],voltages[halflength:],"o-",markersize=5,label="Backwards")
        # ax2.plot(angles[150:224],voltages[150:224],"o-",markersize=5)
        # ax2.plot(angle,ave_voltage,"o-")#Average data
        ax2.set_xticks([0,90,180,270,360])
        ax2.tick_params(axis='both', which='major', labelsize=15)
        plt.legend()
        # ax2.set_ylim(0.1,0.45)
        # ax2.set_ylim(0.10,0.22)
        plt.tight_layout()
        # plt.savefig(savepath+"All Data - "+temp+".png",dpi=300)
        
        plt.close()
    return(popt,angle,ave_voltage)
    

def call_plotter(datapath,savepath):
    temps=os.listdir(datapath)
    amplitude1={}
    wavenumber1={}
    phase1={}
    amplitude2={}
    wavenumber2={}
    phase2={}
    constant={}
    data={}
    for i in temps:
        (p2p,angle,ave_voltage)=plotter(datapath+i,savepath,i)
        start=i.find('K')+3
        stop=i.find('O')
        print()
        amplitude1[float(i[start:stop])]=p2p[0]
        # wavenumber1[float(i[start:stop])]=p2p[1]   
        # phase1[float(i[start:stop])]=p2p[2]
        amplitude2[float(i[start:stop])]=p2p[1]
        # wavenumber2[float(i[start:stop])]=p2p[4]   
        # phase2[float(i[start:stop])]=p2p[5]
        # constant[float(i[start:stop])]=p2p[6]     
        data[float(i[start:stop])]=(angle,ave_voltage)
    return (amplitude1,amplitude2,data)



save=True
path=dirname(os.getcwd())

datapath=path+r"\Data\\"
savepath=path+r"\Graphs\\"

(amplitude1,amplitude2,data)=call_plotter(datapath, savepath)

# plot_field_dependence(amplitude1,wavenumber1,phase1,amplitude2,wavenumber2,phase2,constant,datapath,savepath)

lists = sorted(amplitude1.items())
field,volt1=zip(*lists)
res1=np.abs(np.asfarray(volt1)/1.5E-3)

lists = sorted(amplitude2.items())
field,volt2=zip(*lists)
res2=np.abs(np.asfarray(volt2)/1.5E-3)

field=np.asfarray(field)/10000

fitpoint=8
xdata=np.linspace(field[4],8.2,100)
# popt, pcov = curve_fit(expo1, field[fitpoint:],res1[fitpoint:],maxfev=100000)
# print(popt)
# fig, ax = plt.subplots()
# ax.plot(field,res1,"o:",markersize=8)
# ax.plot(xdata,expo1(xdata,*popt),color="red")
# ax.set_title("2$^{nd}$ Harmonic Amplitude vs. Magnetic Field",size=16)
# ax.set_xlabel("Magnetic Field (T)",size=15)
# ax.set_ylabel("R$_{xy}^{2\omega}$ ($\mu\Omega$)",size=15)
# ax.tick_params(axis='both', which='major', labelsize=15)
# plt.tight_layout()
# # plt.savefig(savepath+"1st Amplitude vs. Field - Fit1.png",dpi=600)
# # plt.savefig(savepath+"1st Amplitude vs. Field - Fit1.pdf",dpi=600)
# plt.close()
# xdata=np.linspace(field[4],8.1,100)
# popt, pcov = curve_fit(expo2, field[fitpoint:],res1[fitpoint:],maxfev=100000)
# print(popt)
# fig, ax = plt.subplots()
# ax.plot(field,res1,"o:",markersize=8)
# ax.plot(xdata,expo2(xdata,*popt),color="red")
# ax.set_title("2$^{nd}$ Harmonic Amplitude vs. Magnetic Field",size=16)
# ax.set_xlabel("Magnetic Field (T)",size=15)
# ax.set_ylabel("R$_{xy}^{2\omega}$ ($\mu\Omega$)",size=15)
# ax.tick_params(axis='both', which='major', labelsize=15)
# plt.tight_layout()
# # plt.savefig(savepath+"1st Amplitude vs. Field - Fit2.png",dpi=600)
# # plt.savefig(savepath+"st Amplitude vs. Field - Fit2.pdf",dpi=600)
# plt.close()


# xdata=np.arange(0.1,0.7,0.05)
# fig, ax = plt.subplots()
# ax.plot(1/field[1:],res1[1:],"o:",markersize=8)
# popt, pcov = curve_fit(linear, 1/field[-8:],res1[-8:],maxfev=100000)
# ax.set_title("2$^{nd}$ Harmonic Amplitude vs. Magnetic Field",size=16)
# ax.plot(xdata,linear(xdata,*popt),color="red")
# ax.set_xlabel("1 / Magnetic Field (T)",size=15)
# ax.set_ylabel("R$_{xy}^{2\omega}$ ($\mu\Omega$)",size=15)
# ax.tick_params(axis='both', which='major', labelsize=15)
# plt.tight_layout()
# print("1/H Popt:"+str(popt))
# # plt.savefig(savepath+"1st Amplitude vs. Field - 1overH.png",dpi=600)
# # plt.savefig(savepath+"1st Amplitude vs. Field - Different Fields.pdf",dpi=600)
# plt.close()


# fig, ax = plt.subplots()
# ax.axhline(y=0,xmin=0,xmax=1,color="Gray",alpha=0.5,linewidth=1)
# ax.plot(field,volt1,"o:",markersize=8,label="FL")
# ax.plot(field,volt2,"o:",markersize=8,label="SSE")
# ax.set_title("2$^{nd}$ Harmonic Amplitude vs. Magnetic Field",size=16)
# ax.set_xlabel("Magnetic Field (T)",size=15)
# ax.set_ylabel("R$_{xy}^{2\omega}$ ($\mu\Omega$)",size=15)
# ax.tick_params(axis='both', which='major', labelsize=15)
# plt.legend()
# plt.tight_layout()
# # plt.savefig(savepath+"1st and 2nd Amplitude vs. Field - Voltages.png",dpi=600)
# # plt.savefig(savepath+"1st and 2nd Amplitude vs. Field - Different Fields.pdf",dpi=600)
# plt.close()




# fig, ax = plt.subplots()
# ax.axhline(y=0,xmin=0,xmax=1,color="Gray",alpha=0.5,linewidth=1)
# ax.plot(field,res1,"o:",markersize=8,label="FL")
# ax.plot(field,res2,"o:",markersize=8,label="SSE")
# ax.set_title("2$^{nd}$ Harmonic Amplitude vs. Magnetic Field",size=15)
# ax.set_xlabel("External Magnetic Field (T)",size=14)
# ax.set_ylabel("R$_{xy}^{2\omega}$ ($\mu\Omega$)",size=14)
# plt.legend()
# plt.tight_layout()
# # plt.savefig(savepath+"1st and 2nd Amplitude vs. Field - Resistances.png",dpi=600)
# # plt.savefig(savepath+"1st and 2nd Amplitude vs. Field - Different Fields.pdf",dpi=600)
# plt.close()

fields=np.sort(np.asfarray(list(data.keys())))

fig, ax1 = plt.subplots()
for (index,ifield) in enumerate(fields):
    plt.plot(data[ifield][0],data[ifield][1]+14+0.3*index,label=str(int(ifield))+" Oe")
    ax1.set_xlabel("Angle $\\alpha$ ($\degree$)",size=17)

ax1.set_ylabel("$V_{xy}^{1\omega}$ ($\mu$V)",size=16)
ax1.set_title("SMR Amplitude vs. Angle",size=16)
ax1.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0,90,180,270,360],[0,90,180,270,360])
plt.legend()
plt.tight_layout()
plt.savefig(savepath+"all fields.png",dpi=600)
plt.savefig(savepath+"all fields.pdf",dpi=600)
    

legendlist=["0 T","0.5 T","0.75 T","1 T","1.5 T","2 T","2.5 T","3 T","4 T","5 T","6 T","7 T","8 T"]
fig, ax1 = plt.subplots()
for (index,ifield) in enumerate(fields[::]):
    plt.plot(data[ifield][0],data[ifield][1],label=legendlist[index],color=colorlist[index])

ax1.set_xlabel("Angle $\\alpha$ ($\degree$)",size=16)

ax1.set_ylabel("$V_{xy}^{2\omega}$ ($\mu$V)",size=16)
# ax1.set_title("2$^{nd}$ Harmonic Voltage vs. Angle",size=14)
plt.xticks([0,90,180,270,360],[0,90,180,270,360])
ax1.tick_params(axis='both', which='major', labelsize=15)
# ax1.set_ylim(0,700)
plt.legend()
plt.tight_layout()
plt.savefig(savepath+"all fields2.png",dpi=600)
plt.savefig(savepath+"all fields2.pdf",dpi=600)
plt.close()


fig, ax1 = plt.subplots()
for (index,ifield) in enumerate(fields[::]):
    plt.plot(data[ifield][0],data[ifield][1]/1.5E-6,label=legendlist[index],color=colorlist[index])

ax1.set_xlabel("Angle $\\alpha$ ($\degree$)",size=16)

# ax1.set_ylabel("$\Delta$$R_{xy}^{2\omega}$ (m$\Omega$)",size=16)
# ax1.set_title("XY Scan - 2$^{nd}$ Harmonic Resistance",size=14)
plt.xticks([0,90,180,270,360],[0,90,180,270,360])
ax1.tick_params(axis='both', which='major', labelsize=15)
# ax1.set_ylim(0,700)
# plt.legend()
plt.tight_layout()
plt.savefig(savepath+"2nd Harmonic Resistance - XY Scan.png",dpi=600)
plt.savefig(savepath+"2nd Harmonic Resistance - XY Scan.pdf",dpi=600)
# plt.savefig(savepath+"all fields2.pdf",dpi=600)
plt.close()


fig, ax = plt.subplots()
# ax.axhline(y=0,xmin=0,xmax=1,color="Gray",alpha=0.5,linewidth=1)
ax.plot(field,10*(res1*field)/160,"o:",markersize=8,label="H$_{FL}$")

ax.set_title("2$^{nd}$ Harmonic Amplitude vs. Magnetic Field",size=16)
ax.set_xlabel("Magnetic Field (T)",size=15)
ax.set_ylabel("H$_{AC}$ (Oe)",size=15)
ax.tick_params(axis='both', which='major', labelsize=15)
# plt.legend()
plt.tight_layout()
plt.savefig(savepath+"Hac.png",dpi=600)
plt.close()

