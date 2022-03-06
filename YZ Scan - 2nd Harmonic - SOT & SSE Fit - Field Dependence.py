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
Scan="YZ"


cmap = cm.get_cmap('viridis', 6)    # PiYG
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

def func4(x,a):
    
    return (a/np.sin(x))

def func5(x,a,b,c,d):
    
    K=np.pi/180
    return ((a*np.cos(K*x))/(b*np.sin(K*x)+c)+d)


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
    # print(len(voltages))
    step_angle=5
    numberofangles=360/step_angle+1
    numberofsweeps=len(voltages)/(numberofangles)

    v=np.split(voltages,numberofsweeps)
    a=np.split(angles,numberofsweeps)
    ave_voltage=[]
    total=np.zeros(int(numberofangles))
    for i in range(int(numberofsweeps)):
        if i%2==0:
            total=total+v[i]
        else:
            total=total+np.flip(v[i])
    ave_voltage=total/(numberofsweeps)
    
    return (ave_voltage,voltages,angles,a[0])

def flatten(ave_voltage):
    start=np.average(ave_voltage[0:3])
    end=np.average(ave_voltage[-3:])
    mult=start/end

    return ave_voltage*mult

def YZScan(thetaHs,SigmaAC,HAC,B):
    #Parameters
    # SigmaAC=5E-5
    # A=0.06

    A=0.04

    H=3
    HE=1800
    DMI=2
    Kz=1E-2
    thetaD=2*(np.pi/180)
    phiD=135*(np.pi/180)
    phiH=90*(np.pi/180)
    
    thetaHs=np.asfarray(thetaHs)*(np.pi/180)
    #Defitions of X, Y and Z unit vectors
    xvector=[1,0,0]
    yvector=[0,1,0]
    zvector=[0,0,1]
    
    FirstHarm=[]
    SecondHarm=[]
    for thetaH in thetaHs:

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
        
        #FirstHarm.append(A*np.cos(thetaH)+B*firstharm)
        SecondHarm.append(A*secondharm+B)
    
    
    return(np.asfarray(SecondHarm))


def plotter(datapath,savepath,temp="Unknown Temp"):
    filename=datapath +"\\V vs Angle - Y.txt"
    
    rawdata=np.asfarray(read(filename))

    (ave_voltage,voltages,angles,angle)=sort_average(rawdata)
    ave_voltage=-1*(ave_voltage-np.average(ave_voltage))
    voltages=-1*(voltages-np.average(voltages))
    ave_resistance=(ave_voltage/1.5E-3)
    xdata = np.append(np.arange(0, 176,5),np.arange(185, 361,5))
    ave_voltage=ave_voltage-np.average(ave_voltage)
    voltages=voltages-np.average(voltages)
    cut=1
    mid=len(angle)//2
    
    angle_f=np.concatenate((angle[cut:mid-cut],angle[mid+cut:-1*cut]))
    ave_resistance_f= np.concatenate((ave_resistance[cut:mid-cut],ave_resistance[mid+cut:-1*cut]))

    # popt, pcov = curve_fit(YZScan, angle_f, ave_resistance_f, [5E-5,5E-3,0],bounds=[[0,1E-3,-1],[1,1E-2,1]],maxfev=100000)
    popt, pcov = curve_fit(YZScan, angle_f, ave_resistance_f, [5E-5,5E-3,0],maxfev=100000)
    # print(pcov)
    np.save(r"C:\Users\Egecan\Desktop\Numerical Fitting\Np Arrays\YZ2ndHarm",ave_resistance)
    print(temp,popt)
    # popt=[0.12,1,1]
    if save==True:
        # fig, ax1 = plt.subplots(figsize=(16/2.6,11.2/2.6))#Plotting Sine FIT and Average Data
        fig, ax1 = plt.subplots()
        # ax1.plot(angle[start:start+step],ave_voltage[start:start+step],"o-",linewidth=0.5,markersize=5)
        ax1.plot(angle,ave_resistance*1E3,"o:",markersize=5.5,label="Data")
        ax1.plot(angle_f[2:], YZScan(angle_f, *popt)[2:]*1E3, 'r-',label="Fit")#Fit
        ax1.set_xlabel("Angle $\\gamma$ ($\degree$)",size=16)
        ax1.set_ylabel("$\\Delta$$R_{xy}^{2\omega}$ (m$\Omega$)",size=16)
        ax1.set_title("YZ Scan - 2$^{nd}$ Harmonic Resistance",size=14)
        # ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        ax1.tick_params(axis="both",labelsize=15)
        ax1.set_xticks([0,90,180,270,360])
        ax1.set_yticks([-0.03,0,0.03,0.06,0.09])
        ax1.ticklabel_format(useOffset=True)
        plt.legend()
        # ax1.set_xticks([0,15,30,45,60,75,90])
        # ax1.set_xticks([0,15,30])
        # ax1.set_ylim(0.05,0.3)
        # ax1.set_ylim(0.13,0.23)
        plt.tight_layout()
        plt.savefig(savepath+"Average - "+temp+".png",dpi=600)
        plt.savefig(savepath+"2nd Harm - YZ Scan - "+temp+".pdf",dpi=600)
        plt.close()
        fig, ax2 = plt.subplots()#Plotting All data and Average Data
        ax2.set_xlabel("Angle $\\gamma$ ($\degree$)",size=16)
        ax2.set_ylabel("$\\Delta$$R_{xy}^{2\omega}$ ($\Omega$)",size=16)
        ax2.set_title("YZ Scan - 2$^{nd}$ Harmonic Resistance",size=14)
        ax2.set_xticks([0,90,180,270,360])
        ax2.tick_params(axis="both",labelsize=15)
        halflength=int(len(voltages)/2)
        ax2.plot(angles[:halflength],voltages[:halflength],"o-",markersize=5,label="Forwards")
        ax2.plot(angles[halflength:],voltages[halflength:],"o-",markersize=5,label="Backwards")
        # ax2.plot(angles[150:224],voltages[150:224],"o-",markersize=5)
        # ax2.plot(angle,ave_voltage,"o-")#Average data
        ax2.set_xticks([0,90,180,270,360])
        # ax2.ticklabel_format(axis="both",style="scientific")
        # ax2.ticklabel_format(useOffset=True)
        plt.legend()
        # ax2.set_ylim(0.1,0.45)
        # ax2.set_ylim(0.10,0.22)
        plt.tight_layout()
        plt.savefig(savepath+"All Data - "+temp+".png",dpi=300)
        
        plt.close()
    return(angle,ave_voltage)
    




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
        (angle,ave_voltage)=plotter(datapath+i,savepath,i)
        start=i.find('K')+3
        stop=i.find('O')
        # print()
        # amplitude1[float(i[start:stop])]=p2p[0]
        # wavenumber1[float(i[start:stop])]=p2p[1]   
        # phase1[float(i[start:stop])]=p2p[2]
        # amplitude2[float(i[start:stop])]=p2p[3]
        # wavenumber2[float(i[start:stop])]=p2p[4]   
        # phase2[float(i[start:stop])]=p2p[5]
        # constant[float(i[start:stop])]=p2p[6]
        data[float(i[start:stop])]=(angle,ave_voltage)
    return (data)



def plot_field_dependence(amplitude1,wavenumber1,phase1,amplitude2,wavenumber2,phase2,constant,data_path,save_path):
    
    # data_path=path+"Data\\"
    # save_path=path+"Graphs\\"

    fields=os.listdir(path+"Data\\")
    amp_dict={}

    data_path=path+"Data\\"
    save_path=path+"Graphs\\"


   
    # amp_dict[Field]=[amp,temp]
    
    
    fig, ax = plt.subplots()
    for Field in fields:
        ax.plot(amp_dict[Field][1],amp_dict[Field][0],"o-.",markersize=5,label=Field)
    ax.set_title("SSE Amplitude vs. Temperature",size=15)
    ax.set_xlabel("Temperature (K)",size=14)
    ax.set_ylabel("Peak-to-Peak Voltage ($\mu$V)",size=14)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path+"SSE Amplitude vs. Temperature - Different Fields.png",dpi=300)
    
    return amp_dict





save=True
path=dirname(os.getcwd())

datapath=path+"\\Data\\"
savepath=path+"\\Graphs\\"

data=call_plotter(datapath, savepath)

# plot_field_dependence(amplitude1,wavenumber1,phase1,amplitude2,wavenumber2,phase2,constant,datapath,savepath)

fields=np.sort(np.asfarray(list(data.keys())))

# fig, ax1 = plt.subplots()
# for (index,field) in enumerate(fields):
#     plt.plot(data[field][0],data[field][1]+14+0.3*index,label=str(int(field))+" Oe")
#     ax1.set_xlabel("Angle $\\alpha$ ($\degree$)",size=17)

# ax1.set_ylabel("$V_{xy}^{1\omega}$ ($\mu$V)",size=16)
# ax1.set_title("SMR Amplitude vs. Angle",size=16)
# plt.xticks([0,90,180,270,360],[0,90,180,270,360])
# plt.legend()
# plt.tight_layout()
# plt.savefig(savepath+"all fields.png",dpi=600)
# plt.savefig(savepath+"all fields.pdf",dpi=600)

legendlist=["0 T","0.5 T","1 T","1.5 T","2 T","2.5 T","3 T","4 T","5 T","6 T","7 T","8 T","9 T"]

fig, ax1 = plt.subplots()
for (index,ifield) in enumerate(fields):
    plt.plot(data[ifield][0],data[ifield][1]/1.5E-6,label=legendlist[index],color=colorlist[index],linewidth=4)

ax1.set_xlabel("Angle $\\gamma$ ($\degree$)",size=16)
# ax1.set_ylabel("$\\Delta$$R_{xy}^{2\omega}$ (m$\Omega$)",size=16)
# ax1.set_title("YZ Scan - 2$^{nd}$ Harmonic Resistance",size=14)
plt.xticks([0,90,180,270,360],[0,90,180,270,360],size=15)
ax1.tick_params(axis="both",labelsize=15)
# ax1.set_ylim(0,700)
# plt.legend()
plt.tight_layout()
plt.savefig(savepath+"2nd Harmonic Resistance - YZ Scan.png",dpi=600)
plt.savefig(savepath+"2nd Harmonic Resistance - YZ Scan.pdf",dpi=600)
# plt.savefig(savepath+"all fields2.pdf",dpi=600)
plt.close()



