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
Scan="XZ"

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

def func4(x,a):
    
    return (a/np.sin(x))

def func5(x,a,b,c,d):
    
    K=np.pi/180
    return ((a*np.cos(K*x))/(b*np.sin(K*x)+c)+d)

def XZScan(thetaHs,SigmaAC,HAC,B):
    #Parameters
    # SigmaAC=5E-5
    # HAC=1E-3
    A=0.06
    
    H=3
    HE=900
    DMI=2
    Kz=1E-2
    thetaD=2*(np.pi/180)
    phiD=135*(np.pi/180)
    phiH=0*(np.pi/180)
    
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
        
        # FirstHarm.append(A*np.cos(thetaH)+B*firstharm+C)
        SecondHarm.append(A*secondharm+B)
    
    
    
    return(np.asfarray(SecondHarm))

def create_datapaths(currentfolder,temps):
    datapaths=[]
    datapath=os.getcwd()+"\\"+currentfolder
    for i in temps:
        datapaths.append(datapath+i)
    return datapaths

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

def plotter(datapath,savepath,temp="Unknown Temp"):
    filename=datapath +"\\V vs Angle - Y.txt"
    
    rawdata=np.asfarray(read(filename))
    (ave_voltage,voltages,angles,angle)=sort_average(rawdata)
    
    ave_voltage=ave_voltage-np.average(ave_voltage)
    voltages=voltages-np.average(voltages)
    ave_resistance=(ave_voltage/1.5E-3)

    cut=2
    mid=len(angle)//2
    
    angle_f=np.concatenate((angle[cut:mid-cut],angle[mid+cut:-1*cut]))
    ave_resistance_f= np.concatenate((ave_resistance[cut:mid-cut],ave_resistance[mid+cut:-1*cut]))
    popt, pcov = curve_fit(XZScan, angle_f, ave_resistance_f, [0.04,0,0],maxfev=10000)

    # np.save(r"C:\Users\Egecan\Desktop\Numerical Fitting\Np Arrays\XZ2ndHarm",ave_resistance)

    # start=1
    # step=19
    # popt, pcov = curve_fit(func5, angle[start:start+step], ave_voltage[start:start+step],initial_guess,bounds=( [-5,-5,-5,-5], [5,5,5,5]),maxfev=1000000)
    
    # print(temp,popt)
    # popt=[0.12,1,1]
    if save==True:
        fig, ax1 = plt.subplots()#Plotting Sine FIT and Average Data
        # ax1.plot(angle[start:start+step],ave_voltage[start:start+step],"o-",linewidth=0.5,markersize=5)
        ax1.plot(angle,ave_resistance*1E3,"o:",linewidth=1,markersize=6,label="Data")
        ax1.plot(angle_f, XZScan(angle_f, *popt)*1E3, 'r-',label="Fit")#Fit
        # ax1.plot(xdata, XZScan(xdata, 0.08,5E-5), 'b-',label="Fit2")
        ax1.set_xlabel("Angle $\\beta$ ($\degree$)",size=16)
        ax1.set_ylabel("$\\Delta$$R_{xy}^{2\omega}$ (m$\Omega$)",size=16)
        ax1.set_title("XZ Scan - 2$^{nd}$ Harmonic Resistance",size=14)
        ax1.tick_params(axis="both",labelsize=15)
        ax1.set_xticks([0,90,180,270,360])
        plt.legend()
        print(temp,np.round(popt,6))
        # ax1.set_xticks([0,15,30,45,60,75,90])
        # ax1.set_xticks([0,15,30])
        # ax1.set_ylim(0.05,0.3)
        # ax1.set_ylim(0.13,0.23)
        plt.tight_layout()
        plt.savefig(savepath+"Average - "+temp+".png",dpi=600)
        plt.savefig(savepath+"2nd Harm - XZ Scan - "+temp+".pdf",dpi=600)
        # plt.savefig(savepath+"Average - "+temp+".pdf",dpi=600)
        plt.close()
        fig, ax2 = plt.subplots()#Plotting All data and Average Data
        ax2.set_xlabel("Angle $\\beta$ ($\degree$)",size=16)
        ax2.set_ylabel("$\\Delta$$V_{xy}^{2\omega}$ ($\mu$V)",size=16)
        ax2.set_title("XZ Scan - 2$^{nd}$ Harmonic Voltage vs. Angle",size=14)
        halflength=int(len(voltages)/2)
        ax2.plot(angles[:halflength],voltages[:halflength],"o-",markersize=5,label="Forwards")
        ax2.plot(angles[halflength:],voltages[halflength:],"o-",markersize=5,label="Backwards")
        # ax2.plot(angles[150:224],voltages[150:224],"o-",markersize=5)
        # ax2.plot(angle,ave_voltage,"o-")#Average data
        ax2.set_xticks([0,90,180,270,360])
        plt.legend()
        # ax2.set_ylim(0.1,0.45)
        # ax2.set_ylim(0.10,0.22)
        plt.tight_layout()
        plt.savefig(savepath+"All Data - "+temp+".png",dpi=600)
        
        plt.close()
    return(popt,angle,ave_voltage)
    



def plot_temp_dependence(amplitude,wavenumber,phase,constant,data_path,save_path):
    lists = sorted(amplitude.items())
    temp,amp=zip(*lists)
    temp=list(temp)
    amp=np.asarray(list(amp))
    if save==True:
        fig, ax3 = plt.subplots()
        ax3.plot(temp,amp,"o-.",markersize=5)
        ax3.set_title("SSE Amplitude vs. Temperature",size=15)
        ax3.set_xlabel("Temperature (K)",size=14)
        ax3.set_ylabel("Peak-to-Peak Voltage ($\mu$V)",size=14)
        plt.tight_layout()
    
        plt.savefig(save_path+"SSE Amplitude vs. Temperature.png",dpi=300)
    plt.close()
    
    lists = sorted(wavenumber.items())
    temp,period=zip(*lists)
    temp=list(temp)
    period=np.asarray(list(period))
    if save==True:
        fig, ax3 = plt.subplots()
        ax3.plot(temp,(2*np.pi)/((period)/(2*np.pi/360)),"o-.",color="sienna",markersize=5)
        ax3.set_title("Period vs. Temperature",size=15)
        ax3.set_xlabel("Temperature (K)",size=14)
        ax3.set_ylabel("Period",size=14)
        plt.yticks([0*np.pi,2*np.pi,4*np.pi],["0","$2\pi$","$4\pi$"])
        plt.tight_layout()
    
        plt.savefig(save_path+"Period vs. Temperature.png",dpi=300)
    plt.close()
    
    lists = sorted(phase.items())
    temp,phi=zip(*lists)
    temp=list(temp)
    phi=np.asarray(list(phi))
    if save==True:
        fig, ax3 = plt.subplots()
        ax3.plot(temp,phi,"o-.",color="darkolivegreen",markersize=5)
        #ax3.yaxis.set_ticks([0,np.pi/2],["0","pi/2"])
        plt.yticks([-np.pi/2,0,np.pi/2],["-$\dfrac{\pi}{2}$","0","$\dfrac{\pi}{2}$"])
        ax3.set_xlabel("Temperature (K)",size=14)
        ax3.set_ylabel("Phase",size=14)
        ax3.set_title("Phase vs. Temperature",size=15)
        plt.tight_layout()
    
        plt.savefig(save_path+"Phase vs. Temperature.png",dpi=300)
    plt.close()
    
    
    lists = sorted(constant.items())
    temp,cons=zip(*lists)
    temp=list(temp)
    cons=list(cons)
    if save==True:
        fig, ax3 = plt.subplots()
        ax3.plot(temp,cons,"o-.",color="red",markersize=5)
        ax3.set_title("Voltage Offset vs. Temperature",size=15)
        ax3.set_xlabel("Temperature (K)",size=14)
        ax3.set_ylabel("Voltage ($\mu$V)",size=14)
        plt.tight_layout()
    
        plt.savefig(save_path+"Offset Voltage vs. Temperature.png",dpi=300)
    plt.close()
    return (amp,temp)


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
        # print()
        amplitude1[float(i[start:stop])]=p2p[0]
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

# data[90000][1][21]=(data[90000][1][20]+data[90000][1][22])/2
# data[90000][1][24]=(data[90000][1][23]+data[90000][1][25])/2
# data[70000][1][21]=(data[70000][1][20]+data[70000][1][22])/2
# data[50000][1][21]=(data[50000][1][20]+data[50000][1][22])/2
# data[30000][1][21]=(data[30000][1][20]+data[30000][1][22])/2


legendlist=["0 T","0.5 T","0.75 T","1 T","1.5 T","2 T","2.5 T","3 T","4 T","5 T","6 T","7 T","8 T","9 T"]
# fig, ax1 = plt.subplots()
# for (index,field) in enumerate(fields[::]):
#     plt.plot(data[field][0],data[field][1]-np.average(data[field][1]),label=legendlist[index],color=colorlist[index])

# ax1.set_xlabel("Angle $\\beta$ ($\degree$)",size=17)

# ax1.set_ylabel("$\\Delta$$V_{xy}^{2\omega}$ ($\mu$V)",size=16)
# ax1.set_title("XZ Scan of 2$^{nd}$ Harmonic Voltage",size=16)
# plt.xticks([0,90,180,270,360],[0,90,180,270,360],size=12)
# ax1.tick_params(axis="both",labelsize=12)
# plt.legend()
# ax1.set_ylim(0,1000)
# plt.tight_layout()
# plt.savefig(savepath+"all fields.png",dpi=600)
# # plt.savefig(savepath+"all fields2.pdf",dpi=600)



fig, ax1 = plt.subplots()
for (index,ifield) in enumerate(fields[::]):
    plt.plot(data[ifield][0],data[ifield][1]/1.5E-6,label=legendlist[index],color=colorlist[index])

ax1.set_xlabel("Angle $\\beta$ ($\degree$)",size=16)

# ax1.set_ylabel("$\Delta$$R_{xy}^{2\omega}$ (m$\Omega$)",size=16)
# ax1.set_title("XZ Scan - 2$^{nd}$ Harmonic Resistance",size=14)
plt.xticks([0,90,180,270,360],[0,90,180,270,360])
ax1.tick_params(axis='both', which='major', labelsize=15)
# ax1.set_ylim(0,700)
# plt.legend()
plt.tight_layout()
plt.savefig(savepath+"2nd Harmonic Resistance - XZ Scan.png",dpi=600)
plt.savefig(savepath+"2nd Harmonic Resistance - XZ Scan.pdf",dpi=600)
# plt.savefig(savepath+"all fields2.pdf",dpi=600)
plt.close()


