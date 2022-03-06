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
from scipy.signal import savgol_filter


from matplotlib import cm
from  matplotlib.colors import rgb2hex



step_angle=10
Scan="YZ"


def linear(x,a,b):
    return a*x+b

def linear0(x,a):
    return a*x

cmap = cm.get_cmap('viridis',6)    # PiYG
colorlist=[]
for i in range(cmap.N):
    rgba = cmap(i)
    # rgb2hex accepts rgb or rgba
    colorlist.append(rgb2hex(rgba))


def YZScan(thetaHs,A,B,C):
    #Parameters
    SigmaAC=4E-5
    HAC=5E-3
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
        
        FirstHarm.append(A*np.cos(thetaH)+B*firstharm+C)
        SecondHarm.append(B*secondharm)
    
    
    
    return(np.asfarray(FirstHarm))

def func2(x, a, b, c, d,e,f,g):
    return a * np.sin(b*x+c) + d*np.sin(e*x+f)+g

def func1(x, a, b, c, g):
    return a * np.sin(b*x+c) + g

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
    
    numberofsweeps=len(voltages)/((360/step_angle)+1)
    
    splitvoltagesv=np.split(voltages,numberofsweeps)
    splitangles=np.split(angles,numberofsweeps)  
    ave_voltage=[]
    total=np.zeros(int(360/step_angle)+1)
    for i in range(int(numberofsweeps)):
        if (i)%2==0:
            total=total+splitvoltagesv[i]
        else:
            total=total+np.flip(splitvoltagesv[i])
    ave_voltage=total/(numberofsweeps)       

#    ave_voltage=np.flip(voltages[int(len(voltages)/2):])   
#    ave_voltage=savgol_filter(ave_voltage,3,1)
    
    return (ave_voltage,voltages,angles,splitangles[0])

def flatten(ave_voltage):
    start=np.average(ave_voltage[0:3])
    end=np.average(ave_voltage[-3:])
    mult=start/end
    
    return ave_voltage*mult

def plotter(datapath,savepath,temp="Unknown Temp"):
    filename=datapath +"\\V vs Angle - R.txt"
    rawdata=np.asfarray(read(filename))
    (ave_voltage,voltages,angles,angle)=sort_average(rawdata)
    xdata = np.linspace(0, 360,361)
    ave_voltage=ave_voltage-np.average(ave_voltage)
    voltages=voltages-np.average(voltages)
    ave_resistance=(ave_voltage/1E-4)
    
    cut=1
    mid=len(angle)//2
    
    angle_f=np.concatenate((angle[cut:mid-cut],angle[mid+cut:-1*cut]))
    ave_resistance_f= np.concatenate((ave_resistance[cut:mid-cut],ave_resistance[mid+cut:-1*cut]))
    
    popt, pcov = curve_fit(YZScan, angle_f, ave_resistance_f,  [-5E-3,-2E-2,-1E-5],maxfev=10000)
    print(popt)
    fig, ax1 = plt.subplots()#Plotting Sine FIT and Average Data

    ax1.plot(angle,ave_resistance*1E3,"o:",markersize=5,color="C0",label="Data")
    ax1.plot(angle_f, YZScan(angle_f, *popt)*1E3,label="Fit",color="red")#Fit


    # ax1.set_xlabel("Angle $\\$\\gamma$$ ($\degree$)\nFit Amplitude="+str(round(popt[0],3)),size=17)
    ax1.set_xlabel("Angle $\\gamma$ ($\degree$)",size=16)
    ax1.set_ylabel("$\Delta R_{xy}^{1\omega}$ (m$\Omega$)",size=16)
    ax1.set_title("YZ Scan - 1$^{st}$ Harmonic Resistance",size=14)
    plt.xticks([0,90,180,270,360],[0,90,180,270,360])
    ax1.tick_params(axis='both', which='major', labelsize=15)
    plt.legend()
    plt.tight_layout()
    if save==True:
        plt.savefig(savepath+"Average - "+temp+".png",dpi=600)
        plt.savefig(savepath+"YZ Scan - Average - "+temp+".pdf",dpi=600)
    if show==False:
        plt.close()

    # fig, ax2 = plt.subplots()#Plotting All data and Average Data
    # ax2.set_xlabel("Angle $\\gamma$ ($\degree$)",size=16)
    # ax2.set_ylabel("$V_{xy}^{1\omega}$ ($\mu$V)",size=16)
    # ax2.set_title("YZ Scan - 1$^{st}$ Harmonic Voltage vs. Angle",size=16)
    # plt.xticks([0,90,180,270,360],[0,90,180,270,360])
    # ax2.tick_params(axis='both', which='major', labelsize=15)
    # halflength=int(len(voltages)/2)
    # ax2.plot(angles[:halflength],voltages[:halflength],"o-",markersize=5,label="Forwards")
    # ax2.plot(angles[halflength:],voltages[halflength:],"o-",markersize=5,label="Backwards")
    # plt.legend()
    # plt.tight_layout()
    # if save==True:
    #     plt.savefig(savepath+"All Data - "+temp+".png",dpi=600)
    # if show==False:        
    #     plt.close()
    return(angle,ave_voltage,popt)
    



def call_plotter(datapath,savepath):
    temps=os.listdir(datapath)
    amplitude={}
    # amplitude2={}
    # wavenumber1={}
    # wavenumber2={}
    # phase1={}
    # phase2={}
    # constant={}
    for i in temps:
        (angles,voltages,popt)=plotter(datapath+i,savepath,i)
        start=i.find('K')+4
        stop=i.find('O')-1
        amplitude[float(i[start:stop])]=(angles,voltages,popt)
        # wavenumber1[float(i[:stop])]=p2p[1]   
        # phase1[float(i[:stop])]=p2p[2]   
        # amplitude2[float(i[:stop])]=p2p[3]
        # wavenumber2[float(i[:stop])]=p2p[4]   
        # phase2[float(i[:stop])]=p2p[5]
        # constant[float(i[:stop])]=p2p[6]           
    return (amplitude)



path=dirname(os.getcwd())
    
data_path=path+"\Data\\"
save_path=path+"\Graphs\\"

save=False
show=False
save2=False

amplitude=call_plotter(data_path,save_path)
# (amp,temp)=plot_temp_dependence(amplitude1,wavenumber1,phase1,constant,amplitude2,wavenumber2,phase2,data_path,save_path)

amplitude2=amplitude
fields=np.sort(np.asfarray(list(amplitude.keys())))

fig, ax1 = plt.subplots()
for (index,field) in enumerate(fields):
    plt.plot(amplitude[field][0],amplitude[field][1],color=colorlist[index],label=str(int(field))+" Oe")
    ax1.set_xlabel("Angle $\\gamma$ ($\degree$)",size=17)

ax1.set_ylabel("$V_{xy}^{1\omega}$ ($\mu$V)",size=16)
ax1.tick_params(axis='both', which='major', labelsize=15)
# ax1.set_title(Scan+" Scan - 1$^{st}$ Harmonic Voltage vs. Angle",size=14)
plt.xticks([0,90,180,270,360],[0,90,180,270,360])
plt.legend()
plt.tight_layout()
plt.savefig(save_path+"Voltages.png",dpi=600)
# plt.savefig(save_path+"all fields.pdf",dpi=600)
    
fig, ax1 = plt.subplots()
# for (index,field) in enumerate([fields[0],fields[1],fields[2],fields[4],fields[8],fields[12]]):
for (index,field) in enumerate(fields):
    plt.plot(amplitude[field][0],(amplitude[field][1]-np.average(amplitude[field][1]))*1E6,color=colorlist[index],label=str(int(field))+" Oe")
ax1.set_aspect(15)
ax1.set_xlabel("Angle $\\gamma$ ($\degree$)",size=17)

ax1.set_ylabel("$\Delta$$V_{xy}^{1\omega}$ ($\mu$V)",size=16)
ax1.set_title(Scan+" Scan - 1$^{st}$ Harmonic Voltage",size=14)
plt.xticks([0,90,180,270,360],[0,90,180,270,360])
ax1.tick_params(axis='both', which='major', labelsize=15)
# plt.legend()
plt.tight_layout()
plt.savefig(save_path+r"\\1st Harmonic YZ Scan - Voltages - All Fields - Delta.png",dpi=600)
plt.savefig(save_path+r"\\1st Harmonic YZ Scan - Voltages - All Fields - Delta.pdf",dpi=600)
# plt.savefig(save_path+"all fields2.pdf",dpi=600)

fig, ax1 = plt.subplots()
# for (index,field) in enumerate([fields[0],fields[1],fields[2],fields[4],fields[8],fields[12]]):
for (index,field) in enumerate(fields):
    plt.plot(amplitude[field][0],(amplitude[field][1]-np.average(amplitude[field][1]))/1E-7,color=colorlist[index],label=str(int(field))+" Oe")

ax1.set_xlabel("Angle $\\gamma$ ($\degree$)",size=17)

# ax1.set_ylabel("$\Delta$ $R_{xy}^{1\omega}$ (m$\Omega$)",size=16)
# ax1.set_title(Scan+" Scan - 1$^{st}$ Harmonic Resistance",size=14)
plt.xticks([0,90,180,270,360],[0,90,180,270,360])
ax1.tick_params(axis='both', which='major', labelsize=15)
# plt.legend()
plt.tight_layout()
plt.savefig(save_path+r"\\1st Harmonic YZ Scan - Resistances - All Fields - Delta.png",dpi=600)
plt.savefig(save_path+r"\\1st Harmonic YZ Scan - Resistances - All Fields - Delta.pdf",dpi=600)
plt.close()


# (amplitude)=call_plotter(data_path,save_path)
# # (amp,temp)=plot_temp_dependence(amplitude1,wavenumber1,phase1,constant,amplitude2,wavenumber2,phase2,data_path,save_path)



# # (amp,temp)=plot_temp_dependence(amplitude1,wavenumber1,phase1,constant,amplitude2,wavenumber2,phase2,data_path,save_path)
# amplitudes=np.array(sorted(list(amplitude.items())),dtype=object)
# X_=np.asfarray(amplitudes[:,0]/10000)
# Y_=np.abs(np.asfarray(amplitudes[:,1]*1000))
# fig, ax = plt.subplots()
# ax.plot(X_,Y_,"o:",color="C0",label="Data")
# popt,pcov=curve_fit(linear0,X_, Y_)
# ax.plot(amplitudes[2:,0]/10000,linear0(amplitudes[2:,0]/10000,*popt),color="C3",label="Linear Fit")
# # ax.set_title("SMR Amplitude vs. Temperature",size=16)
# ax.set_ylabel("AHE Coefficient ($\\times$10$^{-3}$)",size=15)
# ax.set_xlabel("Magnetic Field (T)",size=15)
# ax.tick_params(axis='both', which='major', labelsize=15)
# ax.set_title("XZ & YZ Scans - AHE Coefficient vs. Field",size=16)
# plt.legend()
# plt.tight_layout()

# plt.savefig(save_path+"Coef1 vs Field.png",dpi=600)
# plt.savefig(save_path+"Coef1 vs Field.pdf",dpi=600)

# # # # amplitude2=amplitude
# # fields=np.sort(np.asfarray(list(amplitude.keys())))


# fieldlabels=["0 T","0.5 T","0.75 T","1 T","1.5 T","2 T","2.5 T","3 T","4 T","5 T","6 T","7 T","8 T","9 T"]




# fig, ax1 = plt.subplots()
# for (index,field) in enumerate(fields):
#     plt.plot(amplitude[field][0],(amplitude[field][1]-np.average(amplitude[field][1]))/1E-7,label=fieldlabels[index]+" Tesla",color=colorlist[index])
# ax1.set_xlabel("Angle $\\gamma$ ($\degree$)",size=16)
# ax1.set_ylabel("$\Delta$$R_{xy}^{1\omega}$ (m$\Omega$)",size=16)
# ax1.set_title(Scan+" Scan - 1$^{st}$ Harmonic Resistance",size=14)
# plt.xticks([0,90,180,270,360],[0,90,180,270,360],size=15)
# ax1.tick_params(axis='both', which='major', labelsize=15)
# # plt.legend()
# plt.tight_layout()
# plt.savefig(save_path+"Resistances_AllFields.png",dpi=600)
# plt.savefig(save_path+"Resistances_AllFields.pdf",dpi=600)
# # plt.savefig(save_path+"all fields.pdf",dpi=600)




# # fig, ax1 = plt.subplots()
# # # for (index,field) in enumerate([fields[0],fields[1],fields[2],fields[4],fields[8],fields[12]]):
# # for (index,field) in enumerate(fields):
# #     plt.plot(amplitude[field][0],amplitude[field][1]-np.average(amplitude[field][1]),label=fieldlabels[index],color=colorlist[index])

# # ax1.set_xlabel("Angle $\\gamma$ ($\degree$)",size=17)

# # ax1.set_ylabel("$\Delta$$V_{xy}^{1\omega}$ ($\mu$V)",size=16)
# # ax1.set_title(Scan+" Scan - 1$^{st}$ Harmonic Resistance vs. Angle",size=16)
# # plt.xticks([0,90,180,270,360],[0,90,180,270,360],size=15)
# # # plt.legend()
# # plt.tight_layout()
# # # ax1.set_ylim(0,1000)
# # plt.savefig(save_path+"Voltages - All Fields - Delta.png",dpi=600)
# # # plt.savefig(save_path+"all fields2.pdf",dpi=600)
# # plt.close()




# # fig, ax1 = plt.subplots()
# # # for (index,field) in enumerate([fields[0],fields[1],fields[2],fields[4],fields[8],fields[12]]):
# # for (index,field) in enumerate(fields):
# #     plt.plot(amplitude[field][0],(amplitude[field][1]-np.average(amplitude[field][1]))/0.1,label=fieldlabels[index]+" Tesla",color=colorlist[index])

# # ax1.set_xlabel("Angle $\\gamma$ ($\degree$)",size=17)

# # ax1.set_ylabel("$\Delta$ $R_{xy}^{1\omega}$ (m$\Omega$)",size=16)
# # ax1.set_title(Scan+" Scan - 1$^{st}$ Harmonic Voltage vs. Angle",size=16)
# # plt.xticks([0,90,180,270,360],[0,90,180,270,360])
# # # plt.legend()
# # plt.tight_layout()
# # plt.savefig(save_path+"Resistances_AllFields.png",dpi=600)
# # plt.close()




# # fitamplitudes=[]

# # for (index,field) in enumerate(fields):
# #     fitamplitudes.append(amplitude[field][2][0])
   
    
# # fitamplitudes=np.asfarray(fitamplitudes)/0.1  
# # fig, ax1 = plt.subplots()  
# # plt.plot(fields/10000,fitamplitudes,"o:",markersize=9,label="SMR Amplitude")

# # ax1.set_xlabel("External Field (Tesla)",size=17)

# # ax1.set_ylabel("SMR Amplitude (m$\Omega$)",size=16)
# # ax1.set_title(Scan+" Scan - SMR Amplitude vs. Angle",size=16)
# # plt.legend()
# # plt.tight_layout()
# # plt.savefig(save_path+"SMR Amplitude Field Dependence.png",dpi=600)
# # plt.close()