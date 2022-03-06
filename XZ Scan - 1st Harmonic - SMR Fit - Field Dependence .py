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

step_angle=10
Scan="XZ"




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
    voltages=rawdata[:,1]*1000000
    
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

    # popt, pcov = curve_fit(func1, angle, ave_voltage, [1.5,0.034,0,49])
    # print(temp,popt)
    
    fig, ax1 = plt.subplots()#Plotting Sine FIT and Average Data
    # ax1.plot(angles,voltages,color="C0",linewidth=1)
    ax1.plot(angle,ave_voltage,"o:",markersize=5,color="C0",label="Data")
    # ax1.plot(xdata, func1(xdata, *popt),label="Fit",color="red")#Fit

    popt=np.asfarray([1.5,0.034,0,49])
    # ax1.set_xlabel("Angle $\\beta$ ($\degree$)\nFit Amplitude="+str(round(popt[0],3)),size=17)
    ax1.set_xlabel("Angle $\\beta$ ($\degree$)",size=17)

    ax1.set_ylabel("$V_{xy}^{1\omega}$ ($\mu$V)",size=16)
    ax1.set_title(Scan+" Scan - 1$^{st}$ Harmonic Voltage vs. Angle",size=16)
    plt.xticks([0,90,180,270,360],[0,90,180,270,360])
    ax1.tick_params(axis='both', which='major', labelsize=15)
    plt.legend()
    plt.tight_layout()
    if save==True:
        plt.savefig(savepath+"Average - "+temp+".png",dpi=600)
    if show==False:
        plt.close()

    fig, ax2 = plt.subplots()#Plotting All data and Average Data
    ax2.set_xlabel("Angle ($\degree$)\nFit Amplitude="+str(round(popt[0],3)))
    
    ax2.set_ylabel(r"Voltage ($\mu$V)")
    ax2.set_title("Voltage vs. Angle - All Data - "+temp)
    halflength=int(len(voltages)/2)
    ax2.plot(angles[:halflength],voltages[:halflength],"o-",markersize=5,label="Forwards")
    ax2.plot(angles[halflength:],voltages[halflength:],"o-",markersize=5,label="Backwards")
    plt.legend()
    plt.tight_layout()
    if save==True:
        plt.savefig(savepath+"All Data - "+temp+".png",dpi=600)
    if show==False:        
        plt.close()
    return(angle,ave_voltage,popt)
    

    
def multiplotter(datapaths,savepath):
    fig, ax = plt.subplots()
    for datapath in datapaths:
        filename=datapath+"\\V vs Angle - R.txt"
        a=np.asfarray(read(filename))
        angles=a[:,0]
        voltage=a[:,1]*1000000
    
        split=len(voltage)/37

        v=np.split(voltage,split)
        a=np.split(angles,split)

        ave_voltage=np.mean([v[0],np.flip(v[1])],axis=0)

        xdata = np.linspace(0, 360,360)

        popt, pcov = curve_fit(func2, a[0], ave_voltage, [0.04,0.017,2,2],bounds=([0,0.016,0,0], [2,0.018,1.9,10]))

#        ax.plot(a[0],ave_voltage,"o",markersize=5)
        ax.plot(xdata, func1(xdata, *popt), '-',label=datapath[-5:])#Fit
    ax.set_xlabel("Angle ($\degree$)\nFit Amplitude="+str(round(popt[0],3)))
    ax.set_ylabel("Voltage ($\mu$V)")
    ax.set_title("Voltage vs. Angle")
    plt.tight_layout()
    plt.legend()
#    plt.savefig(savepath+"Average Chosen Temps- "+".png",dpi=600)
#    plt.close()
    # return (voltages,angles)


def plot_temp_dependence(amplitude1,wavenumber1,phase1,constant,amplitude2,wavenumber2,phase2,data_path,save_path):
    lists = sorted(amplitude1.items())
    temp,amp1=zip(*lists)
    temp=np.asarray(list(temp))
    amp1=np.asarray(list(amp1))

    # fig, ax3 = plt.subplots()
    # ax3.plot(temp,amp1,"o-.",markersize=5)
    # ax3.set_title("SSE Amplitude vs. Temperature - Part I",size=16)
    # ax3.set_xlabel("Temperature (K)",size=15)
    # ax3.set_ylabel("SSE Amplitude ($\mu$V)",size=15)
    # ax3.tick_params(axis='both', which='both', labelsize=14)
    # plt.tight_layout()
    # if save2==True:
    #     plt.savefig(save_path+"SSE Amplitude vs. Temperature - Part 1.png",dpi=1200)
    #     # plt.savefig(save_path+"SSE Amplitude vs. Temperature.pdf",dpi=1200)
    # plt.close()
    
    lists = sorted(amplitude2.items())
    temp,amp2=zip(*lists)
    temp=np.asarray(list(temp))
    amp2=np.asarray(list(amp2))

    fig, ax3 = plt.subplots()
    ax3.plot(temp,amp1,"o-.",markersize=5,label="1st Component")
    ax3.plot(temp,amp2,"o-.",markersize=5,label="2nd Component")
    ax3.set_title("SSE Amplitude vs. Temperature",size=16)
    ax3.set_xlabel("Temperature (K)",size=15)
    ax3.set_ylabel("SSE Amplitude ($\mu$V)",size=15)
    ax3.tick_params(axis='both', which='both', labelsize=14)
    plt.legend()
    plt.tight_layout()
    
    if save2==True:
        plt.savefig(save_path+"SSE Amplitude vs. Temperature - Part 1&2.png",dpi=1200)
        # plt.savefig(save_path+"SSE Amplitude vs. Temperature.pdf",dpi=1200)
    plt.close()
    
    
    lists = sorted(wavenumber1.items())
    temp,period1=zip(*lists)
    temp=np.asarray(list(temp))
    period1=np.asarray(list(period1))

    fig, ax3 = plt.subplots()
    ax3.plot(temp,(2*np.pi)/((period1)/(2*np.pi/360)),"o-.",color="sienna",markersize=5)
    ax3.set_title("Period vs. Temperature - Part 1",size=15)
    ax3.set_xlabel("Temperature (K)",size=14)
    ax3.set_ylabel("Period",size=14)
    plt.yticks([0*np.pi,np.pi,2*np.pi],["0","$\pi$","$2\pi$"])
    plt.tight_layout()
    if save2==True:
        plt.savefig(save_path+"Period vs. Temperature - Part 1.png",dpi=600)

    plt.close()
    
    
    lists = sorted(wavenumber2.items())
    temp,period2=zip(*lists)
    temp=np.asarray(list(temp))
    period2=np.asarray(list(period2))

    fig, ax3 = plt.subplots()
    ax3.plot(temp,(2*np.pi)/((period1)/(2*np.pi/360)),"o-.",color="sienna",markersize=5,label="1st Component")
    ax3.plot(temp,(2*np.pi)/((period2)/(2*np.pi/360)),"o-.",color="darksalmon",markersize=5,label="2nd Component")
    ax3.set_title("Period vs. Temperature",size=15)
    ax3.set_xlabel("Temperature (K)",size=14)
    ax3.set_ylabel("Period",size=14)
    plt.yticks([0*np.pi,(1/3)*np.pi,(2/3)*np.pi,np.pi,(4/3)*np.pi,(5/3)*np.pi,2*np.pi],["0","","2$\pi$/3","$\pi$","","","$2\pi$"])
    plt.legend()
    plt.tight_layout()
    if save2==True:
        plt.savefig(save_path+"Period vs. Temperature - Part 1&2.png",dpi=600)

    plt.close()
    
    
    
    lists = sorted(phase1.items())
    temp,phi1=zip(*lists)
    temp=np.asarray(list(temp))
    phi1=np.asarray(list(phi1))

    fig, ax3 = plt.subplots()
    ax3.plot(temp,phi1,"o-.",color="darkolivegreen",markersize=5)
    #ax3.yaxis.set_ticks([0,np.pi/2],["0","pi/2"])
    # plt.yticks([-np.pi/2,0,np.pi/2],["-$\dfrac{\pi}{2}$","0","$\dfrac{\pi}{2}$"])
    plt.yticks([0,np.pi,2*np.pi],["0","$\pi$","$2\pi$"])
    ax3.set_xlabel("Temperature (K)",size=14)
    ax3.set_ylabel("Phase",size=14)
    ax3.set_title("Phase vs. Temperature - Part 1",size=15)
    plt.tight_layout()
    if save2==True:
        plt.savefig(save_path+"Phase vs. Temperature - Part 1.png",dpi=600)
        
    plt.close()
    
    
        
    lists = sorted(phase2.items())
    temp,phi2=zip(*lists)
    temp=np.asarray(list(temp))
    phi2=np.asarray(list(phi2))

    fig, ax3 = plt.subplots()
    ax3.plot(temp,phi1,"o-.",color="darkolivegreen",markersize=5,label="1st Component")
    ax3.plot(temp,phi2,"o-.",color="darkcyan",markersize=5,label="2nd Component")
    #ax3.yaxis.set_ticks([0,np.pi/2],["0","pi/2"])
    plt.yticks([-np.pi/2,0,np.pi/2],["-$\dfrac{\pi}{2}$","0","$\dfrac{\pi}{2}$"])
    plt.yticks([0,np.pi,2*np.pi],["0","$\pi$","$2\pi$"])
    ax3.set_xlabel("Temperature (K)",size=14)
    ax3.set_ylabel("Phase",size=14)
    ax3.set_title("Phase vs. Temperature",size=15)
    plt.legend()
    plt.tight_layout()
    if save2==True:
        plt.savefig(save_path+"Phase vs. Temperature - Part 1&2.png",dpi=600)
        
    plt.close()
    
    
    lists = sorted(constant.items())
    temp,cons=zip(*lists)
    temp=np.asarray(list(temp))
    cons=np.asarray(list(cons))

    fig, ax3 = plt.subplots()
    ax3.plot(temp,cons,"o-.",color="red",markersize=5)
    ax3.set_title("Offset vs. Temperature",size=16)
    ax3.set_xlabel("Temperature (K)",size=16)
    ax3.set_ylabel("$Voltage_{xy}^{2\omega}$ ($\muV)",size=16)
    ax3.set_ylabel("$Voltage$ ($\mu$V)")
    plt.tight_layout()
    if save2==True:
        plt.savefig(save_path+"Offset Voltage vs. Temperature.png",dpi=600)
#    plt.close()
    return (amp1,temp)


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


def plot_field_dependence(path):
    
    fields=os.listdir(path)
    amp_dict={}
    for Field in fields:    
        data_path=path+Field+"\Data\\"
        save_path=path+Field+"\Graphs\\"

        (amplitude,wavenumber,phase,constant)=call_plotter(data_path,save_path)
        (amp,temp)=plot_temp_dependence(amplitude,wavenumber,phase,constant,data_path,save_path)
        amp_dict[Field]=(amp,temp)
    
    
    fig, ax = plt.subplots()
    for Field in fields:
        ax.plot(amp_dict[Field][1],amp_dict[Field][0],"o-.",markersize=5,label=Field)
    ax.set_title("SMR Amplitude vs. Temperature",size=16)
    ax.set_xlabel("Temperature (K)",size=15)
    ax.set_ylabel("SMR Amplitude ($\mu$V)",size=15)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path+"SMR Amplitude vs. Temperature - Different Fields.png",dpi=600)


path=dirname(os.getcwd())
    
data_path=path+"\Data\\"
save_path=path+"\Graphs\\"

save=True
show=False
save2=True

(amplitude)=call_plotter(data_path,save_path)
# (amp,temp)=plot_temp_dependence(amplitude1,wavenumber1,phase1,constant,amplitude2,wavenumber2,phase2,data_path,save_path)

amplitude2=amplitude
fields=np.sort(np.asfarray(list(amplitude.keys())))


fieldlabels=["0","0.5","0.75","1","1.5","2","2.5","3"]

fig, ax1 = plt.subplots()
for (index,field) in enumerate(fields):
    plt.plot(amplitude[field][0],amplitude[field][1],label=fieldlabels[index]+" Tesla")
    ax1.set_xlabel("Angle $\\beta$ ($\degree$)",size=17)

ax1.set_ylabel("$V_{xy}^{1\omega}$ ($\mu$V)",size=16)
ax1.set_title(Scan+" Scan - 1$^{st}$ Harmonic Voltage vs. Angle",size=16)
plt.xticks([0,90,180,270,360],[0,90,180,270,360])
plt.legend()
plt.tight_layout()
plt.savefig(save_path+"Voltages .png",dpi=600)
# plt.savefig(save_path+"all fields.pdf",dpi=600)
    
fig, ax1 = plt.subplots()
# for (index,field) in enumerate([fields[0],fields[1],fields[2],fields[4],fields[8],fields[12]]):
for (index,field) in enumerate(fields):
    plt.plot(amplitude[field][0],amplitude[field][1]-np.average(amplitude[field][1]),label=fieldlabels[index]+" Tesla")

ax1.set_xlabel("Angle $\\beta$ ($\degree$)",size=17)

ax1.set_ylabel("$\Delta$$V_{xy}^{1\omega}$ ($\mu$V)",size=16)
ax1.set_title(Scan+" Scan - 1$^{st}$ Harmonic Voltage vs. Angle",size=16)
plt.xticks([0,90,180,270,360],[0,90,180,270,360])
plt.legend()
plt.tight_layout()
plt.savefig(save_path+"Voltages - All Fields - Delta.png",dpi=600)
# plt.savefig(save_path+"all fields2.pdf",dpi=600)




fig, ax1 = plt.subplots()
# for (index,field) in enumerate([fields[0],fields[1],fields[2],fields[4],fields[8],fields[12]]):
for (index,field) in enumerate(fields):
    plt.plot(amplitude[field][0],(amplitude[field][1]-np.average(amplitude[field][1]))/0.1,label=fieldlabels[index]+" Tesla")

ax1.set_xlabel("Angle $\\beta$ ($\degree$)",size=17)

ax1.set_ylabel("$\Delta$ $R_{xy}^{1\omega}$ (m$\Omega$)",size=16)
ax1.set_title(Scan+" Scan - 1$^{st}$ Harmonic Voltage vs. Angle",size=16)
plt.xticks([0,90,180,270,360],[0,90,180,270,360])
plt.legend()
plt.tight_layout()
plt.savefig(save_path+"Resistances_AllFields.png",dpi=600)
plt.close()




fitamplitudes=[]

for (index,field) in enumerate(fields):
    fitamplitudes.append(amplitude[field][2][0])
   
    
fitamplitudes=np.asfarray(fitamplitudes)/0.1  
fig, ax1 = plt.subplots()  
plt.plot(fields/10000,fitamplitudes,"o:",markersize=9,label="SMR Amplitude")

ax1.set_xlabel("External Field (Tesla)",size=17)

ax1.set_ylabel("SMR Amplitude (m$\Omega$)",size=16)
ax1.set_title(Scan+" Scan - SMR Amplitude vs. Angle",size=16)
plt.legend()
plt.tight_layout()
plt.savefig(save_path+"SMR Amplitude Field Dependence.png",dpi=600)
plt.close()