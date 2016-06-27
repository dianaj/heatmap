
# coding: utf-8





import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate
#from PyAstronomy import pyasl
#import healpy as hp
from Parcel_healpy4 import parcel as p

if __name__ == '__main__':
    colors = ['r','b','k','g']
    linestyles = ['-','--','-.', ':']
    eps = [ 0.1, 2*np.pi/10, 2*np.pi,1000]
    phi0 = [0, np.pi/3, np.pi/2, 3*np.pi/4]
    thetas = [np.pi/2.0, np.pi/4.0, np.pi/8.0]
    N=[10000,1000,600,500]
    
    #ns = np.arange(-2, 3, 0.1)
    #eps1 = 10**ns
    #phis = np.arange(0, np.pi/2,0.01)
    
    
    fig, ax = plt.subplots(figsize = (16,16))
    
    ax0 = plt.subplot2grid((4,2), (0,0), rowspan = 2, colspan=2)
    ax2 = plt.subplot2grid((4,2), (2,0), rowspan = 2, colspan=2)
    #ax2 = plt.subplot2grid((4,2), (2,0), colspan=2)
    #ax4 = plt.subplot2grid((4,2), (3,0), colspan=2)
    ax3 = ax2.twinx()
    #ax5 = ax4.twinx()
    
    'Plot 1'
    labels = []
    names = []
    #gas = p(e =0, epsilon=eps[1])
    #days, T = gas.DE(pmax = 10, NSIDE = 8)
    #print T
    
    for i in range(len(eps)):
        for j in range(1,4):
            gas = p(e =0, epsilon=eps[i])
            
            days, T = gas.DE(pmax = 5, NSIDE = 8)
            ax0.plot(days,T[:, j*150, 2], color = colors[i], ls = linestyles[j])
            #ax0.plot(phi,T, color = colors[i], ls = linestyles[j])
            #"line"+str(i)= ax0.plot(phi,T, color = colors[i], ls = linestyles[j])
            #labels.append("line" + str[i])
            #names.append("eps = "+str(eps[i])+", \Theta = "+ str(thetas[j]) )
            
    #ax0.legend(labels, names, loc  = 4)
            
    ax0.set_xlim(2.5,3.5)
            
            #for tick in ax0.xaxis.get_ticklabels():
                #tick.set_fontsize('large')
            #   tick.set_weight('bold')
                
    #ax0.set_xticks([1.5, 2, 2.5, 3, 3.5])
    #ax0.set_xticklabels(['$3\pi/2$', '$2\pi$', '$5\pi/2$','$3\pi$', '$7\pi/2$'], fontsize  = 14)
    ax0.set_xlabel("Time (days) ", fontsize = 14)
    #ax0.set_xlabel(" $\phi$ ", fontsize = 14)
    ax0.set_ylabel("$T/T0*sin\Theta^{1/4}$",fontsize = 14)
    ax0.set_title('Temperature of gas parcel from dawn till dawn - circular orbit')
    
    
    
    'plot 3'
    
    ee = 0.3
    gas = p(e = ee, epsilon=3, Porb = 3.2)
    #gas.print_stuff()
    days, T = gas.DE(pmax = 20, steps = 300, NSIDE = 8)
    for i in range(8):
        
        #phi= np.linspace(0,2*30.0, num=300*30)
        #temp, = ax2.plot((phi*np.pi)/(gas.wadv*gas.P),T)
        pos = np.where(np.abs(T[0,:, 0]- np.pi/2)<0.01)
        #print pos
        temp, = ax2.plot(days,T[:,328+3*i,2])
    t,r = gas.radius(pmax = 20)
    
    #phi2 = self.wadv*t
    #phi = phi2/np.pi
    
    radius, = ax3.plot(t/gas.P,r/gas.AU, linestyle = '--', color = 'k')
    #ax2.set_xlim(2,8)
    ax2.legend([temp], ["Temperature of gas parcel "], loc  = 2)
    ax3.legend([radius],["Distance of planet from star"], loc = 1)
    ax2.set_ylabel("T/T0", fontsize = 14)
    ax3.set_ylabel("r(t) in AU ",fontsize = 14)
    ax2.set_xlabel('Time ($planet \, days -- P_{rot}$)',fontsize = 14)
    
    ax3.set_title('Temperature of gas parcel - eccentric orbit (e = '+ str(ee)+')')
    
    fig.savefig('Many parcels-hp')
    
    
    # In[9]:
    
    """TO DO: figure out how t make a legend for the first figure  """
    """I changes something here. theta is no longer defines as a property of the gas but has to be 
    specified as an argument of the DE"""
    
    """
    colors = ['r','b','k','g']
    linestyles = ['-','--','-.']
    eps = [ 0.1, 2*np.pi/10, 2*np.pi,1000]
    phi0 = [0, np.pi/2, 3*np.pi/4]
    thetas = [np.pi/2.0, np.pi/4.0, np.pi/8.0]
    N=[10000,1000,600,500]
    
    ns = np.arange(-2, 3, 0.1)
    eps1 = 10**ns
    phis = np.arange(0, np.pi/2,0.01)
    
    
    fig, ax = plt.subplots(figsize = (12,24))
    
    ax0 = plt.subplot2grid((4,2), (0,0), colspan=2)
    ax1 = plt.subplot2grid((4,2), (1,0), colspan=2)
    ax2 = plt.subplot2grid((4,2), (2,0), colspan=2)
    ax4 = plt.subplot2grid((4,2), (3,0), colspan=2)
    ax3 = ax2.twinx()
    ax5 = ax4.twinx()
    
    'Plot 1'
    labels = []
    names = []
    for i in range(len(eps)-1):
        for j in range(2):
            gas = p(e =0, wadv = 1.1, epsilon=eps[i], theta = thetas[j])
            #gas.print_stuff()
            phi, days, T = gas.DE(pmax = 3.0)
            ax0.plot(phi,T, color = colors[i], ls = linestyles[j])
            #"line"+str(i)= ax0.plot(phi,T, color = colors[i], ls = linestyles[j])
            #labels.append("line" + str[i])
            #names.append("eps = "+str(eps[i])+", \Theta = "+ str(thetas[j]) )
            
    #ax0.legend(labels, names, loc  = 4)
            
    ax0.set_xlim(3.0/2,7.0/2)
            
            #for tick in ax0.xaxis.get_ticklabels():
                #tick.set_fontsize('large')
            #   tick.set_weight('bold')
                
    ax0.set_xticks([1.5, 2, 2.5, 3, 3.5])
    ax0.set_xticklabels(['$3\pi/2$', '$2\pi$', '$5\pi/2$','$3\pi$', '$7\pi/2$'], fontsize  = 14)
    ax0.set_xlabel(" $\phi$ ", fontsize = 14)
    ax0.set_ylabel("$T/T0*sin\Theta^{1/4}$",fontsize = 14)
    ax0.set_title('Temperature of gas parcel from dawn till dawn - circular orbit')
    
    
    'Plot 2'
    Tm =[]
    Tda = []
    Tdu = []
    
    w = np.linspace(1, 20, num = len(eps))
    for i in range(len(eps1)):
        
            gas = p(epsilon=eps1[i], wadv = 1.1)
            
            Tmaxi, Tdawni, Tduski  = gas.findT(pmax = 10, steps =300)
            Tm.append(Tmaxi)
            Tda.append(Tdawni)
            Tdu.append(Tduski)
            
            
    maxim, = ax1.semilogx(eps1,Tmax(eps1), color = 'r')
    dawn, = ax1.semilogx(eps1,Tdawn(eps1),color = 'b')
    dusk, = ax1.semilogx(eps1,Tdusk(eps1),color = 'k')
    ax1.semilogx(eps1,Tm, 'ro',linestyle = 'None')
    ax1.semilogx(eps1,Tda, 'bD',linestyle = 'None')
    ax1.semilogx(eps1,Tdu, 'k*',linestyle = 'None')
    ax1.legend([maxim, dusk, dawn], ["$T_{max}$","$T_{dusk}$", "$T_{dawn}$"], loc  = 4)
    ax1.set_ylim(0,1.1)
    ax1.set_xlabel("$\epsilon$", fontsize = 14)
    ax1.set_ylabel("T/T0", fontsize = 14)
    
    ax1.set_title('$T_{max}$, $T_{dusk}$, $T_{dawn}$ as functions of $\epsilon$')
    
    
    'plot 3'
    
    ee = 0.3
    gas = p(e = ee, epsilon=3, wadv = 1.1, theta = np.pi/2, Porb = 3.2)
    #gas.print_stuff()
    for i in range(8):
        phi, days, T = gas.DE(pmax = 30, steps = 300, phi0=(2.0*np.pi/8.0)*i)
        #phi= np.linspace(0,2*30.0, num=300*30)
        #temp, = ax2.plot((phi*np.pi)/(gas.wadv*gas.P),T)
        temp, = ax2.plot(days,T)
    
    t,r = gas.radius(pmax = 30)
    
    #phi2 = self.wadv*t
    #phi = phi2/np.pi
    
    radius, = ax3.plot(t/gas.P,r/gas.AU, linestyle = '--', color = 'k')
    #ax2.set_xlim(0,15)
    ax2.legend([temp], ["Temperature of gas parcel "], loc  = 2)
    ax3.legend([radius],["Distance of planet from star"], loc = 1)
    ax2.set_ylabel("T/T0", fontsize = 14)
    ax3.set_ylabel("r(t) in AU ",fontsize = 14)
    ax2.set_xlabel('Time ($planet \, days -- P_{rot}$)',fontsize = 14)
    
    ax3.set_title('Temperature of gas parcel - eccentric orbit (e = '+ str(ee)+')')
    
    'plot 4'
    
    ee = 0.5
    gas = p(e = ee, epsilon=3, wadv = 1.1, theta = np.pi/2, Porb = 2.5)
    #gas.print_stuff()
    for i in range(5):
        phi, days, T = gas.DE(pmax = 15, steps = 300, phi0=(2.0*np.pi/5.0)*i)
        #phi= np.linspace(0,2*30.0, num=300*30)
        #temp, = ax2.plot((phi*np.pi)/(gas.wadv*gas.P),T)
        temp, = ax4.plot(days,T)
    #phi, days, T = gas.DE(pmax = 12, steps = 500)
    t,r = gas.radius(pmax = 15)
    
    #temp, = ax4.plot(days,T, color ='g')
    #temp, = ax4.plot((phi*np.pi)/(gas.wadv*gas.P),T, color ='g')
    
    
    radius, = ax5.plot(t/gas.P,r/gas.AU, linestyle = '--', color = 'k')
    
    ax4.legend([temp], ["Temperature of gas parcel "], loc  = 2)
    ax5.legend([radius],["Distance of planet from star"], loc = 1)
    ax4.set_ylabel("T/T0", fontsize = 14)
    ax5.set_ylabel("r(t) in AU ",fontsize = 14)
    ax4.set_xlabel('Time ($planet \, days -- P_{rot}$)',fontsize = 14)
    
    ax5.set_title('Temperature of gas parcel - eccentric orbit (e = '+ str(ee)+')')
    
    fig.savefig('gas_parcel')     
    
    
    
            #plt.xlim(0,0.2)
    #plt.ylim(0.1,1.5)
    
    #plt.ylim(0.2,1.4)
    #plt.ylim(0.4,0.8)
    #dawn = (2*np.pi/gas.wadv)/4*np.array([3,7,11,15])
    #dusk = (2*np.pi/gas.wadv)/4 *np.array([1,5,9,13])
    
    #plt.vlines(dawn,0,1)
    #plt.vlines(dusk,0,1,'r')
    
    """
    
    
