
# coding: utf-8





import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate
#from PyAstronomy import pyasl
import healpy as hp
from Parcel_healpy6 import parcel as p

    
     # In[9]:

if __name__ == '__main__':
    colors = ['r','b','k','g']
    linestyles = ['-','--','-.', ':']
    eps = [ 0, 2*np.pi/10, 2*np.pi,100]
    phi0 = [0, np.pi/3, np.pi/2, 3*np.pi/4]
    thetas = [np.pi/2.0, np.pi/4.0, np.pi/8.0]
    N=[10000,1000,600,500]
    
    #ns = np.arange(-2, 3, 0.1)
    #eps1 = 10**ns
    #phis = np.arange(0, np.pi/2,0.01)
    
    
    fig, ax = plt.subplots(figsize = (8,24))
    
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
        gas = p(e =0, epsilon=eps[i],wadv = 2, pmax = 6)
            
        t, T = gas.DE()
        for j in range(0,2):
            spot = hp.ang2pix(gas.NSIDE, np.pi/2 + np.pi/4*j,0)
            #ax0.plot(t*gas.wadv/(2*np.pi),T[:, spot, 2], color = colors[i], ls = linestyles[j])
            ax0.plot(t/gas.P,T[:, spot, 2], color = colors[i], ls = linestyles[j],linewidth = 2, 
                     label = "eps = "+"%2.1f"%(eps[i])+ ", Theta = "+ "%2.1f"%(np.pi/2 + np.pi/4*j))
            #ax0.plot(phi,T, color = colors[i], ls = linestyles[j])
            #"line"+str(i)= ax0.plot(phi,T, color = colors[i], ls = linestyles[j])
            #labels.append("line" + str[i])
            #names.append("eps = "+str(eps[i])+", \Theta = "+ str(thetas[j]) )
            
    ax0.legend(loc  = 4, fontsize = 7)
          
    ax0.set_xlim(4+1.5/2,5+1.5/2)
            
            #for tick in ax0.xaxis.get_ticklabels():
                #tick.set_fontsize('large')
            #   tick.set_weight('bold')
                
    ax0.set_xticks([4.75, 5, 5.25, 5.5, 5.75])
    ax0.set_xticklabels(['$-\pi/2$', '$0$', '$\pi/2$','$\pi$', '$3\pi/2$'], fontsize  = 14)
    #ax0.set_xlabel("Time (days) ", fontsize = 14)
    ax0.set_xlabel(" $\phi$ ", fontsize = 14)
    ax0.set_ylabel("$T/T0*sin\Theta^{1/4}$",fontsize = 14)
    ax0.set_title('Temperature of gas parcel from dawn till dawn - circular orbit')
    
    
    
    'plot 3'
    
    ee = 0.3
    gas = p(e = ee, epsilon=3, Porb = 3.2, pmax = 15)
    #gas.print_stuff()
    t, T = gas.DE()
    for i in range(2):
        
        #phi= np.linspace(0,2*30.0, num=300*30)
        #temp, = ax2.plot((phi*np.pi)/(gas.wadv*gas.P),T)
        #pos = np.where(np.abs(T[0,:, 0]- np.pi/2)<0.01)
        #print pos
        temp, = ax2.plot(t/gas.Porb,T[:,328+3*i,2])
    r = gas.radius
    
    #phi2 = self.wadv*t
    #phi = phi2/np.pi
    
    radius, = ax3.plot(t/gas.Porb,r/gas.AU, linestyle = '--', color = 'k')
    #ax2.set_xlim(2,8)
    ax2.legend([temp], ["Temperature of gas parcel "], loc  = 2, fontsize = 8)
    ax3.legend([radius],["Distance of planet from star"], loc = 1, fontsize =8)
    ax2.set_ylabel("T/T0", fontsize = 14)
    ax3.set_ylabel("r(t) in AU ",fontsize = 14)
    ax2.set_xlabel('Time (orbital periods)',fontsize = 14)
    
    ax3.set_title('Temperature of gas parcel - eccentric orbit (e = '+ str(ee)+')')
    
    plt.tight_layout()
    plt.show()
    #fig.savefig('Many parcels-hp.pdf')
    
    
    # In[9]:
    
'Prsentation circular case figure'
colors = ['r','b','k','g']
linestyles = ['-','--','-.', ':']
eps = [ 0, 2*np.pi/10, 2*np.pi,100]
phi0 = [0, np.pi/3, np.pi/2, 3*np.pi/4]
thetas = [np.pi/2.0, np.pi/4.0, np.pi/8.0]
N=[10000,1000,600,500]

fig, ax0 = plt.subplots(1,1,figsize = (12,6))



'Plot 1'
labels = []
names = []
#gas = p(e =0, epsilon=eps[1])
#days, T = gas.DE(pmax = 10, NSIDE = 8)
#print T

for i in range(len(eps)):
    gas = p(e =0, epsilon=eps[i],wadv = 2, pmax = 6)
        
    t, T = gas.DE()
    for j in range(0,2):
        spot = hp.ang2pix(gas.NSIDE, np.pi/2 + np.pi/4*j,0)
        #ax0.plot(t*gas.wadv/(2*np.pi),T[:, spot, 2], color = colors[i], ls = linestyles[j])
        ax0.plot(t/gas.P,T[:, spot, 2], color = colors[i], ls = linestyles[j],linewidth = 2, 
                 label = "eps = "+"%2.1f"%(eps[i])+ ", Theta = "+ "%2.1f"%(np.pi/2 + np.pi/4*j))
        #ax0.plot(phi,T, color = colors[i], ls = linestyles[j])
        #"line"+str(i)= ax0.plot(phi,T, color = colors[i], ls = linestyles[j])
        #labels.append("line" + str[i])
        #names.append("eps = "+str(eps[i])+", \Theta = "+ str(thetas[j]) )
        
ax0.legend(loc  = 4, fontsize = 7)
      
ax0.set_xlim(4+1.5/2,5+1.5/2)
        
        #for tick in ax0.xaxis.get_ticklabels():
            #tick.set_fontsize('large')
        #   tick.set_weight('bold')
            
ax0.set_xticks([4.75, 5, 5.25, 5.5, 5.75])
ax0.set_xticklabels(['$-\pi/2$ (Dawn)', '$0$ (Noon)', '$\pi/2$ (Dusk)' ,'$\pi$ (Midnight)', 
                     '$3\pi/2$ (Dawn)'], 
                    fontsize  = 12, fontweight = 'bold')
#ax0.set_xlabel("Time (days) ", fontsize = 14)
ax0.set_xlabel(" $\phi$ - Local Stellar Time", fontsize = 14, fontweight = 'bold')
ax0.set_ylabel("$T/T0*sin\Theta^{1/4}$",fontsize = 14,fontweight = 'bold')
ax0.set_title('Temperature of gas parcel from dawn till dawn - circular orbit')

fig.savefig('PRES_de_e=0.pdf')

# In[9]:

'Presentation ecentric case'



colors = ['r','b','k','g']
linestyles = ['-','--','-.', ':']
eps = [ 0, 2*np.pi/10, 2*np.pi,100]
phi0 = [0, np.pi/3, np.pi/2, 3*np.pi/4]
thetas = [np.pi/2.0, np.pi/4.0, np.pi/8.0]
N=[10000,1000,600,500]

#ns = np.arange(-2, 3, 0.1)
#eps1 = 10**ns
#phis = np.arange(0, np.pi/2,0.01)


fig, ax = plt.subplots(figsize = (14,5))

#ax0 = plt.subplot2grid((4,2), (0,0), rowspan = 2, colspan=2)
ax2 = plt.subplot2grid((2,2), (0,0), rowspan = 2, colspan=2)
#ax2 = plt.subplot2grid((4,2), (2,0), colspan=2)
#ax4 = plt.subplot2grid((4,2), (3,0), colspan=2)
ax3 = ax2.twinx()
#ax5 = ax4.twinx()

'plot 3'

ee = 0.3
gas = p(e = ee, epsilon=3, Porb = 3.2, pmax = 15)
#gas.print_stuff()
t, T = gas.DE()
for i in range(6):
    
    #phi= np.linspace(0,2*30.0, num=300*30)
    #temp, = ax2.plot((phi*np.pi)/(gas.wadv*gas.P),T)
    #pos = np.where(np.abs(T[0,:, 0]- np.pi/2)<0.01)
    #print pos
    if i == 0:
        temp, = ax2.plot(t/gas.Porb,T[:,328+3*i,2], linewidth = 2)
        
    else:
        temp, = ax2.plot(t/gas.Porb,T[:,328+3*i,2])
        
r = gas.radius

#phi2 = self.wadv*t
#phi = phi2/np.pi

radius, = ax3.plot(t/gas.Porb,r/gas.AU, linestyle = '--', color = 'k')
ax2.set_ylim(0.2,1.2)
ax2.legend([temp], ["Temperature of gas parcel "], loc  = 2, fontsize = 8)
ax3.legend([radius],["Distance of planet from star"], loc = 1, fontsize =8)
ax2.set_ylabel("T/T0", fontsize = 14, fontweight = 'bold')
ax3.set_ylabel("r(t) in AU ",fontsize = 14, fontweight = 'bold')
ax2.set_xlabel('Time (orbital periods)',fontsize = 14, fontweight = 'bold')

ax3.set_title('Temperature of gas parcel - eccentric orbit (e = '+ str(ee)+')')

plt.tight_layout()
plt.show()
    
fig.savefig('PRES_de_eNOT0.pdf')
    
    

# In[9]:

taurads = np.linspace(0.5,20, num = 20)*2
wadv = np.linspace (-10,10, num =20)



'plot 3'

ee = 0.3

for i in range(len(taurads)):
    fig, ax = plt.subplots(figsize = (6,4))

    #ax0 = plt.subplot2grid((4,2), (0,0), rowspan = 2, colspan=2)
    ax2 = plt.subplot2grid((2,2), (0,0), rowspan = 2, colspan=2)
    #ax2 = plt.subplot2grid((4,2), (2,0), colspan=2)
    #ax4 = plt.subplot2grid((4,2), (3,0), colspan=2)
    ax3 = ax2.twinx()
    #ax5 = ax4.twinx()
    gas = p(e = ee, wadv = 2, tau_rad = taurads[i], Porb = 3.2, pmax = 7)
    #gas.print_stuff()
    t, T = gas.DE()
    
    #phi= np.linspace(0,2*30.0, num=300*30)
    #temp, = ax2.plot((phi*np.pi)/(gas.wadv*gas.P),T)
    #pos = np.where(np.abs(T[0,:, 0]- np.pi/2)<0.01)
    #print pos
    
    temp, = ax2.plot(t/gas.Porb,T[:,368,2], linewidth = 2, color = 'orangered')
        
    r = gas.radius
    
    #phi2 = self.wadv*t
    #phi = phi2/np.pi
    
    radius, = ax3.plot(t/gas.Porb,r/gas.AU, linestyle = '--', color = 'k')
    ax2.set_ylim(0.2,1.2)
    ax2.set_xlim(3,6)
    ax2.legend([temp], ["Temp, tau_rad = "+"%2.1f"%taurads[i]], loc  = 2, fontsize = 12)
    #ax3.legend([radius],["Distance of planet from star"], loc = 1, fontsize =8)
    ax2.set_ylabel("T/T0", fontsize = 14, fontweight = 'bold')
    ax3.set_ylabel("r(t) in AU ",fontsize = 14, fontweight = 'bold')
    ax2.set_xlabel('Time (orbital periods)',fontsize = 14, fontweight = 'bold')
    
    ax3.set_title('Temperature of gas parcel - eccentric orbit (e = '+ str(ee)+')')
    
    plt.tight_layout()
    #plt.show()
        
    fig.savefig('PRES_taurad#'+str(i)+'.png')
    



# In[9]:
taurads = np.linspace(0.25,20, num = 20)*2
wadv = np.linspace (-1,3, num =20)



'plot 3'

ee = 0.3

for i in range(len(wadv)):
    fig, ax = plt.subplots(figsize = (6,4))

    #ax0 = plt.subplot2grid((4,2), (0,0), rowspan = 2, colspan=2)
    ax2 = plt.subplot2grid((2,2), (0,0), rowspan = 2, colspan=2)
    #ax2 = plt.subplot2grid((4,2), (2,0), colspan=2)
    #ax4 = plt.subplot2grid((4,2), (3,0), colspan=2)
    ax3 = ax2.twinx()
    #ax5 = ax4.twinx()
    gas = p(e = ee, wadv =wadv[i] , tau_rad = 0.5, Porb = 3.2, pmax = 7)
    #gas.print_stuff()
    t, T = gas.DE()
    
    #phi= np.linspace(0,2*30.0, num=300*30)
    #temp, = ax2.plot((phi*np.pi)/(gas.wadv*gas.P),T)
    #pos = np.where(np.abs(T[0,:, 0]- np.pi/2)<0.01)
    #print pos
    
    temp, = ax2.plot(t/gas.Porb,T[:,368,2], linewidth = 2, color = 'cornflowerblue')
        
    r = gas.radius
    
    #phi2 = self.wadv*t
    #phi = phi2/np.pi
    
    radius, = ax3.plot(t/gas.Porb,r/gas.AU, linestyle = '--', color = 'k', alpha = 0.6)
    ax2.set_ylim(0.0,1.2)
    ax2.set_xlim(3,6)
    ax2.legend([temp], ["Temp, wadv = "+"%2.1f"%wadv[i]], loc  = 2, fontsize = 12)
    #ax3.legend([radius],["Distance of planet from star"], loc = 1, fontsize =8)
    ax2.set_ylabel("T/T0", fontsize = 14, fontweight = 'bold')
    ax3.set_ylabel("r(t) in AU ",fontsize = 14, fontweight = 'bold')
    ax2.set_xlabel('Time (orbital periods)',fontsize = 14, fontweight = 'bold')
    
    ax3.set_title('Temperature of gas parcel - eccentric orbit (e = '+ str(ee)+')')
    
    plt.tight_layout()
    #plt.show()
        
    fig.savefig('PRES_wadv#'+str(i)+'.png')
    


# In[9]:    
'''TEST DE 2'''


"""TO DO: figure out how t make a legend for the first figure  
Also, make sure this actually works. it's got problems """



colors = ['r','b','k','g']
linestyles = ['-','--','-.']
eps = [ 0.1, 2*np.pi/10, 2*np.pi,1000]
#phi0 = [0, np.pi/2, 3*np.pi/4]
#thetas = [np.pi/2.0, np.pi/4.0, np.pi/8.0]
N=[10000,1000,600,500]

ns = np.arange(-2, 3, 0.5)
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
        gas = p(e =0, wadv = 1, epsilon=eps[i], Porb = 3.0)
        #gas.print_stuff()
        t, d = gas.DE()
        #plt.plot(d[(gas.pmaxi-2.5)*gas.stepsi:(gas.pmaxi-1.5)*gas.stepsi:,338 - 150*j,1],
                   #d[(gas.pmaxi-2.5)*gas.stepsi:(gas.pmaxi-1.5)*gas.stepsi:,338-150*j,2], 
                     #color = colors[i], ls = linestyles[j])
                     
        ax0.plot(t[(gas.pmaxi-2.25)*gas.stepsi:(gas.pmaxi-1.25)*gas.stepsi:]/gas.P,
                   d[(gas.pmaxi-2.25)*gas.stepsi:(gas.pmaxi-1.25)*gas.stepsi:,336-224*j,2], 
                     color = colors[i], ls = linestyles[j])
        #"line"+str(i)= ax0.plot(phi,T, color = colors[i], ls = linestyles[j])
        #labels.append("line" + str[i])
        #names.append("eps = "+str(eps[i])+", \Theta = "+ str(thetas[j]) )
       
#ax0.legend(labels, names, loc  = 4)      
ax0.set_xlim(0.75,1.75)
        
        #for tick in ax0.xaxis.get_ticklabels():
            #tick.set_fontsize('large')
        #   tick.set_weight('bold')
            
ax0.set_xticks([0.75, 1, 1.25, 1.5, 1.75])
ax0.set_xticklabels(['$-\pi/2$', '$0$','$\pi/2$','$\pi$','$3\pi/2$'], fontsize  = 14)
ax0.set_xlabel(" $\phi$ ", fontsize = 14)
ax0.set_ylabel("$T/T0*sin\Theta^{1/4}$",fontsize = 14)
ax0.set_title('Temperature of gas parcel from dawn till dawn - circular orbit')

    
    # In[9]:
'Plot 2'
ns = np.arange(-2, 3, 0.5)
eps1 = 10**ns
Tm =[]
Tda = []
Tdu = []

w = np.linspace(1, 20, num = len(eps))

for i in range(len(eps1)):
    
        gas = p(epsilon=eps1[i], wadv = 11, Porb = 3.0, pmax = 40)
        
        Tmaxi, Tdawni, Tduski  = gas.findT()
        Tm.append(Tmaxi)
        Tda.append(Tdawni)
        Tdu.append(Tduski)
        
        
maxim, = plt.semilogx(eps1,gas.Tmax(eps1), color = 'r')
dawn, = plt.semilogx(eps1,gas.Tdawn(eps1),color = 'b')
dusk, = plt.semilogx(eps1,gas.Tdusk(eps1),color = 'k')
plt.semilogx(eps1,Tm, 'ro',linestyle = 'None')
plt.semilogx(eps1,Tda, 'bD',linestyle = 'None')
plt.semilogx(eps1,Tdu, 'k*',linestyle = 'None')

    
    # In[9]:
ax1.legend([maxim, dusk, dawn], ["$T_{max}$","$T_{dusk}$", "$T_{dawn}$"], loc  = 4)
ax1.set_ylim(0,1.1)
ax1.set_xlabel("$\epsilon$", fontsize = 14)
ax1.set_ylabel("T/T0", fontsize = 14)

ax1.set_title('$T_{max}$, $T_{dusk}$, $T_{dawn}$ as functions of $\epsilon$')

    
    # In[9]:
'plot 3'

ee = 0.3
gas = p(e = ee, epsilon=3, wadv = 1.1, Porb = 3.2, pmax = 10)
#gas.print_stuff()
for i in range(8):
    t, d = gas.DE()
    #phi= np.linspace(0,2*30.0, num=300*30)
    #temp, = ax2.plot((phi*np.pi)/(gas.wadv*gas.P),T)
    temp, = ax2.plot(t,d[:,368-30*i,2])

t,r = gas.t, gas.radius

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
gas = p(e = ee, epsilon=3, wadv = 1.1, Porb = 2.5, pmax = 10)
#gas.print_stuff()
for i in range(5):
    t, d = gas.DE()
    #phi= np.linspace(0,2*30.0, num=300*30)
    #temp, = ax2.plot((phi*np.pi)/(gas.wadv*gas.P),T)
    temp, = ax4.plot(t,d[:,368-30*i,2])
#phi, days, T = gas.DE(pmax = 12, steps = 500)
t,r = gas.t, gas.radius

#temp, = ax4.plot(days,T, color ='g')
#temp, = ax4.plot((phi*np.pi)/(gas.wadv*gas.P),T, color ='g')


radius, = ax5.plot(t,r/gas.AU, linestyle = '--', color = 'k')

ax4.legend([temp], ["Temperature of gas parcel "], loc  = 2)
ax5.legend([radius],["Distance of planet from star"], loc = 1)
ax4.set_ylabel("T/T0", fontsize = 14)
ax5.set_ylabel("r(t) in AU ",fontsize = 14)
ax4.set_xlabel('Time ($planet \, days -- P_{rot}$)',fontsize = 14)

ax5.set_title('Temperature of gas parcel - eccentric orbit (e = '+ str(ee)+')')

plt.tight_layout()
plt.show()
fig.savefig('gas_parcel')     



        #plt.xlim(0,0.2)
#plt.ylim(0.1,1.5)

#plt.ylim(0.2,1.4)
#plt.ylim(0.4,0.8)
#dawn = (2*np.pi/gas.wadv)/4*np.array([3,7,11,15])
#dusk = (2*np.pi/gas.wadv)/4 *np.array([1,5,9,13])

#plt.vlines(dawn,0,1)
#plt.vlines(dusk,0,1,'r')



