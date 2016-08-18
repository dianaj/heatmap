# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 19:01:50 2016

@author: diana
"""

#problem : this is not the damn phase angle where the flux is a maximum, it's
#the phi max approximation 


import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from Parcel_healpy6 import parcel as p


# In[10]:    
gas = p(name = "HD189733b",Teff = 5938, e=0,Porb = 2.21, a = 0.031, wadv = 11.0, 
                          epsilon = 6, argp = 90, Rstar = 0.781, Mstar = 0.846,
                          Rplanet = 0.838, pmax = 10, steps = 300, NSIDE = 8 )      
    
t, d, Fwv = gas.Fobs(BOLO = True)

plt.plot(t,Fwv)

# In[10]: 

if __name__ == '__main__':
    

    epsi = 100.0/(10**(np.linspace(0,3, num = 20)))
    eps = np.concatenate((epsi,[0]), axis = 0)
    maxf = []
    minf = []
    phi = []
    for i in range(len(eps)):
        if eps[i] >= 15:
            gas = p(name = "HD189733b",Teff = 5938, e=0,Porb = 2.21, a = 0.031, wadv = 11.0, 
                          epsilon = eps[i], argp = 90, Rstar = 0.781, Mstar = 0.846,
                          Rplanet = 0.838, pmax = 60, steps = 500, NSIDE = 16 )        
        elif 15 > eps[i] > 10:
            gas = p(name = "HD189733b",Teff = 5938, e=0,Porb = 2.21, a = 0.031, wadv = 11.0, 
                          epsilon = eps[i], argp = 90, Rstar = 0.781, Mstar = 0.846,
                          Rplanet = 0.838, pmax = 30, steps = 500, NSIDE = 16 )
                          
        else : 
            gas = p(name = "HD189733b",Teff = 5938, e=0,Porb = 2.21, a = 0.031, wadv = 11.0, 
                          epsilon = eps[i], argp = 90, Rstar = 0.781, Mstar = 0.846,
                          Rplanet = 0.838, pmax = 12, steps = 800, NSIDE = 16 )
        #Fwv actually means total flux now, not wavelenght.
        t, d, Fwv = gas.Fobs(BOLO = True)
        backend = int((gas.pmaxi-1.5)*gas.stepsi)
        maxf.append(max(Fwv[backend::]))
        minf.append(min(Fwv[backend::]))
        #phi.append(d[np.argmax(Fwv),np.argmax(d[np.argmax(Fwv),:,2]),1])
    
        
        #phi.append(gas.alpha[np.argmax(Fwv[backend::])])
        maxphi = (gas.t[np.argmax(Fwv[backend::])+backend]/gas.P)%1*(2*np.pi)
        phi.append(maxphi)
        
    
    maxf = np.array(maxf)
    minf = np.array(minf)
    phi = np.array(phi).flatten()
   
    phiT = []
    for i in range(len(eps)):    
        phiT.append(gas.phi_max(eps[i]))
    
    phiT = np.array(phiT).flatten()
    print phi
    

# In[10]: 
for i in [5]:
    gas = p(name = "HD189733b",Teff = 4938, e=0,Porb = 2.21, a = 0.031, wadv = 11.0, 
                              epsilon = eps[i], argp = 90, Rstar = 0.781, Mstar = 0.846,
                              Rplanet = 1.138, pmax = 20, steps = 800, NSIDE = 32 )
    t, d, Fwv = gas.Fobs()
    
    
    backend = int((gas.pmaxi-1.5)*gas.stepsi)
    maxf[i]  = (max(Fwv[backend::]))
    minf[i] = (min(Fwv[backend::]))
    phi[i] = (gas.t[np.argmax(Fwv[backend::])+backend]/gas.P)%1*(2*np.pi)



 # In[10]: 
np.savetxt('hotspot(try2)_NSIDE'+str(gas.NSIDE)+'BOLO'+'_e='+str(gas.e)+'.out', 
                np.array(zip(eps,maxf,minf,phi,phiT)).reshape(-1,5), 
                header = "eps,max_flux, min_flux, phi_max, phi_theory")

# In[10]: 

np.savetxt('hotspot(NEW)_NSIDE'+str(gas.NSIDE)+'BOLO'+'_e='+str(gas.e)+'.out', 
                np.array(zip(eps,maxf,minf,phi,phiT)).reshape(-1,5), 
                header = "eps,max_flux, min_flux, phi_max, phi_theory")
                
                

 # In[10]: 'DIAG FIGURE'

import scipy

def f(x, a, b, c, d):
    return a*np.arctan(b-x*c)+d
    
flux_contrast = (maxf-minf)/ maxf
popt, pcov = scipy.optimize.curve_fit(f, flux_contrast,-(phi-2*np.pi)*180/(np.pi)-180)
   
fig, ax = plt.subplots(3,1,figsize = (12,16))
flux_contrast = (maxf-minf)/ maxf
ax[0].plot(eps, flux_contrast, linestyle = '--', color = 'green',marker='o',
     markerfacecolor='blue', markersize=6)
ax[0].set_ylabel('Flux Contrast ')
ax[0].set_xlabel('epsilon')
'(-phi+180)*np.pi/180'
ax[1].plot(eps,-(phi-2*np.pi)*180/(np.pi)-180 , linestyle ='--', color = 'green',marker='o',
     markerfacecolor='blue', markersize=6)
     
ax[1].plot(eps,(phiT)*180/(np.pi) , linestyle ='--', color = 'green',marker='o',
     markerfacecolor='red', markersize=6)
     
ax[1].set_ylabel('phi_max (degrees)')
ax[1].set_xlabel('epsilon')

#ax[2].plot((maxf-minf)/maxf, phi*180/np.pi, color='green', linestyle='dashed', marker='o',
#     markerfacecolor='blue', markersize=8, label = 'Flux contrast vs. phase offset (model)')
ax[2].plot(flux_contrast[::], -(phi-2*np.pi)[::]*180/(np.pi)-180, color='green', linestyle='dashed', marker='o',
     markerfacecolor='blue', markersize=8, label = 'Flux contrast vs. phase offset (model)')     
ax[2].plot(flux_contrast, f(flux_contrast, popt[0],popt[1],popt[2],popt[3]), 'r',
label = '%2.2f'%popt[0]+'*ArcTan('+'%2.2f'%popt[1]+'-x*'+'%2.2f'%popt[2]+')+'+'%2.2f'%popt[3] )

ax[2].plot(flux_contrast[::], (phiT)[::]*180/(np.pi), color='green', linestyle='dashed', marker='o',
     markerfacecolor='red', markersize=6, label = 'Flux contrast vs. phase offset (Theoretical)')

ax[2].set_ylabel('phi_max (degrees)')
ax[2].set_xlabel('Flux Contrast')
ax[2].legend(fontsize = 6)
plt.tight_layout()


fig.savefig('hotspot_flux2.pdf')

# In[10]:

'FITTING FIGURE'

fig, ax = plt.subplots(2,1, figsize = (10,16))


def f(x, a, b, c, d):
    return a*np.arctan(b-x*c)+d
    
flux_contrast = (maxf-minf)/ maxf
popt, pcov = scipy.optimize.curve_fit(f, flux_contrast,-(phi-2*np.pi)*180/(np.pi)-180)
ax[0].plot(flux_contrast[::], -(phi-2*np.pi)[::]*180/(np.pi)-180, color='green', linestyle='dashed', marker='o',
     markerfacecolor='blue', markersize=8, label = 'Flux contrast vs. phase offset (model)')     
ax[0].plot(flux_contrast, f(flux_contrast, popt[0],popt[1],popt[2],popt[3]), 'r',
label = '%2.2f'%popt[0]+'*ArcTan('+'%2.2f'%popt[1]+'-x*'+'%2.2f'%popt[2]+')+'+'%2.2f'%popt[3] )

ax[0].plot(flux_contrast[::], (phiT)[::]*180/(np.pi), color='green', linestyle='dashed', marker='o',
     markerfacecolor='red', markersize=6, label = 'Flux contrast vs. phase offset (Theoretical)')

ax[0].set_ylabel('phi_max (degrees)')
ax[0].set_xlabel('Flux Contrast')
ax[0].legend(fontsize = 6)

def f2(x, a, b, c):
    return a*x**b+c
    

popt2, pcov2 = scipy.optimize.curve_fit(f2, flux_contrast[4::],(-(phi-2*np.pi)*180/(np.pi)-180)[4::], p0 = [-80,2,80])

ax[1].plot(flux_contrast[::]**popt2[1], -(phi-2*np.pi)[::]*180/(np.pi)-180, color='green', linestyle='dashed', marker='o',
     markerfacecolor='blue', markersize=8, label = 'Flux contrast vs. phase offset (model)')     
ax[1].plot(flux_contrast**popt2[1], f2(flux_contrast, popt2[0],popt2[1],popt2[2]), 'r',
label = '%2.2f'%popt2[0]+'* x^'+'%2.2f'%popt2[1]+'+'+'%2.2f'%popt2[2] )

ax[1].plot(flux_contrast[::]**popt2[1], (phiT)[::]*180/(np.pi), color='green', linestyle='dashed', marker='o',
     markerfacecolor='red', markersize=6, label = 'Flux contrast vs. phase offset (Theoretical)')

ax[1].set_ylabel('phi_max (degrees)')

ax[1].set_xlabel('Flux Contrast ^ '+'%1.2f'%popt2[1])
ax[1].legend(fontsize = 6)

fig.savefig('phimax_fitting.pdf')
plt.show()




# In[10]:

'FIGURE: actual flux picture and vlines'

from matplotlib import colors
colorss = np.array((colors.cnames.keys()))[(eps.astype(int))]

fig, ax = plt.subplots(2,1, figsize = (10, 16))

ax[0].plot(flux_contrast[::], -(phi-2*np.pi)[::]*180/(np.pi)-180, color='green', linestyle='dashed', marker='o',
     markerfacecolor='blue', markersize=8, label = 'Flux contrast vs. phase offset (model)')     
#ax[2].plot(flux_contrast, f(flux_contrast, popt[0],popt[1],popt[2],popt[3]), 'r',
#label = '%2.2f'%popt[0]+'*ArcTan('+'%2.2f'%popt[1]+'-x*'+'%2.2f'%popt[2]+')+'+'%2.2f'%popt[3] )

ax[0].plot(flux_contrast[::], (phiT)[::]*180/(np.pi), color='green', linestyle='dashed', marker='o',
     markerfacecolor='red', markersize=6, label = 'Flux contrast vs. phase offset (Theoretical)')

ax[0].set_ylabel('phi_max (degrees)')
ax[0].set_xlabel('Flux Contrast')
ax[0].legend(fontsize = 7)

for i in range(len(eps)):
        gas = p(name = "HD189733b",Teff = 4938, e=0,Porb = 2.21, a = 0.031, wadv = 11.0, 
                      epsilon = eps[i], argp = 90, Rstar = 0.781, Mstar = 0.846,
                      Rplanet = 1.138, pmax = 30, steps = 500, NSIDE = 4 )
        t, d, Fwv = gas.Fobs()
        backend = int((gas.pmaxi-1.5)*gas.stepsi)
        ax[1].plot((t/gas.Porb)[backend::],Fwv[backend::], color = colorss[i], label = 'eps = '+'%3.2f'%gas.epsilon)

        ax[1].vlines(t[np.argmax(Fwv[backend::])+backend]/gas.Porb, 0, max(Fwv[backend::]), color = colorss[i])
        
for i in range(len(eps)):
        gas = p(name = "HD189733b",Teff = 4938, e=0,Porb = 2.21, a = 0.031, wadv = 2.0, 
                      epsilon = eps[i], argp = 90, Rstar = 0.781, Mstar = 0.846,
                      Rplanet = 1.138, pmax = 30, steps = 500, NSIDE = 4 )
        t, d, Fwv = gas.Fobs()
        backend = int((gas.pmaxi-1.5)*gas.stepsi)
        ax[1].plot((t/gas.Porb)[backend::],Fwv[backend::], color = colorss[i], label = 'eps = '+'%3.2f'%gas.epsilon)

        ax[1].vlines(t[np.argmax(Fwv[backend::])+backend]/gas.Porb, 0, max(Fwv[backend::]), color = colorss[i])
ax[1].legend(loc = 2, fontsize = 7)        
plt.tight_layout()

fig.savefig('hotspot_flux3try.pdf')


# In[10]:
'ANALYTIC TEMP STUFF'

gas = p(name = "HD189733b",Teff = 5938, e=0,Porb = 2.21, a = 0.031, wadv = 11.0, 
                          epsilon = 2, argp = 90, Rstar = 0.781, Mstar = 0.846,
                          Rplanet = 0.838, pmax = 60, steps = 500, NSIDE = 16 ) 
phiT1 = []
Tmax = []
Tmin = []
eps1 = eps
for i in range(len(eps1)):    
    phiT1.append(gas.phi_max(eps1[i]))

    Tmax.append(np.cos(gas.phi_max(eps1[i]))**(0.25))
    if eps1[i] == 0:
        Tmin.append(0)
    else:
        Tmin.append((np.pi+ (3*np.pi/eps1[i])**(4/3))**(-0.25))
    

Tmax = np.array(Tmax)
Tmin = np.array(Tmin)

#Amp = (Tmax - Tmin)/Tmax

Amp = (Tmax - Tmin)
phiT1 = np.array(phiT1)




 # In[10]:
import scipy 
'IMPORTANT FIGURE'

#File ='hotspot(try2)_NSIDE16BOLO_e=0.out' 
File ='hotspot(NEW)_NSIDE16BOLO_e=0.out' 

eps = np.loadtxt(File, usecols = [0], skiprows = 1)
#eps = np.delete(eps,[4,5])
max_flux = np.loadtxt(File, usecols = [1], skiprows = 1)
min_flux = np.loadtxt(File, usecols = [2], skiprows = 1)

flux_contrast = (max_flux - min_flux)/max_flux
area = np.pi*gas.Rplanet**2
T_contrast = ((max_flux/(gas.sigmaB*area))**0.25-(min_flux/(gas.sigmaB*area))**0.25)/gas.T0
#flux_contrast = np.delete(flux_contrast,[4,5])
phi = np.loadtxt(File, usecols = [3], skiprows = 1)
#phi[5] = phi[5]-0.05
#phi = np.delete(phi,[4,5])
phiT = np.loadtxt(File, usecols = [4], skiprows = 1)

fig, ax = plt.subplots(figsize =(10,8) )
cm = plt.cm.get_cmap('RdYlBu')

# make some temorary arrays
#X = flux_contrast
X = T_contrast
Y = -(phi-2*np.pi)[::]*180/(np.pi)-180
Y2  = phiT*180/np.pi
Y3 =  (phiT1)[::]*180/(np.pi)
C = np.log10(eps+0.001)
C1 = np.log10(eps1+0.001)


X = np.hstack(X)
Y = np.hstack(Y)
C = np.hstack(C)

s = ax.scatter(X,Y,c=C, cmap = cm )



def f(x, a, b, c, d):
    return a*np.arctan(b-x*c)+d
    

popt, pcov = scipy.optimize.curve_fit(f, T_contrast,-(phi-2*np.pi)*180/(np.pi)-180)







'''
plt.plot(flux_contrast[::], (phiT)[::]*180/(np.pi), color='darkgreen', linestyle='dashed')
s2 = ax.scatter(X,Y2,c=C, cmap = cm, marker=u'D', s = 60, vmin = -1, vmax = 2, 
                label = 'Disc integrated flux' )

'''
plt.plot(Amp, phiT1*180/np.pi, linestyle='dashed', color = 'k',alpha = 0.6)
         
s3 = ax.scatter(Amp,Y3,c=C1, cmap = cm, marker=u'o', s = 60, vmin = -2, vmax = 2, 
         label = 'Parcel of gas' ) 
         

         
         
plt.plot(T_contrast, -(phi-2*np.pi)*180/(np.pi)-180, color='midnightblue', linestyle='dashed')  
s1 = ax.scatter(X,Y,c=C, cmap = cm, marker=u'*', s = 200, vmin = -2, vmax = 2 , 
          label = 'Disc integrated flux')
func, = plt.plot(T_contrast, f(T_contrast, popt[0],popt[1],popt[2],popt[3]), 'r',
                label = 'Analytic approximation', 
                color = 'darkred', linewidth = 2, alpha = 0.6)

#ax.plot(flux_contrast, 151.6*np.arctan(1.06-flux_contrast*0.79)-39.8 )

ax.set_xlim(0,1.05)
ax.set_ylim(-1,85)

plt.ylabel ('Phase offset ($^{\circ}$)', fontsize = 12, fontweight = 'bold')
plt.xlabel ('Normalized Temperature Contrast', fontsize = 12, fontweight = 'bold')

plt.legend([s1,func, s3],['Disc integrated ','Analytic approximation','Parcel of gas'],loc = 1, fontsize = 10)


cb = plt.colorbar(s1, ticks = np.log10(eps+0.01), format = '%2.1f' )

#cb.set_xticks([4.75, 5, 5.25, 5.5, 5.75])

labels = np.round(eps+0.001, decimals =2)
cb.ax.set_yticklabels(labels)
cb.set_label('Value of epsilon')

'lines'
#xx = zip(Amp,T_contrast)
#yy = zip(Y3, Y)

#for i in range(len(Amp)):
#    ax.plot(xx[i], yy[i], 'k-', linewidth = 2, alpha = 0.6)


fig.savefig('contrasts_and_phimax_temp.pdf')

print ('%2.2f'%popt[0]+'*ArcTan('+'%2.2f'%popt[1]+'-x*'+'%2.2f'%popt[2]+')+'+'%2.2f'%popt[3])

# In[10]:


# In[10]:

gas = p(name = "HD189733b",Teff = 4938, e=0,Porb = 2.21, a = 0.031, epsilon =20, wadv = 11,
                      argp = 90, Rstar = 0.781, Mstar = 0.846,
                      Rplanet = 1.138, pmax = 30 )
                      
t, d, Fwv = gas.Fobs()

#gas2 = p(name = "HD189733b",Teff = 4938, e=0.0000001,Porb = 2.21, a = 0.031,  
#                      epsilon=8,  argp = 90, Rstar = 0.781, Mstar = 0.846,
#                      Rplanet = 1.138, pmax = 6 )
                      
#t2, d2, Fwv2 = gas2.Fobs()

plt.plot(t,d[:,368+100,2], label = 'not eccentric')
#plt.plot(t2,Fwv2, label = 'eccentric')

#plt.plot(t/gas.Porb,d[:,368,1], label = 'not eccentric')
#plt.plot(t2/gas2.Porb,d2[:,368,1], label = 'eccentric')

plt.legend()
plt.show()


