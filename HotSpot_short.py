# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 19:44:56 2016

@author: diana
"""



import numpy as np
import matplotlib.pyplot as plt
#import healpy as hp
from Parcel_healpy6 import parcel as p


# In[10]:    
gas = p(name = "HD189733b",Teff = 5938, e=0,Porb = 2.21, a = 0.031, wadv = 11.0, 
                          epsilon = 6, argp = 90, Rstar = 0.781, Mstar = 0.846,
                          Rplanet = 0.838, pmax = 10, steps = 300, NSIDE = 8 )      
    
#t, d, Fwv = gas.Fobs(BOLO = True)

#plt.plot(t,Fwv)

epsi = 100.0/(10**(np.linspace(0,3, num = 20)))
eps = np.concatenate((epsi,[0]), axis = 0)

 # In[10]:

"""START HERE IF NOT MAKING NEW DATA"""

File ='hotspot(try2)_NSIDE16BOLO_e=0.out' 
File2 ='hotspot(TEMP).out' 
 

'''LOAD FLUX STUFF'''
eps = np.loadtxt(File, usecols = [0], skiprows = 1)

max_flux = np.loadtxt(File, usecols = [1], skiprows = 1)
min_flux = np.loadtxt(File, usecols = [2], skiprows = 1)
#maxf = max_flux
#minf = min_flux

flux_contrast = (max_flux - min_flux)/max_flux
area = np.pi*gas.Rplanet**2
T_contrast = ((max_flux/(gas.sigmaB*area))**0.25-(min_flux/(gas.sigmaB*area))**0.25)/gas.T0

phi = np.loadtxt(File, usecols = [3], skiprows = 1)

phiT = np.loadtxt(File, usecols = [4], skiprows = 1)

'''LOAD TEMP STUFF'''
eps1 = np.loadtxt(File2, usecols = [0], skiprows = 1)
Tmax = np.loadtxt(File2, usecols = [1], skiprows = 1)
Tmin = np.loadtxt(File2, usecols = [2], skiprows = 1)
phiT1 = np.loadtxt(File2, usecols = [3], skiprows = 1)

Amp = (Tmax - Tmin)


 # In[10]:

'''IMPORTANT FIGURE'''

#max_flux= maxf
#min_flux = minf

flux_contrast = (max_flux - min_flux)/max_flux
area = np.pi*gas.Rplanet**2
T_contrast = ((max_flux/(gas.sigmaB*area))**0.25-(min_flux/(gas.sigmaB*area))**0.25)/gas.T0

import scipy 


fig, ax = plt.subplots(figsize =(10,8) )
cm = plt.cm.get_cmap('RdYlBu')

X = T_contrast
Y = -(phi-2*np.pi)[::]*180/(np.pi)-180
Y2  = phiT*180/np.pi
Y3 =  (phiT1)[::]*180/(np.pi)

#i don't want my log color scale to start at 0, it messes things up
C = np.log10(eps+0.001)
C1 = np.log10(eps+0.001)

np.savetxt('hotspot(Tcontrast).out', 
                np.array(zip(X,Y)).reshape(-1,2), 
                header = "t-contrast, phase_angle")
X = np.hstack(X)
Y = np.hstack(Y)
C = np.hstack(C)

s = ax.scatter(X,Y,c=C, cmap = cm )



def f(x, a, b, c, d):
    return a*x+b*x**2 + c*x**3+d
    

popt, pcov = scipy.optimize.curve_fit(f, T_contrast,-(phi-2*np.pi)*180/(np.pi)-180)



'lines'
xx = zip(Amp,T_contrast)
yy1 = zip(Y3, Y3)

xx1 = zip(T_contrast,T_contrast)
yy = zip(Y3,Y)

for i in range(len(Amp)):
    ax.plot(xx[i], yy[i], 'k-', linewidth = 2, alpha = 0.6)




#plt.plot(flux_contrast[::], (phiT)[::]*180/(np.pi), color='darkgreen', linestyle='dashed')
#s2 = ax.scatter(X,Y2,c=C, cmap = cm, marker=u'D', s = 60, vmin = -1, vmax = 2, 
#                label = 'Disc integrated flux' )


plt.plot(Amp, phiT1*180/np.pi, linestyle='dashed', color = 'k',alpha = 0.4)
         
s3 = ax.scatter(Amp,Y3,c=C1, cmap = cm, marker=u'o', s = 60, vmin = -2, vmax = 2, 
         label = 'Parcel of gas' ) 
         

         
         
plt.plot(T_contrast, -(phi-2*np.pi)*180/(np.pi)-180, color='midnightblue', linestyle='None')  
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


#ax.plot(xx1[i], yy[i], 'k-', linewidth = 2, alpha = 0.6)


fig.savefig('contrasts_and_phimax_temp_linesTRY.pdf')

print ('%2.2f'%popt[0]+'x+'+'%2.2f'%popt[1]+'x**2+'+'%2.2f'%popt[2]+'x**3'+'+'+'%2.2f'%popt[3])




# In[10]: 

"""
START HERE IF YOU'RE GONNA RUN THE SIMULATION

THIS LOOP FINDS THE MAXIMUM FLUX AND LOCATION. 
RUN ONLY IF YOU GOT TIME .
"""


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
        #Fwv actually means total flux now, not wavelenght dependant.
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
    
    """This little piece finds the theoretical (not disk integrated maximum)"""
   
    phiT = []
    for i in range(len(eps)):    
        phiT.append(gas.phi_max(eps[i]))
    
    phiT = np.array(phiT).flatten()
    print phi
    


# In[10]:

'''TRYING TO MAKE IT WORK WITH SMALLER VALUES OF EPS. IT DOESN'T WORK IT OVERFLOWS'''

eps = np.concatenate((eps,[0.01,0.001]), axis = 0)
maxf = np.concatenate((maxf,[0.01,0.001]), axis = 0)
minf = np.concatenate((minf,[0.01,0.001]), axis = 0)
phi = np.concatenate((phi,[0.01,0.001]), axis = 0)

eps[21] = 0.05
eps[22] = 0.01
for i in [21,22]:
    gas = p(name = "HD189733b",Teff = 4938, e=0,Porb = 2.21, a = 0.031, wadv = 11.0, 
                              epsilon = eps[i], argp = 90, Rstar = 0.781, Mstar = 0.846,
                              Rplanet = 1.138, pmax = 10, steps = 500, NSIDE = 16 )
    t, d, Fwv = gas.Fobs(BOLO = True)
    
    
    backend = int((gas.pmaxi-1.5)*gas.stepsi)
    maxf[i]  = (max(Fwv[backend::]))
    minf[i] = (min(Fwv[backend::]))
    phi[i] = (gas.t[np.argmax(Fwv[backend::])+backend]/gas.P)%1*(2*np.pi)



  # In[10]:

"""IF IT WORKED,  SAVE IT FOR LATER"""
 
np.savetxt('hotspot(try2)_NSIDE'+str(gas.NSIDE)+'BOLO'+'_e='+str(gas.e)+'.out', 
                np.array(zip(eps,maxf,minf,phi,phiT)).reshape(-1,5), 
                header = "eps,max_flux, min_flux, phi_max, phi_theory")

                


# In[10]:
'''ANALYTIC TEMP STUFF - doesn't take time '''

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

"""SAVE THE ANALYTIC STUFF TOO """

np.savetxt('hotspot(TEMP).out', 
                np.array(zip(eps1,Tmax, Tmin, phiT1)).reshape(-1,4), 
                header = "eps,Tmax, Tmin, phi_theory")


 # In[10]: