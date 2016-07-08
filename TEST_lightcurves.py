# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 12:05:22 2016

@author: diana
"""
"""DRAWS A FEW LIGHTCURVES FOR DIFFERENT PARAETERS AND DIFFERENT PLANETS . 
SEEMS TO WORK WELL NOW"""


import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate
#from PyAstronomy import pyasl
import healpy as hp
from Parcel_healpy4 import parcel as p
from timer import Timer



#hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, wadv = 1.0/2, 
#                  tau_rad = 20 , argp = 121.71, Rstar = 1.5, Mstar = 1.275)
                  

# In[11]:

hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, wadv = 1.0/2, 
                  epsilon = 1, argp = 121.71, Rstar = 1.5, Mstar = 1.275)
Fstar, Fstarwv = hd.Fstar()
eps = [0.2, 3.2]
wadv = [ 1, 2]
argp = [61, 121]

fig, ax = plt.subplots(4, 2, sharex='col', sharey='row', figsize = (34, 24))
for i in eps:
    for j in wadv:
        for k in argp: 
            hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, wadv = j, 
                  epsilon = i, argp = k, Rstar = 1.5, Mstar = 1.275)
                  

            t, d, Fwv = hd.Fobs(NSIDE=4)
            



            t, d, Fmap_wv, Fmap_wvpix, Fleavingwv, Ftotal  = hd.Fleaving(NSIDE = 4, MAP =True)
            a = eps.index(i)
            b = wadv.index(j)
            "there's 2 ways of making 1!!! gotta define your subplots better. "
            c=  argp.index(k)
            
            
            ax[b+2*c, a].plot( t/hd.P, Fwv/Fstarwv, label = "Observed Flux")
            ax[b+2*c, a].plot( t/hd.P, Fleavingwv/Fstarwv, label = 'Emitted flux')
            #plt.plot( t/hd.P, Fleavingwvpix/Fstarwv)
            ax[b+2*c, a].set_ylabel('F_planet/ F_star')
            
            ax[b+2*c, a].legend(fontsize = 8)
            ax[b+2*c, a].set_title("tau_rad = "+"%.1f"%(hd.tau_rad/ hd.days * 24)+ 
                    " hours, wadv = "+"%.1f"%(hd.wadv/(2*np.pi/hd.P))+ " x wrot"+
                    " argp = " + "%.1f"%(hd.argp)+ "degrees")
ax[3, 0].set_xlabel('Time (s)')
ax[3, 1].set_xlabel('Time (s)')
plt.subplots_adjust(wspace = 0.2)
fig.savefig("Fluxtries.pdf")



# In[12]:

"Lets' try other planets!!!"

HD189733b = p(name = "HD189733b",Teff = 4938, e=0,Porb = 2.21, a = 0.031, wadv = 1.0/2, 
                  tau_rad = 20, argp = 0, Rstar = 0.781, Mstar = 0.846,
                  Rplanet = 1.138 )

Wasp_17b =p(name = "Wasp_17b",Teff = 6650, e=0.13,Porb = 3.74, a = 0.0515, wadv = 1.0/2, 
                  tau_rad = 20, argp = 290, Rstar = 1.38, Mstar = 1.2,
                  Rplanet = 1.991)

GI436b = p(name = "GI436b",Teff = 3684, e=0.15, Porb = 2.64, a = 0.02887, wadv = 1.0/2, 
                  tau_rad = 20, argp = 351, Rstar = 0.464, Mstar = 0.452, 
                  Rplanet = 0.38)

XO_3b = p(name = "XO-3b",Teff = 6781, e=0.26, Porb = 3.19, a = 0.0454, wadv = 1.0/2, 
                  tau_rad = 20, argp = 346, Rstar = 1.49, Mstar = 1.41, 
                  Rplanet = 1.217)

HAT_P_11b = p(name = "HAT_P_11b",Teff = 4780, e=0.198, Porb = 4.8878, a = 0.053, wadv = 1.0/2, 
                  tau_rad = 20, argp = 355.2, Rstar = 0.75, Mstar = 0.81, 
                  Rplanet = 0.422)

HAT_P_2b = p(name = "HAT_P_2b", Teff = 6414, e=0.5171, Porb =5.63, a = 0.0674, wadv = 1.0/2, 
                  tau_rad = 20, argp = 185.22, Rstar = 1.54, Mstar = 1.34, 
                  Rplanet = 0.951)



# In[12]:




planets = [HD189733b, Wasp_17b, GI436b, XO_3b,HAT_P_11b, HAT_P_2b]


#fig, ax = plt.subplots(2, 2, sharex='col', sharey='row', figsize = (16, 16))

for planet in planets : 

                      
    fig = plt.figure(figsize = (8,8))
    taurads = [1,20]
    wadv = [0.5, 1.5, 2]
    
    #for planet in planets:
    Fstar, Fstarwv = planet.Fstar()
    
    planet.tau_rad = 0.5 * 3600.0
    
    planet.wadv = 1*(2*np.pi/planet.P)
    
    t, d, Fwv = planet.Fobs(pmax=4, NSIDE=4)
                
    plt.plot(t/planet.Porb, Fwv/Fstarwv, 
                         label = ("$\tau_{rad} = $"+str(planet.tau_rad/ 3600.0) +
                         " hrs, $\omega_{rot} = $" + str(planet.wadv/((2*np.pi)/planet.P)) + "$\omega_{max}$") )
    plt.legend()
    plt.xlim(3, 4)
    #plt.ylim(0, 0.003)   
    for j in [20]:
            for k in wadv: 
                #XO_3b.tau_rad = j
                #XO_3b.wadv = k
                      
                planet.tau_rad = j * 3600.0
    
                planet.wadv = k*(2*np.pi/planet.P)
                t, d, Fwv = planet.Fobs(pmax=4, NSIDE=4)
                
                plt.plot(t/planet.Porb, Fwv/Fstarwv, 
                         label = ("$\tau_{rad} = $"+str(planet.tau_rad/ 3600.0) +
                         " hrs, $\omega_{rot} = $" + str(planet.wadv/((2*np.pi)/planet.P)) + "$\omega_{max}$") )
                plt.legend()
                plt.xlim(3, 4)
                #plt.ylim(0, 0.003)
    plt.title(str(planet.name)+": P = "+ str(planet.Porb/planet.days)+
                        ", e = " + str(planet.e)+", w = " + str(planet.argp) +" degrees")
    
    fig.savefig(str(planet.name)+".pdf")
            
