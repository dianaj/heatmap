
# coding: utf-8

# In[1]:


# In[2]:

import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate
from PyAstronomy import pyasl
import healpy as hp
from Parcel_healpy3 import parcel


# In[4]:

HD17156b = parcel(Teff = 6079, e=0.6768, Porb = 21.2, a = 0.1589, P = 0.3, wadv = 2, 
                  epsilon = 3.2*11, argp = 90-121.71, Rstar = 1.5, Mstar = 1.275)


# In[5]:

"""planet integrated properties: diagnostics"""

import time as tme

tic = tme.time()



nps = 240

fig, ax = plt.subplots(6,1, sharex = True, figsize = (15,24))

#ax0 = plt.subplot2grid((4,2), (0,0), colspan=2)
#ax1 = plt.subplot2grid((4,2), (1,0), colspan=2)
#ax2 = plt.subplot2grid((4,2), (2,0), colspan=2)
#ax4 = plt.subplot2grid((4,2), (3,0), colspan=2)
'plot 1- orbital separation'
t,r = HD17156b.radius(pmax = nps)
ax[0].plot(t,r/HD17156b.AU)
ax[0].set_ylabel("r - orbital separation (AU)")


'plot 2 - orbital angular velocity - what is this good for??'

t1,worb = HD17156b.ang_vel(pmax = nps)
ax[1].plot(t, worb/worb[0] )
ax[1].set_ylabel("Orbital angular veocity : $\omega_{orb}/\omega_{periastron}$")


'plot 3 - Eq Temprature at the substellar point at each spot in the orbit '

T0 = HD17156b.T0*((HD17156b.a*(1-HD17156b.e))/r)#**0.5
temp, = ax[2].plot (t,T0)
ax[2].set_ylabel("T (K) - temp of substellar point$")

#m=40
#AllT = np.empty((m,300*nps))

#for i in range(m):
days, d = HD17156b.DE(pmax = nps, NSIDE = 8)
AllT = d[:,:,3]
    #AllT[i]= np.array(T)
    #phi= np.linspace(0,2*30.0, num=300*30)
    #temp, = ax2.plot((phi*np.pi)/(gas.wadv*gas.P),T[0,i,:])
    #ax[2].plot((days*gas.P)[600*16::],(np.array(T)*gas.T0)[600*16::], color = 'r')
Tmax = np.max(AllT, axis = 1)
maxtemp, = ax[2].plot(t,(Tmax*HD17156b.T0)[::], color = 'g', ls ='--')
ax[2].legend([temp, maxtemp], ["Instantateous Temperature","Gas max temp"], loc = 1)

#phi, days, T = gas.DE(pmax = 12, steps = 300)



#highest temperature on the model planet 

#'''----- this is harder. will need temperatures on each spot of the planet (as a function of phi) 
#at each moment in time----'''
 
    
'plot 4 - Total absorbed flux / total emitted flux '

#total absorbed flux

F = HD17156b.Finc_hemi(pmax = nps)
F0 =  (HD17156b.sigmaB * HD17156b.Teff**4*np.pi*HD17156b.Rplanet**2*(HD17156b.Rstar/min(r))**2) 
Fnorm = F/F0
ax[3].plot(t, Fnorm)
ax[3].set_ylabel("$F/F_{max}$ - Total absorbed flux")

# In[1]:

days, d ,Fwv, Fwvpix, Fmap_wv, Fmap_wvpix, Fleavingwv = HD17156b.Fleaving(pmax = nps, 
                                                                                      steps = 300.0, 
                                                                                      NSIDE=4, 
                                                                                      wavelenght= 8, TEST =True)



FnormDE = Fleavingwv/F0
ax[3].plot(t, FnormDE)

# In[1]:


#total emmitted flux ???

#''' --- integral  of sigmaB * T(at each spot on planet)**4 dtheta dphi at each point in time '''


'Plot 5 - planets illuminated fraction'

ax[4].set_ylabel("f - planet's illuminated fraction")
t,alpha, f = HD17156b.f_ratio(pmax = nps)
ax[4].plot(t, f)



'plot 6'

days, d, Fmapwvobs, Fmap_wvobspix, weight, weightpix, Fwv = HD17156b.Fobs(pmax = nps, 
                                                                             steps = 300, 
                                                                             NSIDE = 4)
#days, tt, Fmap, Fmapwv, Fobs, Fobswv = hd.Fobs(pmax = nps, steps = 300, Nthetas = 20, Nphis = 40, wavelenght = 8)


ax[5].set_xlabel("Time (s)")
ax[5].set_ylabel("Fplanet / Fstar")
#fobs, = ax[5].plot(tt, Fobs/ HD17156b.Fstar()[0], color = 'b')
#ftotal, = ax[5].plot(t, FDE/ HD17156b.Fstar()[0], color = 'r')
fobswv, = ax[5].plot(t, Fwv/ HD17156b.Fstar()[1], color = 'g')
ax[5].legend([fobswv], [ "F-obs at 8 um"], loc = 1)
#ax[5].set_ylim(0,0.0001)
plt.subplots_adjust(hspace = 0.1)

fig.savefig('HD_diag')
plt.show()
toc = tme.time()

print (toc-tic)


# In[ ]:



