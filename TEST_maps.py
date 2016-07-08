# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 12:18:33 2016

@author: diana
"""

# In[10]:

import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate
#from PyAstronomy import pyasl
import healpy as hp
from Parcel_healpy4 import parcel as p
#from timer import Timer



hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, wadv = 1.0/2, 
                  tau_rad = 20 , argp = 121.71, Rstar = 1.5, Mstar = 1.275)
                  

# In[11]:
t, d, Fmapwvobs, weight, weightpix, Fwv = hd.Fobs(MAP = True)
t, d, Fmap_wv, Fmap_wvpix, Fleavingwv, Ftotal = hd.Fleaving(MAP = True)
t, d, Fweight = hd.illum()


# In[11]:
d, Fweightpix = hd.shuffle(d = d, quantity = Fweight)

d, dshuffle = hd.shuffle()

"""SHUFFLE ISNT WORKING AS EXPECTED SOME OF THE TIME!!! CHECK IT"""
# In[11]:


"""CUTE FIGURE"""

window = np.zeros(weight.shape)
tstart = 15300
t,r = hd.radius(pmax = 3, steps = 300)
F0 =  (hd.sigmaB * hd.Teff**4*np.pi*hd.Rplanet**2*(hd.Rstar/min(r))**2) 
MAX = 0.41
MIN = 0.25

fig = plt.figure(figsize  =(6,4))
hp.visufunc.mollview(dshuffle[tstart-400,:], title = "Planet Temperature at a Moment in Time", xsize = 600, 
                         hold = False , sub =111, min = MIN, max = MAX,
                         unit = "($T/T0$)")
fig.savefig("sometemp")      

# In[ ]:
"""Figure that shows all the maps"""

MAX = 2*1*10**22
MIN = 1*10**16

fig= plt.figure(figsize = (10,16))
tstart = 15000
window[tstart,np.where((np.array(weightpix[tstart,:]) > 0.0 ) & (np.array(Fweightpix[tstart,:])> 0.0))] = 1
    
hp.visufunc.mollview(Fmap_wvpix[tstart,:], title = "Flux t="+str(tstart), xsize = 600, 
                         sub = (8,5,1),hold = False, min = MAX/4.0, max = MAX*3)#, unit = 'W/m')
                         
hp.visufunc.mollview(Fmapwvobs[tstart,:], title = "Flux-obs t="+str(tstart), xsize = 600, 
                         sub = (8,5,2),hold = False, min = MAX/4.0, max = MAX*3)#, unit = 'Flux')
    
hp.mollview(weightpix[tstart,:],title = "Visibility t="+str(tstart),min = 0, max=1,
                sub = (8,5,3))#, unit = 'Window ')  
                
hp.visufunc.mollview(Fweightpix[tstart,:], title = "Illum t="+str(tstart), xsize = 600, 
                         sub = (8,5,4),hold = False, min = 0, max = 1)#, unit = 'Window')
                         
hp.visufunc.mollview(window[tstart, :], title = "Intersection t="+str(tstart), xsize = 600, 
                         sub = (8,5,5),hold = False, min = 0, max = 1.5)
for i in range(1,8):
    
    time = 15000+i*50
    
    window[time,np.where((np.array(weightpix[time,:]) > 0.0 ) & (np.array(Fweightpix[time,:])> 0.0))] = 1
    
    hp.visufunc.mollview(Fmap_wvpix[time,:], title = "t="+str(time), xsize = 600, 
                         sub = (8,5,int(5*i+1)),hold = False, min = MAX/4.0, max = MAX*3)#, unit = 'Flux')
                         
    hp.visufunc.mollview(Fmapwvobs[time,:], title = "t="+str(time), xsize = 600, 
                         sub = (8,5,int(5*i+2)),hold = False, min = MAX/4.0, max = MAX*3)#, unit = 'Flux')
    
    hp.mollview(weightpix[time,:],title = "t="+str(time),min = 0, max=1,
                sub = (8,5,int(5*i+3)))#, unit = 'Window ')  
                
    hp.visufunc.mollview(Fweightpix[time,:], title = "t="+str(time), xsize = 600, 
                         sub = (8,5,int(5*i+4)),hold = False, min = 0, max = 1)#, unit = 'Window')
                         
    hp.visufunc.mollview(window[time, :], title = "t="+str(time), xsize = 600, 
                         sub = (8,5,int(5*i+5)),hold = False, min = 0, max = 1.5)#, unit = 'Intersection')
                         
fig.savefig('Fobs and window intersection - short time')
#sub = (int('42'+str((2*(i+1)))))





# In[8]:

""" THIS IS AN UNFINISHED FIGURE THAT SHOULD LOOK AT THE CONTOURS """

MAX = 2*10**(22)
MIN = 10**(14)
#SphericalProjAxes
fig= plt.figure(figsize = (10,8))

for i in range(4):
    """Define a mollweide Axes to handle mollweide projection.
    Input:
      - rot=, coord= : define rotation and coordinate system. See rotator.
      - coordprec= : number of digit after floating point for coordinates display.
      - format= : format string for value display.
      
      Other keywords from Axes (see Axes).
    """    
    time = 21.2*3*300+20*i
    #ax = 'ax'+str(i)
    #ax=hp.projaxes.MollweideAxes(fig,coord='G',rot=(0,180,0), coordprec = 4, format = 'format2')
    #                           format=format2,flipconv=flip)
    #fig.add_axes(ax)
    #name = hp.projaxes.MollweideAxes()
    #fig.add_subplot(int('33'+str(i+1))) 
    hp.mollview(Fmap_wvobspix[time,:],title = "t="+str(time),min = MIN, max=MAX,
                sub = (4,4,int(4*i+1)), unit = 'Flux' )
    hp.mollview(weight[time,:],title = "t="+str(time),min = 0, max=1,
                sub = (4,4,int(4*i+2)), unit = 'Window ')  
    hp.visufunc.mollview(Fmap_wvpix[time,:], title = "t="+str(time), xsize = 600, 
                         sub = (4,4,int(4*i+3)),hold = False, min = 4.5*MIN*10**8, max = MAX*1.05, unit = 'Flux')
    hp.visufunc.mollview(Fweight[time,:], title = "t="+str(time), xsize = 600, 
                         sub = (4,4,int(4*i+4)),hold = False, min = 0, max = 1, unit = 'Illum')
    #hp.visufunc.projplot(c_thetas1v[time],c_phis1v[time],'m*',markersize = 12)
    #hp.visufunc.projplot(c_thetas7v[time],c_phis7v[time],'ko',markersize = 6)
    
    #hp.visufunc.projplot(c_thetas1i[time],c_phis1i[time],'m*',markersize = 12)
    #hp.visufunc.projplot(c_thetas7i[time],c_phis7i[time],'go',markersize = 6)
    #hp.mollview(Fweightpix[t,:],title = "t="+str(t), 

fig.savefig('Observed flux compared with window weight and emmitted flux')

#sub = (int('42'+str((2*(i+1)))))
# In[ ]:
""" UNFINISHED CLACULATION OF ILLUMINATED FRACTION... SHOULD BE SOME MOD ARITHMETIC TRICK"""

phimin = np.empty(len(t))                      
for i in range(len(t)):
    window[i,np.where((np.array(weightpix[i,:]) > 0.0 ) & (np.array(Fweightpix[i,:])> 0.0))] = 1

phiwindow = phis*window
thetawindow = thetas*window
phimiddle = np.zeros(window.shape)
phimax = np.zeros(len(t))
phimin = np.zeros(len(t))
for i in range (len(t)):
    #phimiddle[i,366:430] = phiwindow[i,np.where((np.abs(thetawindow[i,:]-np.pi/2)<0.1))]
    a = np.max(np.concatenate((np.where(np.abs(thetawindow[i,:]-np.pi/2)<0.1),np.array(0).reshape(1,1)), axis = 1),axis=1)
    b = np.min(np.concatenate((np.where(np.abs(thetawindow[i,:]-np.pi/2)<0.1),np.array(0).reshape(1,1)), axis = 1),axis=1)
    phimax[i] = phis[i,a]
    phimin[i] = phis[i,b]


    
fraction = ((phimax-phimin))/(2*np.pi)

plt.plot(t[::], fraction[::])                           
                         
                         
#sub = (int('42'+str((2*(i+1)))))
# In[ ]: