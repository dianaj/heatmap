
# coding: utf-8

# In[9]:

get_ipython().magic(u'matplotlib inline')


# In[10]:

import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate
from PyAstronomy import pyasl
import healpy as hp
from Parcel_healpy import parcel as p


# In[11]:

hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, P = 0.3, wadv = 2, 
                  epsilon = 3.2*11, argp = 90-121.71, Rstar = 1.5, Mstar = 1.275)


# In[12]:

days, d, Fmapwvobs, Fmap_wvpix, weight, weightpix, Fwv = hd.Fobs(pmax=200)


# In[5]:

Fstar, Fstarwv = hd.Fstar()

fig = plt.figure(figsize=(16,4))
plt.plot(tt, Fwv/Fstarwv)


# In[6]:

#(t, d2, thetas, phis, c_thetas1, c_phis1, 
# c_thetas7, c_phis7, weight, weightpix)= hd.visibility(pmax = 200, TEST = True)


# In[7]:

( t, thetas, phis, c_thetas1, c_phis1, c_thetas7, c_phis7,F, Fnorm, 
 Fweight, Fweightpix) = hd.illum(pmax = 200, TEST = True)


# In[8]:

MAX = 5*10**(22)
MIN = 10**(16)
#SphericalProjAxes
#fig= plt.figure(figsize = (10,8))

for i in range(9):
    """Define a mollweide Axes to handle mollweide projection.
    Input:
      - rot=, coord= : define rotation and coordinate system. See rotator.
      - coordprec= : number of digit after floating point for coordinates display.
      - format= : format string for value display.
      
      Other keywords from Axes (see Axes).
    """    
    time = 8000+15*i
    #ax = 'ax'+str(i)
    #ax=hp.projaxes.MollweideAxes(fig,coord='G',rot=(0,180,0), coordprec = 4, format = 'format2')
    #                           format=format2,flipconv=flip)
    #fig.add_axes(ax)
    #name = hp.projaxes.MollweideAxes()
    hp.mollview(Fmap_wvpix[time,:],title = "t="+str(time),min = MIN, max=MAX)#, sub = (int('33'+str(i+1))),
                
    hp.visufunc.projplot(c_thetas1[time],c_phis1[time],'m*',markersize = 12)
    hp.visufunc.projplot(c_thetas7[time],c_phis7[time],'ko',markersize = 6)
    
    hp.visufunc.projplot(c_thetas1i[time],c_phis1i[time],'m*',markersize = 12)
    hp.visufunc.projplot(c_thetas7i[time],c_phis7i[time],'go',markersize = 6)
    #hp.mollview(Fweightpix[t,:],title = "t="+str(t), 


# In[ ]:



