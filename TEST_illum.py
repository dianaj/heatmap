
# coding: utf-8

# In[75]:


# In[76]:

import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate
from PyAstronomy import pyasl
import healpy as hp
from Parcel_healpy import parcel as p


# In[83]:

hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, P = 0.3, wadv = 2, 
                  epsilon = 3.2*11, argp = 90-121.71, Rstar = 1.5, Mstar = 1.275)

(t, thetas, phis, c_thetas1, c_phis1, c_thetas7, c_phis7,
 F, Fnorm, Fweight, Fweightpix) = hd.illum(pmax =100, TEST = True)


# In[84]:

#t, zt, SSP, SOP = hd.SSP(pmax = 50, steps = 300)
#days, d = hd.DE(pmax=50)
MAX = 1
MIN = 0
#SphericalProjAxes
#fig= plt.figure(figsize = (10,8))
fig = plt.figure(figsize = (12,12))
for i in range(8):
    """Define a mollweide Axes to handle mollweide projection.
    Input:
      - rot=, coord= : define rotation and coordinate system. See rotator.
      - coordprec= : number of digit after floating point for coordinates display.
      - format= : format string for value display.
      
      Other keywords from Axes (see Axes).
    """    
    t = 20430-i*100
    #ax = 'ax'+str(i)
    #ax=hp.projaxes.MollweideAxes(fig,coord='G',rot=(0,180,0), coordprec = 4, format = 'format2')
    #                           format=format2,flipconv=flip)
    #fig.add_axes(ax)
    #name = hp.projaxes.MollweideAxes()
    ax = fig.add_subplot(int('42'+str(i+1))) 
    hp.mollview(Fweightpix[t,:],title = "t="+str(t),min = MIN, max=MAX, sub = int('42'+str(i+1)))
                
    #ax.plot(hp.visufunc.projplot(c_thetas1[t],c_phis1[t],'m*',markersize = 12))
    #ax.plot(hp.visufunc.projplot(c_thetas7[t],c_phis7[t],'ko',markersize = 6))
    #hp.mollview(Fweightpix[t,:],title = "t="+str(t), 
                #min = MIN, max=MAX,sub = (int('33'+str(i+1))))
    #plt.subplot(int('33'+str(i+1))) 
    #hp.visufunc.projplot(c_thetas1[t],c_phis1[t],'m*',markersize = 12, sub = (int('33'+str(i+1))))
                     
    #hp.visufunc.projplot(d[t,np.where(np.abs(d[t,:,2]-SSP[t])<0.5),1 ],
    #                 d[t,np.where(np.abs(d[t,:,2]-SSP[t])<0.5),2], 
    #                 'k*',markersize = 16)

    #hp.visufunc.projplot(c_thetas7[t],c_phis7[t],'ko',markersize = 6)

#for label in labels:
'''
for ax in fig.get_axes():
    label = ax.get_label()
    print label
    #i=(int(label[2])-1)
    #t=i*15+8000
    hp.visufunc.projplot(c_thetas1[t],c_phis1[t],'m*',markersize = 12)
    hp.visufunc.projplot(c_thetas7[t],c_phis7[t],'ko',markersize = 6)
    
ax.set_title('Illumination map')
fig.savefig("TEST-illumination")'''


# In[ ]:



