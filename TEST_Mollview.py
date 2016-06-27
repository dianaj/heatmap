
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate

#from PyAstronomy import pyasl
import healpy as hp
from Parcel_healpy2 import parcel as p


# In[ ]:

gas = p(e=0.3,argp= 90-121.0)


# In[ ]:
fig, ax = plt.subplots(3,  figsize = (18,24))
#t, ang_vel = gas.ang_vel(pmax =23, steps = 100)
#t, weight = gas.weight(pmax =100)
t, alpha, f = gas.f_ratio(pmax =100)
t, zt, SSP, SSO = gas.SSP(pmax = 100)
ax[0].set_title("~4 orbital periods, 100 rotations, e=0.3")
ax[0].plot(t, SSP, label = 'SubStelar angle')
ax[0].plot(t,SSO, label = 'SubObs angle')
ax[0].set_xlabel('Time(seconds)', fontsize = 16)
ax[0].set_ylabel('Angle ( rads)', fontsize = 16)
ax[0].plot(t,(-alpha+180)/57.29, label = "Alpha the phase angle")
ax[0].legend()

ax[1].set_title("ZOOM")
ax[1].plot(t, SSP, label = 'SubStelar angle')
ax[1].plot(t,SSO, label = 'SubObs angle')
ax[1].set_xlabel('Time(seconds)', fontsize = 16)
ax[1].set_ylabel('Angle ( rads)', fontsize = 16)
#ax[1].plot(t,(-alpha+180)/57.29, label = "Alpha the phase angle")
ax[1].set_xlim(0,4000000)
#ax[1].set_ylim(-6,0)

ax[2].set_title("Difference")
ax[2].plot(t, SSO-SSP, label = 'SubStelar - Sub obs angle')
#ax[2].plot(t,SSO, label = 'SubObs angle')
ax[2].set_xlabel('Time(seconds)', fontsize = 16)
ax[2].set_ylabel('Angle ( rads)', fontsize = 16)
ax[2].plot(t,(-alpha+180)/57.29, label = "Alpha the phase angle")
ax[2].set_xlim(0,4000000)
#plt.plot(t,weight[:,[200,300]])
#plt.xlim(0,2000000)
#plt.ylim(-1,1)
#plt.plot(t, ang_vel)
fig.tight_layout()
fig.savefig("SSP and SSO")
#print gas.Porb/gas.P
#print np.sum(ang_vel[:(len(ang_vel)):])*gas.P/100
#plt.plot(t, zt)




# In[ ]:

hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, P = 0.3, wadv = 2, 
                  epsilon = 3.2*11, argp = 90-121.71, Rstar = 1.5, Mstar = 1.275)
t, zt, SSP, SOP = hd.SSP(pmax = 200, steps = 300)

(days, d , F, Fwv, Fpix, Fwvpix, Fmap, Fmap_pix, Fmap_wv, 
 Fmap_wvpix, Fleaving, Fleavingwv)  = hd.Fleaving(pmax = 20, TEST = True)

# In[42]:


fig = plt.figure(figsize = (10,8)) 
MAX = 10**(23)
MIN = 10**10

for i in range(9):
    t = 4000+10*i

    hp.visufunc.mollview(Fmap_wvpix[t,:],title = "Map at t="+str(t), min = MIN, max=MAX
    , sub = int('33'+str(i+1)))
    
    hp.visufunc.projplot(d[t,np.where(np.abs(d[t,:,2]-SSP[t])<0.5),1 ],
                     d[t,np.where(np.abs(d[t,:,2]-SSP[t])<0.5),2], flip='astro',
                     'k*',markersize = 6, )

    hp.visufunc.projplot(d[t,np.where(np.abs(d[t,:,2]-(SOP[t]))<0.2),1 ],
                     d[t,np.where(np.abs(d[t,:,2]-(SOP[t]))<0.2),2], flip='astro',
                     'r*',markersize = 6)


#hp.mollview(Fmap_wvpix[1000,:],title = "Map at t=1000", min = MIN, max=MAX)

#hp.visufunc.projplot(d[1000,np.where(np.abs(d[1000,:,2]-SSP[1000])<0.1),1 ],
#                     d[1000,np.where(np.abs(d[1000,:,2]-SSP[1000])<0.1),2], 
#                     'k*',markersize = 16)
'THIS ISNT PLOTTING!!! BECAUSE MY ANGLES KEEP GETTING BIGGER BUT THE SOP GOES TO @ P OR SO'
#hp.visufunc.projplot(d[1000,np.where(np.abs(d[1000,:,2]-(SOP[1000]))<0.1),1 ],
#                     d[1000,np.where(np.abs(d[1000,:,2]-(SOP[1000]))<0.1),2], 
#                     'r*',markersize = 16)


# In[12]:

t, zt, SSP, SOP,  thetas, phis, c_thetas1, c_phis1, c_thetas7, c_phis7,F, Fnorm, Fweight, Fweightpix = gas.illum(pmax=10, TEST = True)
#print c_thetas
#t, d, thetas, phis, c_thetas1, c_phis1, c_thetas7, c_phis7, weight, weightpix = gas.visibility(pmax=10)
#dayss, d = gas.DE(pmax=10)


tv, thetasv, phisv, c_thetas1v, c_phis1v, c_thetas7v, c_phis7v, coordsSSP, coordsSOP, weightv, weightpixv = gas.visibility(pmax = 10, TEST = True)

# In[12]:
#plt.plot(c_thetas[400,:], c_phis[400,:], '+')
#plt.plot(t,Fnorm)
#c_t = thetas[:,np.where(np.abs(Fnorm[226,:]-0.7)<0.1)]
#c_f = phis[:,np.where(np.abs(Fnorm[226,:]-0.7)<0.1)]

#plt.plot(c_t[500,:],c_f[500,:], '+', markersize = 16)
plt.plot( c_phis1[400],c_thetas1[400],'ok', markersize = 12 )
plt.plot(c_phis7[400],c_thetas7[400],'ob', markersize = 12)
plt.plot(c_phis1[500],c_thetas1[500],'or', markersize = 12)
plt.plot(c_phis7[500],c_thetas7[500],'ok', markersize = 12)




MIN = 0
MAX = 1
hp.mollview(Fweightpix[0,:], min = MIN, max =MAX)
hp.visufunc.projplot(c_thetas1[0],-c_phis1[0], 'k*',markersize =12)
hp.visufunc.projplot(c_thetas7[0],-c_phis7[0], 'k*', markersize =12)

hp.visufunc.projplot(c_thetas1v[0],c_phis1v[0], 'b*',markersize =12)
hp.visufunc.projplot(c_thetas7v[0],c_phis7v[0], 'k*', markersize =12)


hp.mollview(Fweightpix[20,:], min = MIN, max =MAX)
#hp.mollview(d[20,d[20,:,0].astype(int),2])
hp.visufunc.projplot(c_thetas1[20],-c_phis1[20], coords = 'E', 'k*', markersize =12)
hp.visufunc.projplot(c_thetas7[20],-c_phis7[20], 'm*', markersize =12)
#weightx = np.empty(Fweight.shape)
#weightx[30] = Fweight[30,d[30,:,0].astype(int)]
hp.visufunc.projplot(c_thetas1v[20],c_phis1v[20], 'b*',markersize =12)
hp.visufunc.projplot(c_thetas7v[20],c_phis7v[20], 'k*', markersize =12)

hp.visufunc.mollview(Fweightpix[30], min = MIN, max =MAX)
hp.visufunc.projplot(c_thetas1[30],-c_phis1[30], 'k*', markersize =12)
hp.visufunc.projplot(c_thetas7[30],-c_phis7[30], 'm*', markersize =12)
#hp.visufunc.projplot(thetas[30,np.where(Fweight[30,:]==np.max(Fweight[30,:]))],
#                     -phis[30,np.where(Fweight[30,:]==np.max(Fweight[30,:]))], 'k*', markersize =12)
hp.visufunc.projplot(c_thetas1v[30],c_phis1v[30], 'b*',markersize =12)
hp.visufunc.projplot(c_thetas7v[30],c_phis7v[30], 'k*', markersize =12)

hp.visufunc.mollview(Fweightpix[50], min = MIN, max =MAX)
hp.visufunc.projplot(c_thetas1[50],-c_phis1[50], 'k*', markersize =12)
hp.visufunc.projplot(c_thetas7[50],-c_phis7[50], 'm*', markersize =12)

hp.visufunc.projplot(c_thetas1v[50],c_phis1v[50], 'b*',markersize =12)
hp.visufunc.projplot(c_thetas7v[50],c_phis7v[50], 'k*', markersize =12)



# In[12]:
from pylab import arange, show, cm
cmap1 = cm.Spectral

cmap1.set_under('w')
cmap2 =cm.Accent

cmap2.set_under('w')

MIN =0
MAX =1
fig = plt.subplots(2, figsize = (12,8))
hp.visufunc.mollview(Fweightpix[500], title = 'Illlumination t=500', min = MIN, max =MAX, cmap=cmap1, sub = 211)
hp.visufunc.projplot(c_thetas1[500],c_phis1[500], 'k*', markersize =12)
hp.visufunc.projplot(c_thetas7[500],c_phis7[500], 'm*', markersize =12)


hp.visufunc.mollview(weightpixv[500], title = 'Visibility t= 500', min = MIN, max =MAX, cmap=cmap2,  sub =212)
hp.visufunc.projplot(np.pi/2,np.pi, 'b*',markersize =24)
#hp.visufunc.projplot(c_thetas1v[500],c_phis1v[500], 'b*',markersize =12)
#hp.visufunc.projplot(c_thetas7v[500],c_phis7v[500], 'k*', markersize =12)
#hp.visufunc.projplot(thetas[50,np.where(Fweight[50,:]==np.max(Fweight[50,:]))],
#                     -phis[50,np.where(Fweight[50,:]==np.max(Fweight[50,:]))], 'k*', markersize =12)
'''
fig = plt.subplots(2, figsize = (12,8))
hp.visufunc.mollview(Fweightpix[550], title = 'Illlumination t=550',min = MIN, max =MAX, cmap=cmap1,sub = 211)
hp.visufunc.projplot(c_thetas1[550],c_phis1[550], 'k*', markersize =12)
hp.visufunc.projplot(c_thetas7[550],c_phis7[550], 'm*', markersize =12)


hp.visufunc.mollview(weightpixv[500], title = 'Visibility t= 550', min = MIN, max =MAX, cmap=cmap2, sub =212)
hp.visufunc.projplot(c_thetas1v[500],c_phis1v[500], 'b*',markersize =12)
hp.visufunc.projplot(c_thetas7v[500],c_phis7v[500], 'k*', markersize =12)

fig = plt.subplots(2, figsize = (12,8))
hp.visufunc.mollview(Fweightpix[800], title = 'Illlumination t=800', min = MIN, max =MAX, cmap=cmap1, sub = 211)
hp.visufunc.projplot(c_thetas1[800],c_phis1[800], 'k*', markersize =12)
hp.visufunc.projplot(c_thetas7[800],c_phis7[800], 'm*', markersize =12)


hp.visufunc.mollview(weightpixv[800], title = 'Visibility t= 800', min = MIN, max =MAX, cmap=cmap2,  sub =212)
hp.visufunc.projplot(c_thetas1v[800],c_phis1v[800], 'b*',markersize =12)
hp.visufunc.projplot(c_thetas7v[800],c_phis7v[800], 'k*', markersize =12)'''

