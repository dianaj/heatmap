
# coding: utf-8




# In[10]:

import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate
#from PyAstronomy import pyasl
import healpy as hp
from Parcel_healpy4 import parcel as p


# In[11]:

hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, wadv = 1, 
                  epsilon = 3.2*11, argp = 121.71, Rstar = 1.5, Mstar = 1.275)


days, d = hd.DE(pmax =50)
d[:,:,2][d[:,:,2] == 0] = 0.01
wv = 8*10**(-6)
Fwv_ = (np.log(2*hd.h*hd.c**2)-np.log(wv**5) + np.log(1) - 
        np.log(np.e**((hd.h*hd.c)/(wv*hd.K*np.array(d[:,:,2].copy())*hd.T0))-1))
        
Fwv = np.e**(Fwv_)  

dA = hp.nside2pixarea(8)*hd.Rplanet**2
  
Fmap_wv = (Fwv.copy()*dA)


# In[12]:
d, Fmap_wvpix = hd.shuffle(d, Fmap_wv, pmax = 50)
# In[12]:
'little argument at periastron test'

hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, wadv = 1, 
                  epsilon = 3.2*11, argp = 121.71, Rstar = 1.5, Mstar = 1.275)



t, alpha,  f = hd.f_ratio(pmax = 20)

plt.plot(t,f)
plt.plot(t,alpha*np.pi/180)
plt.plot(t,TA)
# In[12]:                

hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, wadv = 2, 
                  epsilon = 3.2, argp = 121.71, Rstar = 1.5, Mstar = 1.275)
                  

days, d, Fmapwvobs, weight, weightpix, Fwv = hd.Fobs(pmax=20, steps = 1000, NSIDE=4)


# In[5]:

#d, Fmap_wvobspix = hd.shuffle(d, Fwvobs, pmax = 100, NSIDE = 4)
#d, weightpix =hd.shuffle(d, weight, pmax= 100, NSIDE = 4)

# In[5]:
days, d, Fmap_wv, Fmap_wvpix, Fleavingwv, Fleavingwvpix  = hd.Fleaving(pmax = 20, steps = 1000,  NSIDE = 4)
 
 
# In[5]:
fig = plt.figure(figsize=(16,4))

Fstar, Fstarwv = hd.Fstar()

#

plt.plot( days* hd.days, Fwv/Fstarwv, label = "Observed Flux")
plt.plot( days* hd.days, Fleavingwv/Fstarwv, label = 'Emitted flux')
plt.plot( days* hd.days, Fleavingwvpix/Fstarwv)
plt.xlabel('F_planet/ F_star')
plt.ylabel('Time (s)')
plt.legend()
plt.title("tau_rad = "+"%.1f"%(hd.tau_rad/ hd.days * 24)+ " hours, wadv = "+"%.1f"%(hd.wadv/(2*np.pi/hd.P))+ " x wrot" )
fig.savefig("Flux_tau_rad = "
            +"%.1f"%(hd.tau_rad/ hd.days * 24)+ " hours, wadv = "+"%.1f"%(hd.wadv/(2*np.pi/hd.P))+ " x wrot.pdf")


# In[6]:

t, weight= hd.visibility(pmax = 100, d=d, NSIDE =4)
crap, weightpix = hd.shuffle(d, weight.copy(), pmax = 100, NSIDE = 4)


# In[7]:

d,Fweight= hd.illum(pmax = 20, steps = 1000, NSIDE =4)

d, Fweightpix = hd.shuffle(d, Fweight.copy(), pmax = 100, NSIDE = 4)

# In[ ]:

window = np.zeros(weight.shape)

MAX = 2*10**(23)
MIN = 10*10**(14)

fig= plt.figure(figsize = (10,16))
tstart = 15000
window[tstart,np.where((np.array(weightpix[tstart,:]) > 0.0 ) & (np.array(Fweightpix[tstart,:])> 0.0))] = 1
    
hp.visufunc.mollview(Fmap_wvpix[tstart,:], title = "Flux t="+str(tstart), xsize = 600, 
                         sub = (8,5,1),hold = False, min = MAX/4.0, max = MAX*3)#, unit = 'Flux')
                         
hp.visufunc.mollview(Fmapwvobs[tstart,:], title = "Flux-obs t="+str(tstart), xsize = 600, 
                         sub = (8,5,2),hold = False, min = MAX/4.0, max = MAX*3)#, unit = 'Flux')
    
hp.mollview(weightpix[tstart,:],title = "Visibility t="+str(tstart),min = 0, max=1,
                sub = (8,5,3))#, unit = 'Window ')  
                
hp.visufunc.mollview(Fweightpix[tstart,:], title = "Illum t="+str(tstart), xsize = 600, 
                         sub = (8,5,4),hold = False, min = 0, max = 1)#, unit = 'Window')
                         
hp.visufunc.mollview(window[tstart, :], title = "Intersection t="+str(tstart), xsize = 600, 
                         sub = (8,5,5),hold = False, min = 0, max = 1.5)
for i in range(1,8):
    
    time = 15000+i*20
    
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
                         
fig.savefig('Fobs and window intersection - longer time')
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
hp.mollview(Fmap_wvobspix[time,:],title = "t="+str(time),min = MIN, max=MAX,
                 unit = 'Flux' )
                 
plt.gca()

                         
                         
                         
                         
                         
                         
                         
                         
                         