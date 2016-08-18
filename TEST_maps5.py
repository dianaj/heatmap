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
from Parcel_healpy5 import parcel as p
#from timer import Timer



hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, wadv = 0.5, 
                  tau_rad = 20 , argp = 121.71, Rstar = 1.5, Mstar = 1.275)
                  

# In[11]:

t, d, Fmapwvobs, weight, weightpix, Fwv = hd.Fobs(MAP = True)
t, d, Fmap_wv, Fmap_wvpix, Fleavingwv, Ftotal = hd.Fleaving(MAP = True)
d, Fweight = hd.illum()

d, Fweightpix = hd.shuffle(d = d, quantity = Fweight)

d, dshuffle = hd.shuffle()

"""SHUFFLE ISNT WORKING AS EXPECTED SOME OF THE TIME!!! CHECK IT"""
# In[11]:


"""CUTE FIGURE"""

window = np.zeros(weight.shape)
tstart = 15300
t = hd.t
r = hd.radius
F0 =  (hd.sigmaB * hd.Teff**4*np.pi*hd.Rplanet**2*(hd.Rstar/min(r))**2) 
MAX = 0.43
MIN = 0.25

fig = plt.figure(figsize  =(6,4))
hp.visufunc.mollview(dshuffle[tstart,:], title = "Planet Temperature at a Moment in Time", xsize = 600, 
                         hold = False , sub =111, min = MIN, max = MAX,
                         unit = "($T/T0$)")
fig.savefig("sometemp.pdf")      

# In[ ]:
"""Figure that shows all the maps"""
MAXf = 23.5
MINf = 22


MAX = 10**23*0.4
MIN = 10**16

fig= plt.figure(figsize = (10,16))
tstart = 11250
window[tstart,np.where((np.array(weightpix[tstart,:]) > 0.0 ) & (np.array(Fweightpix[tstart,:])> 0.0))] = 1
    
hp.visufunc.mollview(np.log10(Fmap_wvpix[tstart,:]), title = "Log10Flux t="+'%1.2f'%((t[tstart]/hd.Porb)%1)+'Porb', xsize = 600, 
                         sub = (8,5,1),hold = False, min = MINf, max = MAXf)#, unit = 'W/m')
                         
hp.visufunc.mollview((Fmapwvobs[tstart,:]), title = "F-obs t="+'%1.2f'%((t[tstart]/hd.Porb)%1)+'Porb', xsize = 600, 
                         sub = (8,5,2),hold = False, min = MIN, max = MAX*3.5)#, unit = 'Flux')
    
hp.mollview(weightpix[tstart,:],title = "Vis t="+'%1.2f'%((t[tstart]/hd.Porb)%1)+'Porb',min = 0, max=1,
                sub = (8,5,3))#, unit = 'Window ')  
                
hp.visufunc.mollview(Fweightpix[tstart,:], title = "Illum t="+'%1.2f'%((t[tstart]/hd.Porb)%1)+'Porb', xsize = 600, 
                         sub = (8,5,4),hold = False, min = 0, max = 1)#, unit = 'Window')
                         
hp.visufunc.mollview(window[tstart, :], title ="f-illum t="+'%1.2f'%((t[tstart]/hd.Porb)%1)+'Porb', xsize = 600, 
                         sub = (8,5,5),hold = False, min = 0, max = 1.5)
for i in range(1,8):
    
    time = tstart+i*800
    
    window[time,np.where((np.array(weightpix[time,:]) > 0.0 ) & (np.array(Fweightpix[time,:])> 0.0))] = 1
    
    hp.visufunc.mollview(np.log10(Fmap_wvpix[time,:]), title = "t="+'%1.2f'%((t[time]/hd.Porb)%1), xsize = 600, 
                         sub = (8,5,int(5*i+1)),hold = False, min = MINf, max = MAXf, cbar=False)#, unit = 'Flux')
                         
    hp.visufunc.mollview((Fmapwvobs[time,:]), title = "t="+'%1.2f'%((t[time]/hd.Porb)%1), xsize = 600, 
                         sub = (8,5,int(5*i+2)),hold = False, min = MIN, max = MAX*3.5, cbar=False)#, unit = 'Flux')
    
    hp.mollview(weightpix[time,:],title = "t="+'%1.2f'%((t[time]/hd.Porb)%1),min = 0, max=1,
                sub = (8,5,int(5*i+3)),cbar=False)#, unit = 'Window ')  
                
    hp.visufunc.mollview(Fweightpix[time,:], title = "t="+'%1.2f'%((t[time]/hd.Porb)%1), xsize = 600, 
                         sub = (8,5,int(5*i+4)),hold = False, min = 0, max = 1, cbar=False)#, unit = 'Window')
                         
    hp.visufunc.mollview(window[time, :], title = "t="+'%1.2f'%((t[time]/hd.Porb)%1), xsize = 600, 
                         sub = (8,5,int(5*i+5)),hold = False, min = 0, max = 1.5, cbar=False)#, unit = 'Intersection')
                         
fig.savefig('Fobs and window intersection - long time.pdf')
#sub = (int('42'+str((2*(i+1)))))



# In[8]:
'Animation that shows all the maps'


MAX = 2*1*10**22
MIN = 1*10**16



for i in range(72):
    fig= plt.figure(figsize = (40,3))
    tstart = 10300
    
    time = 10300+i*100
    
    window[time,np.where((np.array(weightpix[time,:]) > 0.0 ) & (np.array(Fweightpix[time,:])> 0.0))] = 1
        
    hp.visufunc.mollview(Fmap_wvpix[time,:], title = "Flux t="+str(time), xsize = 600, 
                             sub = (1,5,1),hold = False, min = MAX/4.0, max = MAX*3)#, unit = 'W/m')
                             
    hp.visufunc.mollview(Fmapwvobs[time,:], title = "Flux-obs t="+str(time), xsize = 600, 
                             sub = (1,5,2),hold = False, min = MAX/4.0, max = MAX*3)#, unit = 'Flux')
        
    hp.mollview(weightpix[time,:],title = "Visibility t="+str(time),min = 0, max=1,
                    sub = (1,5,3))#, unit = 'Window ')  
                    
    hp.visufunc.mollview(Fweightpix[time,:], title = "Illum t="+str(time), xsize = 600, 
                             sub = (1,5,4),hold = False, min = 0, max = 1)#, unit = 'Window')
                             
    hp.visufunc.mollview(window[time, :], title = "Intersection t="+str(time), xsize = 600, 
                             sub = (1,5,5),hold = False, min = 0, max = 1.5)
                             
    fig.savefig('Animate - Fobs and window #'+ str(i))

'''
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
'''                         
#fig.savefig('Fobs and window intersection - short time')
#sub = (int('42'+str((2*(i+1)))))

# In[8]:

Fstar, Fstarwv = hd.Fstar()
for i in range(72):
    fig= plt.figure(figsize = (40,3))
    tstart = 10300
    
    time = 10300+i*100
    
    plt.plot((t/hd.Porb)[tstart:time:], Fwv[tstart:time:]/Fstarwv, linewidth = 2, color = 'midnightblue')
    plt.xlim(t[tstart]/hd.Porb,t[tstart+72*100]/hd.Porb)
    plt.ylim(0,0.00075)
    plt.xlabel('Time (orbits)', fontsize = 14, fontweight = 'bold')
    plt.ylabel('Observed flux ($F_{star}$)', fontsize = 14, fontweight = 'bold')
    
    fig.savefig('Animate - Fwv #'+ str(i))
    



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