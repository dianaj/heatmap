# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 10:57:34 2016

@author: diana
"""
""" THIS IS A FILE FOR TESTING THE RADIUS, ANGULAR VELOCITY, AND SSP SOP CALCULATIONS"""

# In[10]:

import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate
#from PyAstronomy import pyasl
import healpy as hp
from Parcel_healpy6 import parcel as p



# In[12]:
gas = p(name = "HD189733b",Teff = 4938, e=0,Porb = 2.21, a = 0.031, epsilon = 1, 
                      argp = 90, Rstar = 0.781, Mstar = 0.846,
                      Rplanet = 1.138, wadv = 2 )
                      
gas2 = p(name = "HD189733b",Teff = 4938, e=0.00001,Porb = 2.21, a = 0.031,  
                      epsilon=1,  argp = 90, Rstar = 0.781, Mstar = 0.846,
                      Rplanet = 1.138 )                      
alpha = gas.alpha
f =  gas.f

t, zt, SSP, SOP = gas.SSP()

d, coordsSOP, weight = gas.visibility(TEST = True)

coordsSSP = d[:,:,1]

# In[12]:

'''THERE's SOME KINKS IN THIS ONE'''

fig, ax = plt.subplots(3,  figsize = (8,24))

ax[0].set_title(gas.name+" $w_{adv} = $"+ '%1.1f'%(gas.wadv/(2*np.pi/gas.P)))
ax[0].plot(t/gas.P, SSP, label = 'SubStelar angle', linewidth = 2)
ax[0].plot(t/gas.P,SOP, label = 'SubObs angle',linewidth = 2)
ax[0].set_xlabel('Time(rotational periods)', fontsize = 14)
ax[0].set_ylabel('Angle ( rads)', fontsize = 14)
ax[0].plot(t/gas.P,(-alpha+180)/57.29, label = "Alpha the phase angle",linewidth = 2, alpha = 0.6)
ax[0].legend(fontsize = 8)

ax[1].set_title("Coordinates")
#ax[1].plot(t, SSP, label = 'SubStelar angle')
#ax[1].plot(t,SSO, label = 'SubObs angle')
ax[1].plot(t/gas.P,d[:,360,1], label = 'parcel location rel. to SSP ',  linewidth = 2)
ax[1].plot(t/gas.P,coordsSOP[:,360], label = 'parcel location rel to SOP', linewidth = 2)
ax[1].set_xlabel('Time(rotational periods)', fontsize = 14)
ax[1].set_ylabel('Angle ( rads)', fontsize = 14)
ax[1].legend(fontsize = 8)
#ax[1].plot(t,(-alpha+180)/57.29, label = "Alpha the phase angle")
#ax[1].set_xlim(10,20)
#ax[1].set_ylim(-6,0)

ax[2].set_title("Difference")
ax[2].plot(t/gas.P, SOP-SSP, label = 'SOP - SSP', linewidth = 2)
ax[2].plot(t/gas.P,-coordsSOP[:,360]+coordsSSP[:,360], label = 'coordsSSP - coordsSOP', linewidth = 2)
#ax[2].plot(t,SSO, label = 'SubObs angle')
ax[2].set_xlabel('Time(rotational periods)', fontsize = 14)
ax[2].set_ylabel('Angle ( rads)', fontsize = 14)
ax[2].plot(t/gas.P,(-alpha+180)/57.29, label = "Alpha the phase angle" , linewidth = 2, alpha = 0.6)
ax[2].plot(t/gas.P,f*5, label = "f_ratio * 5", linestyle = '--', linewidth = 2)
ax[2].legend(fontsize = 8)
#ax[2].set_xlim(0,4000000)
#plt.plot(t,weight[:,[200,300]])
#plt.xlim(0,2000000)
#plt.ylim(-1,1)
#plt.plot(t, ang_vel)
fig.tight_layout()
fig.savefig("SSP and SSO_circ.pdf")


# In[11]:

hd = p(name = "HD17156b", Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, wadv = 0.5, 
                  tau_rad = 20 , argp = 121.71, Rstar = 1.5, Mstar = 1.275)


# In[12]:

alpha = hd.alpha
f =  hd.f

t, zt, SSP, SOP = hd.SSP()

d, coordsSOP, weight = hd.visibility(TEST = True)

coordsSSP = d[:,:,1]

fig, ax = plt.subplots(3,  figsize = (8,24))

ax[0].set_title(hd.name+" $w_{adv} = $"+ '%1.1f'%(hd.wadv/(2*np.pi/hd.P)))#'%1.1f'%gas.wadv)
ax[0].plot(t/hd.P, SSP, label = 'SubStelar angle', linewidth = 2)
ax[0].plot(t/hd.P,SOP, label = 'SubObs angle',linewidth = 2)
ax[0].set_xlabel('Time(rotational periods)', fontsize = 14)
ax[0].set_ylabel('Angle ( rads)', fontsize = 14)
ax[0].plot(t/hd.P,(-alpha+180)/57.29, label = "Alpha the phase angle",linewidth = 2, alpha = 0.6)
ax[0].set_xlim(0,22)
ax[0].legend(fontsize = 8)

ax[1].set_title("Coordinates")
#ax[1].plot(t, SSP, label = 'SubStelar angle')
#ax[1].plot(t,SSO, label = 'SubObs angle')
ax[1].plot(t/hd.P,d[:,360,1], label = 'parcel location rel. to SSP ',  linewidth = 2)
ax[1].plot(t/hd.P,coordsSOP[:,360], label = 'parcel location rel to SOP', linewidth = 2)
ax[1].set_xlabel('Time(rotational periods)', fontsize = 14)
ax[1].set_ylabel('Angle ( rads)', fontsize = 14)
ax[1].set_xlim(0,22)
ax[1].legend(fontsize = 8)
#ax[1].plot(t,(-alpha+180)/57.29, label = "Alpha the phase angle")
#ax[1].set_xlim(10,20)
#ax[1].set_ylim(-6,0)

ax[2].set_title("Difference")
ax[2].plot(t/hd.P, SOP-SSP, label = 'SOP - SSP', linewidth = 2)
ax[2].plot(t/hd.P,-coordsSOP[:,360]+coordsSSP[:,360], label = 'coordsSSP - coordsSOP', linewidth = 2)
#ax[2].plot(t,SSO, label = 'SubObs angle')
ax[2].set_xlabel('Time(rotational periods)', fontsize = 14)
ax[2].set_ylabel('Angle ( rads)', fontsize = 14)
ax[2].plot(t/hd.P,(-alpha+180)/57.29, label = "Alpha the phase angle" , linewidth = 2, alpha = 0.6)
ax[2].plot(t/hd.P,f*5, label = "f_ratio * 5", linestyle = '--', linewidth = 2)
ax[2].set_xlim(0,22)
ax[2].legend(fontsize = 8)
#ax[2].set_xlim(0,4000000)
#plt.plot(t,weight[:,[200,300]])
#plt.xlim(0,2000000)
#plt.ylim(-1,1)
#plt.plot(t, ang_vel)
fig.tight_layout()
fig.savefig("SSP and SSO.pdf")

# In[12]

"""Figure that shows how the radius calculation works, angular velocity calculation works using pyasl """
#from __future__ import print_function, division
#import numpy as np
from PyAstronomy import pyasl
import matplotlib.pylab as plt

# Instantiate a Keplerian elliptical orbit with
# semi-major axis of 1.3 length units,
# a period of 2 time units, eccentricity of 0.5,
# longitude of ascending node of 70 degrees, an inclination
# of 10 deg, and a periapsis argument of 110 deg.

ec = 0.6
ec1= 0.5

Om = 90
Om1 = 180

inc = 87.0
inc1 = 70.0

w0 = 90.0
w1 = 221.0


ke = pyasl.KeplerEllipse(1.3, 2., e=ec, Omega=Om, i=inc, w=w0)
ke1 = pyasl.KeplerEllipse(1.3, 2., e=ec1, Omega=Om1, i=inc1, w=w1)
# Get a time axis
t = np.linspace(0, 4, 200)

# Calculate the orbit position at the given points
# in a Cartesian coordinate system.
'true anomaly'
pos, TA = ke.xyzPos(t,getTA = True)

#print (pos)
pos1, TA1 = ke1.xyzPos(t, getTA = True)
#print("Shape of output array: ", pos.shape)

# x, y, and z coordinates for 50th time point
#print("x, y, z for 50th point: ", pos[20, ::])

# Calculate orbit radius as a function of the
radius = ke.radius(t)
radius1 = ke1.radius(t)


# Calculate velocity on orbit
vel = ke.xyzVel(t)

absvel = np.sqrt(vel[:,0]**2+vel[:,1]**2+vel[:,2]**2)
#print (vel)
vel1 = ke1.xyzVel(t)

# Find the nodes of the orbit (Observer at -z)
ascn, descn = ke.xyzNodes_LOSZ()
ascn1, descn1 = ke1.xyzNodes_LOSZ()

'true anomaly'

#1. get mean anomaly

#MA = ke.meanAnomaly(t)
#EA =[]
#for ma in MA:
#    ea = pyasl.MarkleyKESolver().getE(ma, e)
#    EA.append(ea)
#EA = np.array(EA)    
#true_anomaly = np.arctan(((1+e)*(np.tan(EA/2))**2/(1-e))**(0.5))*2


# Plot x and y coordinates of the orbit
fig, ax = plt.subplots(3,2, figsize = (14,14))

ax[0,0].set_title("Orbital Position: e = "+str(ec)+", Omega = "+str(Om)+", Inc = "+str(inc)+", w = "+str(w0))

ax[0,0].set_ylabel("North ->")
ax[0,0].plot([0], [0], 'k+', markersize=9)
ax[0,0].plot(pos[::,1], pos[::,0], 'bp')
# Point of periapsis
per, = ax[0,0].plot([pos[0,1]], [pos[0,0]], 'md', markersize = 10)
# Nodes of the orbit
asc, = ax[0,0].plot([ascn[1]], [ascn[0]], 'go', markersize=10)
des, = ax[0,0].plot([descn[1]], [descn[0]], 'ro', markersize=10)
ax[0,0].legend([per, asc, des], ["Periapsis","Ascending Node", "Descending Node"])
ax[0,0].set_xlim(-3,3)
ax[0,0].set_ylim(-3,3)

ax[0,1].set_title("Orbital Position: e = "+str(ec1)+", Omega = "+ str(Om1)+", Inc = "+str(inc1)+", w = "+str(w1))
ax[0,1].set_xlabel("East ->")
ax[0,1].set_ylabel("North ->")
ax[0,1].plot([0], [0], 'k+', markersize=9)
ax[0,1].plot(pos1[::,1], pos1[::,0], 'bp')
# Point of periapsis
per, = ax[0,1].plot([pos1[0,1]], [pos1[0,0]], 'md', markersize = 10)
# Nodes of the orbit
asc, = ax[0,1].plot([ascn1[1]], [ascn1[0]], 'go', markersize=10)
des, = ax[0,1].plot([descn1[1]], [descn1[0]], 'ro', markersize=10)
ax[0,1].legend([per, asc, des], ["Periapsis","Ascending Node", "Descending Node"])
ax[0,1].set_xlim(-3,3)
ax[0,1].set_ylim(-3,3)
# Plot RV

ax[1,0].set_xlabel("Time")
ax[1,0].set_ylabel("Radial velocity [length/time]")
ax[1,0].plot(t, absvel, 'r.-')
ax[1,0].set_ylim(-6,6)

ax[1,1].set_xlabel("Time")
ax[1,1].set_ylabel("Radial velocity [length/time]")
ax[1,1].plot(t, vel1[::,2], 'r.-')
ax[1,1].set_ylim(-6,6)

ax[2,0].set_xlabel('Time')
ax[2,0].set_ylabel('Observer angle')
ax[2,0].plot(t, 90-np.arctan(pos[::,0]/pos[::,2]))

ax[2,1].set_xlabel('time')
ax[2,1].set_ylabel('illuminated fraction of planet')
alpha = 180.0-np.array(TA1)*57.2958
ax[2,1].plot(t, pyasl.lambertPhaseFunction(alpha), ls ='None', marker ='o')
ax[2,1].plot(t, 0.5*(1+np.cos(alpha/57.2958)),ls ='None', marker ='o')
plt.tight_layout()
fig.savefig('radius (t)')


# In[12]:


# In[12]:


# In[12]: