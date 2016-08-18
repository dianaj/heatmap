
# coding: utf-8




# In[10]:

import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate
#from PyAstronomy import pyasl
import healpy as hp
from Parcel_healpy6 import parcel as p
from Parcel_healpy6 import fitter as fit
from timer import Timer



# In[11]:

x0 = p(name = "XO-3b",Teff = 6781, e=0.26, Porb = 3.19, a = 0.0454, wadv = 1.0/2, 
                  tau_rad = 20, argp = 346, Rstar = 1.49, Mstar = 1.41, 
                  Rplanet = 1.217, pmax = 3, NSIDE = 8, steps = 100)
                  
x0f = fit(x0, ts = np.linspace(-160000,460000,num = 1500),
          me = "XO-3b",Teff = 6781, e=0.26, Porb = 3.19, a = 0.0454, wadv = 1.0/2, 
                  tau_rad = 20, argp = 346, Rstar = 1.49, Mstar = 1.41, 
                  Rplanet = 1.217, pmax = 3, NSIDE = 8, steps = 100)

# In[11]:

import cProfile

def do_cprofile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()
    return profiled_func



@do_cprofile
def expensive_function():
    
    return x0f.Fobs()

# perform profiling
result = expensive_function()
# In[11]:


# In[11]:

hd = p(Teff = 6079, e=0.6768,Porb = 21.2, a = 0.1589, wadv = 1.0/2, 
                  tau_rad = 20 , argp = 121.71, Rstar = 1.5, Mstar = 1.275, NSIDE = 4)
                  


# In[11]:
with Timer() as t:
    hd.Fleaving()
print "=> elasped Fleaving: %s s" % t.secs

with Timer() as t:
    hd.Fobs()
print "=> elasped Fobs: %s s" % t.secs

with Timer() as t:
    hd.DE()
print "=> elasped DE: %s s" % t.secs


# In[11]:



# In[12]:                

 
 
# In[5]:



# In[6]:


# In[ ]:
                  

# In[ ]:

                         
                         
                         
                         
                         
                         