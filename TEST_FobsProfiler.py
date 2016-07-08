
# coding: utf-8




# In[10]:

import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate
#from PyAstronomy import pyasl
import healpy as hp
from Parcel_healpy5 import parcel as p
from timer import Timer



# In[11]:

x0 = p(name = "XO-3b",Teff = 6781, e=0.26, Porb = 3.19, a = 0.0454, wadv = 1.0/2, 
                  tau_rad = 20, argp = 346, Rstar = 1.49, Mstar = 1.41, 
                  Rplanet = 1.217, pmax = 4, NSIDE = 4)
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
    
    return x0.Fobs()

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

"""
try:
    from line_profiler import LineProfiler

    def do_profile(follow=[]):
        def inner(func):
            def profiled_func(*args, **kwargs):
                try:
                    profiler = LineProfiler()
                    profiler.add_function(func)
                    for f in follow:
                        profiler.add_function(f)
                    profiler.enable_by_count()
                    return func(*args, **kwargs)
                finally:
                    profiler.print_stats()
            return profiled_func
        return inner

except ImportError:
    def do_profile(follow=[]):
        "Helpful if you accidentally leave in production!"
        def inner(func):
            def nothing(*args, **kwargs):
                return func(*args, **kwargs)
            return nothing
        return inner

def get_number():
    for x in xrange(5000000):
        yield x

@do_profile(follow=[])
def expensive_function():
    
    return hd.Fobs()

result = expensive_function()

# In[12]:
class SomeClass(object):

    def __setattr__(self, name, value):
        print(name, value)
        self.__dict__[name] = value*12

    def __init__(self, attr1, attr2):
        self.attr1 = attr1
        self.attr2 = attr2


sc = SomeClass(attr1=1, attr2=2)

sc.attr1 = 3

#plt.plot(t,TA)"""
# In[12]:                

 
 
# In[5]:



# In[6]:


# In[ ]:
                  

# In[ ]:

                         
                         
                         
                         
                         
                         