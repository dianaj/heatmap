# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 17:22:17 2016

@author: diana
"""
""" TO DO:

- make sure this works for negative wadv 
- especially check SOP, SSP, coordinates, illum, visibility and DE

- make visibility a separate function called by the DE 
- Try to get Temperatures to not start at 0. 
        They should start somewhere good... based on taurad and wadv and T0
        
- start documenting!!!

- include TESTs for most functions that are easy to run. 
        maybe have a TEST file and in case of TEST = TRUE have the function call that file and 
        draw what it does. 
        
        
-  NEED A GET ITEM METHOD SO WE CAN CHANGE VALUES IN AN OBJECT WITHOUT CHANGING THE WHOLE THING
"""



import numpy as np
#import matplotlib.pyplot as plt
#from scipy import integrate
from PyAstronomy import pyasl
import healpy as hp
import time 



class parcel:
    
    '''
    A gas parcel on a planet. Default are values kind of made up 
    INPUT FOR ALL FUNCTIONS 
    
    :::  pmax in orbital periods --- default is 3
    :::  steps --- meaning timesteps used to integrate the DE per 24 hour period.
                           --- default is 300
    :::  NSIDE --- healpix parameter that specifies number of pixels that describe location on the planet. 
               --- default is 8 - means 798 equal area subdivisions of planet surface
               --- 4 means 192 pixels and works fine as well (and faster)
    :::  wavelenght ---  used for flux calculations. values should be entered in micrometers
                    ---  default is 8 um                    
    '''                     
    AU = 1.49597870700* 10**11 #m
    sigmaB = 5.670367*10**(-8) #W m−2 K−4) 
    K = 1.38064852*10**(-23) #gas constant #(m^2*kg)/(s^2*K) 
    days = 86400 #seconds
    
    rsun = 695700000 #m
    Msun = 1.989 *10**(30) #kg
    G = 6.67408*10**(-11) #N*m^2/ kg^2
    MJup = 1.898 *10**27 #kg
    RJup  = 69911000 #m
    
    h = 6.626070040*10**(-34) #J*s
    c = 299792458 #m/s
    
    
    def __init__(self, name = 'HotWaterEarth', Teff =6000.0, Rstar = 1.0, Mstar = 1.0, Rplanet = 1.0870, a = 0.1589, 
                 e = 0.0, argp = 0, theta = np.pi/2, 
                 A = 0, ro = 100.0 , cp = 4200.0, H= 5.0, hh = 14121.0, Porb = -1, 
                 wadv = 1.2, tau_rad = 20, epsilon = None):
        self.name = name    # instance variable unique to each instance


        self.Teff = Teff #temp star
        self.Rstar = Rstar * parcel.rsun #radius star
        self.Mstar = 1.0* parcel.Msun #mass of the star in solar masses
        self.Rplanet = Rplanet * parcel.RJup
        self.a = a * parcel.AU #semimajor axis
        self.e = e #eccentricity
        self.argp = argp #angle betwen periatron and transit in degrees 
        self.theta = theta #latitude of gas parcel in angle starting from north pole(0 at north pole, pi/2 at equator)
        self.A = A #Bond albedo
        self.ro = ro #mass density of gas parcel kg/m^3
        self.cp = cp #heat capacity of gas parcel
        self.hh =  hh #scale height of planet 
        self.H = H* self.hh #thickness of parcel (will be a few scale heights)
        if Porb <= 0.0:
            self.Porb = 2*np.pi*(self.a**3/(self.Mstar*parcel.G))**(0.5)
        else:
            self.Porb =  Porb* self.days #orbital period in seconds formula --
            
        self.P =  self.Prot() #self.Prot() P* parcel.days#period of rotation planet in seconds 
        
        self.ch=self.ro*self.cp*self.H # J/m^2 heat capacity of gas column
        self.wadv = (2*np.pi/self.P)*wadv #think this is a param we need to fit for. In units of how much faster 
                                        #the parcel moves compared to the period of the planet. 
                                        #+ moves with rotation of planet, - against 
        self.T0 = Teff*(1-A)**(0.25)*(self.Rstar/(self.a*(1-self.e)))**(0.5) # * np.sin(theta))**0.25 when the parcel 
                                                                        #lives away from the equator
        #self.tau_rad = self.ch/(sigmaB*self.T0**3)
        #self.epsilon = self.wadv*self.tau_rad
        if epsilon is None:        
            
            self.tau_rad = tau_rad * 3600.0 #it was in hours now its in seconds
            self.epsilon = self.tau_rad * np.abs(self.wadv)
            
        else:
            self.tau_rad = np.abs(self.epsilon/self.wadv)
        
        self.rotationsPerOrbit = int(self.Porb/self.P)+1 #used for giving the default time lenght for DE
        self.rotationsPerDay = int(self.P/self.days) #used for giving the default precision for DE
        #self.Teq = (((1-A)*F(t)*np.sin(theta)*np.max(cos(phi(t)),0))/sigmaB)**(0.25)
        
    
    
    def print_stuff(self):
            print("""Name  = {0} (Name of model)
                    Teff = {1} K (Temperature Star)
                    Rstar = {2} R-sun (Radius Star in solar radii)
                    Rplanet = {3} R-sun (Radius Planet in solar radii)
                    a = {4} AU (semimajor axis)
                    e = {5} (eccentricity)
                    argp = {6} #angle betwen periatron and transit in degrees 
                    theta = {7} (latitude of gas parcel - 0 North Pole, pi/2 equator
                    A = {8} (Bond Albedo)
                    ro = {9} kg/m^3 (Mass density of gas parcel)
                    cp = {10} J/kg (heat capacity of gas parcel)
                    h = {11} km (scale height of planet atmosphere)
                    H = {12} scale heights (thickness of parcel)
                    P = {13} days (period of rotation of planet) **might need fixing
                    Porb = {14} days (orbital period of planet)
                    ch = {15} - some type of units??? (heat capacity of gas column) )
                    wadv = {16} - units of angular frequency 2Pi/planetperiod *wadv??
                    T0 = {17} K (Temperature of substellar point at periastron)
                    tau_rad = {18} hrs (radiative time scale)
                    epsilon = {19} dimensionless circulation efficeincy param
                    default rotationsPerOrbit = {20} #used for giving the default time lenght for DE
                    default rotationsPerDay = {21} #used for giving the default precision for DE
                    **** all of these are converted to SI units but they need 
                    to be entered in the units mentioned here ***
                    """.format (self.name, self.Teff, self.Rstar/parcel.rsun, self.Rplanet/parcel.rsun, 
                                self.a/parcel.AU, self.e, self.argp, 
                                self.theta, 
                                self.A, self.ro,
                               self.cp, self.hh/1000, self.H/self.hh,self.P/parcel.days, self.Porb/parcel.days, 
                                self.ch, self.wadv*self.P/ (2*np.pi), 
                                self.T0, self.tau_rad/ 3600.0, 
                                self.epsilon, int(self.Porb/self.P), int(self.P/self.days) ))
                                
     
    """ FUNCTIONS FOR DEFINING THE PLANET CONDITIONS, STELLAR FLUX """
       
    def Prot(self):
        ''' 
        Calculates rotational period :: wmax :: of the planet when given orbital period
        :::  wmax is chosen to match orbital angular velocity of planet at periastron ::::
        '''
        
        Prot = self.Porb*(1-self.e)**(1.5)*(1+self.e)**(-0.5)        
        return Prot
    
    

    
    def Fstar(self, wavelenght = 8.0):
        '''   
        Calculates total flux emmitted by the STAR and wavelenght dependent flux.
        INPUT   ::: wavelenght in micrometers. Default is 8 um.
        RETURNS ::: total flux (F), wavelenght dependent flux (Fwv)
        '''
        wv = wavelenght * 10**(-6) #wavelenght was in micrometers, it's now in meters
        F = self.sigmaB*self.Teff**4*np.pi*self.Rstar**2
        Fwv = (2*self.h*self.c**2/wv**5)*(1/(np.e**((self.h*self.c)/(wv*self.K*self.Teff))-1))*np.pi*self.Rstar**2 
        return F, Fwv
        
        
    def Finc(self, pmax = 3.0, steps = 300.0):
        
        '''
        Calculates flux incident on the substellar point of the planet.
        Used in DE for '"incoming energy" value as Finc/ Finc_max
        
        INPUT   :::  pmax in orbital periods --- default is 3
                :::  steps :: meaning timesteps used to integrate the DE per 24 hour period.
                           :: default is 300
                          
        CALLS   :::  self.radius (to get distance from planet to star as a function of time )
        RETURNS :::  incident flux at substellar point array --- 1 flux value for each time step 
        '''
        
        return self.sigmaB*self.Teff**4*(self.Rstar/np.array(self.radius(pmax,steps)[1]))**2
        
    def Finc_hemi(self, pmax = 3.0, steps = 300.0):
        
        '''
        Calculates flux incident on the planet. Used for normalizing things mostly. 
        
        INPUT   ::: pmax , steps
        
        CALLS   ::: self.Finc
        
        RETURNS ::: array :: flux incident on planet hemisphere at each time step
        '''
        
        #pmaxi = pmax * self.rotationsPerOrbit
        #stepsi = steps * self.rotationsPerDay
        # Einc = Finc* integral(phi: 0->2pi, theta: 0->pi/2) rp^2*cos(theta)sin(theta) dtheta dphi
        # ---- integral overplanet coordinates
        return self.Finc(pmax, steps)*np.pi*self.Rplanet**2
            
    
    
    """FUNCTIONS FOR DEFINING THE COORDINATE SYSTEM AND MOVEMENT AS A FUNCTION OF TIME """
    """::: they all use module pyasl from PyAstronomy ::: """
    
    def radius(self, pmax=3.0, steps = 300):
        '''
        INPUT   ::: pmax, steps
                
        RETURNS ::: t :: a time array in seconds 
                    r :: orbital separation radius array (as a function of time) 
        ::: uses module pyasl from PyAstronomy :::          
        '''
        
        pmaxi = pmax * self.rotationsPerOrbit
        stepsi = steps * self.rotationsPerDay
        
        ke = pyasl.KeplerEllipse(self.a, self.Porb, e=self.e, Omega=180., i=90.0, w=self.argp)

        # Get a time axis
        tmax = self.P*pmaxi
        Nmin = int((pmaxi)*stepsi)
        t = np.linspace(0,tmax,num=Nmin)

        radius = ke.radius(t)
        return t, radius
    
    def ang_vel(self, pmax = 3.0, steps = 300):
        '''
        INPUT   ::: pmax, steps
                
        
        RETURNS ::: t --- a time array in seconds 
                    ang_vel --- angular velocity at each step in the orbit (as a function of time) '''
                  
        "COULD BE CONBINED WITH RADIUS, IT USES THE SAME STUFF"
        
        pmaxi = pmax * self.rotationsPerOrbit
        stepsi = steps * self.rotationsPerDay
        
        ke = pyasl.KeplerEllipse(self.a, self.Porb, e=self.e, Omega=180., i=90.0, w=self.argp)
        tmax = self.P*pmaxi
        Nmin = int((pmaxi)*stepsi)
        t = np.linspace(0,tmax,num=Nmin)
        vel = ke.xyzVel(t)
        absvel = np.sqrt(vel[:,0]**2+vel[:,1]**2+vel[:,2]**2)
        ang_vel = absvel/ self.radius(pmax, steps)[1]
        #print max(ang_vel)*self.a*(1-self.e)
        
        return t, ang_vel
        
    
    def f_ratio(self, pmax = 3.0, steps = 300):
        '''
        Takes : pmax  :: number of rotational periods (float)
                steps :: number of steps per period (int)
        
        Returns : t     :: a time array in seconds 
                  alpha :: the phase angle in degrees = 90 - argp - TrueAnomaly ::THIS IS A BIT WEIRD 
                  f     :: illuminated fraction of planet -- 1/2(1-cos(alpha))'''
                  
        "COULD BE CONBINED WITH RADIUS, IT USES THE SAME STUFF"
        
        pmaxi = pmax * self.rotationsPerOrbit
        stepsi = steps * self.rotationsPerDay
        
        ke = pyasl.KeplerEllipse(self.a, self.Porb, e=self.e, Omega=180., i=90.0, w=self.argp)
        tmax = self.P*pmaxi
        Nmin = int((pmaxi)*stepsi)
        t = np.linspace(0,tmax,num=Nmin)
        pos, TA = ke.xyzPos(t,getTA = True)
        
        # i want this to be 0 at transit. f = 90-w at transit so alpha = 90-w - f
        alpha = (90-self.argp) - np.array(TA)*57.2958
        
        f = 0.5*(1-np.cos(alpha*np.pi/180.0))
        #f = pyasl.lambertPhaseFunction(alpha)
        return t, alpha, f
    
    def SSP(self, pmax =3, steps =300):
        '''
        INPUT   ::: pmax, steps
                
        
        RETURNS ::: t --- a time array in seconds 
                    ang_vel --- angular velocity at each step in the orbit (as a function of time) 
        '''
        '''
        x: int or float.

        returns: True if x is odd, False otherwise
        '''
        
        pmaxi = pmax * self.rotationsPerOrbit
        stepsi = steps * self.rotationsPerDay
        
        t, ang_vel = self.ang_vel(pmax, steps)
        t, alpha, f = self.f_ratio(pmax, steps)
        'z0= orbital position at t = 0'
        #z0=self.argp*np.pi/180
        z0=0
        'orbital phase'
        zt = np.empty(int(pmaxi*stepsi))
        deltat = (self.P*1.0)/(stepsi*1.0)
        zt[0]=z0
        for i in range(1,int(pmaxi*stepsi)):  
            zt[i]= zt[i-1]+ang_vel[i-1]*deltat
        
        SSP =(zt-((self.wadv)*t -((self.wadv*t)/(2*np.pi)).astype(int)*2*np.pi)-
        (t/self.Porb).astype(int)*2*np.pi)
        
        SOP = (((-alpha[0]+180)*np.pi/180)-
                ((self.wadv)*t -((self.wadv*t)/(2*np.pi)).astype(int)*2*np.pi))
        #SOP = SSP + ((-alpha+180)*np.pi/180)
        return t, zt, SSP, SOP
    

    """ILLUMINATION AND VISIBILITY """    
    
    def illum(self, pmax =3, steps =300, NSIDE = 8) :  

        """Class methods are similar to regular functions.

        Note:
            Do not include the `self` parameter in the ``Args`` section.

        Args:
            param1: The first parameter.
            param2: The second parameter.

        Returns:
            True if successful, False otherwise.

        """
        """CREATES COORDINATE MATRIX WRT SUBSTELLAR POINT AND STARTING TEMPERATURES ARRAY.
        REPLACES 0 TEMPERATURES WITH 0.01.
        COORD/ TEMPERATUTE ARRAY D AND WEIGHT FUNCTION """
        
        pmaxi = pmax * self.rotationsPerOrbit
        stepsi = steps * self.rotationsPerDay              
            
        #TIME
        Nmin = int((pmaxi)*stepsi)
        tmax = self.P*pmaxi
        t = np.linspace(0,tmax, num=Nmin)
        
        #in case you want the time in planet days        
        dayss =(t)*(np.abs(self.wadv))/(2*np.pi)  
            
        #POSITIONS
        coords =np.empty(((hp.nside2npix(NSIDE)),2))
        for i in range(hp.nside2npix(NSIDE)):
            coords[i,:]=(np.array(hp.pix2ang(NSIDE, i)))
        phis = coords[:,1]
        thetas = coords[:,0]
            
        #STARTING TEMPERATURES
        T = (np.sin(thetas)**0.25)*(np.cos(phis)+np.abs(np.cos(phis)))/2
        T[np.where(T<0.01)]=0.01 
        
        #3D ARRAY TO CONTAIN EVERYTHING. CALLED IT d   
        c= np.array(zip(thetas,phis,T)).reshape(-1,hp.nside2npix(NSIDE),3)
        d = np.repeat(c,Nmin, axis = 0)
            
        #CREATE THE COORDINATE ARRAY. THERE'S 2 CASES: CIRCULAR ORBIT AND ECCENTRIC ORBIT    
        if self.e == 0.0:
                
            deltaphi = (2*np.pi/stepsi)* self.wadv/(2*np.pi/self.P)
            d[:,:,1]= (phis+deltaphi* np.array(range(0,Nmin)).reshape(-1,1))%(2*np.pi) #-((self.wadv)*t[i] - int(self.wadv*t[i])) #update location of gas parcel
                    
            
        else:
            t, zt, SSP, SOP = self.SSP(pmax, steps)
            d[:,:,1]= (phis.copy())+ (SSP.reshape(-1,1))*(np.sign(self.wadv))


        #ILLUMINATION WEIGHT FUNCTION. WILL GET PASSED TO DE ALONG WITH D THE COORDINATE MATRIX AND THE STUPID TIME MATRIX
        Fweight = ((np.cos(d[:,:,1])+ np.abs(np.cos(d[:,:,1])))/2.0 * np.sin(d[:,:,0]))
                  

        return t, d, Fweight

    
    
        
    def visibility(self, pmax =3, steps =300, NSIDE = 8, d = None, TEST = False):
        '''Calculates the visibility of each gas parcel on the planet..

        Note:
             weight = ((np.cos(coordsSOP)+np.abs(np.cos(coordsSOP)))/2.0)*np.sin(thetas)

        Args:
            pmax (int): The first parameter.  
            steps (int): The second parameter.   
            NSIDE (power of 2):   Will determine how many gas parcels we'll have on a planet.       
            d (3d array) :  time, location:(thetas, phis,) Temperature (at each location).    
            TEST (bool): Usually false   

        Returns:
            if TEST = False
                t      : time array in seconds\
                weight : visibility array to be applied to flux array in coordinates wrt SSP\
                
            if TEST = True
                t      : time array in seconds\
                d      : coordinates wrt SSP and temperature aray\
                coordsSOP : coordinates relative to the suborserver point\
                weight : visibility array to be applied to flux array in coordinates wrt SSP\

        '''
        if d is None:
            t, d, Fweight = self.illum(pmax, steps, NSIDE)
        
        else: 
            
            d = d
                
        #t, zt, SSP, SOP = self.SSP(pmax, steps)
        t, alpha, f = self.f_ratio(pmax, steps)
        
        phis = d[:,:,1] # location of gas parcel on planet relative to SSP
        thetas = d[:,:,0]
        
        #coordsSSP = (phis)
        coordsSOP = phis-(-alpha.reshape(-1,1)+180)*np.pi/180
        
        #2#coordsSOP = phis+(zt%(2*np.pi)).reshape(-1,1)
        
        #coordsSOP = phis+(zt).reshape(-1,1)+(-alpha[0]+180)*np.pi/180
        #coordsSOP = phis - (SOP - SSP).reshape(-1,1)    

        weight = ((np.cos(coordsSOP)+np.abs(np.cos(coordsSOP)))/2.0)*np.sin(thetas)


        if TEST :
            return t, d, coordsSOP, weight
        else:
            return t, weight

        
        
        
    def DE(self, pmax=3.0, steps = 300.0, NSIDE = 8):
            'pmax is in number of periods you want to integrate for.'
            'The number of steps is per period'
            
            """returns: -time (days) and 
                        -d:  pixel number , angular location of parcel(thata,phi), temperature of gas parcel
                        at each time value
            """
            
            tic = time.time()
            print "Starting DE"
            
            pmaxi = pmax * self.rotationsPerOrbit
            stepsi = steps * self.rotationsPerDay
            
            t, d, Fweight  =self.illum(pmax, steps, NSIDE)
            
            
            if self.e == 0.0:
                
                deltaphi = (2*np.pi/stepsi)* self.wadv/(2*np.pi/self.P)
                
                for i in range(1,len(t)):#phis.shape[2]
                    
                    dT =(1.0/self.epsilon*Fweight[i-1]-(d[i-1,:,2])**4 )*deltaphi #calculate change

                    d[i,:,2]= d[i-1,:,2].copy()+ dT #update temperature array
                    
                toc = time.time()
                print ("Time this took is: " , str(toc-tic), "seconds")
                
                return t, d
            
            else:
                
                deltat = self.P/stepsi
                
                deltat_ = deltat/self.tau_rad
                wrot = (2*np.pi/self.P)* self.wadv/(2*np.pi/self.P)                   
                deltaphi = wrot*deltat 
                    
                #parcel.sigmaB*self.Teff**4*(self.Rstar/np.array(self.radius(pmax,steps)[1]))**2    
                'normalized flux -- (minimum radius/ radius(t))**2'   
                Fstar = (self.Finc(pmax, steps).reshape(-1,1)) #*Fweight
             
                F = Fstar/(parcel.sigmaB*self.Teff**4*(self.Rstar/(self.a*(1-self.e)))**2)
                #F = ((self.a*(1-self.e)/self.radius(pmax, steps)[1])**2)
               
                for i in range(1,len(t)):
                        
                        
                        dT =(( F[i-1]*Fweight[i-1] - (d[i-1,:,2])**4 )* (deltat_))
                        
                        d[i,:,2]= d[i-1,:,2].copy()+dT #update temperature array
                        
                        
                toc = time.time()
                print ("Time DE took: "+ str(toc-tic))
                return t, d 

    
    def Fleaving(self, pmax = 3.0, steps = 300.0, NSIDE = 8, wavelenght = 8.0):#, TEST = False):
        "The substellar point is moving in this case"
        #wavelenght was in micrometers
        pmaxi = pmax * self.rotationsPerOrbit
        stepsi = steps * self.rotationsPerDay
        
        wv = wavelenght*10**(-6)        
        
        t, d = self.DE(pmax, steps, NSIDE)
        
        
        'this should be done already'
        #d[:,:,2][d[:,:,2] < 0.01] = 0.01 #replace the zeroes at the beginning so they dont overflow
        
        #print min(d[:,:,2])
        #Fwv_ = ((2*self.h*self.c**2/wavelenght**5)*
        #(1/(np.e**((self.h*self.c*10**6)/(wavelenght*self.K*np.array(d[:,:,2].copy())*self.T0))-1)))
        
        #Fwv = Fwv_*10**30
        
        a = (2*self.h*self.c**2/wv**5)
        
        #print a 
        b = (self.h*self.c*10**6)/(wavelenght*self.K*self.T0)
        
        #print b
        
        #Fwv= a* 1/(np.e**(b/np.array(d[:,:,2].copy()))-1)
        
        #Fwv = np.zeros(d[:,:,2].shape)
        #for i in range(len(days)):
            
        
        Fwv = a* 1/(np.expm1(b/np.array(d[:,:,2].copy())))
            
            
        
            #print min(d[i,:,2])
        
        
        
        #Fwv_ = (np.log(2*self.h*self.c**2)-np.log(wv**5) + np.log(1) - 
        #np.log(np.e**((self.h*self.c)/(wv*self.K*np.array(d[:,:,2].copy())*self.T0))-1))
        
        #Fwv = np.e**(Fwv_)        
        
        "Get the flux"
        dA = hp.nside2pixarea(NSIDE)*self.Rplanet**2
  
        Fmap_wv = (Fwv.copy()*dA)#/Fwvstar
        
        Ftotal_ = (self.sigmaB * (d[:,:,2]*self.T0)**4)*dA
        
        crap, Fmap_wvpix = self.shuffle(d, Fmap_wv, pmax, steps, NSIDE)

        
        """IM INTEGRATING TWICE BECAUSE I WANT TO MAKE SURE LEAVING FLUX AND SHUFFLED
        LEAVING FLUX ARE THE SAME. THEY ARE """
        
        Fleavingwv = np.zeros(int(stepsi * pmaxi))
        Ftotal = np.zeros(int(stepsi * pmaxi))
        for i in range(int(stepsi * pmaxi)):
            
            Fleavingwv[i] = np.sum(Fmap_wv[i,:])
            Ftotal[i] = np.sum(Ftotal_[i,:])
            
        #print Fleavingwv == Fleavingwvpix


            
        #if TEST :
        #    return days, d, Fmap_wv, Fmap_wvpix, Fleavingwv, Fleavingwvpix     
        #else:
        return t, d, Fmap_wv, Fmap_wvpix, Fleavingwv, Ftotal
        

        
        

    
    def Fobs(self, pmax = 3.0, steps = 300.0, NSIDE = 8, wavelenght = 8.0, PRINT = False):
        "surface flux (leaving) from one hemisphere of the planet "
        "the flux leaving the planet will be different"
        
        tic = time.time()
        print ("Starting Fobs")
        
        'call the functions we need'
        t, d, Fmap_wv, Fmap_wvpix,Fleavingwv, Ftotal = self.Fleaving(pmax, steps, NSIDE, wavelenght)
        
        t, weight = self.visibility(pmax, steps, NSIDE, d)

        crap, weightpix = self.shuffle(d, weight, pmax, steps, NSIDE)
        
        
        'start doing what this function does'
        #Fmapwvobs = Fmap_wv*weight
        
        pmaxi = pmax * self.rotationsPerOrbit
        stepsi = steps * self.rotationsPerDay
        
        Nmin = int((pmaxi)*stepsi)
        
        Fmapwvobs = Fmap_wvpix*weightpix
        Fmapwvobs_check = Fmap_wv*weight
        Fwv = np.empty(Nmin)
        Fwv_check = np.empty(Nmin)
        for i in range(Nmin):

            Fwv[i] = np.sum(Fmapwvobs[i,:])
            Fwv_check[i] = np.sum(Fmapwvobs_check[i,:])
            
        print "DO we get the same flux before and after shuffling? " 
        print (Fwv == Fwv_check)

        if PRINT == True:
            tt = self.P*np.array(days).reshape(-1,1)
            #flux = np.array(Fleaving).reshape(-1,1)
            fluxwv = np.array(Fwv).reshape(-1,1)
            np.savetxt('observedflux_planet'+str(steps)+'_steps_per_period_'+str(pmaxi)+'periods_'
                       +str(NSIDE)+'_NSIDE.out', 
                       zip( tt, fluxwv), header = "time(planetdays),time(s), outgoing flux, outgoing flux per wv")
        toc = time.time()
        print ('Done with Fobs')
        print ('Time Fobs took: ' + str(toc-tic) + 'seconds')
        
        if PRINT == False:
            return  t, d, Fmapwvobs, weight, weightpix, Fwv #*weight
    

    
                    
    def shuffle (self, d = None, quantity = None, pmax = 3.0, steps = 300, NSIDE = 8):
        "d is going to be a array of the form [time, [NSIDE, [phis, thetas, other quantities]]]"
        #pmaxi = pmax * self.rotationsPerOrbit
        #stepsi = steps * self.rotationsPerDay
        
        if (d is None) and (quantity is None):        
            t, d = self.DE(pmax, steps, NSIDE)
            dd = np.array(d.copy())
            quantity = np.array(d[:,:,2].copy())
            
        else:
            dd = np.array(d.copy())
        
        timesteps = len(dd[:,0,0])
            #quantity = quantity
        
        #order = (hp.ang2pix(NSIDE, d[:,:,0], d[:,:,1])).astype(int)
        #print order.shape
        #quantity = quantity[order]
        order = np.zeros((timesteps, hp.nside2npix(NSIDE)))
        for i in range(timesteps):        
            order[i,:] = (hp.ang2pix(NSIDE, dd[i,:,0], dd[i,:,1])).astype(int)
                    
            quantity[i,:] = quantity[i,order[i,:].astype(int)]
        
        print("shuffled quantity" )
        return d, quantity
            
    def contours_illum():
                #contours of equal illumination theta
        
        c_thetas1 = []
        c_phis1 = []
        for i in range(int(pmax*steps)):
            c_t = thetas[i,np.where(np.abs(Fweight[i,:]-1)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            c_p = coordsSSP[i,np.where(np.abs(Fweight[i,:]-1)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            c_thetas1.append(c_t)
            c_phis1.append(c_p)
        
        c_thetas7 = []
        c_phis7 = []
        for i in range(int(pmax*steps)):
            c_t = thetas[i,np.where(np.abs(Fweight[i,:]-0.7)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            c_p = coordsSSP[i,np.where(np.abs(Fweight[i,:]-0.7)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            'should we put these angles in substellat point coordinates or what? '
            c_thetas7.append(c_t)
            c_phis7.append(c_p)
            
            
    def contours_vis():
                #contours of equal illumination theta
        
        c_thetas1 = []
        c_phis1 = []
        for i in range(int(pmax*steps)):
            c_t = thetas[i,np.where(np.abs(weight[i,:]-1)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            c_p = coordsSSP[i,np.where(np.abs(weight[i,:]-1)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            c_thetas1.append(c_t)
            c_phis1.append(c_p)
        
        c_thetas7 = []
        c_phis7 = []
        
        for j in range(int(pmax*steps)):
            c_t = thetas[j,np.where(np.abs(weight[j,:]-0.7)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            c_p = coordsSSP[j,np.where(np.abs(weight[j,:]-0.7)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            'should we put these angles in substellar point coordinates or what? '
            c_thetas7.append(c_t)
            c_phis7.append(c_p)                   
    
    def findT (self, pmax=3.0, steps = 300.0):
        
        #tmax = self.P*pmax
        #Nmin = int((pmax)*300)
        #deltat = tmax/Nmin
        pmaxi = pmax * self.rotationsPerOrbit
        stepsi = steps * self.rotationsPerDay
        
        phi, days, T = self.DE(pmax, steps)
        
        
        deltaphi = 2.0/steps
        Tmax = np.max(T[int(stepsi*(pmaxi-2))::])
        
            
        for i in range(int(stepsi*(pmaxi-2)), len(phi)):
            
            
            
                if deltaphi >= np.abs(1.5 - (phi[i]-2*(pmaxi-2))):
                    #print np.abs(1.5 - phi[i])*np.pi, 'dawn difference' 
                    Tdawn = T[i]
                
                if deltaphi >= np.abs(2.5 - (phi[i]-2*(pmaxi-2))):
                    #print np.abs(2.5 - phi[i])*np.pi, 'dusk difference'

                    Tdusk = T[i]
            
        return Tmax, Tdawn, Tdusk
        
    

    """Tmax/ Tmin Plot and analytic approx - for Testing"""

    def phi_max(self,eps):
            'analytic approximation for the angle at which the gas reaches a maximum temperature'
                       
            x0 = 2.9685
            x1 = 7.0623
            x2 = 1.1756
            x3 = -0.2958
            x4 = 0.1846
            f = x0*(1.0+x1*eps**(-x2+(x3/(1.0+x4*eps))))**(-1.0)
            return np.arctan(f)


#plt.semilogx(eps, phi_max(eps))

    def Tmax (self,eps):
            
            return np.cos(self.phi_max(eps))**(0.25)


    def Tdusk (self,eps):
            
            y0 = 0.69073
            y1 = 7.5534
            f = (np.pi**2*(1.0+y0/eps)**(-8.0) + y1*eps**(-8.0/7.0))**(-1.0/8.0)
            return f

    def Tdawn (self,eps):
            
            f = (np.pi + (3*np.pi/eps)**(4.0/3.0))**(-0.25)
            return f    

        

'''SAVING THE MOLLVIEW PLOTS!!'''
"""
import matplotlib.pyplot as plt
gas = parcel (e=0.3)

#days, d , F, Fwv, Fpix, Fwvpix = gas.Fleaving(pmax = 6, steps = 100, NSIDE =8)
days, tt, d, Fmap, Fmap_wv, Fmap_pix , Fmap_wvpix, Fleaving, Fleavingwv=gas.Fleaving_map(pmax = 12, 
                                                                                         steps = 200, 
                                                                                        NSIDE =2)

fig, ax = plt.subplots(figsize = (6,12))
hp.visufunc.mollview(Fmap_pix[900,:], xsize = 100, hold = False, sub = (211))
hp.visufunc.mollview(Fmap_pix[1100,:], xsize = 600, sub = (212))
plt.subplots_adjust(hspace=0.1)
fig.savefig("Moll trial")
#fig, ax = plt.subplots(10)
#for i in range(10):
#    hp.projaxes.MollweideAxes.projplot(Fmap_wvpix[1000+i*20,:])
    #hp.mollview(Fmap_wv[1000+i*20,:])
"""



"""TEST PLANET RADIUS STUF"""
"""
# Plot x and y coordinates of the orbit
fig, ax = plt.subplots(3,2, figsize = (14,14))

ax[0,0].set_title("Orbital Position: e = "+str(ec)+", Omega = "+str(Om)+", Inc = "+str(inc)+", w = "+str(w0))
ax[0,0].set_xlabel("East ->")
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
alpha = 180.0-np.array(TA)*57.2958
ax[2,1].plot(t, pyasl.lambertPhaseFunction(alpha+345), ls ='None', marker ='o')
ax[2,1].plot(t, 0.5*(1+np.cos(alpha/57.2958+345/57.2958)),ls ='None', marker ='o')
plt.tight_layout()
fig.savefig('radius (t)')
"""

