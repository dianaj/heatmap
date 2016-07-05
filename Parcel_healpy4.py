# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 17:22:17 2016

@author: diana

This module allows you to create a planet object and calculate the planetary 
phase curve for an arbitrary number of orbits (default is 3).  


Example
-------

.. _NumPy Documentation HOWTO:
   https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

TO DO
-----

- make sure this works for negative wadv 
- especially check SOP, SSP, coordinates, illum, visibility and DE

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


    """This class allows you to create a planet object and assign it appropriate orbital and planetary parameters.
    It uses class functions to calculate the planetary phase curve (emmitted flux) for an arbitrary number of orbits (default is 3). 

    Note
    ----
        Need to create a _getitem_ method. At the moment, i can create an obeject and give it properties, but 
        i can't change one of the properties without defining it all over again. 
        
    Object Attributes --  see __init__ documentation 
    ------------------------------------------------
    (instance variables unique to each instance???)

    Params that will appear in functions 
    ------------------------------------
        pmax
                int; Number of orbital periods we will integrate for. Default  is 3.
        steps
                int; number of steps PER 24 hours. Default is 300.
            
        NSIDE
                power of 2; healpix parameter that determines 
                number of pixels that will subdivide the planet surface ;
                NPIX = 12* NSIDE**2.
                (ex: 192 pixels --> NPIX = 192, NSIDE = 4; 798 pixels --> NPIX = 798, NSIDE = 8)
                Default is 8 but might work fine with 4.
                
        wavelenght
                Values shoyld be entered in micrometers. 
                Used for inc/ emmitted flux calculations at this particular wavelenght.
                Default is 8.
                
    
                
    """

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
                 e = 0.0, argp = 0, 
                 A = 0, ro = 100.0 , cp = 4200.0, H= 5.0, hh = 14121.0, Porb = -1, 
                 wadv = 1.2, tau_rad = 20, epsilon = None):
        
        """The __init__ method allows to set attributes unique to each parcel instance.
        It takes some parameters in the units specified in the docstring. Some are 
        converted to SI units, some are calculated from a few parameters. It has some default
        values that dont have any particular meaning. 
        
        
        Note
        ----
        Some parameters are not used in calculations. They're just there is case i need them for something in the future. 

        Parameters
        ----------
        
        name (str): can give it a name if you want 

        Teff (float): 
            Temperature of star (K)
        
        Rstar (float): 
            Radius of star (in units of solar radii)
        
        Mstar (float): 
            mass of the star in solar masses
        
        Rplanet (float): 
            Radius of planet in Jupiter Radii
        
        a (float): 
            Semimajor axis in AU 
        
        e (float, 0 to 1):
            eccentricity
        
        argp (float):
            Argument at periastron in degrees - angle betwen periatron and transit in degrees 

        A : Bond albedo (set to 0)
            Model doesn not handle reflected light right now. Setting albedo to a different value 
            will have no effect than to reduce incoming flux by a certain fraction. 
        
        ro : not used 
            mass density of atmospheric gas parcel kg/m^3 
        cp : not used
            heat capacity of gas in J/K 
        hh : not used
            scale height of planet
        H : not used
            H* hh --> thickness of gas parcel in scale heights     
        ch : not used 
            heat capacity of gas column in J/m^2  
        
        Porb (float): orbital period in days 
            will be calculated from Kepler's laws if Porb param set to -1    
            
            
        P : 
            Rotational period of the gas around the planet; calculated by self.Prot(). 
            
  
        wadv : PARAM WE WOULD FIT FOR
            
            multiple of wmax (ex: 2 or 0.5)
            wrot = (2*Pi/P) is chosen to match wmax (orbital angular velocity at periastron );
            wadv is expressed as a multiple of wmax, with ( - ) meaning a rotation in the oposite direction.            
            
        T0 : not used in calculations
            Initial temperature of gas at substellar point at periastron. 
            
            Teff*(1-A)**(0.25)*(self.Rstar/(self.a*(1-self.e)))**(0.5) 

        tau_rad (and epsilon): PARAM WE WOULD FIT FOR
            
            epsilon = tau_rad * wadv

            For eccentric orbit, value for tau_rad should be entered in hours and epsilon
            should be left blank. 

            For a circular orbit, epsilon (efficiency parameter) and wadv should be provided 
            and tau_rad left blank.              

        
        rotationsPerOrbit : int(self.Porb/self.P)+1 
            
            used for giving the default time lenght for DE
        
        rotationsPerDay : int(self.P/self.days) 
            
            used for giving the default precision for DE

        """
        self.name = name    # instance variable unique to each instance


        self.Teff = Teff #temp star
        self.Rstar = Rstar * parcel.rsun #radius star
        self.Mstar = 1.0* parcel.Msun #mass of the star in solar masses
        self.Rplanet = Rplanet * parcel.RJup
        self.a = a * parcel.AU #semimajor axis
        self.e = e #eccentricity
        self.argp = argp #angle betwen periatron and transit in degrees 
        
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
            """Simply prints the planetary characteristics that you assigned to your object,
            as well as the ones that were calculated from the same information; Takes no arguments, 
            returns nothing. """
        
            print("""Name  = {0} (Name of model)
                    Teff = {1} K (Temperature Star)
                    Rstar = {2} R-sun (Radius Star in solar radii)
                    Rplanet = {3} R-sun (Radius Planet in solar radii)
                    a = {4} AU (semimajor axis)
                    e = {5} (eccentricity)
                    argp = {6} #angle betwen periatron and transit in degrees 
                    theta = {7} (latitude of gas parcel - 0 North Pole, pi/2 equator)
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
                                "Nothing", 
                                self.A, self.ro,
                               self.cp, self.hh/1000, self.H/self.hh,self.P/parcel.days, self.Porb/parcel.days, 
                                self.ch, self.wadv*self.P/ (2*np.pi), 
                                self.T0, self.tau_rad/ 3600.0, 
                                self.epsilon, int(self.Porb/self.P), int(self.P/self.days) ))
                                
     
    """ FUNCTIONS FOR DEFINING THE PLANET CONDITIONS, STELLAR FLUX """
       
    def Prot(self):

        """Calculates rotational period :: Prot :: of the planet when given orbital period.
        
        Used only by __init__ . 
        
        Note
        ----
             wrot = 2*Pi/Prot is chosen to match wmax (orbital angular velocity at periastron );
             wadv is expressed as a multiple of wmax, with ( - ) meaning a rotation in the oposite direction. 
             
    
        Parameters 
        ----------
            None. Uses eccentricity and orbital period, from the planet's description. 
        
        Returns
        -------
            
            Prot in seconds.

            """
        Prot = self.Porb*(1-self.e)**(1.5)*(1+self.e)**(-0.5)        
        return Prot
    
    

    
    def Fstar(self, wavelenght = 8.0):

        """Calculates total flux emmitted by the star AND wavelenght dependent flux.
        
        Used to express flux emmitted by the planet as a fraction of stellar flux. 
        
        Note
        ----
             F = sigmaB * Teff**4 * Pi * Rstar**2.
             fwv = BBflux(wv) * Pi * Rstar**2
             
    
        Parameters
        ----------
            wv (float)
                Wavelenght we are interested in micrometers.

            
        Returns
        -------
            
                F (float)
                    Total flux
                
                Fwv (float)
                    Flux emmitted at the specified wavelenght.

            """
        wv = wavelenght * 10**(-6) #wavelenght was in micrometers, it's now in meters
        F = self.sigmaB*self.Teff**4*np.pi*self.Rstar**2
        Fwv = (2*self.h*self.c**2/wv**5)*(1/(np.e**((self.h*self.c)/(wv*self.K*self.Teff))-1))*np.pi*self.Rstar**2 
        return F, Fwv
        
        
    def Finc(self, pmax = 3.0, steps = 300.0):
                
        """Total (all wavelenghts) flux incident on substellar point as 
        a function of time (position in the orbit).
        
        Used in self.DE for '"incoming energy" value as (Finc/ Finc_max)*Fweight;
        Fweight comes from self.illum
        

        Note
        ----
             sigmaB*Teff**4*(Rstar/r(t))**2.
             
    
        Parameters
        ----------
            pmax (int)
                Number of orbital periods we will integrate for.
            steps (int)
                number of steps PER 24 hours.

        Calls
        -------
            
            self.radius(pmax, steps) to get orbital separation.
            
        Returns
        -------
            
                Finc (1D array)
                    Theoretical total flux incident on substellar point at each moment on time;
                    lenght pmax * int(Porb/ 24hrs)*steps. 

            """  

        
        return self.sigmaB*self.Teff**4*(self.Rstar/np.array(self.radius(pmax,steps)[1]))**2
        
    def Finc_hemi(self, pmax = 3.0, steps = 300.0):
        
        """Flux incident on the lighted hemisphere of planet. Used for normalizing things mostly.

        Note
        ----
             Analytic integral of incoming flux over hemisphere.
             
    
        Parameters
        ----------
            pmax (int)
                Number of orbital periods we will integrate for.
            steps (int)
                number of steps PER 24 hours.

        Calls
        -------
            
            self.Finc(pmax, steps) to get flux incident on the substellar point.
            
        Returns
        -------
            
                Finc_hemi (1D array)
                    Theoretical total flux incident on planet at each moment on time;
                    lenght pmax * int(Porb/ 24hrs)*steps. 

            """           

        # Einc = Finc* integral(phi: 0->2pi, theta: 0->pi/2) rp^2*cos(theta)sin(theta) dtheta dphi
        # ---- integral overplanet coordinates
        return self.Finc(pmax, steps)*np.pi*self.Rplanet**2
            
    
    
    """FUNCTIONS FOR DEFINING THE COORDINATE SYSTEM AND MOVEMENT AS A FUNCTION OF TIME """
    """::: they all use module pyasl from PyAstronomy ::: """
    
    def radius(self, pmax=3.0, steps = 300):

        """Calculates orbital separation (between planet and its star) as a function of time.
        Used in calculating incident flux.

        Note
        ----
             Could be combined with the other functions in this section, it uses the same stuff.
             Might want to review the Omega and w arguments that i'm passing to 
             pyasl.KeplerEllipse(self.a, self.Porb, e=self.e, Omega=180., i=90.0, w=self.argp) 
             
    
        Parameters
        ----------
            pmax (int)
                Number of orbital periods we will integrate for.
            steps (int)
                number of steps PER 24 hours.

        Uses
        -------
            
            pyasl from PyAstronomy to calculate true anomaly
            
        Returns
        -------
            
                t (1D array)
                    Time array in seconds, of lenght pmax * int(Porb/ 24hrs)*steps. 
                    
                radius (1D array)
                    Same lenght as t. Orbital separation radius array (as a function of time) 

            """   
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

        """Calculates phase angle (alpha) and illuminated fraction of the planet (f) as a function of time.
        Alpha is used in self.visibility; f is used only for testing.

        Note
        ----
             Could be combined with self.radius, it uses the same stuff.
             Alpha is a bit weird.
             
    
        Parameters
        ----------
            pmax (int)
                Number of orbital periods we will integrate for.
            steps (int)
                number of steps PER 24 hours.

        Uses
        -------
            
            pyasl from PyAstronomy to calculate true anomaly
            
        Returns
        -------
            
                t (1D array)
                    Time array in seconds, of lenght pmax * int(Porb/ 24hrs)*steps. 
                    
                ang_vel (1D array)
                    Same lenght as t. Angular velocity at each step in the orbit (as a function of time) 

            """        

        
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

        """Calculates phase angle (alpha) and illuminated fraction of the planet (f) as a function of time.
        Alpha is used in self.visibility; f is used only for testing.

        Note
        ----
             Could be combined with self.radius, it uses the same stuff.
             Alpha is a bit weird.
             
    
        Parameters
        ----------
            pmax (int)
                Number of orbital periods we will integrate for.
            steps (int)
                number of steps PER 24 hours.

        Uses
        -------
            
            pyasl from PyAstronomy to calculate true anomaly
            
        Returns
        -------
            
                t (1D array)
                    Time array in seconds, of lenght pmax * int(Porb/ 24hrs)*steps. 
                    
                alpha (1D array)
                    Same lenght as t. Phase angle in degrees. 
                    alpha = 90 - argp - TrueAnomaly ::THIS IS A BIT WEIRD 
                    
                f (1D array)
                    Same lenght as t. Illuminated fraction of planet: f = 1/2(1-cos(alpha)).

    
            """
        
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

        """Calculates coordinates of substellar point wrt to the location of the 
        substellar point at periastron (theta = Pi/2, phi =0). 
        
        Also calculates coordinates of subobserver location wrt subobserver location at periastron, 
        which is only used for testing. 
        
        Note
        ----
             Could be combined with self.radius, it uses the same stuff.
             Location of SOP at periastron should be related to argument at periastron. 
             My calculation is a bit iffy but seems to work. 
    
        Parameters
        ----------
            pmax (int)
                Number of orbital periods we will integrate for.
            steps (int)
                number of steps PER 24 hours.

        Uses
        -------
            
            pyasl from PyAstronomy to calculate true anomaly
            
        Returns
        -------
            
                t (1D array)
                    Time array in seconds, of lenght pmax * int(Porb/ 24hrs)*steps. 
                    
                zt (1D array)
                    Same lenght as t. Cumulative orbital angular displacement. 
                    
                SSP (1D array)
                    Same lenght as t. SSP = ((zt mod (2 Pi/ Rotation)) mod (2 Pi/ orbit));
                    Gives coordinate of substellar point relative 
                    to the substellar point at periastron (located at theta = Pi/2, phi = 0).
                    Used to rotate the coordinate array array as a function of time .
                
                SOP (1D array)
                    Coordinates of sub-observer point mod (2Pi/Rotation)
    
            """
        
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

        """Creates coordinate matrix wrt substellar point at each point in time. 
        Creates initial temperature array. 
        Calculates weight function to that will be applied to stellar flux to obtain 
        incident flux a each location on the planet, at each point in time.

        Note
        ----
                In places where initial temperature is 0, we replace T = 0 with T = 0.01 to avoid overflows.
                At t = 0, the planet is at periastron and the substellar point 
            is located at theta = Pi	/2, phi = 0
    
        Parameters
        ----------
            pmax
                int; Number of orbital periods we will integrate for.
            steps
                int; number of steps PER 24 hours.
            
            NSIDE
                power of 2; healpix parameter that determines 
                number of pixels that will subdivide the planet surface ;
                NPIX = 12* NSIDE**2.
                (ex: 192 pixels --> NPIX = 192, NSIDE = 4; 798 pixels --> NPIX = 798, NSIDE = 8)
                
        Calls
        -------
            
            self.SSP(pmax, steps) to get:
                
                t
                    1D time array of lenght pmax * int(Porb/ 24hrs)*steps; 
                    in seconds
                
                SSP
                    1D array, same lenght as t. Gives coordinate of substellar point relative 
                    to the substellar point at periastron (located at theta = Pi/2, phi = 0).
                    Used to rotate the coordinate array array as a function of time . 
            
        Returns
        -------
            
                t
                    1D time array of lenght pmax * int(Porb/ 24hrs)*steps; 
                    in seconds
                
                d 
                    3D position and temperature array;
                    
                    shape = (len(time), NPIX, 3)
                    
                    d[:,:,0] = thetas (latitude -- 0 to Pi, equator at Pi/2)
                    *remains constant as a function of time 
                    
                    d[:,:,1] = phis (longitude -- 0 to 2Pi); phi(t) = phi(0)+SSP(t)
                    
                    
                    d[:,:,2] = starting temperature array 
                    
                Fweight
                    2-D array that represents the weight applied to the stellar flux 
                    at each location on the planet to obtain incident flux at each moment in time.
                    --- shape is (lenght time, NPIX)

    
            """
        
        
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
        
        """Calculates the visibility of each gas parcel on the planet, i.e. how much flux is recieved by an 
        observer from each location on the planet as a function of time.

        Note
        ----
            weight = ((np.cos(coordsSOP)+np.abs(np.cos(coordsSOP)))/2.0)*np.sin(thetas)
            the phase angle alpha is a bit of a mystery
    
        Parameters
        ----------
            pmax (int)
                Number of orbital periods we will integrate for.
            steps (int)
                number of steps PER 24 hours.
            
            NSIDE (power of 2)
                healpix parameter that determines 
                number of pixels that will subdivide the planet surface ;
                NPIX = 12* NSIDE**2.
                (ex: 192 pixels --> NPIX = 192, NSIDE = 4; 798 pixels --> NPIX = 798, NSIDE = 8)
            
            d (3D array, optional)
                position and temperature array; is provided as an argument by self.Fobs; 
                if not provided, self.illum(pmax, steps) will be called and the array d will be calculated 
                    
                    shape = (len(time), NPIX, 3)
                    
                    d[:,:,0] = thetas (latitude -- 0 to Pi, equator at Pi/2)
                    *remains constant as a function of time 
                    
                    d[:,:,1] = phis (longitude -- 0 to 2Pi); phi(t) = phi(0)+SSP(t)
                    
                    d[:,:,2] = starting temperature array
                    
                    d[458,234,2] = position wrt to substellar point of parcel # 234 at timestep 458
                    
            TEST (bool)
                Usually False     
        
        Calls
        -------
            
            self.f_ratio(pmax, steps) to get:
                
                t
                    1D time array of lenght pmax * int(Porb/ 24hrs)*steps; 
                    in seconds
                
                alpha
                    1D array, same lenght as t. The phase angle. Gives angle between substellar point 
                    and subobserver location. 
                    
            
        Returns
        -------
            
            if TEST = False
                t      
                    time array in seconds
                weight 
                    visibility array to be applied to flux array in coordinates wrt SSP
                
            if TEST = True
                t (1D array)      
                    time array in seconds
                
                d (3D array)     
                    coordinates wrt SSP and temperature aray
                
                coordsSOP (2D array)
                    coordinates relative to the suborserver point
                
                weight (2D array )
                    visibility array to be applied to flux array in coordinates wrt SSP

    
            """        
       
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
            
            
            """DE that calculates temperature of each gas parcel as a 
            function of time . Relies on self.illum(pmax, steps, NSIDE) to 
            pass it a time array and coordinates. 

            Note
            ----
            Solves the DE by repeatedly adding dT to previous T value. 
            Might want to change this to a more sophisticated 
            differential equation solver.
    
            Parameters
            ----------
            pmax
                int; Number of orbital periods we will integrate for.
            steps
                int; number of steps PER 24 hours.
            
            NSIDE
                power of 2; healpix parameter that determines 
                number of pixels that will subdivide the planet surface ;
                NPIX = 12* NSIDE**2.
                (ex: 192 pixels --> NSIDE = 4; 798 pixels --> NSIDE = 8)
                
            Calls
            -------
            
            self.illum (pmax, steps, NSIDE) to get:
                t
                    1D time array of lenght pmax * int(Porb/ 24hrs)*steps; 
                    in seconds
                
                d 
                    3D position and temperature array;
                    shape = (len(time), NPIX, 3)
                    
                    3 refers to the 3 columns : 
                    d[:,:,0] = thetas - latitude -- 0 to Pi	 
                    
                    d[:,:,1] = phis - longitude -- 0 to Pi	
                    (at t = 0, the planet is at periastron and the substellar point 
                    is located at theta = Pi	/2, phi = 0)
                    
                    d[:,:,2] = starting temperature array              
            
            
            Returns
            -------
            t
                unchanged
                
            d 
                only change is to replace the starting temperature values
                with values calculated by the DE
                                        
    
    
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
        """Calculates outgoing planetary flux (Total and wavelenght dependant)
        from the temperature values coming from the DE. 

            Note
            ----
            Has an overflow problem sometimes. 
        
    
            Parameters
            ----------
            pmax
                int; Number of orbital periods we will integrate for.
            steps
                int; number of steps PER 24 hours.
            
            NSIDE
                power of 2; healpix parameter that determines 
                number of pixels that will subdivide the planet surface ;
                NPIX = 12* NSIDE**2.
                (ex: 192 pixels --> NSIDE = 4; 798 pixels --> NSIDE = 8)
                
            wavelenght
                in micrometers; wavelenght to calculate the flux at. 
                
            Calls
            -------
            
            self.DE (pmax, steps, NSIDE) to get:
                t
                    1D time array of lenght pmax * int(Porb/ 24hrs)*steps; 
                    in seconds
                
                d 
                    3D position and temperature array;
                    shape = (len(time), NPIX, 3)
                    
                    3 refers to the 3 columns : 
                    d[:,:,0] = thetas - latitude -- 0 to Pi	 
                    
                    d[:,:,1] = phis - longitude -- 0 to Pi	
                    (at t = 0, the planet is at periastron and the substellar point 
                    is located at theta = Pi	/2, phi = 0)
                    
                    d[:,:,2] = surface temperature array
                    
            self.shuffle(d, Fmap_wv, pmax, steps, NSIDE) to get 
                
                Fmap_wvpix
                    2D flux array rearraged so the pixels are drawn at the right spots. 
                    see shuffle()
                    
            
            
            Returns
            -------

            t
                unchanged
                
            d 
                unchanged
                
            Fmap_wv 
                2D array[time, NPIX flux values]. outgoing flux map 
                
            Fmap_wvpix 
                2D array[time, NPIX flux values]. outgoing flux map rearraged 
                for drawing. see shuffle() 
                
            Fleavingwv
                1D array, contains  flux (wavelenght dependant) integrated over planet surface
                at each moment in time. 
                
            Ftotal
                1D array, contains flux (all wavelenghts) integrated over planet surface
                at each moment in time. 
   
            """
        
        pmaxi = pmax * self.rotationsPerOrbit
        stepsi = steps * self.rotationsPerDay
        
        wv = wavelenght*10**(-6)        
        
        t, d = self.DE(pmax, steps, NSIDE)
        
        
        'Sometimes this has an overflow problem'
        
        
        a = (2*self.h*self.c**2/wv**5)
        
        b = (self.h*self.c*10**6)/(wavelenght*self.K*self.T0)

        Fwv = a* 1/(np.expm1(b/np.array(d[:,:,2].copy())))
      
        
        "Get the flux"
        dA = hp.nside2pixarea(NSIDE)*self.Rplanet**2
  
        Fmap_wv = (Fwv.copy() *dA)#/Fwvstar
        
        Ftotal_ = (self.sigmaB * (d[:,:,2]*self.T0)**4)*dA
        
        crap, Fmap_wvpix = self.shuffle(d, Fmap_wv, pmax, steps, NSIDE)

        
        Fleavingwv = np.zeros(int(stepsi * pmaxi))
        Ftotal = np.zeros(int(stepsi * pmaxi))
        for i in range(int(stepsi * pmaxi)):
            
            Fleavingwv[i] = np.sum(Fmap_wv[i,:])
            Ftotal[i] = np.sum(Ftotal_[i,:])

        return t, d, Fmap_wv, Fmap_wvpix, Fleavingwv, Ftotal

    
    def Fobs(self, pmax = 3.0, steps = 300.0, NSIDE = 8, wavelenght = 8.0, PRINT = False):
        """ Calculates outgoing planetary flux as seen by an observer (wavelenght dependant).
        

            Note
            ----
            THIS IS THE FUNCTION THAT WILL GIVE YOU THE LIGHT CURVE. 
        
    
            Parameters
            ----------
            pmax
                int; Number of orbital periods we will integrate for.
            steps
                int; number of steps PER 24 hours.
            
            NSIDE
                power of 2; healpix parameter that determines 
                number of pixels that will subdivide the planet surface ;
                NPIX = 12* NSIDE**2.
                (ex: 192 pixels --> NSIDE = 4; 798 pixels --> NSIDE = 8)
                
            wavelenght
                in micrometers; wavelenght to calculate the flux at. 
                
            PRINT (bool):
                Option to print the results to a text file. 
                
            Calls
            -------
            
            self.Fleaving(pmax, steps, NSIDE, wavelenght) to get :

                Fmap_wv 
                    2D array[time, NPIX flux values]. outgoing flux map 
                    
                Fmap_wvpix 
                    2D array[time, NPIX flux values]. outgoing flux map rearraged 
                    for drawing. see shuffle() 

                
            self.visibility(pmax, steps, NSIDE, d) to get:
                t
                    time array in seconds
                    
                d
                    3D coordinate and temperature array
                    
                weight
                    2D array; [time, position (angle)] ; visibility function


            self.shuffle(d, weight, pmax, steps, NSIDE) to get:
                
                weightpix
                    2D array; [time, position(pixel number)] ; rearranged visibility 
                    function for drawing with hp.mollview(). 

   
            Returns
            -------
            
            t, d, Fmapwvobs, weight, weightpix, Fwv

            t
                unchanged
                
            d 
                unchanged
                
            Fmapwvobs 
                2D array[time, NPIX flux values]; outgoing flux map; 
                rearranged for drawing on a hp.mollview map. 
                *you want this one
                
            weight 
                2D array; [time, position (angle)] ; visibility function
                
            weightpix
                2D array; [time, position(pixel number)] ; rearranged visibility 
                function for drawing with hp.mollview().  
                
            Fwv
                1D array, contains observed flux (wavelenght dependant) integrated over planet surface
                at each moment in time.
                *and this one 
   
            """
        
        tic = time.time()
        print ("Starting Fobs")
        
        'call the functions we need'
        t, d, Fmap_wv, Fmap_wvpix,Fleavingwv, Ftotal = self.Fleaving(pmax, steps, NSIDE, wavelenght)
        
        t, weight = self.visibility(pmax, steps, NSIDE, d)

        crap, weightpix = self.shuffle(d, weight, pmax, steps, NSIDE)
        
        
        'start doing what this function does'

        pmaxi = pmax * self.rotationsPerOrbit
        stepsi = steps * self.rotationsPerDay
        
        Nmin = int((pmaxi)*stepsi)
        
        
        #UNCOMMENT IF YOU WANNA CHECK THAT THE SHUFFLING ISNT DESTROYING ANYTHING
        Fmapwvobs = Fmap_wvpix*weightpix
        #Fmapwvobs_check = Fmap_wv*weight
        Fwv = np.empty(Nmin)
        #Fwv_check = np.empty(Nmin)
        for i in range(Nmin):

            Fwv[i] = np.sum(Fmapwvobs[i,:])
            #Fwv_check[i] = np.sum(Fmapwvobs_check[i,:])
            
        #print "DO we get the same flux before and after shuffling? " 
        #print (Fwv == Fwv_check)

        if PRINT == True:
            
            fluxwv = np.array(Fwv).reshape(-1,1)
            np.savetxt('observedflux_planet'+str(steps)+'_steps_per_period_'+str(pmaxi)+'periods_'
                       +str(NSIDE)+'_NSIDE.out', 
                       zip( t, fluxwv), header = "time(planetdays),time(s), outgoing flux, outgoing flux per wv")
        toc = time.time()
        print ('Done with Fobs')
        print ('Time Fobs took: ' + str(toc-tic) + 'seconds')
        
        if PRINT == False:
            return  t, d, Fmapwvobs, weight, weightpix, Fwv
    

    
                    
    def shuffle (self, d = None, quantity = None, pmax = 3.0, steps = 300, NSIDE = 8):
        """ This function will take an array for a quantity as well as it's coordinates
        and rearrange it so it will correspond to healpy pixel number.
    

        Note
        ----
        Normally you would provide a quantity defined on different surface patches 
        of the planet and the matching coordinate array, and this function will rearrange it to
        correspong to pixel number.
        
        For testing purposes, if these are not provided, shuffle() with get the d array from the DE
        and shuffle the temperature values. 

        Parameters
        ----------
        d
            3D numpy array. see DE()
        quantity
            2D array [time, valuethat needs rearranging] 
        pmax
            int; Number of orbital periods we will integrate for.
        steps
            int; number of steps PER 24 hours.
        
        NSIDE
            power of 2; healpix parameter that determines 
            number of pixels that will subdivide the planet surface ;
            NPIX = 12* NSIDE**2.
            (ex: 192 pixels --> NSIDE = 4; 798 pixels --> NSIDE = 8)
            
        
        Calls
        -------
        
        self.DE(pmax, steps, NSIDE) if d and quantity is not provided.


   
        Returns
        -------
   
        d 
            unchanged
            
        quantity 
            2D array[time, values for pixel]; 
            quantity rearranged for drawing on a hp.mollview map. 

        """

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
            
    
    def findT (self, pmax=3.0, steps = 300.0):
        """ Finds numeric approximation of Max/ Min temperature on the planet.
    

        Note
        ----
        Used for testing. Supposed to compare to the analytic approximations in 
        the functions phi-max, Tmax, Tdusk, Tdawn, to
        check that the DE is working well. Or to check that the analytic approx. 
        is working well. 
        
        Only works for circular orbits.
        
        THIS WORKS WITH THE OLD VERSION. HAVE TO UPDATE COORDINATES BEFORE 
        TRYING TO USE IT AGAIN. 

        Parameters
        ----------

        pmax
            int; Number of orbital periods we will integrate for.
        steps
            int; number of steps PER 24 hours.

        Calls
        -------
        
        self.DE(pmax, steps, NSIDE), the 0 eccentricity branch.


   
        Returns
        -------
   
        Tmax (float)
            Maximum temperature on the planet in T/T0
            
        Tdawn 
            Dawn temperature on the planet in T/T0
        
        Tdusk
            Dusk temperature on the planet in T/T0

        """
        
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
        


    def phi_max(self,eps):
        """ Finds analytic approximation for location of  Max temperature on the planet
        for a circular orbit.
    

        Note
        ----
        Used for testing. Only works for circular orbits.
        

        Parameters
        ----------

        eps (float) 
            efficiency parameter

        Returns
        -------
   
        phi (rads)
            Longitude at which gas reaches maximum temperature on planet.


        """
                       
        x0 = 2.9685
        x1 = 7.0623
        x2 = 1.1756
        x3 = -0.2958
        x4 = 0.1846
        f = x0*(1.0+x1*eps**(-x2+(x3/(1.0+x4*eps))))**(-1.0)
        return np.arctan(f)



    def Tmax (self,eps):
        """ Finds analytic approximation for Max temperature on the planet
        for a circular orbit.
    

        Note
        ----
        Used for testing. Only works for circular orbits.
        

        Parameters
        ----------

        eps (float) 
            efficiency parameter

        Returns
        -------
   
        Tmax
            Theoretical max temperature on planet in T/T0.

        """
            
        return np.cos(self.phi_max(eps))**(0.25)


    def Tdusk (self,eps):
        """ Finds analytic approximation for temperature at dusk on the planet
        for a circular orbit.
    

        Note
        ----
        Used for testing. Only works for circular orbits.
        

        Parameters
        ----------

        eps (float) 
            efficiency parameter

        Returns
        -------
   
        Tdusk
            Theoretical temperature at dusk on planet (T/T0.)

        """
            
        y0 = 0.69073
        y1 = 7.5534
        f = (np.pi**2*(1.0+y0/eps)**(-8.0) + y1*eps**(-8.0/7.0))**(-1.0/8.0)
        return f

    def Tdawn (self,eps):
        """ Finds analytic approximation for temperature at dawn on the planet
        for a circular orbit.
    

        Note
        ----
        Used for testing. Only works for circular orbits.
        

        Parameters
        ----------

        eps (float) 
            efficiency parameter

        Returns
        -------
   
        Tdawn
            Theoretical temperature at dawn on planet (T/T0.)

        """
            
        f = (np.pi + (3*np.pi/eps)**(4.0/3.0))**(-0.25)
        return f    

        
"""NOT FINISHED


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
    
"""