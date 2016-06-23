# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 17:22:17 2016

@author: diana
"""

# Problems : 
# 

# 


import numpy as np
#import matplotlib.pyplot as plt
#from scipy import integrate
from PyAstronomy import pyasl
import healpy as hp



class parcel:
    'A gas parcel on a planet. Default are values kind of made up '
    
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
                 A = 0, ro = 100.0 , cp = 4200.0, H= 5.0, hh = 14121.0, P = 1.0, Porb = -1, wadv = 1.2, epsilon = 6):
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
            self.Porb =  Porb*parcel.days #orbital period in seconds formula --
            
        self.P =  P* parcel.days #self.Prot() #self.Prot() P* parcel.days#period of rotation planet in seconds 
        
        self.ch=self.ro*self.cp*self.H # J/m^2 heat capacity of gas column
        self.wadv = (2*np.pi/self.P)*wadv #think this is a param we need to fit for. In units of how much faster 
                                        #the parcel moves compared to the period of the planet. 
                                        #+ moves with rotation of planet, - against 
        self.T0 = Teff*(1-A)**(0.25)*(self.Rstar/(self.a*(1-self.e)))**(0.5) # * np.sin(theta))**0.25 when the parcel 
                                                                        #lives away from the equator
        #self.tau_rad = self.ch/(sigmaB*self.T0**3)
        #self.epsilon = self.wadv*self.tau_rad
        self.epsilon = epsilon
        self.tau_rad = self.epsilon/self.wadv
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
                    tau_rad = {18} s (radiative time scale)
                    epsilon = {19} dimensionless circulation efficeincy param
                    
                    **** all of these are converted to SI units but they need 
                    to be entered in the units mentioned here ***
                    """.format (self.name, self.Teff, self.Rstar/parcel.rsun, self.Rplanet/parcel.rsun, 
                                self.a/parcel.AU, self.e, self.argp, 
                                self.theta, 
                                self.A, self.ro,
                               self.cp, self.hh/1000, self.H/self.hh,self.P/parcel.days, self.Porb/parcel.days, 
                                self.ch, self.wadv/(2*np.pi/self.P), 
                                self.T0, self.tau_rad, 
                                self.epsilon ))
            
    
    def radius(self, pmax=2.0, steps = 300):
        
        
        ke = pyasl.KeplerEllipse(self.a, self.Porb, e=self.e, Omega=90., i=90.0, w=90.0)

        # Get a time axis
        tmax = self.P*pmax
        Nmin = int((pmax)*steps)
        t = np.linspace(0,tmax,num=Nmin)

        radius = ke.radius(t)
        return t, radius
    
    def ang_vel(self, pmax = 2.0, steps = 300):
        
        ke = pyasl.KeplerEllipse(self.a, self.Porb, e=self.e, Omega=90., i=90.0, w=90.0)
        tmax = self.P*pmax
        Nmin = int((pmax)*steps)
        t = np.linspace(0,tmax,num=Nmin)
        vel = ke.xyzVel(t)
        absvel = np.sqrt(vel[:,0]**2+vel[:,1]**2+vel[:,2]**2)
        ang_vel = absvel/ self.radius(pmax, steps)[1]
        #print max(ang_vel)*self.a*(1-self.e)
        
        return t, ang_vel
    
    def SSP(self, pmax =2, steps =300):
        t, ang_vel = self.ang_vel(pmax, steps)
        t, alpha, f = self.f_ratio(pmax, steps)
        'z0= orbital position at t = 0'
        #z0=self.argp*np.pi/180
        z0=0
        'orbital phase'
        zt = np.empty(int(pmax*steps))
        deltat = (self.P*1.0)/(steps*1.0)
        zt[0]=z0
        for i in range(1,int(pmax*steps)):  
            zt[i]= zt[i-1]+ang_vel[i-1]*deltat
        
        SSP =(zt-((self.wadv)*t -((self.wadv*t)/(2*np.pi)).astype(int)*2*np.pi)-
        (t/self.Porb).astype(int)*2*np.pi)
        SOP = (((-alpha[0]+180)*np.pi/180)-
                ((self.wadv)*t -((self.wadv*t)/(2*np.pi)).astype(int)*2*np.pi))
        #SOP = SSP + ((-alpha+180)*np.pi/180)
        return t, zt, SSP, SOP
    
 
           
    def Prot(self, pmax=2.0, steps =300):
        k = self.Mstar* parcel.G
        vmax = ((k*(1+self.e))/(self.a*(1-self.e)))**(0.5)
        #print vmax
        Prot = (2*np.pi*self.Rplanet)/(vmax)
        return Prot
    
    
    def f_ratio(self, pmax = 2.0, steps = 300):
        'returns the illuminated fraction of the planet'
        ke = pyasl.KeplerEllipse(self.a, self.Porb, e=self.e, Omega=90., i=90.0, w=90.0)
        tmax = self.P*pmax
        Nmin = int((pmax)*steps)
        t = np.linspace(0,tmax,num=Nmin)
        pos, TA = ke.xyzPos(t,getTA = True)

        alpha = 180 - np.array(TA)*57.2958 + self.argp
        #plt.plot(t,alpha)
        f = pyasl.lambertPhaseFunction(alpha)
        return t, alpha, f
    
    def Fstar(self, wavelenght = 8.0):
        wv = wavelenght * 10**(-6) #wavelenght was in micrometers, it's now in meters
        F = parcel.sigmaB*self.Teff**4*np.pi*self.Rstar**2
        Fwv = (2*self.h*self.c**2/wv**5)*(1/(np.e**((self.h*self.c)/(wv*self.K*self.Teff))-1))*np.pi*self.Rstar**2 
        return F, Fwv
        
        
    def Finc(self, pmax = 2.0, steps = 300.0):
        "flux incident on a point of planet - the substellar point "
        
        return parcel.sigmaB*self.Teff**4*(self.Rstar/np.array(self.radius(pmax,steps)[1]))**2
        
    
    def illum(self, pmax =2, steps =300, NSIDE = 8, days = None, d = None, TEST = False) :  
        
        if days == None and d == None:
            days, d = self.DE(pmax, steps, NSIDE)
        
        else: 
            days  = days
            d = d
        
        t, zt, SSP, SOP = self.SSP(pmax, steps)
        

                
        phis = d[:,:,2]
        thetas = d[:,:,1]
        'i think these are phis now'
        coordsSSP = (phis)
        
        
        
        """!!!!!!"""
        
        """OH shit! coords SSP are already sorted!!!"""          
        Fweight = ((np.cos(coordsSSP)+
                   np.abs(np.cos(coordsSSP)))/2.0 * np.sin(thetas))
                   
                    
        F = (self.Finc(pmax, steps).reshape(-1,1)* Fweight)
             
        Fnorm = F/(parcel.sigmaB*self.Teff**4*(self.Rstar/(self.a*(1-self.e)))**2)
        
        Fweightpix = np.empty(Fweight.shape)
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
            Fweightpix[i,:]= Fweight[i,hp.ang2pix(NSIDE, thetas[i,:],coordsSSP[i,:])] #phis[i,:])]
        
        if TEST == True:
            return t, zt, SSP, SOP, thetas, phis, c_thetas1, c_phis1, c_thetas7, c_phis7,F, Fnorm, Fweight, Fweightpix
        
        if TEST == False:
            return Fweight, Fweightpix
            
        else:
            print "Do you want the contours or not? If yes, set TEST = True. Else, set TEST = False"
            
    def visibility(self, pmax =2, steps =300, NSIDE = 8, days = None, d = None, TEST = False):
        
        if days == None and d == None:
            days, d = self.DE(pmax, steps, NSIDE)
        
        else: 
            days  = days
            d = d
                
        t, zt, SSP, SOP = self.SSP(pmax, steps)
        phis = d[:,:,2]#-(angles.reshape(-1,1)) # location of gas parcel on planet relative to where it started
        thetas = d[:,:,1]
        'i think these are phis now'
        coordsSSP = (phis) #-
        
        #coordsSOP = phis+(zt%(2*np.pi)).reshape(-1,1)
        
        coordsSOP = phis+(zt).reshape(-1,1)
        #coordsSOP = phis - (SOP - SSP).reshape(-1,1)    
        """coords SOP are already sorted too!!!"""
        weight = ((np.cos(coordsSOP)+np.abs(np.cos(coordsSOP)))/2)*np.sin(thetas)
                
        weightpix = np.empty(weight.shape)
        #contours of equal illumination theta
        
        c_thetas1 = []
        c_phis1 = []
        for i in range(int(pmax*steps)):
            c_t = thetas[i,np.where(np.abs(weightpix[i,:]-1)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            c_p = coordsSSP[i,np.where(np.abs(weightpix[i,:]-1)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            c_thetas1.append(c_t)
            c_phis1.append(c_p)
        
        c_thetas7 = []
        c_phis7 = []
        
        for j in range(int(pmax*steps)):
            c_t = thetas[j,np.where(np.abs(weightpix[j,:]-0.7)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            c_p = coordsSSP[j,np.where(np.abs(weightpix[j,:]-0.7)< 0.02)]# for c in (1, 0.7, 0.4, 0.1))]
            'should we put these angles in substellar point coordinates or what? '
            c_thetas7.append(c_t)
            c_phis7.append(c_p)
            weightpix[j]= weight[j,hp.ang2pix(NSIDE, thetas[j,:], coordsSSP[j,:])] #phis[j,:])]

        if TEST == True:
            return t, thetas, phis, c_thetas1, c_phis1, c_thetas7, c_phis7, coordsSSP, coordsSOP, weight, weightpix
        
        if TEST == False:
            return t, zt, coordsSSP, coordsSOP, weight, weightpix
            
        else:
            print "Do you want the contours or not? If yes, set TEST = True. Else, set TEST = False"
        
        
        
    def Finc_hemi(self, pmax = 2.0, steps = 300.0):
        "total energy (flux) incident on a hemisphere of planet "
        "depends on temperature of star!!!"
        # Einc = Finc* integral(phi: 0->2pi, theta: 0->pi/2) rp^2*cos(theta)sin(theta) dtheta dphi
        # ---- integral overplanet coordinates
        return self.Finc(pmax, steps)*np.pi*self.Rplanet**2
    
    def Fleaving(self, pmax = 2.0, steps = 300.0, NSIDE = 8, wavelenght = 8.0, TEST = False):
        "The substellar point is moving in this case"
        wv = wavelenght*10**(-6) #wavelenght was in micrometers, it's now in meters
        
        days, d = self.DE(pmax, steps, NSIDE)
        #these 2 will match with the coordinates of parcel
        "F = (parcel.sigmaB*(np.array(d[:,:,3])*self.T0)**4)"
        
        "flux planet/ flux star here"
        d[:,:,3][d[:,:,3] == 0] = 0.0001 #replace the zeroes at the beginning so they dont overflow
        Fwv = ((2*self.h*self.c**2/wavelenght**5)*10**(30)*
        (1/(np.e**((self.h*self.c*10**6)/(wavelenght*self.K*np.array(d[:,:,3])*self.T0))-1)))
        
        "Fpix  = (parcel.sigmaB*(np.array(d[:,:,4])*self.T0)**4)"
        d[:,:,4][d[:,:,4] == 0] = 0.0001 #replace the zeroes at the beginning so they dont overflow
        Fwvpix = ((2*self.h*self.c**2/wavelenght**5)*10**(30)*
        (1/(np.e**((self.h*self.c*10**6)/(wavelenght*self.K*np.array(d[:,:,4])*self.T0))-1)))
        #weight[j,hp.ang2pix(NSIDE, thetas[j,:], coordsSSP[j,:])]
        #print F.shape
        #print Fwv.shape
        
        "Get the flux"
        dA = hp.nside2pixarea(NSIDE)*self.Rplanet**2
        
        'Fmap = (F*dA)#/Fstar'
        Fmap_wv = (Fwv*dA)#/Fwvstar
        'Fmap_pix = (Fpix*dA)#/Fstar'
        Fmap_wvpix = (Fwvpix*dA)#/Fwvstar
        
        'Fleaving = np.zeros(int(steps * pmax))'
        Fleavingwv = np.zeros(int(steps * pmax))
        for i in range(int(steps * pmax)):
            'Fleaving[i] = np.sum(Fmap[i,:])'
            Fleavingwv[i] = np.sum(Fmap_wv[i,:])
            Fleavingwvpix = np.sum(Fmap_wvpix[i,:])
            #if np.abs(Fleavingwvpix-Fleavingwv[i]) > 1000:
                #print (Fleavingwvpix-Fleavingwv[i])
                #print ("WOW, WHATS WRONG!" + str(i))
            
        if TEST==True:
            return days, d ,Fwv, Fwvpix, Fmap_wv, Fmap_wvpix, Fleavingwv     
        if TEST == False:
            return days, d, Fmap_wv, Fmap_wvpix, Fleavingwv
        
        else:
            print "you want all the crap or not?? If yes, set Test = TRUE. Else, set Test = False"
        
        

    
    def Fobs(self, pmax = 2.0, steps = 300.0, NSIDE = 8, wavelenght = 8.0, PRINT = False):
        "surface flux (leaving) from one hemisphere of the planet "
        "the flux leaving the planet will be different"
        import time
        tic = time.time()
        print ("Starting Fobs")
        Nmin = int((pmax)*steps)
        days, d, Fmap_wv, Fmap_wvpix, Fleavingwv = self.Fleaving(pmax, steps, NSIDE, wavelenght)
        t, zt, coordsSSP, coordsSOP, weight, weightpix = self.visibility(pmax, steps, NSIDE, days, d, TEST = False)
        #t, alpha, f = self.f_ratio(pmax, steps)
        thetas = d[:,:,1]
        
        
        'used to just use phis = d[:,:,2], same thetas'

        #"Fmapobs = Fmap*weight"
        
        Fmapwvobs = Fmap_wv*weight
        Fmap_wvobspix = Fmap_wvpix*weight

        
        #"Fleaving = np.empty(Nmin)"
        Fwv = np.empty(Nmin)
        for i in range(Nmin):
            
            
            #"Fmap_pix[i,:] = Fmapobs[i,d[i,:,0].astype(int)]"
            #Fmap_wvpix[i,:] = Fmapwvobs[i,hp.ang2pix(NSIDE, thetas[i], coordsSSP[i])]
            'used to have'
            #Fmap_wvpix[i,:] = Fmapwvobs[i,d[i,:,0].astype(int)]
            #Fleaving[i] = np.sum(Fmapobs[i,:]) 
            Fwv[i] = np.sum(Fmap_wvobspix[i,:]) 
            
        
        """
        Fleaving = np.empty(int(steps * pmax))
        phistart = np.empty(int(steps * pmax))
        phiend = np.empty(int(steps * pmax))
        
        phistart[np.where(alpha>=0)] = ((np.pi/2 - f*np.pi)+2*np.pi*days.astype(int))[np.where(alpha>=0)] #array for all t
        phiend[np.where(alpha>=0)] = phistart[np.where(alpha>=0)] + np.pi
        
        phiend[np.where(alpha<0)] = ((3*np.pi/2 + f*np.pi)+2*np.pi*days.astype(int))[np.where(alpha<0)] #array for all t
        phistart[np.where(alpha<0)] = phiend[np.where(alpha<0)] - np.pi
        #Fmapwvobs = Fmap_wv*weight    
        for i in range(int(steps * pmax)):

             
            angles = np.linspace(0, 2*np.pi, num =Nphis)+2*np.pi*days.astype(int)[i]
            chunk = np.where((phistart[i] <= angles)  & (angles <= phiend[i]))
            Fleaving[i] = np.sum(Fmap[i, chunk ])"""
        if PRINT == True:
            tt = self.P*np.array(days).reshape(-1,1)
            #flux = np.array(Fleaving).reshape(-1,1)
            fluxwv = np.array(Fwv).reshape(-1,1)
            np.savetxt('observedflux_planet'+str(steps)+'_steps_per_period_'+str(pmax)+'periods_'
                       +str(NSIDE)+'_NSIDE.out', 
                       zip( tt, fluxwv), header = "time(planetdays),time(s), outgoing flux, outgoing flux per wv")
        toc = time.time()
        print ('Done with Fobs')
        print ('Time Fobs took: ' + str(toc-tic) + 'seconds')
        
        if PRINT == False:
            return  days, d, Fmapwvobs, Fmap_wvobspix, weight, weightpix, Fwv #*weight
    

        
    def DE(self, pmax=2.0, steps = 300.0, NSIDE = 8):
            'pmax is in number of periods you want to integrate for.'
            'The number of steps is per period'
            
            """returns: -time (days) and 
                        -d:  pixel number , angular location of parcel(thata,phi), temperature of gas parcel
                        at each time value
                        """
            #theta = self.theta               
            import time 
            tic = time.time()
            print "Starting DE"
            
            Nmin = int((pmax)*steps)
            tmax = self.P*pmax
            
            #time
            
            t = np.linspace(0,tmax, num=Nmin)
            dayss =(t)*(self.wadv)/(2*np.pi)
            
            
            
            #position
            coords =np.empty(((hp.nside2npix(NSIDE)),2))
            for i in range(hp.nside2npix(NSIDE)):
                coords[i,:]=(np.array(hp.pix2ang(NSIDE, i)))
                
            phis = coords[:,1]
            thetas = coords[:,0]
            pixnum = np.arange(len(phis))
            
            #starting Temperatures
            T = (np.sin(thetas)**0.25)*(np.cos(phis)+np.abs(np.cos(phis)))/2
            Tpix = np.empty(hp.nside2npix(NSIDE))
            #3D array to contain everything    
            c= np.array(zip(pixnum, thetas,phis,T, Tpix)).reshape(-1,hp.nside2npix(NSIDE),5)
            d = np.repeat(c,Nmin, axis = 0)
            
            
            
            if self.e == 0.0:
                
                deltaphi = (2*np.pi/steps)* self.wadv/(2*np.pi/self.P)
                
                for i in range(1,Nmin):#phis.shape[2]
                    phis = d[i-1,:,2]
                    dT =(1.0/self.epsilon *(np.cos(phis)+np.abs(np.cos(phis)))/2
                         *np.sin(thetas) - (d[i-1,:,3])**4 )*deltaphi #calculate change

                    temp = d[i-1,:,3]+ dT 

                    d[i,:,3]= temp #update temperature array
                    "coordinates WRT substellar point!!"
                    d[i,:,2]= (phis+deltaphi)%(2*np.pi) #-((self.wadv)*t[i] - int(self.wadv*t[i])) #update location of gas parcel
                    d[i,:,0]= hp.ang2pix(NSIDE, d[i,:,1], d[i,:,2]) #update pixel number to draw gas parcel at
                    d[i,:,4]= d[i,d[i,:,0].astype(int),3] # shuffle the temperatures to put them at the right pixel number


                toc = time.time()
                print ("Time this took is: " , str(toc-tic), "seconds")
                
                return dayss, d
            
            else:
                t, zt, SSP, SOP = self.SSP(pmax, steps)
                #t_=t/self.tau_rad
                deltat = self.P/steps
                deltat_ = deltat/self.tau_rad
                wrot = (2*np.pi/self.P)* self.wadv/(2*np.pi/self.P)                   
                deltaphi = wrot*deltat 
                    
                    
                F = ((self.a*(1-self.e)/self.radius(pmax, steps)[1])**2)
                phis0 = d[0,:,2]
                for i in range(1,len(t)):
                        phis = d[i-1,:,2]+SSP[0]
                        
                        dT =(( F[i]*(np.cos(phis)+np.abs(np.cos(phis)))/2*np.sin(thetas) - 
                              (d[i-1,:,3])**4 )* (deltat_))
                        
                        d[i,:,3]= d[i-1,:,3]+dT #update temperature array
                        'tricky part. Trying to take into acount how teh SSP moves'
                        d[i,:,2]= ((phis0 )- SSP[i])
                                #deltaphi)) #%(2*np.pi)) #deltaphi)%(2*np.pi)-
                                #(zt-(t/self.Porb).astype(int)*2*np.pi)[i]) #update location of gas parcel
                        d[i,:,0]= hp.ang2pix(NSIDE, d[i,:,1], d[i,:,2]) #update pixel number to draw gas parcel at
                        d[i,:,4]= d[i,d[i,:,0].astype(int),3] # shuffle the temperatures to put them at the right pixel number
                        
                toc = time.time()
                print ("Time DE took: "+ str(toc-tic))
                return dayss, d 
            
            
                        
    def findT (self, pmax=2.0, steps = 300.0):
        
        #tmax = self.P*pmax
        #Nmin = int((pmax)*300)
        #deltat = tmax/Nmin
        
        phi, days, T = self.DE(pmax, steps)
        
        
        deltaphi = 2.0/steps
        Tmax = np.max(T[int(steps*(pmax-2))::])
        
            
        for i in range(int(steps*(pmax-2)), len(phi)):
            
            
            
                if deltaphi >= np.abs(1.5 - (phi[i]-2*(pmax-2))):
                    #print np.abs(1.5 - phi[i])*np.pi, 'dawn difference' 
                    Tdawn = T[i]
                
                if deltaphi >= np.abs(2.5 - (phi[i]-2*(pmax-2))):
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

