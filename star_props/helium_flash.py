import subprocess 

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy.interpolate as interp
import scipy.integrate as integrate
import pandas as pd

import star_props.integrate as zint
from star_props import HFDIR        

GN = 6.674e-8


class Star():
    def __init__(self):#,stellar_mass, stellar_composition= (0,1,0),):# extend to more general compositions than x,y,z?
        #gen_profiles()
        'do nothing'
        self.rsol=6.957e10 #solar radius in cm
        self.msol=1.98844e33 #solar mass in grams
        self.rcore=0.021*self.rsol
        self.rmax=4.811*self.rsol
      
    # def gen_profiles(self, stellar_mass, stellar_composition):
    #       'interpolate data tables to generate density, temperature, etc. as function of time, radius'
    #       'self.rho = array of time, mass'
    #       rho

    def M(self, r): #stupid stepwise M
        
        if 0<r<self.rcore:
            return 4/3*np.pi*(r/self.rcore)**3*self.rho(r)
        elif self.rmax>r>self.rcore:
            return 0.21*self.msol + 4/3*np.pi*self.rho(r)*(r**3-self.rcore**3)
        else: 
            return self.msol

    def rho(self,r): #stupid stepwise rho

        if r>self.rmax:
            return 0
        elif r>self.rcore:
            return 1e-2
        else:
            return pow(10,4.5)

    def update(self):
        'to implement for scanning purposes'

class DarkMatter():
      
    def __init__(self, mass, sigma, profile='NFW',massunit='GeV'):
        self.gev2g = 1.783e-24    
        self.sigma = sigma
        if massunit== 'GeV':
            self.mgev = mass
            self.mg =  self.mgev*self.gev2g
        elif massunit.lower()=='g' or massunit.lower()=='gram':
            self.mg= mass
            self.mgev =  self.mgev/gev2g
        else:
            print("Error: improper massunit declaration.")          

    def velocity_profile(self):
        'create later'

    def update(self, mass = None, massunit='GeV'):
        'to implement for scanning mass, sigma'
        

class HeliumFlash():
    rsol=6.957e10 #solar radius in cm
    c = 2.9979e10
    matplotlib.rcParams['figure.figsize'] = [15.0, 10.0]
    matplotlib.rcParams.update({'font.size': 22})


    def __init__(self):
        'nothing'
        self.default_m=1e42      
        self.default_sigma=1e0
        self.dm = DarkMatter(self.default_m ,self.default_sigma  , massunit='GeV')
        self.star = Star()

    def update( self, mass = None): 
        'implement later'

    def f(self,t, yt): #yt= [r,vr,theta,omega]          
        r, vr, theta, omega = yt          
        Fr, Ftheta = -self.star.rho(r)*self.dm.sigma*np.sqrt(vr**2 + r**2 * omega**2)*np.array([vr, r*omega])
        
        #if abs(Fr)!=0: print(np.array([Fr, Ftheta])/(GN*self.star.M(r) / r**2)) 
        f1 = vr
        f2 = Fr/self.dm.mg-GN*self.star.M(r) / r**2 + r*omega**2
        f3 = omega
        f4 = Ftheta/(self.dm.mg*r)-2*vr*omega/r
        return np.array([f1,f2,f3,f4])

    def trajectory(self, gamma, y0 ): #gamma is incoming DM angle in DM coordinates
        inter = zint.rk4(self.f)
        r0,vr0,theta0,omega0 = y0
        inter.initialize_values(y0=y0,t0=0)
        
        np.set_printoptions(precision=2)
        dt=pow(10,3)
        t_data=[]
        y_data=[]
        while inter.y[0]<r0+1:# and inter.y[1]<0:
            t_data.append(inter.t)
            newy=list(inter.integrate(inter.t + dt)  )
            y_data.append(newy)
            #print(np.array(newy))

        t_data=np.array(t_data)
        y_data=np.array(y_data)

        return [t_data, y_data]


    def plot_traj(self,data):
        t_data, y_data = data

        plot_ydata = y_data[::,0]*np.cos(y_data[::,2])
        plot_xdata = y_data[::,0]*np.sin(y_data[::,2])

        fig, ax = plt.subplots()
        plt.plot(plot_xdata,plot_ydata)
        print(self.star.rmax,self.rsol*200/self.star.rmax)
        outer_star=ax.add_patch(plt.Circle([0,0],radius=self.star.rmax))        
        outer_star.set_alpha(0.1)
        core=ax.add_patch(plt.Circle([0,0],radius=self.star.rcore))
        core.set_alpha(0.3)
        ax.axis('equal')
        ax.set_xlim(-2*self.star.rmax,2*self.star.rmax)
        ax.set_ylim(-2*self.star.rmax,2*self.star.rmax)
        plt.show()    


class HeliumDeflagration():
    matplotlib.rcParams.update({'font.size': 20})
      
    def __init__(self, rho, tcold, tcrit=None, thermal_width = 'linear', k=2):
        self.tcold=tcold
        self.rho=rho
        self.sigmab = 5.67e-5 #stefan-boltzmann constant
        self.c = 2.99e10 #light in cm/s
        self.k=k
        
        if tcrit==None:            
            self.tcrit=HeliumDeflagration.get_tcrit(self.rho,self.tcold)
        else:
            self.tcrit=tcrit
        
        width_dict={'linear':self.thermal_width_linear ,'simple':self.thermal_width_simple,
                    'timmes':self.thermal_width_timmes }           
        try: 
            width_dict[thermal_width]()
            self.trigger_mass()
        except KeyError as e:
            print(e)
            print("options are ", width_dict.keys())
    
    
    @staticmethod        
    def get_tcrit(rho,tcold):
        tcrit_data=pd.read_csv(HFDIR+"/star_props/Tcrit.dat")
        tcrit_data.columns = tcrit_data.columns.str.replace(' ','')
        tinterp=interp.LinearNDInterpolator(tcrit_data[['rho','Tcold']], tcrit_data['Tcrit'])        
        return tinterp(rho,tcold)        
    
    def get_properties(self,temp):
        out= subprocess.check_output([HFDIR+"/star_props/timmes/eosfxt.so", str(temp),str(self.rho)])
        out= map(float, out.strip().split())
        pep, eta, xne, ener = out

        out= subprocess.check_output([HFDIR+"/star_props/opac.o",str(temp), str(self.rho), str(pep), str(eta), str(xne)])
        opac, orad, ocond, conductivity=  map(float,out.strip().split())
        #print("Testing if conductivity, {0} is related to opacity, {1}".format(conductivity, (16/3.0)*self.sigmab*temp**3/ (self.rho*opac) ))
        return opac, orad, ocond, conductivity, ener


    def s3alpha(self,temp):
        t9=temp/10**9
        rate=5.1e8*self.rho**2*(1/t9)**3*np.exp(-4.4027/t9)
        return rate

    def linear_tprofile(self, ksi, a=5/3., **kwargs): #ksi = r/lambda_t, so distance from center of thermal region to edge
        tmax = a*self.tcrit             
        temp = tmax - ksi*(tmax - self.tcold)
        return temp


    def thermal_width_simple(self): #lambda_t
        opac = self.get_properties( self.tcrit )[0]
        self.lambda_t = np.sqrt((16*self.sigmab/3.0)*self.tcrit**4/( self.s3alpha(self.tcrit)*opac*self.rho**2))
        return 


    def thermal_width_linear(self, tprof =linear_tprofile): #k is one for cylinder, two for sphere
        conductivity = self.get_properties(self.tcold)[3] #Not clear where to evaluate, probably at edge of heated region 
        u_der = abs((tprof(self,0.01) - tprof(self,-0.01)) /0.01)
        prof_integral, error = integrate.quad(lambda ksi:  (ksi**self.k)*self.s3alpha(tprof(self,ksi)) , 0, 1)
        self.lambda_t = np.sqrt(  u_der*conductivity / (self.rho*prof_integral ) )
        return 

    def thermal_width_timmes(self):
        #props = self.get_properties(self.tcrit)
        opac  = self.get_properties(self.tcrit)[0]
        eperg = self.get_properties(self.tcrit)[4]-self.get_properties(self.tcold)[4]
        self.lambda_t = np.sqrt(self.c * eperg/(opac*self.rho*self.s3alpha(self.tcrit) )  )
        return 

    def trigger_mass(self, width_func=thermal_width_linear,**kwargs):    
        self.mtrig = 4/3*np.pi*(self.lambda_t**3)*self.rho       
        return 
    
    def dedxtrigger(self):
        'fill in'
    
    def tau_diffusion(self): #only for a sphere
        T=self.tcrit
        props= self.get_properties(T)
        k = props[3]
        #dski =1 #ksi is r/lambda_t, or the normalized distance from center of deflagration
        #dT =T -tcold
        dksi=0.01 
        dT = self.linear_tprofile(1-dksi)-self.linear_tprofile(1)
        Edot = abs(4*np.pi*self.lambda_t*k*( dT/dksi ))
        E, error = 4*np.pi*self.rho*self.lambda_t**3*np.array(integrate.quad(
            lambda ksi:  (ksi**2)*(self.get_properties(self.linear_tprofile(ksi))[4]-self.get_properties(self.Tcold)[4]) , 0, 1))
        tau = E/Edot
        return tau
    
    def tau_nuclear(self):
        T=self.tcrit
        props= self.get_properties(T)
        Edot = #np.array(integrate.quad(lambda ksi: ))
        
        E, error = 4*np.pi*self.rho*self.lambda_t**3*np.array(integrate.quad(
            lambda ksi:  (ksi**2)*(self.get_properties(self.linear_tprofile(ksi))[4]-self.get_properties(self.Tcold)[4]) , 0, 1))
        tau = E/Edot
        return tau
    
    
