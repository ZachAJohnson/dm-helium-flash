import subprocess 

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy
import scipy.interpolate as interp
import scipy.integrate as integrate
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
import pandas as pd

import star_props.integrate as zint
from star_props import HFDIR        

GN = 6.674e-8


class Star():
    def __init__(self, mass=1.3):#,stellar_mass, stellar_composition= (0,1,0),):# extend to more general compositions than x,y,z?
        #gen_profiles()
        'do nothing'
        self.rsol=6.957e10 #solar radius in cm
        self.msol=1.98844e33 #solar mass in grams
        #self.rcore=0.021*self.rsol
        #self.rmax=4.811*self.rsol
        self.folder='../data/1.3_solar_mass/'
        self.gen_profiles(self.folder)
        self.age = 0 #Changed for every new trajectory
        self.age_index = 0 #Changed for every new trajectory
      
    def gen_profiles(self, folder): #Generates profile of some specific star from MESA data
        files = {'T':'temperature_K.npy', 'age': 'stellar_age_yrs.npy', 'r': 'radius_cm.npy',
         's3':'nuclear_triple_alpha_erg_s.npy','cs': 'sound_speed_cm_s.npy',
         'eta':'degeneracy.npy','rho':'density_g_cm3.npy', 'cell':'cell_width_cm.npy',
         'Y':'helium_fraction_Y.npy', 'X': 'hydrogen_fraction_X.npy'}

        datas = {}
        for key, val in files.items():
            datas[key] = np.load(self.folder + val, allow_pickle=True, encoding='latin1') 
        
        rmaxs=[prof[0] for prof in datas['r']]
        s3max=[max(prof) for prof in datas['s3'] ]
        eta = [max(prof) for prof in datas['eta']]
        flashindex=int(np.where(max(s3max)==s3max  )[0])
        flashtime=datas['age'][flashindex]

        start_eta = 4
        start = np.where(np.array(eta)>start_eta)[0][2]
        end = flashindex
        starttime=datas['age'][start]
        endtime=datas['age'][end]
        print("Valid times between where eta > {0}, at time {1} Gyr, and where helium flash occurs at {2} Gyr".format(start_eta,starttime/1e9,endtime/1e9))
        #print(datas['cell'][0][-3:-1])
        #print(datas['r'][0][-3:-1])
        #print(np.array([sum(datas['cell'][]) -max(datas['r'][1]) ]))
        #print((datas['r'][0][-i-2]-datas['r'][0][-i-1])/datas['cell'][0][-i-2])
        #print("max error is:", max([abs((datas['r'][0][-i-2]-datas['r'][0][-i-1])/datas['cell'][0][-i-2]-1) for i in range(len(datas['cell'][0])-1) ]))
        #print(datas['cell'][0][-10:-1],'\n',datas['r'][0][-10:-1])


        

        self.frmax = interp1d(datas['age'],rmaxs )

        for key in files.keys():  #shortens time scale to only what is relevant
            datas[key] = datas[key][start : end]
        points=[]
        #yrarray=[]
        for ti in range(end-start):
            t=datas['age'][ti]
            for ri in range(len(datas['r'][ti])): 
                r=datas['r'][ti][ri]
                #yrarray.append(t)
                points.append([t,r])

        self.rdata=datas['r']
        self.agedata=datas['age']
        self.rhodata=datas['rho']
        self.etadata = datas['eta']
        self.Tdata = datas['T']
        self.Xdata = datas['X']
        self.Ydata = datas['Y']


        #Now making mass as a function of radius, density
        self.mdata=[]
        for it in range(len(self.rdata)):
            rev_rho = np.flip(datas['rho'][it]) #reversing because mesa gives backwards for god knows why
            rev_r = np.flip(datas['r'][it])
            rev_m=[4/3*np.pi*rev_rho[0]*rev_r[0]**3]
            for ir in range(1,len(rev_rho)):
                rev_m.append( rev_m[ir-1] + 4/3*np.pi*rev_rho[ir]*(rev_r[ir]**3 - rev_r[ir-1]**3)  )

            rev_m[0]=0 #Crucial! Gotta avoid that 1/r^2 divergence....
            self.mdata.append(np.flip(rev_m))

        #self.yrarray=yrarray
        self.rarray=np.concatenate(datas['r'])
        #self.marray=np.concatenate(self.mdata)
        self.points=points
        self.Tarray=np.concatenate(datas['T'])
        self.rhoarray=np.concatenate(datas['rho'])        

        self.set_age(starttime)

        print("Loaded all arrays")

        # rlin=np.geomspace(1e3,max(rarray),num=1000)
        # yrlin=np.geomspace(starttime,endtime,num=1000)

        # self.LOGR, self.LOGYR = np.meshgrid(np.log10(rlin),np.log10(yrlin))
        # INSTAR = self.LOGR<np.log10(rmax(self.LOGYR))

        # T = griddata((np.log10(rarray),np.log10(yrarray)),Tarray,(self.LOGR,self.
        #     LOGYR) ,method='nearest',rescale=True)
        # RHO = griddata((np.log10(rarray),np.log10(yrarray)),rhoarray,(self.LOGR,self.LOGYR) ,method='nearest',rescale=True)
        
        # self.T = np.where(INSTAR,T,0) #Kills values outside star
        # self.RHO = np.where(INSTAR,RHO,0) #Kills values outside star

        return 
    
    def set_age(self,age):
        self.age= age
        self.age_index=(abs(np.array(self.agedata)-age)).argmin()
        self.now_rmax = self.frmax(self.age)
        #self.set_trigger(deflag)

        try:
            min_eta=4
            core_index = np.where(self.etadata[self.age_index] > min_eta)[0][0]
            self.now_rcore=self.rdata[self.age_index][core_index]
        except Exception as e:
            print("No valid place where eta>{0}, no core".format(min_eta))
            print(e)

        return None

    def M(self, r): #stupid stepwise M       
        nearest_index=np.array(abs(self.rdata[self.age_index]-r)).argmin() 
        m_nearest=self.mdata[self.age_index][nearest_index]
        return m_nearest

    def rho(self, r):
        nearest_index=np.array(abs(self.rdata[self.age_index]-r)).argmin() 
        rho_nearest=self.rhodata[self.age_index][nearest_index]
        return rho_nearest

    def temp(self, r):
        nearest_index = np.array(abs(self.rdata[self.age_index]-r)).argmin() 
        T_nearest = self.Tdata[self.age_index][nearest_index]
        return T_nearest

    def comp(self,r):
        nearest_index = np.array(abs(self.rdata[self.age_index]-r)).argmin() 
        X_nearest = self.Xdata[self.age_index][nearest_index]
        Y_nearest = self.Ydata[self.age_index][nearest_index]
        return X_nearest, Y_nearest

    def set_trigger(self):
        #def min_trig(self,yt):
        #r, vr, theta, omega = yt
        #nearest_index = np.array(abs(self.rdata[self.age_index]-r)).argmin()     
        #rho = self.rho(r)  
        #temp = self.temp(r)    
        Temps = np.arange(1.0e8,1.5e9,3.0e8)  
        Temps = np.array([7.0e8])      

        self.mintrigdata=[]
        for ir in range(len(self.rdata[self.age_index])) :            
            if self.rdata[self.age_index][ir] > self.now_rcore:
                self.mintrigdata.append(1e100)
            else:
                rho = self.rhodata[self.age_index][ir]
                tcold = self.Tdata[self.age_index][ir]
                trig = lambda T: HeliumDeflagration(rho,tcold,tcrit = T).trigger()
                trigs = [trig(T) for T in Temps]
                self.mintrigdata.append(min(trigs))

    def mintrig(self,r):
        nearest_index = np.array(abs(self.rdata[self.age_index]-r)).argmin() 
        mintrig_nearest = self.mintrigdata[nearest_index]
        return mintrig_nearest
        

    def update(self):
        'to implement for scanning purposes'

class DarkMatter():
      
    def __init__(self, mass, sigman, profile='NFW',massunit='GeV'):
        self.gev2g = 1.783e-24    
        self.sigman = sigman        
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

    def sigmaA(self, A):
        return self.sigman*0.1*A**4
        

class HeliumFlash():
    c = 2.9979e10
    Rsol2cm = 6.957e10
    Msol2g = 1.989e33

    matplotlib.rcParams['figure.figsize'] = [15.0, 10.0]
    matplotlib.rcParams.update({'font.size': 22})
    np.set_printoptions(precision=2,suppress=False)

    def __init__(self, dm, star):
        'nothing'
        self.dm = dm
        self.star = star
        self.incore=False
        self.instar=False
        self.av_A = lambda r: pow(sum(np.array(self.star.comp(r)) * np.array([1,4])**4), 0.25)

    def update( self, mass = None): 
        'implement later'

    def f(self,t, yt): #yt= [r,vr,theta,omega]          
        r, vr, theta, omega = yt          
        Fr, Ftheta = -self.star.rho(r)*self.dm.sigmaA(self.av_A(r))*np.sqrt(vr**2 + r**2 * omega**2)*np.array([vr, r*omega])

        if self.star.rho(r)<1e3:
            self.incore=False
        elif self.incore == False: 
            #print("HIT CORE!!")
            self.incore=True
        
        if Fr ==0:
            self.instar=False
        elif self.instar==False: 
            #print("Hit Star!! , density-", self.star.rho(r))
            self.instar=True

        f1 = vr
        f2 = Fr/self.dm.mg-GN*self.star.M(r) / r**2 + r*omega**2
        f3 = omega
        f4 = Ftheta/(self.dm.mg*r)-2*vr*omega/r
        
        return np.array([f1,f2,f3,f4])


    def trajectory(self, gamma, y0 ): #gamma is incoming DM angle in DM coordinates
        #np.set_printoptions(precision=3,suppress=True)        
        r0,vr0,theta0,omega0 = y0
        v0=np.sqrt(vr0**2 + r0**2*omega0**2)
        y = y0
        size=len(y)
        t = 0
   
        dt=pow(10,2)
        t_data=[t]
        y_data=[y]
        not_stopped=True
        hit_core=False
        n=0
        while y[0]<r0+1 and not_stopped and n<1e5:     
            k1 = dt*self.f(t,      y)
            k2 = dt*self.f(t+dt/2, y + k1/2 )
            k3 = dt*self.f(t+dt/2, y + k2/2 )
            k4 = dt*self.f(t+dt,   y + k3 )

            dy = (k1 + 2*(k1+k3) + k4)/6.0

            v=np.sqrt(y[1]**2 + y[0]**2*y[3]**2)
            print("r is {0}, and v is {1}".format(y[0]/self.star.now_rcore,v/self.c),end='\r')

            too_big= [abs(dy[i])>0.001*abs(y[i]) and abs(y[i]>1e-7) for i in range(0,size)]
            too_small= [abs(dy[i])<1e-3*abs(y[i]) for i in range(size)]

            if any(too_big):
                dt=dt/2
                #print("Decreasing stepsize- ", too_big, n, y[0])
            elif all(too_small):
                dt=2*dt
                #print("Increasing stepsize- ", too_small)
            else:                                  
                y += dy
                t += dt
                t_data.append(t)                                                        
                y_data.append(list(y))
                n+=1
            if v**2/v0**2 < 0.001: #Dm roughly stopped if KE/KE_0 <1%
                not_stopped=False

            #rho_new = self.star.rho(y[0])
            if y[0]<self.star.now_rcore:# and rho_new != rho_old:
                #hit_core=True
                #return True
                #self.t_data=np.array(t_data)
                #self.y_data=np.array(y_data)
    
                
                min_trig = self.star.mintrig(y[0])
                print("v = {0}, KE/KE_0 = {1}".format(v/self.c, v**2/v0**2 ), "dEdx = ",self.dEdx(y)," and minimum needed is: ", min_trig,'\n',end='\r')
                if self.dEdx(y) > min_trig:
                    return hit_core, True
            rho_old = self.star.rho(y[0])

        self.t_data=np.array(t_data)
        self.y_data=np.array(y_data)
        return hit_core, False

    def plot_traj(self,t_data, y_data, r1,r2):
        plot_ydata = y_data[:,0]*np.cos(y_data[::,2])
        plot_xdata = y_data[:,0]*np.sin(y_data[::,2])
        
        x1range=(-r1,r1)
        y1range=(-r1,r1)
        x2range=(-r2,r2)
        y2range=(-r2,r2)
        
        fig, (ax1,ax2) = plt.subplots(2)
        ax1.plot(plot_xdata,plot_ydata,'ro')
        ax2.plot(plot_xdata,plot_ydata,'ro')

        outer_star=ax1.add_patch(plt.Circle([0,0],radius=self.star.now_rmax))        
        outer_star.set_alpha(0.1)
        outer_star=ax2.add_patch(plt.Circle([0,0],radius=self.star.now_rmax))        
        outer_star.set_alpha(0.1)
        
        core=ax2.add_patch(plt.Circle([0,0],radius=self.star.now_rcore))
        core.set_alpha(0.3)
        ax1.axis('equal')
        ax2.axis('equal')
        
        ax1.set_xlim(x1range[0],x1range[1])
        ax1.set_ylim(y1range[0],y1range[1])
        ax2.set_xlim(x2range)
        ax2.set_ylim(y2range)
        plt.show()  
        


    def dEdx(self,yt):
        r, vr, theta, omega = yt          
        return self.star.rho(r)*self.dm.sigmaA(self.av_A(r))*(vr**2 + omega**2*r**2)

    #def explode(self,yt):
    #    rho = self.star.rho(r)
    #    return 

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
        tmp=tinterp(rho,tcold)        
        return tmp   
    
    def get_properties(self,temp):
        out= subprocess.check_output([HFDIR+"/star_props/timmes/eosfxt.so", str(float(temp)),str(self.rho)])
        out= map(float, out.strip().split())
        pep, eta, xne, ener = out

        out= subprocess.check_output([HFDIR+"/star_props/opac.o",str(float(temp)), str(self.rho), str(pep), str(eta), str(xne)])
        opac, orad, ocond, conductivity=  map(float,out.strip().split())
        #print("Testing if conductivity, {0} is related to opacity, {1}".format(conductivity, (16/3.0)*self.sigmab*temp**3/ (self.rho*opac) ))
        return opac, orad, ocond, conductivity, ener, eta


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
        #print("\nu'= ",np.array([u_der]),", k = ", np.array([conductivity]), ", Sdot int = ", np.array([prof_integral]))
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
            lambda ksi:  (ksi**2)*(self.get_properties(self.linear_tprofile(ksi))[4]-self.get_properties(self.tcold)[4]) , 0, 1))
        tau = E/Edot
        return tau
    
    def tau_nuclear(self):
        #T=self.tcrit
        #props= self.get_properties(T)
        T = lambda ksi: self.linear_tprofile(ksi)
        Eperg = lambda T: self.get_properties(T)[4]
        EpergCold = self.get_properties(self.tcold)[4]
        Edot, error = (4*np.pi)*(self.lambda_t**3)*self.rho*np.array(integrate.quad(
            lambda ksi: (ksi**2)*self.s3alpha(T(ksi)),0,1))
        
        print("At Tcrit: ",np.array([self.tcrit]),'\n E/s: ', np.array([Edot]))      
        print("Using M_trigger = ", np.array([self.mtrig]) )  
        E, error = (4*np.pi)*self.rho*(self.lambda_t**3)*np.array(integrate.quad(
            lambda ksi:  (ksi**2)*(Eperg(T(ksi))-EpergCold) , 0, 1))
        print("E: ", E)
        print("Average E/g/s: ", np.array([E/self.mtrig/self.tau_diffusion()]))
        print("lambda = ", self.lambda_t, " cm")
        tau = E/Edot
        return tau

    def E_trigger(self):
        T = lambda ksi: self.linear_tprofile(ksi)
        Eperg = lambda T: self.get_properties(T)[4]
        EpergCold = self.get_properties(self.tcold)[4]
        E, error = (4*np.pi)*self.rho*(self.lambda_t**3)*np.array(integrate.quad(
            lambda ksi:  (ksi**2)*(Eperg(T(ksi))-EpergCold) , 0, 1))
        self.Etrig = E
        return E
    
    def trigger(self):
        print("for rho={0}, t_cold = {1}, t_crit = {2} ".format(self.rho, self.tcold, self.tcrit))
        print("Etrig = ",np.array([self.E_trigger()]), ", lambda = ", np.array([self.lambda_t]),", M = ",np.array([self.mtrig]))
        return self.E_trigger()/self.lambda_t
    
