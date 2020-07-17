import numpy as np

class rk4():
      def __init__(self, f): #f is array of functions for derivative
            self.f = f 
            

      def initialize_values(self, y0, t0=0 ):
            self.y=y0
            self.size=len(self.y)                          
            self.t=t0

      def integrate(self,tf,nsteps=1000):
            dt = (tf-self.t)/nsteps
            
            n=0
            while self.t <= tf:                      
                  k1 = dt*self.f(self.t,      self.y)
                  k2 = dt*self.f(self.t+dt/2, self.y + k1/2 )
                  k3 = dt*self.f(self.t+dt/2, self.y + k2/2 )
                  k4 = dt*self.f(self.t+dt,   self.y + k3 )
                  
                  dy = (k1 + 2*(k1+k3) + k4)/6.0

                  too_big= [abs(dy[i])>0.1*abs(self.y[i]) and abs(self.y[i])>1e-4 for i in range(self.size)]
                  too_small= [abs(dy[i])<1e-3*abs(self.y[i]) for i in range(self.size)]
                  # if any(too_big):
                  #       dt=dt/2
                  #       print("Decreasing stepsize- ", too_big)
                  # elif all(too_small):
                  #       dt=2*dt
                  #       print("Increasing stepsize- ", too_small)
                  # else:
                  n+=1
                  self.y+= dy
                  self.t += dt


            return self.y

