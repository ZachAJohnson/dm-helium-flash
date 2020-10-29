import numpy as np
import sys
sys.path.append('../')
import star_props.helium_flash
from star_props.helium_flash import HeliumDeflagration as deflag
from star_props.helium_flash import HeliumFlash 
from star_props.helium_flash import Star
from star_props.helium_flash import DarkMatter 

#Helpful units
Rsol2cm = 6.957e10
Msol2g = 1.989e33
gev2g = 1.7826619e-24
GN = 6.674e-8
rsol=6.957e10 #solar radius in cm
c = 2.9979e10


star_folder='/home/zajohns/dm-helium-flash/data/0.790_solar_mass/'

star=Star(star_folder)
rho_core_data=[rho_array[-1] for rho_array in star.rhodata]
rho_sample = np.geomspace(min(rho_core_data),max(rho_core_data),num=10)
index_sample=[np.argmin(abs(rho-rho_core_data)) for rho in rho_sample]
t_sample=star.agedata[index_sample]

# Use largest(oldest) star core size to set impact parameter values
star.set_age(t_sample[-1])


b_array    = star.now_rcore*np.geomspace(1e-2,20,num=10)
m_array    = np.geomspace(1e14,5e22,num=20) # in grams
t_array    = t_sample
sigma_array= np.geomspace(1,1e9, num=20)
v_array    = np.geomspace(1e-5,5e-3,num=7)
theta0=0
FullScan=np.array([[[[[ [b,m,t,sigma,v]  for b in b_array] for m in m_array] for t in t_array] \
            for sigma in sigma_array] for v in v_array])
FullScan=FullScan.reshape(int(FullScan.size/5),5)
scan_size = len(FullScan)

#Take input which breaks FullScan into PartScan
n_pieces, piece = [int(arg) for arg in sys.argv[1:]]
piece-=1
PartScan = np.array_split(FullScan,n_pieces)[piece]


trajectory_output=[];old_age=0
for info in PartScan: 
	b, m, age, sigma, v = info
	print("# b   m   age  sigma  v")
	print(info)
	#Setting age of star and it's size then
	if age!=old_age:
		star.set_age(age)

	dm = DarkMatter(m,sigma,massunit='g')
	flash = HeliumFlash(dm,star)

	#Now initial DM position/velocity
	r0 = 10*star.now_rmax
	gamma = np.arcsin(b/r0)
	vr0 = -np.cos(gamma)*c*v
	phi0=  np.sin(gamma)*c*v/r0
	#Calculate trajectory, outputs hit_core, ignited
	traj=flash.trajectory(gamma, [r0,vr0,theta0,phi0],calc_flash=True) 
	trajectory_output.append( [info, traj])
	old_age=age
	print(traj,'\n')

savefile='/home/zajohns/dm-helium-flash/data/0.790_solar_mass/flash_scan-{0}.npy'.format(piece)
np.save(savefile,trajectory_output)

