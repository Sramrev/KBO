part = 0

import rebound
import pandas as pd
import numpy as np
from astropy import constants as c
from astropy import units as u

orbital_elements = [
 {'a': 5.196315723518543,
  'e': 0.04925135211374762,
  'inc': 0.02275332972970286,
  'Omega': 1.7540383003717264,
  'omega': 4.7850414455522134,
  'f': 2.7234113181814785},
 {'a': 9.550594418717841,
  'e': 0.05483484925416918,
  'inc': 0.043404524059630224,
  'Omega': 1.982964238068382,
  'omega': 5.924705414335117,
  'f': 0.21417049814055383},
 {'a': 19.18754079397335,
  'e': 0.04734412922165621,
  'inc': 0.013477515455083794,
  'Omega': 1.2914266480630847,
  'omega': 1.6936180377756733,
  'f': 2.841023323085272},
 {'a': 30.07166723909742,
  'e': 0.008705460262078923,
  'inc': 0.03089645254168813,
  'Omega': 2.3000452024492146,
  'omega': 4.765373307711288,
  'f': 4.6900275931363815},
 {'a': 39.48809423216381,
  'e': 0.2489902520038915,
  'inc': 0.2991600055654475,
  'Omega': 1.9251177520188532,
  'omega': 1.9857950718906283,
  'f': 0.6209330744283078}]


# Get L7 Model 
#L7_res = pd.read_csv('upload/L7Model/resonant.dat',skiprows=28,sep='\s+',header=0)
data = pd.read_csv('upload/modified_L7_32.csv')
#data = L7_res[L7_res['n1']/L7_res['n2']==3/2] # Asteroid:Neptune resonances



def segments(part):
    min = len(data)//80
    rem = len(data)-80*min

    if part <= rem:
        index = (part-1)*(min+1)
    else:
        index = rem*(min+1) + (part-rem-1)*(min)

    return index

data = data[segments(part):segments(part+1)]


n_out = 10000  # Number of timesteps 
p_len = len(data)
phi_list = np.zeros((p_len,n_out)) # Array with shape [particle_index][timestep]

for it in range(p_len):

    sim = rebound.Simulation()
    sim.G = c.G.value
    sim.units = ('yr', 'AU', 'Msun')
    sim.integrator = "WHFAST"
    sim.dt = 2

    #############################################
    #### Select planets ####
    # For all planets use 'JSUN_', to remove a planet change letter to _
    used_planets = 'JSUN_'
    mask = [char!='_' for char in used_planets]
    #############################################

    sim.add(m=1+5.98e-6)
    masses = [0.00095465,0.00028558,0.00004344,0.00005149,0]
    for i in np.array([0,1,2,3,4])[mask]:
        mass = masses[i]
        elem = orbital_elements[i]
        sim.add(m=mass, a=elem["a"], e=elem["e"], inc=elem["inc"], omega=elem["omega"], Omega=elem["Omega"], f=elem["f"])
    sim.move_to_com()

    # add 3:2 object
    sim.add(m=0,a=list(data['a'])[it],e=list(data['e'])[it],Omega=np.deg2rad(list(data['node'])[it]),inc=np.deg2rad(list(data['i'])[it]),omega=np.deg2rad(list(data['peri'])[it]),M=np.deg2rad(list(data['M'])[it]))
    # add planet 9
    #sim.add(m=1.8621634e-5,a=382.4,inc=np.deg2rad(15.6),e=0.2,pomega=np.deg2rad(247.6),Omega=np.deg2rad(97.5))
    particles = sim.particles

    # Integrate the simulation, saving relevant orbital elements

    for i,t in enumerate(np.linspace(0, 1e9,n_out)):
        sim.integrate(t)

        ln              = particles[4].l
        l_particle      = particles[5].l
        pomega_particle = particles[5].pomega

        phi_list[it][i]     = np.rad2deg(3*l_particle-2*ln-pomega_particle)%360

np.save(f'lists/phi_list_32_{part}',phi_list)
