part = 0

import rebound
import pandas as pd
import numpy as np
from astropy import constants as c
from astropy import units as u

orbital_params = [pd.read_csv(f'upload/commands/script_files/planets_parameters/{planet}.csv') for planet in ['Jupiter','Saturn','Uranus','Neptune']] # List of JSUNP, with orbital elements for each epoch
data = list(np.load('upload/commands/script_files/21_MPC.npy',allow_pickle=True)) #KBO data

# Separate the data according to code part
def segments(part):
    min = len(data)//70
    rem = len(data)-70*min

    if part <= rem:
        index = (part-1)*(min+1)
    else:
        index = rem*(min+1) + (part-rem-1)*(min)

    return index
data = data[segments(part):segments(part+1)]


n_out = 10000  # Number of timesteps 
p_len = len(data)
phi_list = np.zeros((p_len,n_out)) # Array with shape [particle_index][timestep]
a_list = np.zeros((p_len,n_out)) # Array with shape [particle_index][timestep]
e_list = np.zeros((p_len,n_out)) # Array with shape [particle_index][timestep]
inc_list = np.zeros((p_len,n_out)) # Array with shape [particle_index][timestep]


for it in range(p_len):

    sim = rebound.Simulation()
    sim.G = c.G.value
    sim.units = ('yr', 'AU', 'Msun')
    sim.integrator = "WHFAST"
    sim.dt = 2

    sim.add(m=1+5.98e-6)
    masses = [0.00095465,0.00028558,0.00004344,0.00005149]
    epoch = data[it]['Epoch'] # Epoch of asteroid data
    orbital_elements = [orbital_params[i][orbital_params[i]['JDTDB']==epoch] for i in range(4)] #List of four dictionaries with orbital parameters

    # Add each planet
    for planet_index in range(4):
        mass = masses[planet_index]
        elem = orbital_elements[planet_index]
        
        # Add each giant planet
        a = float(elem["A"])
        e = float(elem["EC"])
        inc = float(np.deg2rad(elem["IN"]))
        omega=float(np.deg2rad(elem["W"]))
        Omega=float(np.deg2rad(elem["OM"]))
        M = float(np.deg2rad(elem["MA"]))
        sim.add(m=mass, a=a, e=e, inc=inc, omega=omega, Omega=Omega, M=M)
    sim.move_to_com()

    
    # Add asteroid
    a = data[it]['a']
    e = data[it]['e']
    inc = np.deg2rad(data[it]['i'])
    Omega = np.deg2rad(data[it]['Node'])
    omega = np.deg2rad(data[it]['Peri'])
    M = np.deg2rad(data[it]['M'])
    sim.add(a=a,e=e,Omega=Omega,inc=inc,omega=omega,M=M)
    sim.move_to_com()
    particles = sim.particles
    

    # Integrate the simulation, saving relevant orbital elements
    for i,t in enumerate(np.linspace(0, 4e9,n_out)):
        sim.integrate(t)

        ln              = particles[4].l
        l_particle      = particles[5].l
        pomega_particle = particles[5].pomega

        a_list[it][i]          = particles[5].a
        e_list[it][i]          = particles[5].e
        inc_list[it][i]        = np.rad2deg(particles[5].inc)
        phi_list[it][i]     = np.rad2deg(2*l_particle-1*ln-pomega_particle)%360

res = '21'
np.save(f'lists/phi_list/phi_list_{res}_{part}',phi_list)
np.save(f'lists/inc_list/inc_list_{res}_{part}',inc_list)
np.save(f'lists/a_list/a_list_{res}_{part}',a_list)
np.save(f'lists/e_list/e_list_{res}_{part}',e_list)

