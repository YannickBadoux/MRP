from initial_conditions import critical_velocity, generate_initial_conditions, add_encounter
from run_simulation import run_simulation
from analyse_result import bound

from amuse.lab import units, constants
from amuse.io import write_set_to_file

import numpy as np
from argparse import ArgumentParser
import os
import time
from tqdm.auto import tqdm

def max_impact_parameter(v, a, C=4, e=0):
    '''Calculates the maximum impact parameter for a given velocity, 
    semi-major axis and eccentricity.'''
    D = 0.6*(1+e)
    return (C/v + D) * a

def v20_from_vinf(v_inf, a, M):
    '''Calculates the velocity of the field star at 20 a_pl'''
    return np.sqrt(v_inf**2 + 4*constants.G*M/(20*a))

if __name__ == '__main__':
    parser = ArgumentParser(description='Run a Monte Carlo simulation of an encounter between an equal mass binary and a field star.')
    parser.add_argument('--velocity', type=float, default=1, help='Velocity of the field star in units of the critical velocity')
    parser.add_argument('--n_sim', type=int, default=1000, help='Number of simulations per impact parameter')
    parser.add_argument('--output', type=str, default='automatic_cs_output', help='Output directory for the results')
    # parser.add_argument('--interested_state', type=int, default=1, help='State of the system to be interested in')
    args = parser.parse_args()

    v = args.velocity
    n_sim = args.n_sim
    save_path = args.output
    # interested_state = args.interested_state

    tqdm.write(f"Running simulation with velocity {v} and {n_sim} simulations per impact parameter annulus")

    start_time = time.time()

    #create output directory
    os.makedirs(save_path, exist_ok=True) #to save the results array
    os.makedirs(f'{save_path}/output_velocity_{v}', exist_ok=True) #to save simulation results

    #reproduce the results of Hut & Bahcall 1983. Equal masses and circular orbit
    m1 = 1 | units.Msun
    m2 = 1 | units.Msun
    m3 = 1 | units.Msun
    a = 1 | units.AU

    #calculate the critical velocity
    v_crit = critical_velocity(m1=m1, m2=m2, m3=m3, a=a)

    #calculate the initial velocity of the field star
    v_inf = v * v_crit
    v20 = v20_from_vinf(v_inf, a, m3)

    #initialize results array
    dtype = [('v_inf', 'f8'), ('v20','f8'), ('b','f8'), ('phi', 'f8'), ('theta', 'f8'), ('psi', 'f8'), ('f', 'f8'), ('end_time', 'f8'), ('state', 'u1'), ('index', 'u4')]
    results = np.zeros((0,), dtype=dtype)

    #iterate over impact parameters
    index = 0
    bmax = 1 | units.AU
    bmin = 0 | units.AU
    while True:
        #reset the counter for the interested state
        # interested_state_counter = 0
        exchange_counter = 0
        ionization_counter = 0
        temp_results = np.zeros((n_sim,), dtype=dtype)

        tqdm.write(f"Impact parameter range: {bmin.in_(units.AU)} to {bmax.in_(units.AU)}")

        #pick n_sim combinations of angles and impact parameters
        impact_parameters = np.random.uniform(bmin.value_in(units.AU), bmax.value_in(units.AU)**2, n_sim)
        impact_parameters = np.sqrt(impact_parameters) | units.AU
        phis = np.random.uniform(0, 2*np.pi, n_sim)
        thetas = np.arccos(np.random.uniform(0, 1, n_sim))
        psis = np.random.uniform(0, 2*np.pi, n_sim)
        true_anomalies = np.random.uniform(0, 2*np.pi, n_sim)

        i_sim = 0
        for b, phi, theta, psi, f in tqdm(zip(impact_parameters, phis, thetas, psis, true_anomalies), total=n_sim):
            #generate host star, planet and moon
            bodies = generate_initial_conditions(m1, m2, a, f_pl=f, 
                                                n_moons=0, save_path=None)

            #add field star to encounter
            bodies = add_encounter(bodies, m3, b, v20, phi, theta, psi)

            #rename particles
            bodies.name = ['star1', 'star2', 'star3']

            #run the simulation until the energy error is small enough
            timestep_parameter = 0.03 #default value for Huayno
            i=0
            while i < 4:
                evolved, _, end_time, stop_code = run_simulation(bodies, integrator='huayno',
                                                       timestep_parameter=timestep_parameter)
                #decrease the timestep if the energy error is too high
                if stop_code == 2:
                    timestep_parameter *= 0.5
                    i+=1
                    tqdm.write(f'Rerunning simulation with {timestep_parameter}')
                else:
                    break

            #check the state of the evolved system
            star1, star2, star3 = evolved
            if bound(star1,star2):
                state = 0 #flyby
            elif bound(star1,star3) or bound(star2,star3):
                state = 1 #exchange
                exchange_counter +=1
            elif not bound(star1,star2) and not bound(star1,star3) and not bound(star2,star3):
                state = 2 #ionization
                ionization_counter +=1
            else:
                state = 3 #other, should not happen, but could be resonance if velocity is below the critical velocity

            #save evolved system to file and save initial conditions and state to array
            # write_set_to_file(evolved, save_path+f'/output_velocity_{v}/{index}_{state}.amuse', 'amuse', overwrite_file=True)
            temp_results[i_sim] = (v, v20.value_in(units.kms), b.value_in(units.AU), phi, theta, psi, f, end_time.value_in(units.yr), state, index)
            i_sim += 1
            index += 1
        
        results = np.concatenate((results, temp_results))

        tqdm.write(f"Finished, {ionization_counter} ionizations and {exchange_counter} exchanges found.")
        if ionization_counter == 0 and exchange_counter == 0:
            tqdm.write(f'Only flybys occured. Stopping simulation.')
            break

        #increase the impact parameter range so the total area is constant
        bmin = bmax
        bmax = np.sqrt(bmax.value_in(units.AU)**2 + 1) | units.AU
        
    
    #save the results array
    np.save(save_path+f'/results_velocity_{v}.npy', results)

    tqdm.write(f"Simulation took {time.time()-start_time} seconds")
    tqdm.write(f"Results saved to {save_path}")