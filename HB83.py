from initial_conditions import critical_velocity, generate_initial_conditions, add_encounter
from run_simulation import run_simulation
from analyse_result import bound

from amuse.lab import units
from amuse.io import write_set_to_file

import numpy as np
from argparse import ArgumentParser
import os

def max_impact_parameter(v, a, C=4, e=0):
    '''Calculates the maximum impact parameter for a given velocity, 
    semi-major axis and eccentricity.'''
    D = 0.6*(1+e)
    return (C/v + D) * a

if __name__ == '__main__':
    parser = ArgumentParser(description='Run a Monte Carlo simulation of an encounter between an equal mass binary and a field star.')
    parser.add_argument('--velocity', type=float, default=1, help='Velocity of the field star in units of the critical velocity')
    parser.add_argument('--n_b', type=int, default=100, help='Number of impact parameters to simulate')
    parser.add_argument('--sim_per_b', type=int, default=100, help='Number of simulations per impact parameter')
    parser.add_argument('--output', type=str, default='output_HB83', help='Output directory for the results')
    args = parser.parse_args()

    v = args.velocity
    n_b = args.n_b
    n_sim = args.sim_per_b
    save_path = args.output

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

    v_inf = v * v_crit
    b_max = max_impact_parameter(v, a)
    print(f"v={v} v_crit, b_max={b_max.in_(units.AU)}")

    #impact parameters with probability homogeneous in b**2
    impact_parameters = np.sqrt(np.linspace(0, b_max.value_in(units.AU)**2, n_b)) | units.AU

    #initialize results array
    dtype = [('v_inf', 'f8'), ('b','f8'), ('phi', 'f8'), ('theta', 'f8'), ('psi', 'f8'), ('end_time', 'f8'), ('state', 'u1'), ('index', 'u4')]
    results = np.zeros((len(impact_parameters) * n_sim,), dtype=dtype)

    #iterate over impact parameters
    index = 0
    for b in impact_parameters:
        print(f"Starting with: b={b.in_(units.AU)}")
        #pick 1000 combinations of angles
        phis = np.random.uniform(0, 2*np.pi, n_sim)
        thetas = np.arccos(np.random.uniform(0, 1, n_sim))
        psis = np.random.uniform(0, 2*np.pi, n_sim)
        true_anomalies = np.random.uniform(0, 2*np.pi, n_sim)

        for phi, theta, psi, f in zip(phis, thetas, psis, true_anomalies):
            #generate host star, planet and moon
            bodies = generate_initial_conditions(m1, m2, a, f_pl=f, 
                                                n_moons=0, save_path=None)

            #add field star to encounter
            bodies = add_encounter(bodies, m3, b, v_inf, phi, theta, psi)

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
                    print(f'Rerunning simulation with {timestep_parameter}')
                else:
                    break

            #check the state of the evolved system
            star1, star2, star3 = evolved
            if bound(star1,star2):
                state = 0 #flyby
            elif bound(star1,star3) or bound(star2,star3):
                state = 1 #exchange
            elif not bound(star1,star2) and not bound(star1,star3) and not bound(star2,star3):
                state = 2 #ionization
            else:
                state = 3 #other, should not happen, but could be resonance if velocity is below the critical velocity

            #save evolved system to file and save initial conditions and state to array
            write_set_to_file(evolved, save_path+f'/output_velocity_{v}/{index}_{state}.amuse', 'amuse', overwrite_file=True)
            results[index] = (v, b.value_in(units.AU), phi, theta, psi, end_time.value_in(units.yr), state, index)

            index += 1
    
    #save the results array
    np.save(save_path+f'/results_velocity_{v}.npy', results)

