from initial_conditions import critical_velocity, generate_initial_conditions, add_encounter
from run_simulation import run_simulation
from analyse_result import bound
import function_file as ff

from amuse.lab import units
from amuse.io import write_set_to_file

import numpy as np
from argparse import ArgumentParser
import os
import time
from tqdm.auto import tqdm

if __name__ == '__main__':
    parser = ArgumentParser(description='Run a Monte Carlo simulation of an encounter between an equal mass binary and a field star.')
    # parser.add_argument('--velocity', type=float, default=1, help='Velocity of the field star in units of the critical velocity')
    parser.add_argument('--a_sp', type=float, default=1, help='Semi-major axis of the planet in AU')
    parser.add_argument('--density', type=int, default=3, help='Density of points per square AU')
    parser.add_argument('--far_away_distance', type=float, default=150, help='Distance at which the field star is considered far away in AU')
    parser.add_argument('--output', type=str, default='automatic_cs_output', help='Output directory for the results')
    args = parser.parse_args()

    np.random.seed(42)

    save_path = args.output
    point_density = args.density | units.AU**-2 #number of simulations per square AU
    far_away_distance = args.far_away_distance | units.AU

    tqdm.write(f"Running simulation with a_sp {args.a_sp} and {point_density} simulation density")

    start_time = time.time()

    #create output directory
    os.makedirs(save_path, exist_ok=True) #to save the results array
    os.makedirs(f'{save_path}/output_far_away_dist{args.far_away_distance}', exist_ok=True) #to save simulation snapshsots

    #define masses of all bodies
    m_host = 1 | units.Msun
    m_pl = 1 | units.Mjupiter
    m_field = 1 | units.Msun

    a_sp = args.a_sp | units.AU

    #calculate critical velocity of the system, approximate planet and moon as one body
    v_crit = critical_velocity(m1=m_host, m2=m_pl, m3=m_field, a=a_sp)
    print(f'Critical velocity: {v_crit.in_(units.kms)}')

    #calculate the initial velocity of the field star
    v_inf = 3 | units.kms # roughly the 3D velocity dispersion of Orion (Wei et al., 2025)
    v20 = v_inf
    v20 = ff.vinit_from_vinf(v_inf, 20*a_sp, m_host+m_pl) #TODO: does 20*a_sp work for close in planets?
    print(f'Initial velocity of the field star at 20 a_pl: {v20.in_(units.kms)}')

    #initialize results array
    dtype = [('a_sp', 'f8'), ('v20','f8'), ('b','f8'), ('phi', 'f8'), ('theta', 'f8'), ('psi', 'f8'), ('f_pl', 'f8'), ('end_time', 'f8'), ('state', 'i1'), ('index', 'u4')]
    results = np.zeros((0,), dtype=dtype)

    #iterate over impact parameters
    index = 0
    bmax = a_sp
    bmin = 0 | units.AU
    area = np.pi*(bmax**2 - bmin**2)
    
    #calculate initial number of simulations
    n_sim = int(point_density * area)
    step_size = a_sp
    while True:
        #reset the counter for the interested states
        ffp_counter = 0
        temp_results = np.zeros((n_sim,), dtype=dtype)

        tqdm.write(f"Impact parameter range: {bmin.in_(units.AU)} to {bmax.in_(units.AU)}, Nsim={n_sim}")

        #pick n_sim combinations of angles and impact parameters
        impact_parameters = np.random.uniform(bmin.value_in(units.AU)**2, bmax.value_in(units.AU)**2, n_sim)
        impact_parameters = np.sqrt(impact_parameters) | units.AU
        phis = np.random.uniform(0, 2*np.pi, n_sim)
        thetas = np.arccos(np.random.uniform(0, 1, n_sim))
        psis = np.random.uniform(0, 2*np.pi, n_sim)
        f_pls = np.random.uniform(0, 2*np.pi, n_sim)

        i_sim = 0
        for b, phi, theta, psi, f_pl in tqdm(zip(impact_parameters, phis, thetas, psis, f_pls), total=n_sim):
            #generate host star, planet and moon
            bodies = generate_initial_conditions(M = m_host,
                                                 m_pl = m_pl,
                                                 a_pl = a_sp,
                                                 f_pl = f_pl,
                                                 n_moons = 0,
                                                 radii = [1|units.Rsun, 1|units.Rjupiter])

            #add field star to encounter
            bodies = add_encounter(bodies, m_field, b, v20, phi, theta, psi, radius = 1|units.Rsun)

            state = 0 #other state

            #run the simulation until the energy error is small enough
            timestep_parameter = 0.03 #0.03 is default value for Huayno and Hermite
            i=0
            while i < 4:
                evolved, _, end_time, stop_code = run_simulation(bodies, integrator='hermite',
                                                       timestep_parameter=timestep_parameter,
                                                       far_away_distance=far_away_distance,
                                                       stop_on_collision=True)
                #decrease the timestep if the energy error is too high
                if stop_code == 2:
                    timestep_parameter *= 0.5
                    i+=1
                    tqdm.write(f'Rerunning simulation with {timestep_parameter}')
                    if i == 4:
                        tqdm.write('Simulation failed, stopping.')
                        state = -1 #simulation failed
                        break
                elif stop_code == 1:
                    tqdm.write('Collision detected, stopping.')
                    state = -2 #collision detected
                    break
                elif stop_code == 3:
                    tqdm.write('Simulation took too long, stopping.')
                    state = -3 #max time reached
                    break
                else:
                    break

            if stop_code == 0:
                #check the state of the evolved system
                host_star, planet, field_star = evolved
                
                #check if planet is ejected
                if not bound(planet, host_star) and not bound(planet, field_star):
                    state = 1 #planet ejected
                    ffp_counter += 1
                elif bound(planet, field_star):
                    state = 2 #planet exchanged
                elif bound(planet, host_star):
                    state = 3 #original state


            #save evolved system to file and save initial conditions and state to array
            write_set_to_file(evolved, save_path+f'/output_far_away_dist{args.far_away_distance}/idx{index}_{state}.amuse', 'amuse', overwrite_file=True)
            temp_results[i_sim] = (a_sp.value_in(units.au), v20.value_in(units.kms), b.value_in(units.AU), phi, theta, psi, f_pl, end_time.value_in(units.yr), state, index)
            i_sim += 1
            index += 1
        
        results = np.concatenate((results, temp_results))

        tqdm.write(f"Finished, {ffp_counter} ({ffp_counter/n_sim*100:.1f}%) ffps.")
        if ffp_counter == 0:
            tqdm.write(f'No planet ejections occured. Stopping simulation.')
            break

        #increase the impact parameter linearly but keep the surface density constant
        bmin = bmax
        bmax += step_size
        area = np.pi*(bmax**2 - bmin**2)
        n_sim = int(point_density * area)
        
    
    #save the results array
    np.save(save_path+f'/results_far_away_dist_{far_away_distance}.npy', results)

    tqdm.write(f"Simulation took {time.time()-start_time} seconds")
    tqdm.write(f"Results saved to {save_path}")