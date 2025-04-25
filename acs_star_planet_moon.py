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


if __name__ == '__main__':
    parser = ArgumentParser(description='Run a Monte Carlo simulation of an encounter between an equal mass binary and a field star.')
    parser.add_argument('--a_sp', type=float, default=1, help='Semi-major axis of the planet in AU')
    parser.add_argument('--density', type=int, default=10, help='Number of scattering experiments per pi*a_sp^2')
    parser.add_argument('--output', type=str, default='automatic_cs_output', help='Output directory for the results')
    parser.add_argument('--time_limit', type=float, default=0, help='Time limit for the simulation in hours, set to the same as slurm time limit. Set to 0 for no limit.')
    parser.add_argument('--b_init', type=float, default=0, help='Initial impact parameter in AU, use if you want to start from a specific impact parameter')
    parser.add_argument('--n_init', type=int, default=0, help='Initial number of simulations, use if you want to start from a specific number of simulations')
    parser.add_argument('--integrator', type=str, default='hermite', help='Integrator to use for the simulation, options are hermite, huayno, smalln and ph4')
    args = parser.parse_args()

    # np.random.seed(2208) #set random seed for reproducibility

    save_path = args.output
    a_sp = args.a_sp | units.AU
    point_density = args.density / (np.pi * a_sp**2) #number of simulations per planet orbit area
    slurm_time_limit = args.time_limit * 3600 #convert to seconds
    b_init = args.b_init | units.AU
    n_init = args.n_init

    print(f"Running simulation with a_sp {args.a_sp} and {args.density} simulation density")

    start_time = time.time()
    if slurm_time_limit > 0:
        max_time = start_time + slurm_time_limit - 5 * 60 # leave 5 minutes

    #create output directory
    os.makedirs(save_path, exist_ok=True) #to save the results array
    os.makedirs(f'{save_path}/output_sm_axis{args.a_sp}', exist_ok=True) #to save simulation snapshsots

    #define masses of all bodies
    m_host = 1 | units.Msun
    m_pl = 1 | units.Mjupiter
    m_moon = 1.4815e23 | units.kg #mass of Ganymede, use 8.93e22 kg for Io
    m_field = 1 | units.Msun

    #calculate critical velocity of the system, approximate planet and moon as one body
    v_crit = critical_velocity(m1=m_host, m2=m_pl+m_moon, m3=m_field, a=a_sp)
    print(f'Critical velocity: {v_crit.in_(units.kms)}')

    #calculate the initial velocity of the field star
    v_inf = 3 | units.kms # roughly the 3D velocity dispersion of Orion (Wei et al., 2025)
    v20 = ff.vinit_from_vinf(v_inf, 20*a_sp, m_host+m_pl+m_moon) #TODO: 20*a_sp does not work for close in planets?

    print(f'Initial velocity of the field star at 20 a_pl: {v20.in_(units.kms)}')

    #initialize results array #TODO: include moon parameters
    dtype = [('a_sp', 'f8'), ('v20','f8'), ('b','f8'), ('phi', 'f8'), ('theta', 'f8'), ('psi', 'f8'), ('f_pl', 'f8'), ('f_moon', 'f8'), ('end_time', 'f8'), ('state', 'i1'), ('index', 'u4')]
    results = np.zeros((0,), dtype=dtype)

    #iterate over impact parameters
    index = 0
    bmax = a_sp
    bmin = b_init
    area = np.pi*(bmax**2 - bmin**2)
    
    if n_init == 0:
        #calculate initial number of simulations based on the area and point density
        n_sim = int(point_density * area)
    else:
        #use the given number of simulations if previous simulations were run
        n_sim = n_init

    step_size = a_sp
    while True:
        start_time_b_range = time.time()
        #reset the counter for the interested states
        ffpm_counter = 0
        ffpnm_counter = 0
        temp_results = np.zeros((n_sim,), dtype=dtype)

        print(f"Impact parameter range: {bmin.in_(units.AU)} to {bmax.in_(units.AU)}, Nsim={n_sim}")

        #pick n_sim combinations of angles and impact parameters
        impact_parameters = np.random.uniform(bmin.value_in(units.AU)**2, bmax.value_in(units.AU)**2, n_sim)
        impact_parameters = np.sqrt(impact_parameters) | units.AU
        phis = np.random.uniform(0, 2*np.pi, n_sim)
        thetas = np.arccos(np.random.uniform(0, 1, n_sim))
        psis = np.random.uniform(0, 2*np.pi, n_sim)
        f_pls = np.random.uniform(0, 2*np.pi, n_sim)
        f_moons = np.random.uniform(0, 2*np.pi, n_sim)

        i_sim = 0
        for b, phi, theta, psi, f_pl, f_moon in zip(impact_parameters, phis, thetas, psis, f_pls, f_moons):
            #generate host star, planet and moon
            bodies = generate_initial_conditions(M = m_host,
                                                 m_pl = m_pl,
                                                 a_pl = a_sp,
                                                 m_moon = m_moon,
                                                 i_moon = 0|units.deg,
                                                 f_pl = f_pl,
                                                 f_moon = f_moon,
                                                 n_moons = 1)

            #add field star to encounter
            bodies = add_encounter(bodies, m_field, b, v20, phi, theta, psi)

            state = 0 #other state
            planet_ejected = False
            moon_ejected = False

            #run the simulation until the energy error is small enough
            timestep_parameter = 0.03 #0.03 is default value for Huayno and Hermite
            i=0
            while i < 4:
                start_time_single = time.time()
                evolved, _, end_time, stop_code = run_simulation(bodies, integrator=args.integrator,
                                                       timestep_parameter=timestep_parameter,
                                                       far_away_distance=70 * a_sp,
                                                       stop_on_collision=True)
                # print(f"Simulation took {time.time()-start_time_single} seconds")
                #decrease the timestep if the energy error is too high
                if stop_code == 2:
                    timestep_parameter *= 0.5
                    i+=1
                    print(f'Rerunning simulation with dt_param={timestep_parameter}')
                    if i == 4:
                        print('Simulation failed, stopping.')
                        state = -1 #simulation failed
                        break
                elif stop_code == 1:
                    print('Collision detected, stopping.')
                    state = -2 #collision
                    break
                elif stop_code == 3:
                    print('Simulation took too long, stopping.')
                    state = -3
                else:
                    break

            if state == 0:
                #check the state of the evolved system
                host_star, planet, moon, field_star = evolved
                
                #check if planet is ejected
                if not bound(host_star, planet) and not bound(field_star, planet):
                    planet_ejected = True
                    #check if moon is bound to planet
                    if bound(planet, moon):
                        moon_ejected = False
                        state = 1 # free foating planet with moon
                        ffpm_counter += 1
                    else:
                        moon_ejected = True
                        state = 2 # free floating planet no moon
                        ffpnm_counter += 1

                #check if moon is ejected while planet is bound to a star
                if not bound(host_star, moon) and not bound(field_star, moon) and not bound(planet, moon) and not planet_ejected:
                    moon_ejected = True
                    state = 3 # free floating moon bound planet


                #save evolved system to file and save initial conditions and state to array
                write_set_to_file(evolved, save_path+f'/output_sm_axis{args.a_sp}/{index}_{state}.amuse', 'amuse', overwrite_file=True)
            temp_results[i_sim] = (a_sp.value_in(units.au), v20.value_in(units.kms), b.value_in(units.AU), phi, theta, psi, f_pl, f_moon, end_time.value_in(units.yr), state, index)
            i_sim += 1
            index += 1

            #check time
            if slurm_time_limit > 0:
                if time.time() > max_time:
                    print(f"SLURM time limit reached. Saving results and stopping simulation.")
                    print(f'Stopped at b_min = {bmin.in_(units.AU)} and b_max = {bmax.in_(units.AU)} at sim {i_sim} of {n_sim}.')
                    break
        
        print(f"Average time per simulation: {(time.time()-start_time_b_range)/n_sim:.2f} seconds")
        
        results = np.concatenate((results, temp_results))

        #check if the simulation was stopped due to time limit
        if slurm_time_limit > 0:
            if time.time() > max_time:
                break

        print(f"Finished, {ffpm_counter} ({ffpm_counter/n_sim*100:.1f}%) ffps with moons and {ffpnm_counter} ({ffpnm_counter/n_sim*100:.1f}%) ffps without moon found.")
        if ffpm_counter == 0 and ffpnm_counter == 0:
            print(f'No planet ejections occured. Stopping simulation.')
            break

        #increase the impact parameter linearly but keep the surface density constant
        bmin = bmax
        bmax += step_size
        area = np.pi*(bmax**2 - bmin**2)
        n_sim = int(point_density * area)
        
    
    #save the results array
    np.save(save_path+f'/results_semi-major_{a_sp}.npy', results)

    print(f"Simulation took {time.time()-start_time} seconds")
    print(f"Results saved to {save_path}")