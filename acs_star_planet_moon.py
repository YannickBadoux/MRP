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
    parser.add_argument('--output', type=str, default='automatic_cs_output', help='Output directory for the results')
    args = parser.parse_args()

    # n_sim = args.n_sim
    save_path = args.output
    point_density = args.density | units.AU**-2 #number of simulations per square AU

    tqdm.write(f"Running simulation with a_sp {args.a_sp} and {point_density} simulation density")

    start_time = time.time()

    #create output directory
    os.makedirs(save_path, exist_ok=True) #to save the results array
    os.makedirs(f'{save_path}/output_sm_axis{args.a_sp}', exist_ok=True) #to save simulation snapshsots

    #define masses of all bodies
    m_host = 1 | units.Msun
    m_pl = 1 | units.Mjupiter
    m_moon = 1.4815e23 | units.kg #mass of Ganymede, use 8.93e22 kg for Io
    m_field = 1 | units.Msun

    a_sp = args.a_sp | units.AU

    #calculate critical velocity of the system, approximate planet and moon as one body
    v_crit = critical_velocity(m1=m_host, m2=m_pl+m_moon, m3=m_field, a=a_sp)
    print(f'Critical velocity: {v_crit.in_(units.kms)}')

    #calculate the initial velocity of the field star
    v_inf = 3 | units.kms # roughly the 3D velocity dispersion of Orion (Wei et al., 2025)
    v20 = v_inf
    # v20 = v20_from_vinf(v_inf, a_sp, m_field) #TODO: update the calculation of v20
    print(f'Initial velocity of the field star at 20 a_pl: {v20.in_(units.kms)}')

    #initialize results array #TODO: include moon parameters
    dtype = [('a_sp', 'f8'), ('v20','f8'), ('b','f8'), ('phi', 'f8'), ('theta', 'f8'), ('psi', 'f8'), ('f_pl', 'f8'), ('f_moon', 'f8'), ('end_time', 'f8'), ('state', 'u1'), ('index', 'u4')]
    results = np.zeros((0,), dtype=dtype)

    #iterate over impact parameters
    index = 0
    bmax = a_sp
    bmin = 0 | units.AU
    area = np.pi*(bmax**2 - bmin**2)
    
    #calculate initail point density
    n_sim = int(point_density * area)
    step_size = a_sp
    while True:
        #reset the counter for the interested states
        ffpm_counter = 0
        ffpnm_counter = 0
        temp_results = np.zeros((n_sim,), dtype=dtype)

        tqdm.write(f"Impact parameter range: {bmin.in_(units.AU)} to {bmax.in_(units.AU)}, Nsim={n_sim}")

        #pick n_sim combinations of angles and impact parameters
        impact_parameters = np.random.uniform(bmin.value_in(units.AU)**2, bmax.value_in(units.AU)**2, n_sim)
        impact_parameters = np.sqrt(impact_parameters) | units.AU
        phis = np.random.uniform(0, 2*np.pi, n_sim)
        thetas = np.arccos(np.random.uniform(0, 1, n_sim))
        psis = np.random.uniform(0, 2*np.pi, n_sim)
        f_pls = np.random.uniform(0, 2*np.pi, n_sim)
        f_moons = np.random.uniform(0, 2*np.pi, n_sim)

        i_sim = 0
        for b, phi, theta, psi, f_pl, f_moon in tqdm(zip(impact_parameters, phis, thetas, psis, f_pls, f_moons), total=n_sim):
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

            #run the simulation until the energy error is small enough
            timestep_parameter = 0.03 #0.03 is default value for Huayno and Hermite
            i=0
            while i < 4:
                evolved, _, end_time, stop_code = run_simulation(bodies, integrator='hermite',
                                                       timestep_parameter=timestep_parameter)
                #decrease the timestep if the energy error is too high
                if stop_code == 2:
                    timestep_parameter *= 0.5
                    i+=1
                    tqdm.write(f'Rerunning simulation with {timestep_parameter}')
                    if i == 4:
                        tqdm.write('Simulation failed, stopping.')
                        state = -1 #simulation failed
                        break
                else:
                    break

            if not state == 3:
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

                # #check if moon is ejected
                # if not bound(host_star, moon) and not bound(field_star, moon) and not bound(planet, moon) and not planet_ejected:
                #     moon_ejected = True
                #     state = 'ffmbp' # free floating moon bound planet


                #save evolved system to file and save initial conditions and state to array
                write_set_to_file(evolved, save_path+f'/output_sm_axis{args.a_sp}/{index}_{state}.amuse', 'amuse', overwrite_file=True)
            temp_results[i_sim] = (a_sp.value_in(units.au), v20.value_in(units.kms), b.value_in(units.AU), phi, theta, psi, f_pl, f_moon, end_time.value_in(units.yr), state, index)
            i_sim += 1
            index += 1
        
        results = np.concatenate((results, temp_results))

        tqdm.write(f"Finished, {ffpm_counter} ({ffpm_counter/n_sim*100:.1f}%) ffps with moons and {ffpnm_counter} ({ffpnm_counter/n_sim*100:.1f}%) ffps without moon found.")
        if ffpm_counter == 0 and ffpnm_counter == 0:
            tqdm.write(f'No planet ejections occured. Stopping simulation.')
            break
        # if time.time()-start_time > 3600:
        #     tqdm.write(f"Simulation took too long. Stopping.")
        #     break

        #increase the impact parameter linearly but keep the surface density constant
        bmin = bmax
        bmax += step_size
        area = np.pi*(bmax**2 - bmin**2)
        n_sim = int(point_density * area)
        
    
    #save the results array
    np.save(save_path+f'/results_semi-major_{a_sp}.npy', results)

    tqdm.write(f"Simulation took {time.time()-start_time} seconds")
    tqdm.write(f"Results saved to {save_path}")