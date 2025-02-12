import initial_conditions as ic
from run_simulation import run_simulation
from analyse_result import analyse_result

import numpy as np

from amuse.units import units

if __name__ == '__main__':
    velocities = [1,5] | units.kms
    phis = np.linspace(0, 2*np.pi, 10)
    thetas = np.linspace(0, np.pi, 10)
    psis = np.linspace(0, 2*np.pi, 10)

    #initialize results array
    dtype = [('v_inf', 'f8'), ('phi', 'f8'), ('theta', 'f8'), ('psi', 'f8'), ('state', int)]
    results = np.zeros((len(velocities) * len(phis) * len(thetas) * len(psis),), dtype=dtype)
    index = 0

    #iterate over angles and velocities
    for v_inf in velocities:
        for phi in phis:
            for theta in thetas:
                print(index)
                for psi in psis:
                        #generate host star, planet and moon
                        bodies = ic.generate_initial_conditions(1|units.Msun,
                                                                1|units.Mjupiter,
                                                                1|units.AU,
                                                                0.01495|units.Mearth,
                                                                plot=False)
                        
                        #add field star to encounter
                        bodies = ic.add_encounter(bodies,
                                                1|units.Msun,
                                                1|units.AU,
                                                v_inf,
                                                phi, theta, psi)

                        #run the simulation
                        bodies, energy_error = run_simulation(bodies,
                                                              plot=False,
                                                              integrator='hermite')
                        
                        state = analyse_result(bodies)
                        results[index] = (v_inf.value_in(units.kms), phi, theta, psi, state)
                        index += 1
    
    #save the results
    np.save('test_results2.npy', results)