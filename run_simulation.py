from amuse.units import units, nbody_system
from amuse.lab import Particles, Particle
from amuse.ext.orbital_elements import generate_binaries, orbital_elements
from amuse.io import write_set_to_file, read_set_from_file
from amuse.community.huayno.interface import Huayno
from amuse.community.hermite.interface import Hermite
from amuse.community.smalln.interface import SmallN

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

def mean_square_distance(bodies):
    'Returns the mean squared distance of the particle set'
    square_distances = []
    for i in range(len(bodies)):
        for j in range(i+1, len(bodies)):
            distance = (bodies[i].position - bodies[j].position).length()
            square_distances.append(distance**2)
    return np.mean(square_distances)

def system_check(bodies):
    '''
    Check if the simulation needs to be stopped according to the procedure described in Appendix A 
    of Hut & Bahcall 1983. 
    Returns a boolean indicating if the simulation needs to be stopped and a code indicating the reason.
    '''
    
    return stop, code

def run_simulation(bodies, plot=False, integrator='hermite', save_path=None):
    converter = nbody_system.nbody_to_si(bodies.mass.sum(), bodies[1].position.length())

    if integrator == 'hermite':
        gravity = Hermite(converter)
    elif integrator == 'huayno':
        gravity = Huayno(converter)
    elif integrator == 'smalln':
        gravity = SmallN(converter)
    else:
        raise ValueError("Integrator must be one of 'hermite', 'huayno', or 'smalln'.")

    # gravity.parameters.timestep_parameter = 0.005
    gravity.particles.add_particles(bodies)

    channel = gravity.particles.new_channel_to(bodies)

    end_time = 2*(bodies[0].position - bodies[3].position).length() / bodies[3].velocity.length()
    times = np.linspace(0, end_time.value_in(units.yr), 1000) | units.yr

    initial_energy = gravity.kinetic_energy + gravity.potential_energy
    energy_error = []
    min_square_distance = min_square_distance(bodies)

    for time in times:
        gravity.evolve_model(time)

        #check and save energy error
        error = np.abs((gravity.kinetic_energy + gravity.potential_energy - initial_energy) / initial_energy)
        energy_error.append(error)
        
        if error > 1e-2:
            print(f"STOP! Energy error too high: {error}, Time: {int(time.value_in(units.yr))}")
            break

        channel.copy()

        min_square_distance = min([min_square_distance, mean_square_distance(bodies)])

        #check if the simulation needs to be stopped
        stop, code = system_check(bodies)
        if stop:
            print(f"STOP! Code {code}, Time: {int(time.value_in(units.yr))}")
            break

        if plot:
            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(111, projection='3d')
            for body in bodies:
                ax.scatter(body.x.value_in(units.AU), body.y.value_in(units.AU), body.z.value_in(units.AU),
                        label=body.name)
                
            ax.legend(loc='upper left')
            ax.set_xlabel('x [AU]')
            ax.set_ylabel('y [AU]')
            ax.set_zlabel('z [AU]')

            ax.set_xlim(-10,10)
            ax.set_ylim(-10,10)
            ax.set_zlim(-10,10)
            fig.savefig(f'../movie/simulation{int(time.value_in(units.day))}.png', dpi=200)
            plt.close()

    gravity.stop()

    if save_path is not None:
        write_set_to_file(bodies, save_path, 'amuse', overwrite_file=True)

    return bodies, energy_error
