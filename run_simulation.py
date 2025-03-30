from amuse.units import units, nbody_system, constants
from amuse.lab import Particles, Particle
from amuse.ext.orbital_elements import generate_binaries, orbital_elements
from amuse.io import write_set_to_file, read_set_from_file
from amuse.community.huayno.interface import Huayno
from amuse.community.hermite.interface import Hermite
from amuse.community.smalln.interface import SmallN

import matplotlib.pyplot as plt
import numpy as np
import os
from tqdm.auto import tqdm

def kinetic_energy(particle):
    'Returns the kinetic energy of a particle'
    return 0.5 * particle.mass * particle.velocity.length()**2

def potential_energy(particle, bodies):
    'Returns the potential energy of a particle'
    potential = 0 | units.J
    for body in bodies:
        if body != particle:
            distance = (particle.position - body.position).length()
            potential += -constants.G * particle.mass * body.mass / distance
    return potential

def mean_square_distance(bodies):
    'Returns the mean squared distance of the particle set'
    square_distances = []
    for i in range(len(bodies)):
        for j in range(i+1, len(bodies)):
            distance = (bodies[i].position - bodies[j].position).length()
            square_distances.append(distance**2)
    return np.mean(square_distances)

def system_check(bodies, far_away_distance=70|units.AU):
    '''
    Check if the simulation needs to be stopped. 
    Returns a boolean indicating if the simulation needs to be stopped and a code indicating the reason.
    '''
    stop = False

    unbound = False
    far_away = False
    moving_away = False

    star_idx = np.where(bodies.mass > 0.1 | units.MSun)[0]

    #get center of mass position and velocity
    com_position = bodies.center_of_mass()
    com_velocity = bodies.center_of_mass_velocity()
    for i in star_idx:
        body = bodies[i]
        #check if the particle is unbound, i.e. has positive total energy
        total_energy = kinetic_energy(body) + potential_energy(body, bodies)
        if total_energy > 0 | units.J:
            unbound = True

        #check if the particle is far away from the COM
        if body.position.length() > far_away_distance:
            far_away = True
        
        #check if the particle is moving away from the COM
        position = body.position.value_in(units.AU) - com_position.value_in(units.AU)
        velocity = body.velocity.value_in(units.km/units.s) - com_velocity.value_in(units.km/units.s)
        if np.dot(position, velocity) > 0:
            moving_away = True

        # if any particle is unbound, far away, and moving away, stop the simulation
        if unbound and far_away and moving_away:
            stop = True
            break
    
    return stop

def run_simulation(bodies, plot=False, integrator='hermite', save_path=None, timestep_parameter=0.03, far_away_distance=70|units.AU, stop_on_collision=False):
    if plot:
        #remove all files in the movie folder
        files = os.listdir('../movie')
        for file in files:
            os.remove(f'../movie/{file}')

    converter = nbody_system.nbody_to_si(bodies.mass.sum(), bodies[1].position.length())
    
    #initialize the integrator
    if integrator == 'hermite':
        gravity = Hermite(converter)
        gravity.parameters.dt_param = timestep_parameter
    elif integrator == 'huayno':
        gravity = Huayno(converter)
        gravity.parameters.timestep_parameter = timestep_parameter
    elif integrator == 'smalln':
        gravity = SmallN(converter)
    else:
        raise ValueError("Integrator must be one of 'hermite', 'huayno', or 'smalln'.")

    gravity.particles.add_particles(bodies)
    channel = gravity.particles.new_channel_to(bodies)

    if stop_on_collision:
        collision_detection = gravity.stopping_conditions.collision_detection
        collision_detection.enable()
    
    field_star = bodies[bodies.name == 'field_star'][0]
    max_time = 5*(far_away_distance + field_star.position.length()) / field_star.velocity.length()

    #calculate initial energy
    initial_energy = gravity.kinetic_energy + gravity.potential_energy
    energy_error = []

    time = 0 | units.yr
    stop_code = None
    while True:
        gravity.evolve_model(time)

        #check and save energy error
        error = np.abs((gravity.kinetic_energy + gravity.potential_energy - initial_energy) / initial_energy)
        energy_error.append(error)
        if error > 1e-2: #is this the correct threshold?
            tqdm.write(f"STOP! Energy error too high.")
            bodies, energy_error, stop_code = None, None, 2
            break

        #update bodies
        channel.copy()
        
        #plot the system and save the figure
        if plot:
            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(111, projection='3d')
            COM = bodies.center_of_mass()
            for body in bodies:
                ax.scatter(body.x.value_in(units.AU)-COM.x.value_in(units.AU),
                           body.y.value_in(units.AU)-COM.y.value_in(units.AU),
                           body.z.value_in(units.AU)-COM.z.value_in(units.AU),
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

        #check if the simulation needs to be stopped
        if system_check(bodies, far_away_distance=far_away_distance):
            stop_code = 0 # simulation finished without problems
            break
        elif stop_on_collision and collision_detection.is_set():
            stop_code = 1
            break
        elif time > max_time:
            stop_code = 3
            break

        time += 0.1 | units.yr


    gravity.stop()

    if save_path is not None:
        write_set_to_file(bodies, save_path, 'amuse', overwrite_file=True)

    return bodies, energy_error, time, stop_code