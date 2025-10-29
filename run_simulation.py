from amuse.units import units, nbody_system, constants
from amuse.lab import Particles
from amuse.io import write_set_to_file
from amuse.community.huayno.interface import Huayno
from amuse.community.hermite.interface import Hermite
from amuse.community.smalln.interface import SmallN
from amuse.community.ph4.interface import ph4

import matplotlib.pyplot as plt
import numpy as np
import os

from analyse_result import find_bound_particles, bound

def kinetic_energy(particle, bodies):
    'Returns the kinetic energy of a particle relative to the center of mass'
    vcom = bodies.center_of_mass_velocity()
    relative_velocity = particle.velocity - vcom
    return 0.5 * particle.mass * relative_velocity.length()**2

def potential_energy(particle, bodies):
    'Returns the potential energy of a particle relative to the center of mass'
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
        body_pos = body.position - com_position
        #check if the particle is unbound, i.e. has positive total energy
        total_energy = kinetic_energy(body, bodies) + potential_energy(body, bodies)
        if total_energy > 0 | units.J:
            unbound = True

        #check if the particle is far away from the COM
        if body_pos.length() > far_away_distance:
            far_away = True
        
        #check if the particle is moving away from the COM
        position = body.position.value_in(units.AU) - com_position.value_in(units.AU)
        velocity = body.velocity.value_in(units.km/units.s) - com_velocity.value_in(units.km/units.s)
        if np.dot(position, velocity) > 0:
            moving_away = True

        # if any star is unbound, far away and moving away, dynamical encounter is over
        if unbound and far_away and moving_away:
            stop = True
            break
    
    if not stop:
        host_star, planet, moon, field_star = bodies
        #check if planet is unbound from both stars and stop simulation
        if not bound(host_star, planet) and not bound(field_star, planet):
            # check if moon is bound to planet and split accordingly
            if bound(planet, moon):
                ffp = Particles()
                ffp.add_particle(planet)
                ffp.add_particle(moon)

                other_bodies = Particles()
                other_bodies.add_particle(host_star)
                other_bodies.add_particle(field_star)
            else:
                ffp = Particles()
                ffp.add_particle(planet)

                other_bodies = Particles()
                other_bodies.add_particle(host_star)
                other_bodies.add_particle(field_star)
                other_bodies.add_particle(moon)

            ffp_com = ffp.center_of_mass()
            ffp_vcom = ffp.center_of_mass_velocity()

            ffp_rel_pos = ffp_com - bodies.center_of_mass()
            ffp_rel_vel = ffp_vcom - bodies.center_of_mass_velocity()
            #check if the ffp is far away and moving away from the COM
            if (ffp_rel_pos.length() > far_away_distance) and (np.dot(ffp_rel_pos.value_in(units.AU), ffp_rel_vel.value_in(units.km/units.s)) > 0):
                stop = True

    return stop

def run_simulation(bodies, plot=False, integrator='hermite', save_path=None, timestep_parameter=0.03, far_away_distance=70|units.AU, stop_on_collision=False):
    '''
    Run single scattering experiment with the given bodies.
    
    Parameters:
    bodies: AMUSE particle set containing the bodies to simulate
    plot: flag to plot the simulation at each diagnostic time
    integrator: integrator to use, must be one of 'hermite' (default), 'huayno', 'smalln' or 'ph4'
    save_path: path to save the simulation results, set to None to not save
    timestep_parameter: parameter for the integrator, default is 0.03
    far_away_distance: distance from the center of mass to consider a particle far away, default is 70 AU
    stop_on_collision: flag to stop the simulation on collision, default is False. bodies must have radius attribute
    
    Returns:
    bodies: AMUSE particle set containing the bodies after the simulation
    energy_error: list of energy errors at each diagnostic time
    time: time at which the simulation stopped
    stop_code: code indicating the reason for stopping the simulation
    '''

    if stop_on_collision:
        #check if bodies have radii
        try:
            radii = bodies.radius
        except AttributeError:
            raise ValueError("Bodies must have radii to use collision detection.")

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
    elif integrator == 'ph4':
        timestep_parameter = 0.014 * timestep_parameter / 0.03 # default value for ph4 is 0.014 instead of 0.03
        gravity = ph4(converter)
        gravity.parameters.timestep_parameter = timestep_parameter
    else:
        raise ValueError("Integrator must be one of 'hermite', 'huayno', or 'smalln'.")

    gravity.particles.add_particles(bodies)
    channel = gravity.particles.new_channel_to(bodies)

    if stop_on_collision:
        collision_detection = gravity.stopping_conditions.collision_detection
        collision_detection.enable()
    
    #max model time, taken from Wang et al. 2024
    field_star = bodies[bodies.name == 'field_star'][0]
    max_time = 5*(far_away_distance + field_star.position.length()) / field_star.velocity.length()
    diag_time = max_time / 5000

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
            print(f"STOP! Energy error too high.")
            bodies, energy_error, stop_code = bodies, None, 2
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
            for ci in range(len(collision_detection.particles(0))):
                collision = Particles(particles=[collision_detection.particles(0)[ci], collision_detection.particles(1)[ci]])
                collision = collision.get_intersecting_subset_in(bodies)
                print(f"Collision between {collision[0].name} and {collision[1].name}, at time {time.value_in(units.yr):.2f} yr.", end=' ')
            break
        elif time > max_time:
            stop_code = 3
            break

        time += diag_time


    gravity.stop()

    if save_path is not None:
        write_set_to_file(bodies, save_path, 'amuse', overwrite_file=True)

    return bodies, energy_error, time, stop_code