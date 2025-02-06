#code to generate the initial conditions for the scattering experiment
from amuse.units import units, constants
from amuse.lab import Particles, Particle
from amuse.ext.orbital_elements import generate_binaries, orbital_elements
from amuse.io import write_set_to_file
from amuse.units.optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np


def hill_radius(M, m, a):
    '''
    Calculate the Hill radius of a planet orbiting a star
    Parameters:
    M: mass of the star
    m: mass of the planet
    a: semi-major axis of the planet
    Returns:
    R_H: Hill radius of the planet
    '''
    return a * (m / (3 * M))**(1/3)

def generate_initial_conditions(M, m_pl, a_pl, m_moon, plot=False, save_path='initial_conditions.amuse'):
    #initialize host star and planet particles
    host_star, planet = generate_binaries(M,
                                          m_pl,
                                          a_pl)
    host_star.name = "host_star"
    planet.name = "planet"

    #place outermost moon at 1/3 of the Hill radius
    a_moon = 1/3*hill_radius(M, m_pl, a_pl)
    _, moon = generate_binaries(planet.mass,
                                 m_moon,
                                 a_moon,
                                 inclination=0|units.deg)
    
    moon.name = "moon0"
    moon.position += planet.position
    moon.velocity += planet.velocity

    #TODO: add code to add more moons in a resonant chain

    #add all to a single particle set
    bodies = Particles()
    bodies.add_particle(host_star)
    bodies.add_particle(planet)
    bodies.add_particle(moon)

    #plot the initial conditions if plot=True
    if plot:
        fig, ax = plt.subplots()
        ax.scatter(bodies.x.value_in(units.AU), bodies.y.value_in(units.AU))
        ax.set_aspect('equal')
        plt.show()

    if save_path is not None:
        write_set_to_file(bodies, save_path, 'amuse', overwrite_file=True)

    return bodies

def add_encounter(bodies, M, impact_parameter, v_inf, phi, theta, psi):
    #initialize the field star
    field_star = Particle()
    field_star.mass = M
    field_star.name = 'field_star'

    #get planet semi-major axis
    planet = bodies[bodies.name == 'planet'][0]
    a_pl = orbital_elements(bodies[0], planet)[2].value_in(units.AU)

    initial_distance = 20 * a_pl
    field_star.position = [initial_distance * np.cos(phi) * np.sin(theta),
                           initial_distance * np.sin(phi) * np.sin(theta),
                           initial_distance * np.cos(theta)] | units.AU
    field_star.position += [impact_parameter.value_in(units.AU) * np.cos(psi),
                            impact_parameter.value_in(units.AU) * np.sin(psi),
                            0] | units.AU
    
    bodies.add_particle(field_star)
    return bodies


if __name__ == '__main__':
    #parse command line arguments
    parser = OptionParser()
    parser.add_option("-M", unit=units.MSun, type="float", default=1.0,
                      help="Mass of the host star in solar masses")
    parser.add_option("--m_pl", unit=units.MJupiter, type="float", default=1.0,
                        help="Mass of the planet in Jupiter masses")
    parser.add_option("--a_pl", unit=units.AU, type="float", default=1.0,
                        help="Semi-major axis of the planet in AU")
    parser.add_option("--m_moon", unit=units.MEarth, type="float", default=0.01495,
                        help="Mass of the moon(s) in Earth masses, default is the mass of the Io (0.01495 MEarth)")
    parser.add_option("--seed", type="int", default=2208,
                        help="Random seed for the initial conditions")
    parser.add_option("--plot", action="store_true",
                        help="Flag to plot the initial conditions")
    parser.add_option("--save_path", type="string", default='initial_conditions.amuse',
                        help="Path to save the initial conditions")
    options, arguments = parser.parse_args()

    #set random seed
    np.random.seed(options.seed)

    #generate the initial conditions
    bodies = generate_initial_conditions(options.M,
                                         options.m_pl,
                                         options.a_pl,
                                         options.m_moon,
                                         options.plot,
                                         options.save_path)

