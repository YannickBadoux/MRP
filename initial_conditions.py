#code to generate the initial conditions for the scattering experiment
from amuse.units import units, constants
from amuse.lab import Particles, Particle
from amuse.ext.orbital_elements import generate_binaries, orbital_elements
from amuse.io import write_set_to_file, read_set_from_file
from amuse.units.optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np

from function_file import hill_radius

def critical_velocity(bodies=None, m1=None, m2=None, m3=None, a=None):
    '''Calculate the critical velocity for a given system.'''
    if bodies is not None:
        try:
            host_star = bodies[bodies.name == 'host_star'][0]
            planet = bodies[bodies.name == 'planet'][0]
            field_star = bodies[bodies.name == 'field_star'][0]

            m1 = host_star.mass
            m2 = planet.mass
            m3 = field_star.mass
            a = orbital_elements(host_star, planet)[2]
        except:
            raise ValueError("Bodies does not contain the correct particles")
        
    elif bodies is None:
        if m1 is None or m2 is None or m3 is None or a is None:
            raise ValueError("If bodies is None, m1, m2, m3, and a must be provided")
    
    vc2 = constants.G * (m1*m2*(m1+m2+m3)) / (m3*(m1+m2)*a)

    return vc2.sqrt()

def generate_initial_conditions(M, m_pl, a_pl, a_pm=None, m_moon=0.01495|units.MEarth, i_moon=0|units.deg, f_pl=0|units.deg, f_moon=0|units.deg,
                                plot=False, save_path=None, n_moons=1, radii=None):
    '''
    Generate the initial conditions for the scattering experiment.
    Parameters:
    M: mass of the host star
    m_pl: mass of the planet
    a_pl: semi-major axis of the planet
    a_pm: semi-major axis of the moon(s), leave as None to use 1/3 of the Hill radius
    m_moon: mass of the moon(s), default is the mass of Io
    i_moon: inclination of the moon(s)
    f_pl: true anomaly of the planet
    f_moon: true anomaly of the moon(s)
    plot: flag to plot the initial conditions
    save_path: path to save the initial conditions, set to None to not save
    n_moons: number of moons to add to the planet
    radii: radii of the bodies, set to None for point particles
    Returns:
    bodies: particle set containing the host star, planet and moon(s)
    '''
    #initialize host star and planet particles
    host_star, planet = generate_binaries(M,
                                          m_pl,
                                          a_pl,
                                          true_anomaly=f_pl)
    host_star.name = "host_star"
    planet.name = "planet"
    
    if n_moons > 0:
        #place outermost moon at 1/3 of the Hill radius
        if a_pm is None:
            a_moon = 1/3*hill_radius(a_pl, M, m_pl)
        else:
            a_moon = a_pm

        _, moon = generate_binaries(planet.mass,
                                    m_moon,
                                    a_moon,
                                    inclination=i_moon,
                                    true_anomaly=f_moon)
        
        moon.name = "moon0"
        moon.position += planet.position
        moon.velocity += planet.velocity

        #TODO: add code to add more moons in a resonant chain
        if n_moons > 1:
            raise NotImplementedError("Only one moon is currently supported")

    #add all to a single particle set
    bodies = Particles()
    bodies.add_particle(host_star)
    bodies.add_particle(planet)
    if n_moons != 0:
        bodies.add_particle(moon)

    if radii is not None:
        bodies.radius = radii

    #put host_star at the origin
    bodies.position -= host_star.position
    # bodies.velocity -= host_star.velocity

    #plot the initial conditions if plot=True
    if plot:
        fig, ax = plt.subplots()
        ax.scatter(bodies.x.value_in(units.AU), bodies.y.value_in(units.AU))
        ax.set_aspect('equal')
        plt.show()

    if save_path is not None:
        write_set_to_file(bodies, save_path, 'amuse', overwrite_file=True)

    return bodies

def add_encounter(bodies, M, impact_parameter, v_inf, phi, theta, psi, radius=None):
    #initialize the field star
    field_star = Particle()
    field_star.mass = M
    field_star.name = 'field_star'

    if radius is not None:
        field_star.radius = radius

    #get planet semi-major axis
    planet = bodies[bodies.name == 'planet'][0]
    a_pl = orbital_elements(bodies[0], planet)[2].value_in(units.AU)

    #place field star at the correct distance
    initial_distance = 20 * a_pl
    field_star.position = [initial_distance * np.cos(phi) * np.sin(theta),
                           initial_distance * np.sin(phi) * np.sin(theta),
                           initial_distance * np.cos(theta)] | units.AU
    
    #change the position of the field star so the impact parameter is correct
    impact_parameter = impact_parameter.value_in(units.AU)
    position_change = np.array([impact_parameter * np.cos(psi),
                                impact_parameter * np.sin(psi),
                                0])
    rotation_matrix = np.array([[np.cos(phi) * np.cos(theta), -np.sin(phi), np.cos(phi) * np.sin(theta)],
                                [np.sin(phi) * np.cos(theta), np.cos(phi), np.sin(phi) * np.sin(theta)],
                                [-np.sin(theta), 0, np.cos(theta)]])
    field_star.position += np.dot(rotation_matrix, position_change) | units.AU
    
    #set the velocity of the field star
    field_star.velocity = [-v_inf.value_in(units.kms) * np.cos(phi) * np.sin(theta),
                           -v_inf.value_in(units.kms) * np.sin(phi) * np.sin(theta),
                           -v_inf.value_in(units.kms) * np.cos(theta)] | units.kms
    
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

    #generate the initial conditions
    bodies = generate_initial_conditions(options.M,
                                         options.m_pl,
                                         options.a_pl,
                                         options.m_moon,
                                         options.plot,
                                         options.save_path)