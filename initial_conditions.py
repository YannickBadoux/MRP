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

def generate_initial_conditions(M, m_pl, a_pl, a_pm=None, m_moon=0.01495|units.MEarth, i_moon=0|units.deg, f_pl=0|units.deg, f_moon=0|units.deg, e_pl=0, e_pm=0,
                                plot=False, save_path=None, n_moons=1, radii=None, MMR=None, path_to_ic=None):
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
    MMR: flag to put moons in mean motion resonance
    path_to_ic: path to the initial conditions file, if None, the initial conditions will be generated from scratch
    
    Returns:
    bodies: particle set containing the host star, planet and moon(s)
    '''
    if path_to_ic is not None:
        #load the initial conditions from a file
        bodies = read_set_from_file(path_to_ic, 'amuse')
        return bodies
    
    #initialize host star and planet particles
    host_star, planet = generate_binaries(M,
                                          m_pl,
                                          a_pl,
                                          true_anomaly=f_pl,
                                          eccentricity=e_pl)
    host_star.name = "host_star"
    host_star.type = "star"
    planet.name = "planet"
    planet.type = "planet"
    
    # add no moons if n_moons == 0
    if n_moons > 0:
        #generate 1 moon
        if n_moons == 1:               
            #place moon at 1/3 of the Hill radius if no semi-major axis is given
            if a_pm is None:
                a_moon = 1/3*hill_radius(a_pl, M, m_pl)
            else:
                a_moon = a_pm

            _, moon = generate_binaries(planet.mass,
                                        m_moon,
                                        a_moon,
                                        inclination=i_moon,
                                        true_anomaly=f_moon,
                                        eccentricity=e_pm)
            
            moon.name = "moon0"
            moon.type = "moon"
            moon.position += planet.position
            moon.velocity += planet.velocity
        
        elif n_moons > 1:
            #check if m_moon is a single value or an array
            try:
                _ = len(m_moon)
            except TypeError:
                # if m_moon is a single value, make it an array of the same size as n_moons
                m_moon = [m_moon] * n_moons

            #check if m_moon is the same size as n_moons
            if len(m_moon) != n_moons:
                raise ValueError("m_moon must be the same size as n_moons or a single scalar")
            


    #add all to a single particle set
    bodies = Particles()
    bodies.add_particle(host_star)
    bodies.add_particle(planet)
    if n_moons != 0:
        bodies.add_particles(moon)

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
    field_star.type = 'star'

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

def generate_planetary_system(M, m_pl, a_pl, e_pl, i_pl, f_pl, lan_pl, aop_pl, m_moon=1.4815e23 | units.kg):
    '''Generate a planetary system with a host star, planet(s), and moon(s).
    M: mass of the host star
    m_pl: mass of the planet(s)
    a_pl: semi-major axis of the planet(s)
    e_pl: eccentricity of the planet(s)
    i_pl: inclination of the planet(s)
    f_pl: true anomaly of the planet(s)
    lan_pl: longitude of ascending node of the planet(s)
    aop_pl: argument of periapsis of the planet(s)
    m_moon: mass of the moon(s), default is the mass of Ganymede
    Returns: Particles object containing the host star, planet(s), and moon(s)
    '''
    bodies = Particles()

    for i in range(len(M)):
        host_star, planet = generate_binaries(M,
                                              m_pl[i],
                                              a_pl[i],
                                              true_anomaly=f_pl[i],
                                              eccentricity=e_pl[i],
                                              inclination=i_pl[i],
                                              longitude_of_ascending_node=lan_pl[i],
                                              argument_of_periapsis=aop_pl[i])
        host_star.name = 'host_star'
        host_star.type = 'star'
        planet.name = f'planet{i}'
        planet.type = 'planet'

        _, moon = generate_binaries(planet.mass,
                                    m_moon,
                                    1/3 * hill_radius(a_pl[i], M, m_pl[i]))
        moon.name = f'moon{i}'
        moon.type = 'moon'
        moon.position += planet.position
        moon.velocity += planet.velocity

        if i == 0: #only add the host star once
            bodies.add_particle(host_star)
        bodies.add_particle(planet)
        bodies.add_particle(moon)

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