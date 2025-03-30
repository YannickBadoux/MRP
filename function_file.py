from amuse.lab import units, constants

import numpy as np

def max_impact_parameter(v, a, C=4, e=0):
    '''Calculates the maximum impact parameter for a given velocity, 
    semi-major axis and eccentricity.'''
    D = 0.6*(1+e)
    return (C/v + D) * a

def hill_radius(a, m1, m2, e=0):
    '''Calculates the Hill radius
    a: semi-major axis
    m1: mass of the most massive body
    m2: mass of the second body
    e: eccentricity of the system
    '''
    return a*(1-e)*((m2)/(3*(m1+m2)))**(1/3)

def vinit_from_vinf(v_inf, d_init, m_system):
    '''
    Calculates the initial velocity of the field star for a given velocity at infinity.
    v_inf: velocity of the field star at infinity
    d_init: initial distance between the field star and the system
    m_system: total mass of the system without the field star
    '''
    return np.sqrt(v_inf**2 + 2*constants.G*m_system/d_init)

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

def total_energy(particle, bodies):
    'Returns the total energy of a particle'
    return kinetic_energy(particle) + potential_energy(particle, bodies)