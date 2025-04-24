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

def cross_section(data, state, normalize=False, a_sp=None):
    state_mask = data['state']==state
    und_mask = (data['state']==-1) | (data['state']==-3) #Failed or Timeout

    b_max = np.max(data['b'][state_mask]) | units.AU

    cs = np.pi * b_max**2 * np.sum(state_mask) / len(state_mask)
    stat_error = 1/np.sqrt(np.sum(state_mask)) * cs
    system_error = np.pi * b_max**2 * np.sum(und_mask) / len(state_mask)

    if normalize and a_sp is not None:
        cs /= np.pi * a_sp**2
        stat_error /= np.pi * a_sp**2
        system_error /= np.pi * a_sp**2
    if normalize and a_sp is None:
        raise ValueError("a_sp must be provided for normalization")

    return cs, stat_error, system_error

def state_dict():
    state_dict = {0: ('Other', 'Other'),
              1: ('FFPM', 'Free Floating Planet Moon Pair'),
              2: ('FFPWM', 'Free Floating Planet Without Moon'),
              3: ('FFMBP', 'Free Floating Moon, Bound Planet')}
    return state_dict