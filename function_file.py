from amuse.lab import units, constants
from amuse.ext.orbital_elements import generate_binaries, orbital_elements

import numpy as np
import matplotlib.pyplot as plt

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

def cross_section(data, state, normalize=False):
    a_sp = data['a_sp'][0] | units.AU

    #calculate the impact parameter ranges in the data
    max_in_data = np.max(data['b']) 
    #round up to the nearest multiple of a_sp
    max_b_range = np.ceil(max_in_data/a_sp.value_in(units.AU)) * a_sp.value_in(units.AU)

    b_range_mins = np.arange(0, max_b_range, a_sp.value_in(units.AU))
    b_range_maxes = np.append(b_range_mins[1:], max_b_range)

    bin_index = 0
    #find the last bin where the interested state is present
    for i in range(len(b_range_maxes)):
        b_range = (b_range_mins[i], b_range_maxes[i])
        mask = (data['b'] >= b_range[0]) & (data['b'] < b_range[1]) & (data['state'] == state)
        if np.sum(mask) > 0:
            bin_index = max(bin_index, i)

    b_max = b_range_maxes[bin_index] | units.AU # maximum possible generated impact parameter for state

    #only consider the impact parameters that are less than the maximum
    n_sim = np.sum(data['b'] < b_max.value_in(units.AU))
    state_mask = (data['state']==state) & (data['b'] < b_max.value_in(units.AU))
    und_mask = (data['state']==-1) | (data['state']==-3) & (data['b'] < b_max.value_in(units.AU)) #Failed or Timeout

    cs = np.pi * b_max**2 * np.sum(state_mask) / n_sim
    stat_error = 1/np.sqrt(np.sum(state_mask)) * cs
    system_error = np.pi * b_max**2 * np.sum(und_mask) / n_sim

    if normalize:
        cs /= np.pi * a_sp**2
        stat_error /= np.pi * a_sp**2
        system_error /= np.pi * a_sp**2

    return cs, stat_error, system_error

def state_dict():
    state_dict = {0: ('Other', 'Other'),
              1: ('FFPM', 'Free Floating Planet Moon Pair'),
              2: ('FFPWM', 'Free Floating Planet Without Moon'),
              3: ('FFMBP', 'Free Floating Moon, Bound Planet')}
    return state_dict

def all_state_dict():
    state_dict = {0: ('Other', 'Other'),
              1: ('FFPM', 'Free Floating Planet Moon Pair'),
              2: ('FFPWM', 'Free Floating Planet Without Moon'),
              3: ('FFMBP', 'Free Floating Moon, Bound Planet'),
              -1: ('Failed', 'Too high energy error'),
              -2: ('Collision', 'Collision'),
              -3: ('Timeout', 'Timeout'),}
    return state_dict

def resonance_semi_major_axis(a1, M, m1, m2, res):
    '''
    Calculates the semi-major axis of a body in resonance with another body.
    a1: semi-major axis of the first body
    M: mass of the central body
    m1: mass of the first body
    m2: mass of the second body, the one in resonance
    res: resonance ratio, >1 for larger orbit, <1 for smaller orbit
    '''
    period1 = (a1**3 * 4*np.pi**2 / (constants.G * (M+m1)))**(1/2)
    period2 = period1 * res
    a2 = (constants.G*(M+m2) * period2**2 / (4*np.pi**2))**(1/3)
    return a2

def plot_orbit(primary, secondary):
    """Calculate and plot the orbit of a binary system.
    Args:
        primary: body object representing the primary object
        secondary: body object representing the secondary object
    Returns:
        orbit_plot: line object representing the plotted orbit
        prograde: boolean indicating if the orbit is prograde or retrograde
    """
    #calculate the orbital elements
    m1, m2, a, e, _, i, lan, ap = orbital_elements(primary, secondary)

    prograde = True if i < 90 | units.deg else False
    
    true_anomalies = np.linspace(0, 2 * np.pi, 100)
    xs = []
    ys = []
    for f in true_anomalies:
        _, sec = generate_binaries(m1, m2, a, e,f, i, lan, ap)
        xs.append(sec.x.value_in(units.AU))
        ys.append(sec.y.value_in(units.AU))

    return plt.plot(xs, ys, lw=1), prograde

def plot_bodies(bodies):
    """Plot the x,y positions of the bodies in the simulation.
    Args:
        bodies: particle set containing the bodies
    """
    markers = {'star': 's', 'planet': 'o', 'moon': '.'}

    for body in bodies:
        plt.scatter(body.x.value_in(units.AU), body.y.value_in(units.AU), label=body.name, marker=markers[body.type])
    plt.xlabel('x [AU]')
    plt.ylabel('y [AU]')
    plt.legend()
    # plt.show()

def angular_momentum_deficit(planet, moon, normalize=False):
    '''Calculates the angular momentum deficit of two particles'''
    mu = constants.G * (planet.mass)
    _,_,a,e,_,i,_,_ = orbital_elements(planet, moon)

    AMD = moon.mass*mu.sqrt()*a.sqrt()*(1 - np.sqrt(1 - e**2)*np.cos(i.value_in(units.rad)))
    if normalize:
        AMD = AMD / (moon.mass * mu.sqrt() * a.sqrt())
    return AMD