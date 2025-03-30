import numpy as np
import matplotlib.pyplot as plt

from initial_conditions import generate_initial_conditions, add_encounter
from run_simulation import run_simulation
from acs_star_planet_moon import hill_radius

from amuse.units import units, nbody_system
from amuse.community.huayno.interface import Huayno
from amuse.community.hermite.interface import Hermite

#generate host star, planet and moon
m_host = 1 | units.Msun
m_pl = 1 | units.Mjupiter
m_moon = 1.4815e23 | units.kg #mass of Ganymede, use 8.93e22 kg for Io
m_field = 1 | units.Msun

a_sp = 5 | units.AU
a_pm = hill_radius(a_sp, m_host, m_pl) / 3 # 1/3 of the Hill radius

vinf = 3 | units.kms
v20 = 3 | units.kms
b = 2.5 | units.AU

bodies = generate_initial_conditions(m_host, m_pl, a_sp, m_moon)
bodies = add_encounter(bodies, m_field, b, v20, 0, 1, 1)

bodies.color = ['red', 'blue', 'green', 'k']

converter = nbody_system.nbody_to_si(bodies.mass.sum(), a_sp)

# gravity = Huayno(converter)
gravity = Hermite(converter)
gravity.particles.add_particles(bodies)
channel = gravity.particles.new_channel_to(bodies)

initial_energy = gravity.kinetic_energy + gravity.potential_energy
energy_error = [0]

times = np.linspace(0, 1000, 5000) | units.yr
plt.ion()
fig, ax = plt.subplots(2,2, figsize=(15,15))
ax[0,0].set_yscale('log')
ax[0,0].set_xlabel('Time [yr]')
ax[0,0].set_ylabel('Energy error')




for time in times:
    fig.suptitle(f'Time: {time.in_(units.yr)}')
    gravity.evolve_model(time)
    channel.copy()
    energy_error.append(np.abs((gravity.kinetic_energy + gravity.potential_energy - initial_energy) / initial_energy))
    
    ax[0,0].plot(gravity.model_time.value_in(units.yr), energy_error[-1], '.', c='tab:blue')

    positions = bodies.position #- bodies[2].position
    # com = bodies[2].position
    com = [0,0,0] | units.AU
    ax[0,0].set_ylim(1e-10, 1)

    # ax[0,1].set_xlim(com.x.value_in(units.AU)-100, com.x.value_in(units.AU)+100)
    # ax[0,1].set_ylim(com.y.value_in(units.AU)-100, com.y.value_in(units.AU)+100)

    # ax[1,0].set_xlim(com.x.value_in(units.AU)-100, com.x.value_in(units.AU)+100)
    # ax[1,0].set_ylim(com.z.value_in(units.AU)-100, com.z.value_in(units.AU)+100)

    # ax[1,1].set_xlim(com.y.value_in(units.AU)-100, com.y.value_in(units.AU)+100)
    # ax[1,1].set_ylim(com.z.value_in(units.AU)-100, com.z.value_in(units.AU)+100)

    ax[0,1].scatter(positions.x.value_in(units.AU), positions.y.value_in(units.AU), c=bodies.color)

    ax[0,1].set_xlabel('x [AU]')
    ax[0,1].set_ylabel('y [AU]')

    ax[1,0].scatter(positions.x.value_in(units.AU), positions.z.value_in(units.AU), c=bodies.color)
    ax[1,0].set_xlabel('x [AU]')
    ax[1,0].set_ylabel('z [AU]')

    ax[1,1].scatter(positions.y.value_in(units.AU), positions.z.value_in(units.AU), c=bodies.color)
    ax[1,1].set_xlabel('y [AU]')
    ax[1,1].set_ylabel('z [AU]')



    plt.pause(0.001)
    plt.draw()
    ax[0,1].cla()
    ax[1,0].cla()
    ax[1,1].cla()

plt.ioff()
plt.show()
plt.close()
gravity.stop()


