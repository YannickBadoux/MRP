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

def run_simulation(bodies, plot=False, integrator='hermite'):
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

    times = np.linspace(0, 10, 100) | units.yr

    initial_energy = gravity.kinetic_energy + gravity.potential_energy
    energy_error = []

    for time in tqdm(times):
        gravity.evolve_model(time)
        error = np.abs((gravity.kinetic_energy + gravity.potential_energy - initial_energy) / initial_energy)
        # print(error)
        energy_error.append(error)
        channel.copy()

        if plot:
            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(111, projection='3d')
            for body in bodies:
                ax.scatter(body.x.value_in(units.AU), body.y.value_in(units.AU), body.z.value_in(units.AU),
                        label=body.name)
                # ax.quiver(body.x.value_in(units.AU), body.y.value_in(units.AU), body.z.value_in(units.AU),
                #             body.vx.value_in(units.kms)/10, body.vy.value_in(units.kms)/10, body.vz.value_in(units.kms)/10,
                #             color=ax._get_lines.get_next_color())
            ax.legend(loc='upper left')
            ax.set_xlabel('x [AU]')
            ax.set_ylabel('y [AU]')
            ax.set_zlabel('z [AU]')

            ax.set_xlim(-10,10)
            ax.set_ylim(-10,10)
            ax.set_zlim(-1,5)
            fig.savefig(f'../movie/simulation{int(time.value_in(units.day))}.png', dpi=200)
            plt.close()

    gravity.stop()

    return bodies, energy_error
