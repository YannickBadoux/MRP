from amuse.ext.orbital_elements import orbital_elements

def bound(particle1, particle2):
    'Check if two particles are in a bound orbit'
    eccentricity = orbital_elements(particle1, particle2)[3]
    return eccentricity < 1

def find_bound_particles(bodies):
    'Returns a dictionary with the indices of particles that are bound to each particle'
    bound_particles = {}
    for i in range(len(bodies)):
        bound_particles[i] = []
        for j in range(len(bodies)):
            if bound(bodies[i], bodies[j]):
                bound_particles[i].append(j)
    return bound_particles

def find_bound_pairs(bodies):
    'Find all bound pairs in a set of particles, returns a list of tuples of indices'
    bound_pairs = []
    for i in range(len(bodies)):
        for j in range(i+1, len(bodies)):
            if bound(bodies[i], bodies[j]):
                bound_pairs.append((i,j))
    return bound_pairs

def analyse_result(bodies):
    # get how many moons are in the system
    moon_idx = [i for i, body in enumerate(bodies) if 'moon' in body.name]
    moon_count = len(moon_idx)
    #host star: 0, planet: 1, field star: -1, moons everything else

    #get bound particles
    bound_particles = find_bound_particles(bodies)

    if moon_count == 1:
        moon_idx = moon_idx[0]
        if len(bound_particles[1]) == 1 and moon_idx in bound_particles[1]:
            state = 1 # free floating planet with moon
        elif len(bound_particles[1]) == 0:
            state = 0 # free floating planet without moon
        else:
            state = 2 # other
        #TODO: are we interested in the other possible outcomes?

    else:
        #TODO: implement state determination for multiple moons
        raise NotImplementedError("Multiple moons not implemented yet")
    
    return state