from podio import root_io
import ROOT
import numpy as np
import math

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)
events = podio_reader.get("events")

layer_radii = [14, 23, 34.5, 141, 316]
disk_z = [303, 635, 945, -303, -635, -945]

def radius(hit):
    """
    Calculates polar radius of particle.
    Inputs: hit, SimTrackerHit object.
    Output: r, number representing polar radius in mm.
    """
    true_radius = np.sqrt(hit.getPosition().x**2 + hit.getPosition().y**2)
    for r in layer_radii:
        if abs(true_radius-r) < 4:
            return r
    raise ValueError(f"Not close enough to any of the layers {np.sqrt(true_radius)}")

def z_coord(hit):
    """
    Calculates z coordinate of hit.
    Inputs: hit, SimTrackerHit object.
    Output: z, int representing z coordinate of hit in mm.
    """
    true_z = hit.getPosition().z
    for z in disk_z:
        if abs(true_z-z) < 30:
            return z
    raise ValueError(f"Not close enough to any of the disks {true_z}")

def delta_squared(particle, other_particle):
    delta_t = particle.getTime() - other_particle.getTime()
    delta_x = particle.getPosition().x - other_particle.getPosition().x
    delta_y = particle.getPosition().y - other_particle.getPosition().y
    delta_z = particle.getPosition().z - other_particle.getPosition().z
    return delta_x**2 + delta_y**2 + delta_z**2

def phi(x,y):
    """
    Calculates phi of particle.
    Inputs: x,y floats.
    Output: phi, float representing angle in radians from 0 to 2 pi.
    """
    phi = math.atan(y/x)
    if x < 0:
        phi +=  math.pi
    elif y < 0:
        phi += 2*math.pi
    return phi

def theta(x,y,z):
    """
    Calculates theta of particle.
    Inputs: x,y,z floats.
    Output: theta, float representing angle in radians from 0 to pi.
    """
    return math.acos(z/np.sqrt(x**2 + y**2 + z**2))

# def delta_squared(hit1, hit2):
#     """
#     Calculates delta of a particle-hit pair.
#     Inputs:
#         hit, SimTracker Hit object.
#         theta_i, phi_i, floats representing direction of mc particle.
#         xi, yi, zi, floats representing position of mc particle.
#     Outputs: delta, float.
#     """

#     xi = hit1.getPosition().x
#     yi = hit1.getPosition().y
#     zi = hit1.getPosition().z

#     xf = hit2.getPosition().x
#     yf = hit2.getPosition().y
#     zf = hit2.getPosition().z

#     delta_x = xf - xi
#     delta_y = yf - yi
#     delta_z = zf - zi

#     theta_i = theta(hit1.getMomentum().x, hit1.getMomentum().y, hit1.getMomentum().z)
#     phi_i = phi(hit1.getMomentum().x, hit1.getMomentum().y)

#     theta_f = theta(delta_x, delta_y, delta_z)
#     phi_f = phi(delta_x, delta_y)

#     delta_theta = theta_f - theta_i
#     delta_phi = min(abs(phi_f - phi_i), 2*math.pi - abs(phi_f - phi_i))

#     return delta_theta**2 + delta_phi**2

for i in range(5):
    event = events[i]
    hits = {coord: [] for coord in layer_radii + disk_z}

    # categorizes all hits by layer
    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):

            if collection != "VTXDCollection":
                hits[radius(hit)].append(hit)
            else:
                hits[z_coord(hit)].append(hit)

    num_matches = 0

    for initial_hit in hits[14]:
        mc = initial_hit.getMCParticle()
        match_hits = []
        for other_hit in hits[14]:
            if other_hit == initial_hit:
                continue
            other_mc = other_hit.getMCParticle()
            if mc == other_mc:
                if delta_squared(initial_hit, other_hit) < 1:
                    match_hits.append(other_hit)
                    num_matches += 1
        if match_hits:
            print(f"Initial hit {initial_hit} matched with {len(match_hits)} hits in the next layer")
            break
    print(f"{num_matches} matches for event {i + 1}")
        # for hit in hits[23]:
        #     if hit.getMCParticle() == mc:
        #         if found_pair:
        #             found_multiple = True
        #         found_pair = True
        #         positions.append((hit.getPosition().x, hit.getPosition().y, hit.getPosition().z))
        # if found_multiple:
        #     print("Found multiple pairs with coordinates:")
        #     for pos in positions:
        #         print(pos)
        #     print(f"From initial hit {initial_hit.getPosition().x, initial_hit.getPosition().y, initial_hit.getPosition().z} \n \
        #           with momentum {initial_hit.getMomentum().x, initial_hit.getMomentum().y, initial_hit.getMomentum().z}")
