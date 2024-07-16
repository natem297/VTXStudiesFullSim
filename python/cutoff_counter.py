from podio import root_io
import ROOT
import math
import numpy as np
import csv

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/electronGun/egun_10GeV_10k_Geant4TrackerAction_edep0.root"
# input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
layer_radii = [14, 23, 34.5, 141, 316]
disk_z = [-945, -635, -303, 303, 635, 945]

def radius(hit):
    """
    Calculates polar radius of particle.
    Inputs: hit, SimTrackerHit object.
    Output: r, float representing polar radius in mm.
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

def phi(x,y):
    """
    Calculates phi of particle.
    Inputs: x,y floats.
    Output: phi, float representing angle in radians from 0 to 2 pi.
    """
    phi = math.atan(y/x)
    if x < 0:
        phi += math.pi
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

def delta_squared(theta_i, phi_i, xi, yi, zi, hit):
    """
    Calculates delta of a particle-hit pair.
    Inputs:
        hit, SimTracker Hit object.
        theta_i, phi_i, floats representing direction of mc particle.
        xi, yi, zi, floats representing position of mc particle.
    Outputs: delta, float.
    """

    xf = hit.getPosition().x
    yf = hit.getPosition().y
    zf = hit.getPosition().z

    delta_x = xf - xi
    delta_y = yf - yi
    delta_z = zf - zi

    theta_f = theta(delta_x, delta_y, delta_z)
    phi_f = phi(delta_x, delta_y)

    delta_theta = theta_f - theta_i
    delta_phi = min(abs(phi_f - phi_i), 2*math.pi - abs(phi_f - phi_i))

    return delta_theta**2 + delta_phi**2

num_explored = 0
all_matches = []
no_hits = 0

for i in range(10000):
    event = events[i]
    layers = {coord: [] for coord in layer_radii + disk_z}
    # sorts hits by layer
    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):
            if collection != "VTXDCollection":
                layers[radius(hit)].append(hit)
            else:
                layers[z_coord(hit)].append(hit)

    mc = event.get("MCParticles")
    if len(mc) > 100:
        continue
    num_matches = 0
    num_explored += 1

    particle = event.get("MCParticles")[0]

    xi = particle.getVertex().x
    yi = particle.getVertex().y
    zi = particle.getVertex().z

    px = particle.getMomentum().x
    py = particle.getMomentum().y
    pz = particle.getMomentum().z

    theta_i = theta(px, py, pz)
    phi_i =  phi(px, py)

    for coord in layers:
        for hit in layers[coord]:

            delta = delta_squared(theta_i, phi_i, xi, yi, zi, hit)
            if delta < 0.001:
                num_matches += 1
                break

    all_matches.append(num_matches)

matches_at_least_1 = 0
for match_count in all_matches:
    if match_count >= 1:
        matches_at_least_1 += 1

print(f"Total number of mc particles explored: {num_explored}")
print(f"Percentage of events with hits with match found: {100 * (matches_at_least_1/num_explored)}")
print(f"Total number of matches found: {sum(all_matches)}")
print(f"Ratio of number of matches found to number of initial hits: {sum(all_matches)/num_explored}")
