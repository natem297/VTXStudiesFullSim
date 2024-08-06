from podio import root_io
import ROOT
import math
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)
events = podio_reader.get("events")

layer_radii = [14, 23, 34.5, 141, 316]

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

cells = {i + j: [0]*100 for i in range(1,81922,16384) for j in range(0,3585,256)}
r = 13.7
for i in range(100):
    print(f"starting event {i}")
    event = events[i]
    for hit in event.get("VTXIBCollection"):
        if radius(hit) == 14:
            cells[hit.getCellID()][i] += 1

pix_per_mod = 430080
module_hit_averages = [np.mean(cells[mod]) for mod in cells.keys()]
max_hits = max(module_hit_averages)
avg_hits = np.mean(module_hit_averages)

print(f"Maximum occupancy: {max_hits*15/pix_per_mod}")
print(f"Average occupancy: {avg_hits*15/pix_per_mod}")
