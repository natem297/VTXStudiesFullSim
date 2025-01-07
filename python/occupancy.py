from podio import root_io
import os
import math
import numpy as np

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
    # raise ValueError(f"Not close enough to any of the layers {np.sqrt(true_radius)}")

folder = "/eos/experiment/fcc/users/j/jaeyserm/VTXStudiesFullSim/IDEA_guineaPig_andrea_June2024_v23"
files = os.listdir(folder)
file_count = len(files)

cells = {i + j: [0]*file_count for i in range(1, 130, 128) for j in range(0, 122881, 8192)}
secondary_edep = []

for file_num in range(20):

    filename = files[file_num]
    print(f"starting {filename}")
    input_file_path = os.path.join(folder, filename)
    podio_reader = root_io.Reader(input_file_path)

    events = podio_reader.get("events")
    for event in events:

        for hit in event.get("VertexBarrelCollection"):
            if radius(hit) == 14:
                if hit.isProducedBySecondary():
                    secondary_edep.append(hit.getEDep())
            # if hit.isProducedBySecondary() or radius(hit) != 14:
            #     continue
            # cells[hit.getCellID()][file_num] += 1

# pix_per_mod = 430080
# module_hit_averages = [np.mean(cells[mod]) for mod in cells.keys()]
# max_hits = max(module_hit_averages)
# avg_hits = np.mean(module_hit_averages)

# print(f"Maximum occupancy: {max_hits*15}")
# print(f"Average occupancy: {avg_hits*15}")

print(f"All secondary edep: {secondary_edep}")
print(f"Average secondary edep: {np.mean(secondary_edep)}")
print(f"Max edep: {max(secondary_edep)}")
print(f"Min Edep: {min(secondary_edep)}")
