from podio import root_io
import os
import math
import numpy as np

##########################################################################################
#  this file is for calculating the average and max occupancies in the first layer
##########################################################################################

layer_radii = [14, 23, 34.5, 141, 316] # IDEA approximate layer radii
# layer_radii = [14, 36, 58] # CLD approximate layer radii

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
    raise ValueError(f"Not close enough to any of the layers {true_radius}")

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/"
files = os.listdir(folder)
file_count = len(files)

cells = {i + j: [0]*file_count for i in range(1, 130, 128) for j in range(0, 122881, 8192)}

# categorizes hits by module id
for file_num in range(file_count):

    filename = files[file_num]
    print(f"starting {filename}")
    input_file_path = os.path.join(folder, filename)
    podio_reader = root_io.Reader(input_file_path)

    events = podio_reader.get("events")
    for event in events:

        for hit in event.get("VertexBarrelCollection"):
            # mc particle not tracked or hit not in first layer
            if hit.isProducedBySecondary() or radius(hit) != 14:
                continue
            cells[hit.getCellID()][file_num] += 1

pix_per_mod = 430080 # IDEA

# averages number of hits in each module over all events
module_hit_averages = [np.mean(cells[mod]) for mod in cells.keys()]
max_hits = max(module_hit_averages)
avg_hits = np.mean(module_hit_averages)

print(f"Maximum occupancy: {max_hits*15}")
print(f"Average occupancy: {avg_hits*15}")
