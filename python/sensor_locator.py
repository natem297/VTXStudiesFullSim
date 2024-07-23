from podio import root_io
import ROOT
import math
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)
events = podio_reader.get("events")

layer_radii = [14, 23, 34.5, 141, 316]
disk_z = [303, 635, 945, -303, -635, -945]

ids = {}

def phi(hit):
    """
    Calculates phi of particle.
    Inputs: hit, SimTrackerHit object.
    Output: phi, float representing angle in radians from 0 to 2 pi.
    """
    x = hit.getPosition().x
    y = hit.getPosition().y

    if x == 0:
        if y > 0:
            return math.pi/2
        elif y < 0:
            return -math.pi/2

    phi = math.atan(y/x)

    if x < 0:
        phi +=  math.pi
    elif y < 0:
        phi += 2*math.pi
    return phi

def theta(hit):
    """
    Calculates theta of particle.
    Inputs: hit, SimTrackerHit object.
    Output: theta, float representing angle in radians from 0 to pi.
    """
    x = hit.getPosition().x
    y = hit.getPosition().y
    z = hit.getPosition().z
    return math.acos(z/np.sqrt(x**2 + y**2 + z**2))

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

hits = {coord: [] for coord in layer_radii + disk_z}

for i in range(1):
    event = events[i]
    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):

            if collection != "VTXDCollection":
                hits[radius(hit)].append(hit)
            else:
                hits[z_coord(hit)].append(hit)

for hit in hits[34.5]:
    sensor = hit.getCellID()
    if sensor not in ids:
        ids[sensor] = []
    ids[sensor].append(hit)

sensors = sorted(list(ids.keys()))
for sensor in sensors:
    print(sensor)
    for hit in ids[sensor]:
        print("\t", "z:", hit.getPosition().z, "phi:", phi(hit))
