import uproot
import numpy as np
import math
import ROOT

input_file_path = "/eos/user/j/jaeyserm/public/mumuH_rec_16478_29.root"
file = uproot.open(input_file_path)

events = file["events"]
vtx_barrel = events[209]

radii = [14, 36, 58]

x_data = vtx_barrel["VertexBarrelCollection.position.x"].array()
y_data = vtx_barrel["VertexBarrelCollection.position.y"].array()
z_data = vtx_barrel["VertexBarrelCollection.position.z"].array()

px_data = vtx_barrel["VertexBarrelCollection.momentum.x"].array()
py_data = vtx_barrel["VertexBarrelCollection.momentum.y"].array()
pz_data = vtx_barrel["VertexBarrelCollection.momentum.z"].array()

cellid_data = vtx_barrel["VertexBarrelCollection.cellID"].array()
edep_data = vtx_barrel["VertexBarrelCollection.EDep"].array()
time_data = vtx_barrel["VertexBarrelCollection.time"].array()
pathlength_data = vtx_barrel["VertexBarrelCollection.pathLength"].array()

def radius(r):
    for true_r in radii:
        if abs(r - true_r) < 5:
            return true_r
    raise ValueError(f"Not close enough to any of the layers {np.sqrt(true_r)}")

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

hit_map = {r: {z: {azimuthal: [0]*1000 for azimuthal in range(0, 360, 3)} for z in range(-110, 110, 2)} for r in radii}
hits = {r: [] for r in radii}
for i in range(1000):

    hits = {coord: [] for coord in radii}
    event_x = x_data[i]
    event_y = y_data[i]
    event_z = z_data[i]

    for j in range(len(event_x)):
        x = event_x[j]
        y = event_y[j]
        r = np.sqrt(x**2 + y**2)
        hits[radius(r)].append(j)
