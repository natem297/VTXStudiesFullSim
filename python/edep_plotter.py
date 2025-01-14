from podio import root_io
import ROOT
import numpy as np
import math
import os

##########################################################################################
# this file is for plotting the number of hits over energy deposited (dE/dx) in each layer
##########################################################################################

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/"
files = os.listdir(folder)

# layer_radii = [14, 23, 34.5, 141, 316] # IDEA approximate layer radii
layer_radii = [14, 36, 58] # CLD approximate layer radii
hits = {r: [] for r in layer_radii}

def radius(hit):
    """
    Calculates polar radius of particle.
    Inputs: hit, SimTrackerHit object.
    Output: r, int representing polar radius in mm.
    """
    true_radius = np.sqrt(hit.getPosition().x**2 + hit.getPosition().y**2)
    for r in layer_radii:
        if abs(true_radius-r) < 4:
            return r
    raise ValueError(f"Not close enough to any of the layers {true_radius}")

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

# categorizes barrel hits by radius
for filename in files:

    print(f"starting {filename}")
    input_file_path = os.path.join(folder, filename)
    podio_reader = root_io.Reader(input_file_path)

    events = podio_reader.get("events")
    for event in events:
        for hit in event.get("VertexBarrelCollection"):
            if hit.isProducedBySecondary(): # mc particle not tracked
                continue
            hits[radius(hit)].append(hit)

ROOT.gStyle.SetPalette(ROOT.kRainBow)

# creates histogram for each layer
for layer_index in range(len(layer_radii)):
    hist = ROOT.TH1F("energy", f"Layer {layer_index + 1} dE/dx", 70, 0, 700)

    for hit in hits[layer_radii[layer_index]]:
        edep = 1000000*hit.getEDep() # convert to keV
        path_length = hit.getPathLength() # mm
        mc = hit.getMCParticle()

        if mc.getGeneratorStatus() != 1:
            continue # particles not input into geant

        hist.Fill(edep/path_length)

    hist.SetStats(0)
    hist.GetXaxis().SetTitle("dE/dx (keV/mm)")
    hist.GetYaxis().SetTitle("Number of Events")
    canvas = ROOT.TCanvas("energy", f"Layer {layer_index + 1} dE/dx")
    hist.Draw()
    canvas.Update()
    canvas.SaveAs(f"../plots/cld/edep/layer{layer_index + 1}_edep.png")
