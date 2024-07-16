from podio import root_io
import ROOT
import math
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/electronGun/egun_10GeV_10k_Geant4TrackerAction_edep0.root"
# input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
e0 = events[0]
layer_radii = [14, 23, 34.5, 141, 316]
layers = {r: [] for r in layer_radii}
z_ranges = {0: (96,4), 1: (162,6), 2: (261,9), 3: (162,6), 4: (324,12)}

cell_hits = {}

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

def phi(hit):
    """
    Calculates phi of particle.
    Inputs: hit, SimTrackerHit object.
    Output: phi, float representing angle in radians from 0 to 2 pi.
    """
    x = hit.getPosition().x
    y = hit.getPosition().y
    phi = math.atan(y/x)
    if x < 0:
        phi +=  math.pi
    elif y < 0:
        phi += 2*math.pi
    return phi * (180/math.pi)

# sorts hits by layer
events = [events[i] for i in range(10000)]
for event in events:
    for collection in ["VTXIBCollection", "VTXOBCollection"]:
        for hit in event.get(collection):
            layers[radius(hit)].append(hit)

ROOT.gStyle.SetPalette(ROOT.kRainbow)

for layer_index in range(5):

    max_z = None
    z_range = z_ranges[layer_index][0]
    divider = z_ranges[layer_index][1]
    coords = [(deg, z) for deg in range(120) for z in range((2*z_range)//divider)]
    edep_map = {coord: 0 for coord in coords}

    for hit in layers[layer_radii[layer_index]]:

        azimuthal = int(phi(hit) // 3)
        z = int((hit.getPosition().z + z_range) // divider)

        # finds max z if one is out of bounds
        if z < 0 or z >= (2*z_range)//divider:

            if not max_z or abs(hit.getPosition().z) > max_z:
                max_z = abs(hit.getPosition().z)
            continue
        # adds energy deposited to coordinate
        edep_map[(azimuthal, z)] += 1000*hit.getEDep()

    if max_z:
        print(f"Z outside bounds in layer {layer_index + 1}: {max_z}")

    hist = ROOT.TH2F("Energy Deposited", f"Layer {layer_index + 1} Energy Deposited", 120, 0, 360, (2*z_range)//divider, -z_range, z_range)
    hist.SetTitle(f"Electron Gun Layer {layer_index + 1} Energy Deposited (MeV);Phi (radians);Z (mm)")
    for coord in coords:
        hist.SetBinContent(coord[0]+1, coord[1]+1, edep_map[coord])
    hist.SetStats(0)
    canvas = ROOT.TCanvas("Energy Deposited", f"Layer {layer_index + 1} Energy Deposited")
    canvas.SetRightMargin(0.15)
    hist.Draw("colz")
    canvas.Update()
    canvas.SaveAs(f"../plots/energy_deposited/electron_gun/edep_maps/layer{layer_index+1}_edep_map_eg.png")
