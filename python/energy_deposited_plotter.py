from podio import root_io
import ROOT
import numpy as np
import math

# input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/electronGun/egun_10GeV_10k_Geant4TrackerAction_edep0.root"
input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
layer_radii = [14, 23, 34.5, 141, 316]
layers = {r: [] for r in layer_radii}

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

# sorts hits by layer
events = [events[i] for i in range(100)]
for event in events:
    for collection in ["VTXIBCollection", "VTXOBCollection"]:
        for hit in event.get(collection):
            layers[radius(hit)].append(hit)

ROOT.gStyle.SetPalette(ROOT.kRainBow)

# creates histogram for each layer
for layer_index in range(5):
    cos_map_hist = ROOT.TH2F("events", f"Guinea Pig Layer {layer_index + 1} Energy Deposited and Cosine Theta", 100, 0, 1, 100, 0, 100)
    hist = ROOT.TH1F("energy", f"Guinea Pig Layer {layer_index + 1} dE/dx", 70, 0, 700)
    phis, zs = [], []

    for hit in layers[layer_radii[layer_index]]:
        edep = 1000000*hit.getEDep()
        path_length = hit.getPathLength()

        # if hit.isProducedBySecondary():
        #     continue
        # tracks location of cells with very low energy deposited
        if edep < 2:
            phis.append(phi(hit.getPosition().x, hit.getPosition().y))
            zs.append(hit.getPosition().z)

        # else:
        hist.Fill(edep/path_length)

        cos_theta = math.cos(theta(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z))
        cos_map_hist.Fill(cos_theta, edep)

    # cos_map_hist.SetTitle(f"Guinea Pig Layer {layer_index + 1} Energy Deposited and Cosine Theta;Cosine Theta;Energy (keV)")
    # cos_map_hist.SetStats(0)
    # cos_map_canvas = ROOT.TCanvas("events", f"Layer {layer_index + 1} Energy Deposited and Cosine Theta")
    # cos_map_hist.Draw("colz")
    # cos_map_canvas.Update()
    # cos_map_canvas.SaveAs(f"../plots/energy_deposited/regular/total_edep/layer{layer_index+1}_edep_and_cos.png")

    hist.SetStats(0)
    hist.GetXaxis().SetTitle("dE/dx (keV/mm)")
    hist.GetYaxis().SetTitle("Number of Events")
    canvas = ROOT.TCanvas("energy", f"Guinea Pig Layer {layer_index + 1} dE/dx")
    hist.Draw()
    canvas.Update()
    canvas.SaveAs(f"../plots/energy_deposited/guinea_pig/total_edep/gp_edep_layer{layer_index + 1}_secondary.png")

    # hist = ROOT.TH2F("locations", "Electron Gun Locations of 0 Energy Deposited", 120, 0, 360, 150, -150, 150)
    # hist.SetTitle(f"Layer {layer_index + 1} Locations of Hits with 0 Energy Deposited;Phi (radians);Z (mm)")
    # for i in range(len(phis)):
    #     phi_deg = phis[i] * (180/math.pi)
    #     hist.Fill(phi_deg, zs[i])
    # hist.SetStats(0)
    # canvas = ROOT.TCanvas("locations", "Locaitons of 0 Energy Deposited")
    # hist.Draw("colz")
    # canvas.Update()
    # canvas.SaveAs(f"../plots/energy_deposited/electron_gun/edep_maps/layer{layer_index+1}_no_edep_locations_eg.png")

hist = ROOT.TH1F("energy", "Guinea Pig Disks dE/dx", 70, 0, 700)
for event in events:
    for hit in event.get("VTXDCollection"):
        edep = 1000000*hit.getEDep()
        path_length = hit.getPathLength()

        # if hit.isProducedBySecondary():
        #     continue
        # if edep < 2:
        #     continue
        # else:
        hist.Fill(edep/path_length)

hist.SetStats(0)
hist.GetXaxis().SetTitle("dE/dx (keV/mm)")
hist.GetYaxis().SetTitle("Number of Events")
canvas = ROOT.TCanvas("energy", "Disk Energy Deposited")
hist.Draw()
canvas.Update()
canvas.SaveAs("../plots/energy_deposited/guinea_pig/total_edep/gp_edep_disks_secondary.png")
