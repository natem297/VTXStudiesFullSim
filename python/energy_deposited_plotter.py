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

no_edep_pdgs = {11: 0, -11: 0, 22: 0}
good_cell_ids = set()
num_no_edep_in_good_cells = 0

ROOT.gStyle.SetPalette(ROOT.kRainBow)

# creates histogram for each layer
for layer_index in range(5):
    cos_map_hist = ROOT.TH2F("events", f"Guinea Pig Layer {layer_index + 1} Energy Deposited and Cosine Theta", 100, 0, 1, 100, 0, 100)
    lepton_hist = ROOT.TH1F("lepton energy", f"Guinea Pig Layer {layer_index + 1} Lepton Energy Deposited", 50, 0, 50)
    non_lepton_hist = ROOT.TH1F("non lepton energy", f"Guinea Pig Layer {layer_index + 1} Non Lepton Energy Deposited", 50, 0, 50)
    phis, zs = [], []
    cos_thetas = [[] for _ in range(100)]

    for particle in layers[layer_radii[layer_index]]:
        edep = 1000000*particle.getEDep()
        mc = particle.getMCParticle()
        pdg = mc.getPDG()

        # tracks location of cells with very low energy deposited
        if edep < 2:
            mc = particle.getMCParticle()
            phis.append(phi(particle.getPosition().x, particle.getPosition().y))
            zs.append(particle.getPosition().z)

        # else:
        if pdg == 11 or pdg == -11:
            lepton_hist.Fill(edep)
        else:
            non_lepton_hist.Fill(edep)
        cos_theta = math.cos(theta(particle.getPosition().x, particle.getPosition().y, particle.getPosition().z))
        cos_map_hist.Fill(cos_theta, edep)
        cos_thetas[int(cos_theta//0.01)].append(edep)

    # cos_map_hist.SetTitle(f"Guinea Pig Layer {layer_index + 1} Energy Deposited and Cosine Theta;Cosine Theta;Energy (keV)")
    # cos_map_hist.SetStats(0)
    # cos_map_canvas = ROOT.TCanvas("events", f"Layer {layer_index + 1} Energy Deposited and Cosine Theta")
    # cos_map_hist.Draw("colz")
    # cos_map_canvas.Update()
    # cos_map_canvas.SaveAs(f"../plots/energy_deposited/regular/total_edep/layer{layer_index+1}_edep_and_cos.png")

    # cos_hist = ROOT.TH1F("energy", f"Electron Gun Layer {layer_index + 1} Energy Deposited over Cosine Theta", 100, 0, 1)
    # for j in range(100):
    #     if not cos_thetas[j]:
    #         cos_hist.SetBinContent(j + 1, 0)
    #     else:
    #         cos_hist.SetBinContent(j + 1, np.mean(cos_thetas[j]))
    # cos_hist.GetXaxis().SetTitle("Cosine Theta")
    # cos_hist.GetYaxis().SetTitle("Energy (keV)")
    # cos_hist.SetStats(0)
    # cos_canvas = ROOT.TCanvas("energy", f"Layer {layer_index + 1} Energy Deposited over Cosine Theta")
    # cos_hist.Draw()
    # cos_canvas.Update()
    # cos_canvas.SaveAs(f"../plots/energy_deposited/regular/total_edep/edep_layer{layer_index + 1}_cos.png")

    lepton_hist.SetStats(0)
    lepton_hist.GetXaxis().SetTitle("Energy (keV)")
    lepton_hist.GetYaxis().SetTitle("Number of Events")
    canvas = ROOT.TCanvas("energy", f"Layer {layer_index + 1} Energy Deposited")
    lepton_hist.Draw()
    canvas.Update()
    canvas.SaveAs(f"../plots/energy_deposited/guinea_pig/total_edep/lepton_energy_deposited_layer{layer_index + 1}_eg.png")

    non_lepton_hist.SetStats(0)
    non_lepton_hist.GetXaxis().SetTitle("Energy (keV)")
    non_lepton_hist.GetYaxis().SetTitle("Number of Events")
    canvas = ROOT.TCanvas("energy", f"Layer {layer_index + 1} Energy Deposited")
    non_lepton_hist.Draw()
    canvas.Update()
    canvas.SaveAs(f"../plots/energy_deposited/guinea_pig/total_edep/nonlepton_energy_deposited_layer{layer_index + 1}_eg.png")

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

# hist = ROOT.TH1F("energy", "Electron Gun Disks Energy Deposited", 50, 0, 50)
# for event in events:
#     for hit in event.get("VTXDCollection"):
#         edep = 1000000*hit.getEDep()
#         if edep < 0.000002:
#             continue
#         else:
#             hist.Fill(edep)

# hist.SetStats(0)
# hist.GetXaxis().SetTitle("Energy (keV)")
# hist.GetYaxis().SetTitle("Number of Events")
# canvas = ROOT.TCanvas("energy", "Disk Energy Deposited")
# hist.Draw()
# canvas.Update()
# canvas.SaveAs("../plots/energy_deposited/electron_gun/total_edep/energy_deposited_disks_zoom_eg.png")

print(no_edep_pdgs)
