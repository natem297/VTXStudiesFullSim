from podio import root_io
import ROOT
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
layer_radii = [14, 23, 35, 141, 316]

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
    raise ValueError(f"Not close enough to any of the layers {np.sqrt(true_radius)}")

hist = ROOT.TH1F("momentum", "Guinea Pig MC Particles with Layer 1 Match", 150, 0, 30)
for e in range(100):
    print(f"starting event {e}")
    visited_mc = []
    event = events[e]
    for hit in event.get("VTXIBCollection"):
        if radius(hit) == 14:
            mc = hit.getMCParticle()
            if mc in visited_mc:
                continue
            px = mc.getMomentum().x
            py = mc.getMomentum().y
            pz = mc.getMomentum().z
            hist.Fill(1000*np.sqrt(px**2 + py**2 + pz**2))

hist.GetXaxis().SetTitle("Momentum (MeV)")
hist.GetYaxis().SetTitle("Number of Particles")
canvas = ROOT.TCanvas("momentum", "MC Particle Momentum")
hist.Draw()
canvas.Update()
canvas.SaveAs("../plots/monte_carlo/gp_mc_momentum_layer1_match.png")
