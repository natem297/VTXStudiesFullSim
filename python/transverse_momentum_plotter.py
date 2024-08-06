from podio import root_io
import ROOT
import numpy as np
import math

# input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/electronGun/egun_10GeV_10k_Geant4TrackerAction_edep0.root"
input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
collection = ["MCParticles", "Monte Carlo", "mc"]
# collection = ["VTXIBCollection", "Inner Barrel", "ib"]
# collection = ["VTXOBCollection", "Outer Barrel", "ob"]
# collection = ["VTXDCollection", "Disk", "disk"]

pt_gt_5_count = 0

pt = ROOT.TH1F("pt", f"Transvere Momentum {collection[1]}", 100, 0, 10)
for event in events:
    # if collection[2] == "mc":
    #     particles = [event.get(collection[0])[0]]
    # else:
    particles = event.get(collection[0])

    # if len(event.get("MCParticles")) > 100:
    #     continue

    for particle in particles:
        if particle.getGeneratorStatus() != 1:
            continue
        p = np.sqrt(particle.getMomentum().x**2 + particle.getMomentum().y**2)

        # if p < 0.1: # eliminate non electron hits
        #     continue
        if p > 5:
            pt_gt_5_count += 1

        pt.Fill(1000*p)

pt.SetXTitle("Transverse Momentum (MeV)")
pt.SetYTitle("Number of Particles")
# pt.SetStats(0)
canvas = ROOT.TCanvas("pt", f"Transverse Momentum ({collection[1]})")
pt.Draw("hist")
canvas.Update()
canvas.SaveAs(f"../plots/transverse_momentum/regular/pt_{collection[2]}.png")

print(f"Total number of particles with pt > 5 GeV: {pt_gt_5_count}")
