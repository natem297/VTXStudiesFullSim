from podio import root_io
import ROOT
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)
events = podio_reader.get("events")

energies = []

for event in events:
    for particle in event.get("MCParticles"):
        if particle.getGeneratorStatus() != 1:
            continue
        energies.append(1000*particle.getEnergy())

energy_hist = ROOT.TH1F("energy", "Monte Carlo All Particles Energy", 300, 0, 15)
for e in energies:
    energy_hist.Fill(e)
energy_hist.GetXaxis().SetTitle("Energy (MeV)")
energy_hist.GetYaxis().SetTitle("Number of Particles")
energy_canvas = ROOT.TCanvas("energy", "Monte Carlo All Particles Energy")
energy_hist.Draw()
energy_canvas.Update()
energy_canvas.SaveAs("../plots/monte_carlo/guinea_pig_all_energy.png")
