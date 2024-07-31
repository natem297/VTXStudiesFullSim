from podio import root_io
import ROOT
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)
events = podio_reader.get("events")

lepton_pdgs = [11, -11, 12, 13, 14, -14]
particle_energies = {"lepton": [], "hadron": [], "photon": []}
particle_momentums = {"lepton": [], "hadron": [], "photon": []}

for event in events:
    for particle in event.get("MCParticles"):
        pdg = particle.getPDG()
        if pdg in lepton_pdgs:
            p_type = "lepton"
        elif pdg == 22:
            p_type = "photon"
        else:
            p_type = "hadron"
        particle_energies[p_type].append(1000*particle.getEnergy())
        particle_momentums[p_type].append(1000*np.sqrt(particle.getMomentum().x**2 + particle.getMomentum().y**2 + particle.getMomentum().z**2))

# energy_hist = ROOT.TH1F("energy", "Monte Carlo All Particles Energy", 500, 0, 50)
# for energy_list in particle_energies.values():
#     for e in energy_list:
#         energy_hist.Fill(e)
# energy_hist.GetXaxis().SetTitle("Energy (MeV)")
# energy_hist.GetYaxis().SetTitle("Number of Particles")
# energy_canvas = ROOT.TCanvas("energy", "Monte Carlo All Particles Energy")
# energy_hist.Draw()
# energy_canvas.Update()
# energy_canvas.SaveAs("../plots/monte_carlo/guinea_pig_all_energy.png")

# momentum_hist = ROOT.TH1F("momentum", "Monte Carlo All Particles Momentum", 500, 0, 50)
# for momentum_list in particle_momentums.values():
#     for p in momentum_list:
#         momentum_hist.Fill(p)
# momentum_hist.GetXaxis().SetTitle("Momentum (MeV)")
# momentum_hist.GetYaxis().SetTitle("Number of Particles")
# momentum_canvas = ROOT.TCanvas("momentum", "Monte Carlo All Particles Momentum")
# momentum_hist.Draw()
# momentum_canvas.Update()
# momentum_canvas.SaveAs("../plots/monte_carlo/guinea_pig_all_momentum.png")

photon_energy_hist = ROOT.TH1F("Photon Energy", "Monte Carlo Particle Energy (Normalized)", 5000, 0, 1000)
lepton_energy_hist = ROOT.TH1F("Lepton Energy", "Monte Carlo Particle Energy (Normalized)", 5000, 0, 1000)
hadron_energy_hist = ROOT.TH1F("Hadron Energy", "Monte Carlo Particle Energy (Normalized)", 5000, 0, 1000)

for p in particle_energies["photon"]:
    photon_energy_hist.Fill(p)
for p in particle_energies["lepton"]:
    lepton_energy_hist.Fill(p)
for p in particle_energies["hadron"]:
    hadron_energy_hist.Fill(p)

energy_legend = ROOT.TLegend(0.65, 0.75, 0.8, 0.88)
energy_legend.AddEntry(photon_energy_hist, "Photons", "l")
energy_legend.AddEntry(lepton_energy_hist, "Leptons", "l")
energy_legend.AddEntry(hadron_energy_hist, "Hadrons", "l")

photon_energy_hist.Scale(1.0 / photon_energy_hist.Integral())
lepton_energy_hist.Scale(1.0 / lepton_energy_hist.Integral())
hadron_energy_hist.Scale(1.0 / hadron_energy_hist.Integral())

photon_energy_hist.SetLineColor(ROOT.kRed)
lepton_energy_hist.SetLineColor(ROOT.kBlue)
hadron_energy_hist.SetLineColor(ROOT.kGreen)

photon_energy_hist.GetXaxis().SetTitle("Energy (MeV)")
photon_energy_hist.GetYaxis().SetTitle("Proportion of Particles")
photon_energy_hist.SetStats(0)
photon_energy_hist.SetMaximum(0.16)

canvas = ROOT.TCanvas("Energy", "Monte Carlo Particle Energies")
photon_energy_hist.Draw("hist l")
lepton_energy_hist.Draw("hist l same")
hadron_energy_hist.Draw("hist l same")
energy_legend.Draw()
canvas.Update()
canvas.SaveAs("../plots/monte_carlo/guinea_pig_energy_combined.png")

photon_momentum_hist = ROOT.TH1F("Photon Momentum", "Monte Carlo Particle Momentum (Normalized)", 1000, 0, 100)
lepton_momentum_hist = ROOT.TH1F("Lepton Momentum", "Monte Carlo Particle Momentum (Normalized)", 1000, 0, 100)
hadron_momentum_hist = ROOT.TH1F("Hadron Momentum", "Monte Carlo Particle Momentum (Normalized)", 1000, 0, 100)

for p in particle_momentums["photon"]:
    photon_momentum_hist.Fill(p)
for p in particle_momentums["lepton"]:
    lepton_momentum_hist.Fill(p)
for p in particle_momentums["hadron"]:
    hadron_momentum_hist.Fill(p)

momentum_legend = ROOT.TLegend(0.65, 0.75, 0.8, 0.88)
momentum_legend.AddEntry(photon_momentum_hist, "Photons", "l")
momentum_legend.AddEntry(lepton_momentum_hist, "Leptons", "l")
momentum_legend.AddEntry(hadron_momentum_hist, "Hadrons", "l")

photon_momentum_hist.Scale(1.0 / photon_momentum_hist.Integral())
lepton_momentum_hist.Scale(1.0 / lepton_momentum_hist.Integral())
hadron_momentum_hist.Scale(1.0 / hadron_momentum_hist.Integral())

photon_momentum_hist.SetLineColor(ROOT.kRed)
lepton_momentum_hist.SetLineColor(ROOT.kBlue)
hadron_momentum_hist.SetLineColor(ROOT.kGreen)

photon_momentum_hist.GetXaxis().SetTitle("Momentum (MeV)")
photon_momentum_hist.GetYaxis().SetTitle("Proportion of Particles")
photon_momentum_hist.SetStats(0)

canvas = ROOT.TCanvas("Momentum", "Monte Carlo Particle Momentums")
photon_momentum_hist.Draw("hist l")
lepton_momentum_hist.Draw("hist l same")
hadron_momentum_hist.Draw("hist l same")
momentum_legend.Draw()
canvas.Update()
canvas.SaveAs("../plots/monte_carlo/guinea_pig_momentum_combined.png")

# for p_type in ["lepton", "hadron", "photon"]:

#     if p_type == "lepton" or p_type == "photon":
#         energy_hist = ROOT.TH1F("energy", F"Monte Carlo {p_type.capitalize()} Energy", 350, 0, 35)
#     else:
#         energy_hist = ROOT.TH1F("energy", F"Monte Carlo {p_type.capitalize()} Energy", 400, 930, 970)
#     for e in particle_energies[p_type]:
#         energy_hist.Fill(e)
#     energy_hist.GetXaxis().SetTitle("Energy (MeV)")
#     energy_hist.GetYaxis().SetTitle("Number of Particles")
#     energy_canvas = ROOT.TCanvas("energy", f"Monte Carlo {p_type.capitalize()} Energy")
#     energy_hist.Draw()
#     energy_canvas.Update()
#     energy_canvas.SaveAs(f"../plots/monte_carlo/guinea_pig_{p_type}_energy.png")

#     if p_type == "lepton" or p_type == "photon":
#         momentum_hist = ROOT.TH1F("momentum", f"Monte Carlo {p_type.capitalize()} Momentum", 350, 0, 35)
#     else:
#         momentum_hist = ROOT.TH1F("momentum", f"Monte Carlo {p_type.capitalize()} Momentum", 400, 0, 200)
#     for p in particle_momentums[p_type]:
#         momentum_hist.Fill(p)
#     momentum_hist.GetXaxis().SetTitle("Momentum (MeV)")
#     momentum_hist.GetYaxis().SetTitle("Number of Particles")
#     momentum_canvas = ROOT.TCanvas("momentum", f"Monte Carlo {p_type.capitalize()} Momentum")
#     momentum_hist.Draw()
#     momentum_canvas.Update()
#     momentum_canvas.SaveAs(f"../plots/monte_carlo/guinea_pig_{p_type}_momentum.png")
