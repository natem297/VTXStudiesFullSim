from podio import root_io
import ROOT
import math
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/electronGun/egun_10GeV_10k_Geant4TrackerAction_edep0.root"
# input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

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

theta_data = ROOT.TH1F("Theta", "Electron Initial Theta", 60, 0, 180)
theta_data.GetXaxis().SetTitle("Theta (deg)")
theta_data.GetYaxis().SetTitle("Number of Events")
theta_data.SetMinimum(0)

phi_data = ROOT.TH1F("Phi", "Electron Initial Phi", 120, 0, 360)
phi_data.GetXaxis().SetTitle("Phi (deg)")
phi_data.GetYaxis().SetTitle("Number of Events")
phi_data.SetMinimum(0)

for event in podio_reader.get("events"):
    if len(event.get("MCParticles")) > 100:
        continue
    particle = event.get("MCParticles")[0]

    px = particle.getMomentum().x
    py = particle.getMomentum().y
    pz = particle.getMomentum().z
    # calculates theta and phi based on momentum data
    theta_i = theta(px, py, pz)*(180/math.pi)
    phi_i = phi(px, py)*(180/math.pi)

    theta_data.Fill(theta_i)
    phi_data.Fill(phi_i)

canvas_theta = ROOT.TCanvas("Theta", "Electron Initial Theta")
theta_data.Draw()
canvas_theta.SaveAs("../plots/monte_carlo/mc_initial_theta_eg.png")

canvas_phi = ROOT.TCanvas("Phi", "Electron Initial Phi")
phi_data.Draw()
canvas_phi.SaveAs("../plots/monte_carlo/mc_initial_phi_eq.png")
