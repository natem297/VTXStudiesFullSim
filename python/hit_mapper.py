from podio import root_io
import ROOT
import math
import numpy as np
import os

##########################################################################################
# this file is for plotting the number of hits in a 2D map of phi and z and purely as a
# function of phi and theta
##########################################################################################

# layer_radii = [14, 23, 34.5, 141, 316] # IDEA approximate layer radii
# max_z = 96 # IDEA first layer

layer_radii = [14, 36, 58] # CLD approximate layer radii
max_z = 110 # CLD first layer

z_step = 2

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

def radius(hit):
    """
    Calculates polar radius of particle.
    Inputs: hit, SimTrackerHit object.
    Output: r, number representing polar radius in mm.
    """
    true_radius = np.sqrt(hit.getPosition().x**2 + hit.getPosition().y**2)
    for r in layer_radii:
        if abs(true_radius-r) < 3:
            return r
    raise ValueError(f"Not close enough to any of the layers {true_radius}")

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/IDEA_guineaPig_andrea_June2024_v23/"
files = os.listdir(folder)
event_count = 100 * len(files)

hit_map = {z: {azimuthal: 0 for azimuthal in range(0, 360, 3)} \
                for z in range(-max_z, max_z, z_step)}
hit_map_phi = {azimuthal: 0 for azimuthal in range(0, 360, 3)}
hit_map_theta = {polar: 0 for polar in range(0, 180, 3)}

for filename in files:

    print(f"starting {filename}")
    input_file_path = os.path.join(folder, filename)
    podio_reader = root_io.Reader(input_file_path)

    events = podio_reader.get("events")
    for event in events:
        # iterates through hits and adds position to each map
        for hit in event.get("VertexBarrelCollection"):
            # mc particle not tracked or hit not in first layer
            if hit.isProducedBySecondary() or radius(hit) != 14:
                continue

            ph = phi(hit.getPosition().x, hit.getPosition().y) * (180 / math.pi)
            th = theta(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z) * (180 / math.pi)
            z = hit.getPosition().z

            if abs(z) > max_z: # hit out of map range
                continue

            hit_map[int((z // z_step) * z_step)][int((ph // 3) * 3)] += 1
            hit_map_phi[int(ph // 3) * 3] += 1
            hit_map_theta[int(th // 3) * 3] += 1

# change by event type
save_location = "../plots/cld/z_nominal/cld_z_nominal"
title = "CLD Nominal Z -> Hadrons Layer 1 Hits/BX"

# 2D hit map
hist = ROOT.TH2F("hit map", f"Layer 1 Hits", (2 * max_z) // z_step, -max_z, max_z, 120, 0, 360)
hist.SetTitle(f"{title};z (mm); Azimuthal Angle (deg)")

for z in range(-max_z, max_z, z_step):
    for azimuthal in range(0, 360, 3):
        hits = hit_map[z][azimuthal] / event_count
        hist.SetBinContent(((z + max_z) // z_step) + 1, (azimuthal // 3) + 1, hits)

hist.SetStats(0)
canvas = ROOT.TCanvas("hit map", "Layer 1 Hits")
canvas.SetRightMargin(0.12)
hist.Draw("colz")
canvas.Update()
canvas.SaveAs(f"{save_location}_hit_map.png")

# phi histogram
hist_phi = ROOT.TH1F("hit map", f"{title}", 120, 0, 360)
hist_phi.GetXaxis().SetTitle("Phi (deg)")
hist_phi.GetYaxis().SetTitle("Average Number of Hits")
hist_phi.SetMinimum(0)

for ph in range(0, 360, 3):
    hits = hit_map_phi[ph] / event_count
    hist_phi.SetBinContent((ph // 3) + 1, hits)

hist_phi.SetStats(0)
canvas = ROOT.TCanvas("hit map", "Layer 1 Hits")
hist_phi.Draw()
canvas.Update()
canvas.SaveAs(f"{save_location}_phi.png")

# theta histogram
hist_theta = ROOT.TH1F("hit map", f"{title}", 60, 0, 180)
hist_theta.GetXaxis().SetTitle("Theta (deg)")
hist_theta.GetYaxis().SetTitle("Average Number of Hits")
hist_theta.SetMinimum(0)

for th in range(0, 180, 3):
    hits = hit_map_theta[th] / event_count
    hist_theta.SetBinContent((th // 3) + 1, hits)

hist_theta.SetStats(0)
canvas = ROOT.TCanvas("hit map", "Layer 1 Hits")
hist_theta.Draw()
canvas.Update()
canvas.SaveAs(f"{save_location}_theta.png")
