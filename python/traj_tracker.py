from podio import root_io
import ROOT
import numpy as np
import math
import os

# IDEA
# layer_radii = [14, 23, 34.5, 141, 316] # approximate layer radii

# comp_index_dict = {0: "ib", 1: "ib", 2: "ib", 3: "ob", 4: "ob", \
#                     5: "ld", 6: "ld", 7: "ld", 8: "rd", 9: "rd", 10: "rd"}

# CLD
layer_radii = [14, 36, 58] # approximate layer radii
disk_z = [-300, -230, -160, 160, 230, 300] # approximate z values of disks
comp_index_dict = {0: "ib", 1: "ib", 2: "ib", 3: "ld", 4: "ld", 5: "ld", \
                    6: "rd", 7: "rd", 8: "rd"}

component_coords = layer_radii + disk_z

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

def z_coord(hit):
    """
    Calculates z coordinate of hit.
    Inputs: hit, SimTrackerHit object.
    Output: z, int representing z coordinate of hit in mm.
    """
    true_z = hit.getPosition().z
    for z in disk_z:
        if abs(true_z-z) < 40:
            return z
    raise ValueError(f"Not close enough to any of the disks {true_z}")

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

def trajectory_length(start_comp, mc_particle):
    """
    Finds trajectory length of a given particle by finding hits with the same MC particle.
    Inputs:
        start_comp: int, representing index of component where starting hit is located.
        mc_particle: MCParticle, particle being tracked.
    Outputs:
        traj_lengths: dict, mapping each subdetector (ib, ob, ld, rd) to the number of hits
            in that subdetector.
    """
    # traj_lengths = {"total": 0, "ld": 1, "rd": 2, "ib": 3, "ob": 4} # IDEA
    traj_lengths = {"total": 0, "ld": 1, "rd": 2, "ib": 3} # CLD
    traj_lengths[comp_index_dict[start_comp]] += 1 # adds 1 to starting subdetector
    # searches through remaining components for mc match
    for comp_index in range(start_comp + 1, len(component_coords)):
        for hit in hits[component_coords[comp_index]]:
            if hit.getMCParticle() == mc_particle:
                traj_lengths[comp_index_dict[comp_index]] += 1
                traj_lengths["total"] += 1
                break

    return traj_lengths

thetas = {i: [] for i in range(0,180,3)}
phis = {j: [] for j in range(0,360,3)}

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/"
files = os.listdir(folder)
file_count = len(files)

for filename in files:

    print(f"starting {filename}")
    input_file_path = os.path.join(folder, filename)
    podio_reader = root_io.Reader(input_file_path)
    events = podio_reader.get("events")

    # sorts hits by layer and disk
    for event in events:
        hits = {coord: [] for coord in component_coords}

        for collection in ["VertexBarrelCollection", "VertexEndcapCollection"]:
            for hit in event.get(collection):
                if hit.isProducedBySecondary(): # mc particle not tracked
                    continue

                if collection != "VertexEndcapCollection":
                    hits[radius(hit)].append(hit)
                else:
                    hits[z_coord(hit)].append(hit)

        visited_mc = []

        for comp_index in range(len(component_coords)):
            for hit in hits[component_coords[comp_index]]:

                mc = hit.getMCParticle()
                if mc in visited_mc or mc.getGeneratorStatus() != 1:
                    continue # particle already tracked or not input into geant
                visited_mc.append(mc)

                polar = theta(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z) * (180 / math.pi)
                polar = int(polar // 3) * 3 # rounds down to nearest multiple of 3

                azimuthal = phi(hit.getPosition().x, hit.getPosition().y) * (180 / math.pi)
                azimuthal = int(azimuthal // 3) * 3 # rounds down to nearest multiple of 3

                traj_lengths = trajectory_length(comp_index, mc)

                thetas[polar].append(list(traj_lengths.values()))
                phis[azimuthal].append(list(traj_lengths.values()))

# title = "IDEA Components Crossed"
title = "CLD Components Crossed"

# creates theta hist objects for each subdetector and plots them
thetas_hist_total = ROOT.TH1F("Total", title, 60, 0, 180)
thetas_hist_ld = ROOT.TH1F("Left Disks", title, 60, 0, 180)
thetas_hist_rd = ROOT.TH1F("Right Disks", title, 60, 0, 180)
thetas_hist_ib = ROOT.TH1F("Inner Barrel", title, 60, 0, 180)
# thetas_hist_ob = ROOT.TH1F("Outer Barrel", title, 60, 0, 180)

for i in range(0,180,3):
    if not thetas[i]: # for theta values with no hits
        thetas_hist_total.SetBinContent((i//3) + 1, 0)
        thetas_hist_ld.SetBinContent((i//3) + 1, 0)
        thetas_hist_rd.SetBinContent((i//3) + 1, 0)
        thetas_hist_ib.SetBinContent((i//3) + 1, 0)
        # thetas_hist_ob.SetBinContent((i//3) + 1, 0)
    else:
        thetas_hist_total.SetBinContent((i//3) + 1, np.mean([lengths[0] for lengths in thetas[i]]))
        thetas_hist_ld.SetBinContent((i//3) + 1, np.mean([lengths[1] for lengths in thetas[i]]))
        thetas_hist_rd.SetBinContent((i//3) + 1, np.mean([lengths[2] for lengths in thetas[i]]))
        thetas_hist_ib.SetBinContent((i//3) + 1, np.mean([lengths[3] for lengths in thetas[i]]))
        # thetas_hist_ob.SetBinContent((i//3) + 1, np.mean([lengths[4] for lengths in thetas[i]]))

thetas_hist_total.SetXTitle("Polar Angle (deg)")
thetas_hist_total.SetYTitle("Average Number of Components Crossed")
thetas_hist_total.SetStats(0)

thetas_hist_total.SetLineColor(ROOT.kBlack)
thetas_hist_ld.SetLineColor(ROOT.kGreen)
thetas_hist_rd.SetLineColor(ROOT.kYellow)
thetas_hist_ib.SetLineColor(ROOT.kRed)
# thetas_hist_ob.SetLineColor(ROOT.kBlue)

thetas_legend = ROOT.TLegend(0.4, 0.62, 0.6, 0.75)
thetas_legend.AddEntry(thetas_hist_total, "Total", "l")
thetas_legend.AddEntry(thetas_hist_ld, "Left Disks", "l")
thetas_legend.AddEntry(thetas_hist_rd, "Right Disks", "l")
thetas_legend.AddEntry(thetas_hist_ib, "Inner Barrel", "l")
# thetas_legend.AddEntry(thetas_hist_ob, "Outer Barrel", "l")

thetas_canvas = ROOT.TCanvas("Theta Hits", title)
thetas_hist_total.Draw("hist")
thetas_hist_ld.Draw("same")
thetas_hist_rd.Draw("same")
thetas_hist_ib.Draw("same")
# thetas_hist_ob.Draw("same")
thetas_legend.Draw()
thetas_canvas.Update()
thetas_canvas.SaveAs("../plots/cld/cld_theta_hits_gp_combined.png")

# creates phi hist objects for each subdetector and plots them
phis_hist_total = ROOT.TH1F("Total", title, 120, 0, 360)
phis_hist_ld = ROOT.TH1F("Left Disks", title, 120, 0, 360)
phis_hist_rd = ROOT.TH1F("Right Disks", title, 120, 0, 360)
phis_hist_ib = ROOT.TH1F("Inner Barrel", title, 120, 0, 360)
# phis_hist_ob = ROOT.TH1F("Outer Barrel", title, 120, 0, 360)

for i in range(0,360,3):
    phis_hist_total.SetBinContent((i//3) + 1, np.mean([lengths[0] for lengths in phis[i]]))
    phis_hist_ld.SetBinContent((i//3) + 1, np.mean([lengths[1] for lengths in phis[i]]))
    phis_hist_rd.SetBinContent((i//3) + 1, np.mean([lengths[2] for lengths in phis[i]]))
    phis_hist_ib.SetBinContent((i//3) + 1, np.mean([lengths[3] for lengths in phis[i]]))
    # phis_hist_ob.SetBinContent((i//3) + 1, np.mean([lengths[4] for lengths in phis[i]]))

phis_hist_total.SetXTitle("Azimuthal Angle (deg)")
phis_hist_total.SetYTitle("Average Number of Components Crossed")
phis_hist_total.SetMinimum(0)
phis_hist_total.SetStats(0)

phis_hist_total.SetLineColor(ROOT.kBlack)
phis_hist_ld.SetLineColor(ROOT.kGreen)
phis_hist_rd.SetLineColor(ROOT.kYellow)
phis_hist_ib.SetLineColor(ROOT.kRed)
# phis_hist_ob.SetLineColor(ROOT.kBlue)

phis_legend = ROOT.TLegend(0.4, 0.57, 0.6, 0.7)
phis_legend.AddEntry(phis_hist_total, "Total", "l")
phis_legend.AddEntry(phis_hist_ld, "Left Disks", "l")
phis_legend.AddEntry(phis_hist_rd, "Right Disks", "l")
phis_legend.AddEntry(phis_hist_ib, "Inner Barrel", "l")
# phis_legend.AddEntry(phis_hist_ob, "Outer Barrel", "l")

phis_canvas = ROOT.TCanvas("Phi Hits", title)
phis_hist_total.Draw("hist")
phis_hist_ld.Draw("same")
phis_hist_rd.Draw("same")
phis_hist_ib.Draw("same")
# phis_hist_ob.Draw("same")
phis_legend.Draw()
phis_canvas.Update()
phis_canvas.SaveAs("../plots/cld/cld_phi_hits_gp_combined.png")
