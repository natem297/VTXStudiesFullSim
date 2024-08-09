from podio import root_io
import ROOT
import numpy as np
import math

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/electronGun/egun_10GeV_10k_Geant4TrackerAction_edep0.root"
# input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
layer_radii = [14, 23, 34.5, 141, 316] # approximate r values of barrels
disk_z = [-303, -635, -945, 303, 635, 945] # approximate z values of disks
component_coords = layer_radii + disk_z
# maps subdetector index to sub detector (inner barrel, outer barrel, left disks, right disks)
comp_index_dict = {0: "ib", 1: "ib", 2: "ib", 3: "ob", 4: "ob", \
                    5: "ld", 6: "ld", 7: "ld", 8: "rd", 9: "rd", 10: "rd"}

def radius(hit):
    """
    Calculates polar radius of particle.
    Inputs: hit, SimTrackerHit object.
    Output: r, number representing polar radius in mm.
    """
    true_radius = np.sqrt(hit.getPosition().x**2 + hit.getPosition().y**2)
    for r in layer_radii:
        if abs(true_radius-r) < 4:
            return r
    raise ValueError(f"Not close enough to any of the layers {np.sqrt(true_radius)}")

def z_coord(hit):
    """
    Calculates z coordinate of hit.
    Inputs: hit, SimTrackerHit object.
    Output: z, int representing z coordinate of hit in mm.
    """
    true_z = hit.getPosition().z
    for z in disk_z:
        if abs(true_z-z) < 30:
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
        start_subdet: int, representing index of subdetector where starting hit is located.
    Outputs:
        traj_lengths: dict, mapping each subdetector (ib, ob, ld, rd) to the number of hits
            in that subdetector.
    """
    traj_lengths = {"ib": 0, "ob": 0, "ld": 0, "rd": 0, "total": 1}
    traj_lengths[comp_index_dict[start_comp]] += 1 # adds 1 to starting subdetector

    for comp_index in range(start_comp + 1, 11):
        for hit in hits[component_coords[comp_index]]:
            if hit.getMCParticle() == mc_particle:
                traj_lengths[comp_index_dict[comp_index]] += 1
                traj_lengths["total"] += 1
                break

    return traj_lengths

thetas = {i: [] for i in range(0,180,3)}
phis = {j: [] for j in range(0,360,3)}

for e in range(10000):
    event = events[e]
    hits = {coord: [] for coord in component_coords}

    # categorizes barrel hits by radius and disk hits by z coord
    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):

            if collection != "VTXDCollection":
                hits[radius(hit)].append(hit)
            else:
                hits[z_coord(hit)].append(hit)

    visited_mc = []

    for comp_index in range(11):
        for hit in hits[component_coords[comp_index]]:

            mc = hit.getMCParticle()
            if mc in visited_mc or mc.getGeneratorStatus() != 1:
                continue
            visited_mc.append(mc)

            polar = theta(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z) * (180 / math.pi)
            polar = int(polar // 3) * 3 # rounds down to nearest multiple of 3

            azimuthal = phi(hit.getPosition().x, hit.getPosition().y) * (180 / math.pi)
            azimuthal = int(azimuthal // 3) * 3 # rounds down to nearest multiple of 3

            traj_lengths = trajectory_length(comp_index, mc)

            thetas[polar].append(list(traj_lengths.values()))
            phis[azimuthal].append(list(traj_lengths.values()))

# creates theta hist objects for each subdetector and plots them
thetas_hist_ib = ROOT.TH1F("Inner Barrel", "FCC-ee IDEA Vertex Detector", 60, 0, 180)
thetas_hist_ob = ROOT.TH1F("Outer Barrel", "FCC-ee IDEA Vertex Detector", 60, 0, 180)
thetas_hist_ld = ROOT.TH1F("Left Disks", "FCC-ee IDEA Vertex Detector", 60, 0, 180)
thetas_hist_rd = ROOT.TH1F("Right Disks", "FCC-ee IDEA Vertex Detector", 60, 0, 180)
thetas_hist_total = ROOT.TH1F("Total", "FCC-ee IDEA Vertex Detector", 60, 0, 180)

for i in range(0,180,3):
    if not thetas[i]: # for theta values with no hits
        thetas_hist_ib.SetBinContent((i//3) + 1, 0)
        thetas_hist_ob.SetBinContent((i//3) + 1, 0)
        thetas_hist_ld.SetBinContent((i//3) + 1, 0)
        thetas_hist_rd.SetBinContent((i//3) + 1, 0)
        thetas_hist_total.SetBinContent((i//3) + 1, 0)
    else:
        thetas_hist_ib.SetBinContent((i//3) + 1, np.mean([lengths[0] for lengths in thetas[i]]))
        thetas_hist_ob.SetBinContent((i//3) + 1, np.mean([lengths[1] for lengths in thetas[i]]))
        thetas_hist_ld.SetBinContent((i//3) + 1, np.mean([lengths[2] for lengths in thetas[i]]))
        thetas_hist_rd.SetBinContent((i//3) + 1, np.mean([lengths[3] for lengths in thetas[i]]))
        thetas_hist_total.SetBinContent((i//3) + 1, np.mean([lengths[4] for lengths in thetas[i]]))

thetas_hist_ld.SetXTitle("Polar Angle (deg)")
thetas_hist_ld.SetYTitle("Average Number of Components Crossed")
thetas_hist_ld.SetStats(0)

thetas_hist_ib.SetLineColor(ROOT.kRed)
thetas_hist_ob.SetLineColor(ROOT.kBlue)
thetas_hist_ld.SetLineColor(ROOT.kGreen)
thetas_hist_rd.SetLineColor(ROOT.kYellow)
thetas_hist_total.SetLineColor(ROOT.kBlack)

thetas_legend = ROOT.TLegend(0.4, 0.75, 0.6, 0.88)
thetas_legend.AddEntry(thetas_hist_total, "Total", "l")
thetas_legend.AddEntry(thetas_hist_ib, "Inner Barrel", "l")
thetas_legend.AddEntry(thetas_hist_ob, "Outer Barrel", "l")
thetas_legend.AddEntry(thetas_hist_ld, "Left Disks", "l")
thetas_legend.AddEntry(thetas_hist_rd, "Right Disks", "l")

thetas_canvas = ROOT.TCanvas("Theta Hits", "FCC-ee IDEA Vertex Detector")
thetas_hist_total.Draw("hist")
thetas_hist_ib.Draw("same")
thetas_hist_ob.Draw("same")
thetas_hist_ld.Draw("same")
thetas_hist_rd.Draw("same")
thetas_legend.Draw()
thetas_canvas.Update()
thetas_canvas.SaveAs("../plots/angle_hits/theta_hits_eg_combined.png")

# creates phi hist objects for each subdetector and plots them
phis_hist_ib = ROOT.TH1F("Inner Barrel", "FCC-ee IDEA Vertex Detector", 120, 0, 360)
phis_hist_ob = ROOT.TH1F("Outer Barrel", "FCC-ee IDEA Vertex Detector", 120, 0, 360)
phis_hist_ld = ROOT.TH1F("Left Disks", "FCC-ee IDEA Vertex Detector", 120, 0, 360)
phis_hist_rd = ROOT.TH1F("Right Disks", "FCC-ee IDEA Vertex Detector", 120, 0, 360)
phis_hist_total = ROOT.TH1F("Total", "FCC-ee IDEA Vertex Detector", 120, 0, 360)

for i in range(0,360,3):
    phis_hist_ib.SetBinContent((i//3) + 1, np.mean([lengths[0] for lengths in phis[i]]))
    phis_hist_ob.SetBinContent((i//3) + 1, np.mean([lengths[1] for lengths in phis[i]]))
    phis_hist_ld.SetBinContent((i//3) + 1, np.mean([lengths[2] for lengths in phis[i]]))
    phis_hist_rd.SetBinContent((i//3) + 1, np.mean([lengths[3] for lengths in phis[i]]))
    phis_hist_total.SetBinContent((i//3) + 1, np.mean([lengths[4] for lengths in phis[i]]))

phis_hist_total.SetXTitle("Azimuthal Angle (deg)")
phis_hist_total.SetYTitle("Average Number of Components Crossed")
phis_hist_total.SetMinimum(0)
phis_hist_total.SetStats(0)

phis_hist_ib.SetLineColor(ROOT.kRed)
phis_hist_ob.SetLineColor(ROOT.kBlue)
phis_hist_ld.SetLineColor(ROOT.kGreen)
phis_hist_rd.SetLineColor(ROOT.kYellow)
phis_hist_total.SetLineColor(ROOT.kBlack)

phis_legend = ROOT.TLegend(0.4, 0.6, 0.6, 0.73)
phis_legend.AddEntry(phis_hist_total, "Total", "l")
phis_legend.AddEntry(phis_hist_ib, "Inner Barrel", "l")
phis_legend.AddEntry(phis_hist_ob, "Outer Barrel", "l")
phis_legend.AddEntry(phis_hist_ld, "Left Disks", "l")
phis_legend.AddEntry(phis_hist_rd, "Right Disks", "l")

phis_canvas = ROOT.TCanvas("Phi Hits", "FCC-ee IDEA Vertex Detector")
phis_hist_total.Draw("hist")
phis_hist_ib.Draw("same")
phis_hist_ob.Draw("same")
phis_hist_ld.Draw("same")
phis_hist_rd.Draw("same")
phis_legend.Draw()
phis_canvas.Update()
phis_canvas.SaveAs("../plots/angle_hits/phi_hits_eg_combined.png")
