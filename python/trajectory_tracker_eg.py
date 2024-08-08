from podio import root_io
import ROOT
import numpy as np
import math

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/electronGun/egun_10GeV_10k_Geant4TrackerAction_edep0.root"
# input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
layer_radii = [14, 23, 34.5, 141, 316]
disk_z = [303, 635, 945, -303, -635, -945]
coords = layer_radii + disk_z

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

def delta_squared(theta_i, phi_i, xi, yi, zi, hit):
    """
    Calculates delta of a particle-hit pair.
    Inputs:
        hit: SimTracker Hit object.
        theta_i, phi_i: floats representing direction of mc particle.
        xi, yi, zi: floats representing position of mc particle.
    Outputs: delta, float.
    """

    xf = hit.getPosition().x
    yf = hit.getPosition().y
    zf = hit.getPosition().z

    delta_x = xf - xi
    delta_y = yf - yi
    delta_z = zf - zi

    theta_f = theta(delta_x, delta_y, delta_z)
    phi_f = phi(delta_x, delta_y)

    delta_theta = theta_f - theta_i
    delta_phi = min(abs(phi_f - phi_i), 2*math.pi - abs(phi_f - phi_i))

    return delta_theta**2 + delta_phi**2

def trajectory_tracker(particle, trajectory, detector_index, detector, monte_carlo = False):
    """
    Finds trajectory of a given particle by finding hits with a delta squared below a certain
    threshold.  Then recursively calls using the last hit to find matching hits in further detectors.
    Inputs:
        particle: MCParticle or SimTrackerHit representing most recent object in trajectory.
        trajectory: list starting with MCParticle followed by SimTrackerHit objects.
        detector_index: int, index of detector to be explored next.
        detector: string representing which detector is being explored.
        monte_carlo: boolean representing whether the particle is an MCParticle.
    Outputs:
        trajectory: list starting with MCParticle followed by Sim TrackerHit objects.
    """
    if monte_carlo:
        xi = particle.getVertex().x
        yi = particle.getVertex().y
        zi = particle.getVertex().z
    else:
        xi = particle.getPosition().x
        yi = particle.getPosition().y
        zi = particle.getPosition().z

    px = particle.getMomentum().x
    py = particle.getMomentum().y
    pz = particle.getMomentum().z

    theta_i = theta(px, py, pz)
    phi_i = phi(px, py)

    if detector == "ib":
        coord = layer_radii[detector_index]
    elif detector == "ob":
        coord = layer_radii[detector_index]
    else:
        coord = disk_z[detector_index]

    for hit in layers[coord]:
        if delta_squared(theta_i, phi_i, xi, yi, zi, hit) < 0.001:
            trajectory.append(hit)
            break

    if detector == "ib" and detector_index == 2:
        return trajectory
    elif detector == "ob" and detector_index == 4:
        return trajectory
    elif detector == "disk" and detector_index == 5:
        return trajectory
    elif not trajectory:
        return trajectory_tracker(particle, trajectory, detector_index + 1, detector, True)
    else:
        return trajectory_tracker(trajectory[-1], trajectory, detector_index + 1, detector)

thetas = {i: [] for i in range(0,180,3)}
phis = {j: [] for j in range(0,360,3)}

for e in range(10000):
    # print(f"starting event {e}")
    event = events[e]
    hits = {coord: [] for coord in layer_radii + disk_z}

    # categorizes all hits by layer
    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):

            if collection != "VTXDCollection":
                hits[radius(hit)].append(hit)
            else:
                hits[z_coord(hit)].append(hit)

    visited_mc = []
    visited_hits = []

    for i in range(11):
        for hit in hits[coords[i]]:

            if hit in visited_hits:
                continue
            mc = hit.getMCParticle()
            if mc in visited_mc:
                continue
            visited_mc.append(mc)

            polar = theta(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z) * (180 / math.pi)
            polar = int(polar // 3) * 3

            azimuthal = phi(hit.getPosition().x, hit.getPosition().y) * (180 / math.pi)
            azimuthal = int(azimuthal // 3) * 3

            traj_length = 1

            for j in range(i+1,11):
                hit_found = False
                for hit2 in hits[coords[j]]:

                    if hit2.getMCParticle() == mc:
                        if not hit_found:
                            traj_length += 1
                            hit_found = True
                        visited_hits.append(hit2)


            thetas[polar].append(traj_length)
            phis[azimuthal].append(traj_length)

thetas_hist_total = ROOT.TH1F("Total", "Detector Hits vs Theta", 60, 0, 180)

for i in range(0,180,3):
    if not thetas[i]:
        thetas_hist_total.SetBinContent((i//3) + 1, 0)
    else:
        thetas_hist_total.SetBinContent((i//3) + 1, np.mean(thetas[i]))

thetas_hist_total.SetXTitle("Polar Angle (deg)")
thetas_hist_total.SetYTitle("Average Number of Hits")
thetas_hist_total.SetStats(0)

thetas_canvas = ROOT.TCanvas("Theta Hits", "Detector Hits vs Theta")
thetas_hist_total.Draw("hist")
thetas_canvas.Update()
thetas_canvas.SaveAs("../plots/angle_hits/theta_hits_electron_gun_test.png")

phis_hist_total = ROOT.TH1F("Total", "Detector Hits vs Phi", 120, 0, 360)

for i in range(0,360,3):
    if not phis[i]:
        phis_hist_total.SetBinContent((i//3) + 1, 0)
    else:
        phis_hist_total.SetBinContent((i//3) + 1, np.mean(phis[i]))

phis_hist_total.SetXTitle("Azimuthal Angle (deg)")
phis_hist_total.SetYTitle("Average Number of Hits")
phis_hist_total.SetMinimum(0)
phis_hist_total.SetStats(0)

phis_canvas = ROOT.TCanvas("Phi Hits", "Detector Hits vs Phi")
phis_hist_total.Draw("hist")
phis_canvas.Update()
phis_canvas.SaveAs("../plots/angle_hits/phi_hits_electron_gun_test.png")

# # categorizes all hits by layer
# test_events = [events[i] for i in range(10000)]
# for event in test_events:
#     layers = {coord: [] for coord in layer_radii + disk_z}

#     for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
#         for hit in event.get(collection):

#             if collection != "VTXDCollection":
#                 layers[radius(hit)].append(hit)
#             else:
#                 layers[z_coord(hit)].append(hit)

#     mc = event.get("MCParticles")
#     if len(mc) > 100:
#         continue
#     particle = mc[0]

#     ib_traj = trajectory_tracker(particle, [], 0, "ib", True)
#     ob_traj = trajectory_tracker(particle, [], 0, "ob", True)
#     disk_traj = trajectory_tracker(particle, [], 0, "disk", True)

#     ib_length = len(ib_traj)
#     ob_length = len(ob_traj) - ib_length
#     disk_length = len(disk_traj)
#     total_length = ib_length + ob_length + disk_length

#     polar = theta(particle.getMomentum().x, particle.getMomentum().y, particle.getMomentum().z)
#     polar = int(((360 / (2 * math.pi)) * polar) // 3) * 3
#     thetas[polar].append((ib_length, ob_length, disk_length, total_length))

#     azimuthal = phi(particle.getMomentum().x, particle.getMomentum().y)
#     azimuthal = int(((360 / (2 * math.pi)) * azimuthal) // 3) * 3
#     phis[azimuthal].append((ib_length, ob_length, disk_length, total_length))

# thetas_hist_ib = ROOT.TH1F("Inner Barrel", "FCC-ee IDEA Vertex Detector", 60, 0, 180)
# thetas_hist_ob = ROOT.TH1F("Outer Barrel", "FCC-ee IDEA Vertex Detector", 60, 0, 180)
# thetas_hist_disk = ROOT.TH1F("Disks", "FCC-ee IDEA Vertex Detector", 60, 0, 180)
# thetas_hist_total = ROOT.TH1F("Total", "FCC-ee IDEA Vertex Detector", 60, 0, 180)

# for i in range(0,180,3):
#     if not thetas[i]:
#         thetas_hist_ib.SetBinContent((i//3) + 1, 0)
#         thetas_hist_ob.SetBinContent((i//3) + 1, 0)
#         thetas_hist_disk.SetBinContent((i//3) + 1, 0)
#         thetas_hist_total.SetBinContent((i//3) + 1, 0)
#     else:
#         thetas_hist_ib.SetBinContent((i//3) + 1, np.mean([lengths[0] for lengths in thetas[i]]))
#         thetas_hist_ob.SetBinContent((i//3) + 1, np.mean([lengths[1] for lengths in thetas[i]]))
#         thetas_hist_disk.SetBinContent((i//3) + 1, np.mean([lengths[2] for lengths in thetas[i]]))
#         thetas_hist_total.SetBinContent((i//3) + 1, np.mean([lengths[3] for lengths in thetas[i]]))

# thetas_hist_total.SetXTitle("Polar Angle (deg)")
# thetas_hist_total.SetYTitle("Average Number of Hits")
# thetas_hist_total.SetStats(0)

# thetas_hist_ib.SetLineColor(ROOT.kRed)
# thetas_hist_ob.SetLineColor(ROOT.kBlue)
# thetas_hist_disk.SetLineColor(ROOT.kGreen)
# thetas_hist_total.SetLineColor(ROOT.kBlack)

# thetas_legend = ROOT.TLegend(0.4, 0.75, 0.6, 0.88)
# thetas_legend.AddEntry(thetas_hist_total, "Total", "l")
# thetas_legend.AddEntry(thetas_hist_ib, "Inner Barrel", "l")
# thetas_legend.AddEntry(thetas_hist_ob, "Outer Barrel", "l")
# thetas_legend.AddEntry(thetas_hist_disk, "Disks", "l")

# thetas_canvas = ROOT.TCanvas("Theta Hits", "FCC-ee IDEA Vertex Detector")
# thetas_hist_total.Draw("hist")
# thetas_hist_ib.Draw("same")
# thetas_hist_ob.Draw("same")
# thetas_hist_disk.Draw("same")
# thetas_legend.Draw()
# thetas_canvas.Update()
# thetas_canvas.SaveAs("../plots/angle_hits/theta_hits_eg_combined.png")

# phis_hist_ib = ROOT.TH1F("Inner Barrel", "FCC-ee IDEA Vertex Detector", 120, 0, 360)
# phis_hist_ob = ROOT.TH1F("Outer Barrel", "FCC-ee IDEA Vertex Detector", 120, 0, 360)
# phis_hist_disk = ROOT.TH1F("Disks", "FCC-ee IDEA Vertex Detector", 120, 0, 360)
# phis_hist_total = ROOT.TH1F("Total", "FCC-ee IDEA Vertex Detector", 120, 0, 360)

# for i in range(0,360,3):
#     if not phis[i]:
#         phis_hist_ib.SetBinContent((i//3) + 1, 0)
#         phis_hist_ob.SetBinContent((i//3) + 1, 0)
#         phis_hist_disk.SetBinContent((i//3) + 1, 0)
#         phis_hist_total.SetBinContent((i//3) + 1, 0)
#     else:
#         phis_hist_ib.SetBinContent((i//3) + 1, np.mean([lengths[0] for lengths in phis[i]]))
#         phis_hist_ob.SetBinContent((i//3) + 1, np.mean([lengths[1] for lengths in phis[i]]))
#         phis_hist_disk.SetBinContent((i//3) + 1, np.mean([lengths[2] for lengths in phis[i]]))
#         phis_hist_total.SetBinContent((i//3) + 1, np.mean([lengths[3] for lengths in phis[i]]))

# phis_hist_total.SetXTitle("Azimuthal Angle (deg)")
# phis_hist_total.SetYTitle("Average Number of Hits")
# phis_hist_total.SetMinimum(0)
# phis_hist_total.SetStats(0)

# phis_hist_ib.SetLineColor(ROOT.kRed)
# phis_hist_ob.SetLineColor(ROOT.kBlue)
# phis_hist_disk.SetLineColor(ROOT.kGreen)
# phis_hist_total.SetLineColor(ROOT.kBlack)

# phis_legend = ROOT.TLegend(0.4, 0.6, 0.6, 0.73)
# phis_legend.AddEntry(phis_hist_total, "Total", "l")
# phis_legend.AddEntry(phis_hist_ib, "Inner Barrel", "l")
# phis_legend.AddEntry(phis_hist_ob, "Outer Barrel", "l")
# phis_legend.AddEntry(phis_hist_disk, "Disks", "l")

# phis_canvas = ROOT.TCanvas("Phi Hits", "FCC-ee IDEA Vertex Detector")
# phis_hist_total.Draw("hist")
# phis_hist_ib.Draw("same")
# phis_hist_ob.Draw("same")
# phis_hist_disk.Draw("same")
# phis_legend.Draw()
# phis_canvas.Update()
# phis_canvas.SaveAs("../plots/angle_hits/phi_hits_eg_combined.png")
