from podio import root_io
import ROOT
import numpy as np
import math

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
layer_radii = [14, 23, 34.5, 141, 316]
disk_z = [303, 635, 945, -303, -635, -945]

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
        hit, SimTracker Hit object.
        theta_i, phi_i, floats representing direction of mc particle.
        xi, yi, zi, floats representing position of mc particle.
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

def trajectory_tracker(particle, trajectory):

    xi = particle.getPosition().x
    yi = particle.getPosition().y
    zi = particle.getPosition().z

    px = particle.getMomentum().x
    py = particle.getMomentum().y
    pz = particle.getMomentum().z

    theta_i = theta(px, py, pz)
    phi_i = phi(px, py)

    for coord in layer_radii + disk_z:
        if coord == 14:
            continue

        for hit in layers[coord]:
            if delta_squared(theta_i, phi_i, xi, yi, zi, hit) < 0.01:
                trajectory.append(hit)
                continue

    return trajectory

thetas = {i: [] for i in range(0,180,3)}
phis = {j: [] for j in range(0,360,3)}

for e in range(100):
    print(f"starting event {e}")
    event = events[e]
    layers = {coord: [] for coord in layer_radii + disk_z}

    # categorizes all hits by layer
    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):

            if collection != "VTXDCollection":
                layers[radius(hit)].append(hit)
            else:
                layers[z_coord(hit)].append(hit)

    for initial_hit in layers[14]:
        traj = trajectory_tracker(initial_hit, [initial_hit])

        polar = theta(initial_hit.getPosition().x, initial_hit.getPosition().y, initial_hit.getPosition().z)
        polar = int(((360 / (2 * math.pi)) * polar) // 3) * 3
        thetas[polar].append(len(traj))

        azimuthal = phi(initial_hit.getPosition().x, initial_hit.getPosition().y)
        azimuthal = int(((360 / (2 * math.pi)) * azimuthal) // 3) * 3
        phis[azimuthal].append(len(traj))

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
thetas_canvas.SaveAs("../plots/angle_hits/theta_hits_guinea_pig_test.png")

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
phis_canvas.SaveAs("../plots/angle_hits/phi_hits_guinea_pig_test.png")
