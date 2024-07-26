from podio import root_io
import ROOT
import numpy as np
import math

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/electronGun/egun_10GeV_10k_Geant4TrackerAction_edep0.root"
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

def match_counter(particle, matches, detector_index, detector, monte_carlo = False):
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
        trajectory: list starting with MCParticle followed by SimTrackerHit objects.
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

    if detector == "barrel":
        coord = layer_radii[detector_index]
    else:
        coord = disk_z[detector_index]

    new_hit = None
    for hit in hits[coord]:
        if delta_squared(theta_i, phi_i, xi, yi, zi, hit) < 0.001:
            matches[coord] += 1
            new_hit = hit

    if detector == "barrel" and detector_index == 4:
        return matches
    elif detector == "disk" and detector_index == 5:
        return matches
    elif new_hit is None:
        return match_counter(particle, matches, detector_index + 1, detector, monte_carlo)
    else:
        return match_counter(new_hit, matches, detector_index + 1, detector)

all_matches = {coord: {ang: [] for ang in range(0, 180, 3)} for coord in layer_radii + disk_z}
for i in range(10000):
    event = events[i]
    hits = {coord: [] for coord in layer_radii + disk_z}

    mc = event.get("MCParticles")
    if len(mc) > 100:
        continue

    # categorizes all hits by layer
    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):

            if collection != "VTXDCollection":
                hits[radius(hit)].append(hit)
            else:
                hits[z_coord(hit)].append(hit)

    particle = mc[0]
    th = theta(particle.getMomentum().x, particle.getMomentum().y, particle.getMomentum().z) * (180 / math.pi)

    barrel_matches = match_counter(particle, {r: 0 for r in layer_radii}, 0, "barrel", True)
    for r in barrel_matches:
        all_matches[r][int((th // 3) * 3)].append(barrel_matches[r])

    disk_matches = match_counter(particle, {z: 0 for z in disk_z}, 0, "disk", True)
    for z in disk_matches:
        all_matches[z][int((th //3) * 3)].append(disk_matches[z])

for layer_index in range(5):
    hist = ROOT.TH1F("size", f"Electron Gun Layer {layer_index + 1} Average Number of Matches", 60, 0, 180)
    for ang in range(0, 180, 3):
        if not all_matches[layer_radii[layer_index]][ang]:
            hist.SetBinContent((ang//3) + 1, 0)
        else:
            hist.SetBinContent((ang//3) + 1, np.mean(all_matches[layer_radii[layer_index]][ang]))
    hist.GetXaxis().SetTitle("Polar Angle (Degrees)")
    hist.GetYaxis().SetTitle("Average Number of SimTrackerHit Matches")
    hist.SetStats(0)
    canvas = ROOT.TCanvas("size", f"Electron Gun Layer {layer_index + 1} Average Number of Matches")
    hist.Draw()
    canvas.Update()
    canvas.SaveAs(f"../plots/hit_matches/electron_gun/eg_layer{layer_index + 1}_match_counter.png")

for disk_index in range(6):
    hist = ROOT.TH1F("size", f"Electron Gun Disk {disk_index + 1} Average Number of Matches", 60, 0, 180)
    for ang in range(0, 180, 3):
        hist.SetBinContent((ang//3) + 1, np.mean(all_matches[disk_z[disk_index]][ang]))
    hist.GetXaxis().SetTitle("Polar Angle (Degrees)")
    hist.GetYaxis().SetTitle("Average Number of SimTrackerHit Matches")
    hist.SetStats(0)
    canvas = ROOT.TCanvas("size", f"Electron Gun Disk {disk_index + 1} Average Number of Matches")
    hist.Draw()
    canvas.Update()
    canvas.SaveAs(f"../plots/hit_matches/electron_gun/eg_disk{disk_index + 1}_match_counter.png")
