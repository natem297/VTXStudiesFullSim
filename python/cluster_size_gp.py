from podio import root_io
import ROOT
import numpy as np
import math

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
layer_radii = [14, 23, 34.5, 141, 316]
disk_z = [-303, -635, -945, 303, 635, 945]

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

def delta_squared(hit1, hit2):
    """
    Calculates delta of 2 hits.
    Inputs:
        hit1, hit2, SimTracker Hit object.
    Outputs:
        delta, float.
    """

    ti = hit1.getTime()
    xi = hit1.getPosition().x
    yi = hit1.getPosition().y
    zi = hit1.getPosition().z

    tf = hit2.getTime()
    xf = hit2.getPosition().x
    yf = hit2.getPosition().y
    zf = hit2.getPosition().z

    delta_t = tf - ti
    delta_x = xf - xi
    delta_y = yf - yi
    delta_z = zf - zi

    return delta_t**2 + delta_x**2 + delta_y**2 + delta_z**2

def cluster_size(particle, clusters, detector_index, detector, first_hit = None, monte_carlo = False):
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

    # if first_hit:
    #     theta_i = theta(xi, yi, zi)
    #     phi_i = phi(xi, yi)
    # else:
    theta_i = theta(px, py, pz)
    phi_i = phi(px, py)

    if detector == "barrel":
        coord = layer_radii[detector_index]
    else:
        coord = disk_z[detector_index]

    # if first_hit:
    #     new_hit = particle
    # else:
    new_hit = None

    for hit in hits[coord]:

        # if first_hit:
        #     theta_i_p = theta(px, py, pz)
        #     phi_i_p = phi(px, py)
        #     theta_f_p = theta(hit.getMomentum().x, hit.getMomentum().y, hit.getMomentum().z)
        #     phi_f_p = phi(hit.getMomentum().x, hit.getMomentum().y)
        #     delta_p = (theta_f_p - theta_i_p)**2 + (min(abs(phi_f_p - phi_i_p), 2*math.pi - abs(phi_f_p - phi_i_p)))**2

        #     theta_f = theta(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z)
        #     phi_f = phi(hit.getPosition().x, hit.getPosition().y)
        #     delta_pos = (theta_f - theta_i)**2 + (min(abs(phi_f - phi_i), 2*math.pi - abs(phi_f - phi_i)))**2

        #     delta = delta_p + delta_pos
        # else:
        delta = delta_squared(theta_i, phi_i, xi, yi, zi, hit)

        if delta < 0.01:
            clusters[coord] += 1
            new_hit = hit
            if first_hit is None:
                first_hit = hit

    if detector == "barrel" and detector_index == 4:
        return clusters, first_hit
    elif detector == "disk" and detector_index == 5:
        return clusters, first_hit
    # elif detector == "disk_r" and detector_index == 5:
    #     return clusters
    elif new_hit is None:
        return cluster_size(particle, clusters, detector_index + 1, detector, first_hit, monte_carlo)
    else:
        return cluster_size(new_hit, clusters, detector_index + 1, detector, first_hit)

all_clusters = {coord: {ang: [] for ang in range(0, 180, 3)} for coord in layer_radii + disk_z}
all_deltas = {coord: {ang: [] for ang in range(0, 180, 3)} for coord in layer_radii + disk_z}
for i in range(100):
    print(f"starting event {i}")
    event = events[i]
    hits = {coord: [] for coord in layer_radii + disk_z}

    # categorizes all hits by layer
    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):

            if collection != "VTXDCollection":
                hits[radius(hit)].append(hit)
            else:
                hits[z_coord(hit)].append(hit)

    for coord in layer_radii + disk_z:
        visited_hits = set()

        for hit1 in hits[coord]:

            if hit1 in visited_hits:
                continue

            visited_hits.add(hit1)
            cluster_size = 1
            polar = theta(hit1.getPosition().x, hit1.getPosition().y, hit1.getPosition().z) * (180 / math.pi)

            for hit2 in hits[coord]:

                if hit2 in visited_hits:
                    continue
                elif hit1.getMCParticle() == hit2.getMCParticle():
                    cluster_size += 1
                    visited_hits.add(hit2)
                    all_deltas[coord][int((polar // 3) * 3)].append(np.sqrt(delta_squared(hit1, hit2)))

            all_clusters[coord][int((polar // 3) * 3)].append(cluster_size)

for layer_index in range(5):
    # print(f"Cluster size for layer {layer_index + 1}: {all_clusters[layer_radii[layer_index]]}")
    hist = ROOT.TH1F("size", f"Guinea Pig Layer {layer_index + 1} Average Cluster Size", 60, 0, 180)
    for ang in range(0, 180, 3):
        if not all_clusters[layer_radii[layer_index]][ang]:
            hist.SetBinContent((ang//3) + 1, 0)
        else:
            hist.SetBinContent((ang//3) + 1, np.mean(all_clusters[layer_radii[layer_index]][ang]))
    hist.GetXaxis().SetTitle("Polar Angle (Degrees)")
    hist.GetYaxis().SetTitle("Average Cluster Size")
    hist.SetStats(0)
    canvas = ROOT.TCanvas("size", f"Guinea Pig Layer {layer_index + 1} Average Cluster Size")
    hist.Draw()
    canvas.Update()
    canvas.SaveAs(f"../plots/cluster_sizes/guinea_pig/gp_layer{layer_index + 1}_cluster_size_test.png")

for disk_index in range(6):
    hist = ROOT.TH1F("size", f"Guinea Pig Disk {disk_index + 1} Average Cluster Size", 60, 0, 180)
    for ang in range(0, 180, 3):
        if not all_clusters[disk_z[disk_index]][ang]:
            hist.SetBinContent((ang//3) + 1, 0)
        else:
            hist.SetBinContent((ang//3) + 1, np.mean(all_clusters[disk_z[disk_index]][ang]))
    hist.GetXaxis().SetTitle("Polar Angle (Degrees)")
    hist.GetYaxis().SetTitle("Average Cluster Size")
    hist.SetStats(0)
    canvas = ROOT.TCanvas("size", f"Guinea Pig Disk {disk_index + 1} Average Cluster Size")
    hist.Draw()
    canvas.Update()
    canvas.SaveAs(f"../plots/cluster_sizes/guinea_pig/gp_disk{disk_index + 1}_cluster_size_test.png")

for ang in all_deltas[14]:
    print(f"Mean delta for layer 1, angle {ang}: {np.mean(all_deltas[14][ang])}")
    print(f"Standard deviation: {np.std(all_deltas[14][ang])}")
