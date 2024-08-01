from podio import root_io
import ROOT
import math
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)
events = podio_reader.get("events")

layer_radii = [14, 23, 34.5, 141, 316]
true_radii = [1.37, 2.37, 3.4, 13, 31.5]
disk_z = [303, 635, 945, -303, -635, -945]
z_ranges = {14: (96,4), 23: (162,6), 34.5: (261,9), 141: (162,6), 316: (324,12)}
particles = {11: "electron", -11: "positron"}

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

hit_map = {coord: {z: {azimuthal: {particle: [0] * 100 for particle in particles.values()} \
                for azimuthal in range(0, 360, 3)} \
                    for z in range(-z_ranges[coord][0], z_ranges[coord][0], z_ranges[coord][1])} \
                        for coord in layer_radii}

# hit_map = {coord: {z: {azimuthal: [0] * 100 for azimuthal in range(0, 360, 3)} \
#                 for z in range(-z_ranges[coord][0], z_ranges[coord][0], z_ranges[coord][1])} \
#                     for coord in layer_radii}

for i in range(100):
    event = events[i]
    hits = {coord: [] for coord in layer_radii + disk_z}

    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):

            if collection != "VTXDCollection":
                hits[radius(hit)].append(hit)
            else:
                hits[z_coord(hit)].append(hit)

    for coord in layer_radii:
        for hit in hits[coord]:
            if hit.getEDep() < 0.000002:
                continue

            mc = hit.getMCParticle()
            pdg = mc.getPDG()

            if pdg in particles:
                particle = particles[pdg]
            else:
                continue

            ph = phi(hit.getPosition().x, hit.getPosition().y) * (180 / math.pi)
            z = hit.getPosition().z
            if abs(z) > z_ranges[coord][0]:
                continue
            hit_map[coord][int((z // z_ranges[coord][1]) * z_ranges[coord][1])][int((ph // 3) * 3)][particle][i] += 1
            # hit_map[coord][int((z // z_ranges[coord][1]) * z_ranges[coord][1])][int((ph // 3) * 3)][i] += 1

for layer_index in range(5):
    for particle in particles.values():

        r = true_radii[layer_index]
        z_range = z_ranges[layer_radii[layer_index]]
        area = z_range[1] * r * math.pi/60
        hist = ROOT.TH2F("hit map", f"Guinea Pig Layer {layer_index + 1} Hits", (2 * z_range[0]) // z_range[1], -z_range[0], z_range[0], 120, 0, 360)
        # hist.SetTitle(f"Guinea Pig Layer {layer_index + 1} {particle.capitalize()} Hits Per Area (cm ^-2);Cosine Theta;Azimuthal Angle (deg)")
        hist.SetTitle(f"Guinea Pig Layer {layer_index + 1} Hits Per Area (cm^-2);z (mm); Azimuthal Angle (deg)")

        for z in range(-z_range[0], z_range[0], z_range[1]):
            for azimuthal in range(0, 360, 3):
                hits = np.mean(hit_map[layer_radii[layer_index]][z][azimuthal][particle])
                # hits = np.mean(hit_map[layer_radii[layer_index]][z][azimuthal])
                hist.SetBinContent(((z + z_range[0]) // z_range[1]) + 1, (azimuthal // 3) + 1, hits / area)

        hist.SetStats(0)
        canvas = ROOT.TCanvas("hit map", f"Guinea Pig Layer {layer_index + 1} Hits")
        canvas.SetRightMargin(0.12)
        hist.Draw("colz")
        canvas.Update()
        canvas.SaveAs(f"../plots/hit_rates/{particle}/gp_layer{layer_index + 1}_{particle}_hit_rate.png")
        # canvas.SaveAs(f"../plots/hit_rates/per_area/gp_layer{layer_index + 1}_hit_rate.png")
