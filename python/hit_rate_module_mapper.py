from podio import root_io
import ROOT
import math
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)
events = podio_reader.get("events")

layer_radii = [14, 23, 34.5, 141, 316]
true_radii = [13.7, 23.7, 34, ]
disk_z = [303, 635, 945, -303, -635, -945]
particles = {11: "electron", -11: "positron"}
z_ranges = {14: 96, 23: 160, 34.5: 256}
phi_bins = {14: 24, 23: 15, 34.5: 10}

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

hit_map = {coord: {z: {azimuthal: [0] * 100 for azimuthal in range(0, 360, phi_bins[coord])} \
            for z in range(-z_ranges[coord], z_ranges[coord], 32)} for coord in layer_radii[:3]}

# for i in range(100):
#     event = events[i]
#     hits = {coord: [] for coord in layer_radii + disk_z}

#     for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
#         for hit in event.get(collection):

#             if collection != "VTXDCollection":
#                 hits[radius(hit)].append(hit)
#             else:
#                 hits[z_coord(hit)].append(hit)

#     for coord in layer_radii[:3]:
#         for hit in hits[coord]:
#             if hit.getEDep() < 0.000002 or abs(hit.getPosition().z) > z_ranges[coord]:
#                 continue

#             z = hit.getPosition().z
#             ph = phi(hit.getPosition().x, hit.getPosition().y) * (180 / math.pi)
#             hit_map[coord][int((z // 32) * 32)][int((ph // phi_bins[coord]) * phi_bins[coord])][i] += 1

# for layer_index in range(3):

#     r = layer_radii[layer_index]
#     hist = ROOT.TH2F("hit map", f"Guinea Pig Layer {layer_index + 1} Module Hits", \
#                     (2 * z_ranges[r]) // 32, -z_ranges[r], z_ranges[r], 360 // phi_bins[r], 0, 360)
#     hist.SetTitle(f"Guinea Pig Layer {layer_index + 1} Module Hits per Bunch Crossing;z (mm);Azimuthal Angle (deg)")

#     for z in range(-z_ranges[r], z_ranges[r], 32):
#         for azimuthal in range(0, 360, phi_bins[r]):
#             hits = np.mean(hit_map[r][z][azimuthal])
#             hist.SetBinContent(((z + z_ranges[r]) // 32) + 1, (azimuthal // phi_bins[r]) + 1, hits)

#     hist.SetStats(0)
#     canvas = ROOT.TCanvas("hit map", f"Guinea Pig Layer {layer_index + 1} Module Hits")
#     hist.Draw("colz")
#     canvas.Update()
#     canvas.SaveAs(f"../plots/hit_rates/by_module/gp_layer{layer_index + 1}_module_hit_rate.png")

# plots for ultra light vertex detector concept
for layer_index in range(1):
    r = true_radii[layer_index]
    layer_r = layer_radii[layer_index]
    hist = ROOT.TH2F(f"hit map {layer_index + 1}", f"Guinea Pig Layer {layer_index + 1} Module Hits", \
                    10, -108.35, 108.35, 8, 0, 360)
    for event in events:
        for particle in event.get("MCParticles"):

            x1 = particle.getVertex().x
            y1 = particle.getVertex().y
            z1 = particle.getVertex().z

            if x1**2 + y1**2 > r**2:
                continue

            x2 = particle.getEndpoint().x
            y2 = particle.getEndpoint().y
            z2 = particle.getEndpoint().z

            if x2**2 + y2**2 < r**2:
                continue

            slope = (x2 - x1) / (y2 - y1)

            # solves for x and y in trajectory that match r
            a = 1 + (slope**2)
            b = 2*x1*slope - 2*y1*(slope**2)
            c = (x1**2) - 2*x1*y1*slope + (y1**2)*(slope**2) - (r**2)

            d = (b**2) - 4*a*c

            y = (-b + np.sqrt(d))/(2*a)
            if y < min(y1, y2) or y > max(y1, y2):
                y = (-b - np.sqrt(d))/(2*a)

            x = x1 + slope*(y - y1)

            # interpolates z value of hit
            z = z1 + ((z2 - z1) / (y2 - y1)) * (y - y1)
            if abs(z) > 108.35:
                continue

            ph = phi(x, y) * (180 / math.pi)

            hist.Fill(z, ph)

    hist.Scale(0.1)
    hist.SetTitle(f"Guinea Pig Layer {layer_index + 1} Ultra Light Module Hits per Bunch Crossing;z (mm);Azimuthal Angle (deg)")
    hist.SetStats(0)
    canvas = ROOT.TCanvas(f"hit map {layer_index + 1}", f"Guinea Pig Layer {layer_index + 1} Module Hits")
    hist.Draw("colz")
    canvas.Update()
    canvas.SaveAs(f"../plots/hit_rates/by_module/gp_layer{layer_index + 1}_module_hit_rate_test.png")
