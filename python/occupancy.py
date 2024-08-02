from podio import root_io
import ROOT
import math
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)
events = podio_reader.get("events")

layer_radii = [14, 23, 34.5, 141, 316]
true_radii = [13.7, 23.7, 34]
disk_z = [303, 635, 945, -303, -635, -945]
pixel_size = 0.025

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

pixel_hits = {(s_bin, z_bin): [0]*100 for s_bin in range(87) for z_bin in range(192)}
for i in range(100):
    print(f"starting event {i}")
    event = events[i]
    hits = {coord: [] for coord in layer_radii + disk_z}

    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):

            if collection != "VTXDCollection":
                hits[radius(hit)].append(hit)
            else:
                hits[z_coord(hit)].append(hit)

    for hit in hits[14]:
        r = true_radii[0]
        azimuthal = phi(hit.getPosition().x, hit.getPosition().y)
        s = r * azimuthal
        z = hit.getPosition().z
        if abs(z) > 96:
            continue
        s_pixel = int(s // 1)
        z_pixel = int((z + 96) // 1)
        pixel_hits[(s_pixel, z_pixel)][i] += 1

pixel_hit_averages = {pix: np.mean(pixel_hits[pix]) for pix in pixel_hits.keys()}
max_hits = max(pixel_hit_averages.values())
avg_hits = np.mean(list(pixel_hit_averages.values()))

hist = ROOT.TH2F("hits", "Guinea Pig Hits Per Area (mm^-2)", 192, -96, 96, 86, 0, 86)
for s, z in pixel_hit_averages.keys():
    hist.SetBinContent(z + 1, s + 1, pixel_hit_averages[(s, z)])

hist.SetStats(0)
hist.SetTitle("Guinea Pig Layer 1 Hits Per Area (mm^-2);z (mm);Arc Distance (mm)")
canvas = ROOT.TCanvas("hits", "Guinea Pig Hits Per Area (mm^-2)")
hist.Draw("colz")
canvas.Update()
canvas.SaveAs("../plots/hit_rates/per_area/occupancy.png")

print(f"Maximum occupancy: {max_hits*25*25*15/1000000}")
print(f"Average occupancy: {avg_hits*25*25*15/1000000}")
