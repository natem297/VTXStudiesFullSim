from podio import root_io
import ROOT
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
e0 = events[0]
layer_radii = [14, 23, 34.5, 141, 316]
disk_z = [-945, -635, -303, 303, 635, 945]
hits = {coord: [] for coord in layer_radii + disk_z}
hits["mc"] = []

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

for particle in e0.get("MCParticles"):
    hits["mc"].append(particle)

for i in range(100):
    event = events[i]
    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):

            if collection != "VTXDCollection":
                hits[radius(hit)].append(hit)
            else:
                hits[z_coord(hit)].append(hit)

graph = ROOT.TGraph2D()
for i in range(len(hits["mc"])):
    graph.SetPoint(i, hits["mc"][i].getVertex().z, hits["mc"][i].getVertex().x, hits["mc"][i].getVertex().y)
graph.SetTitle("Monte Carlo Vertex Positions;z;x;y")
graph.SetMarkerColor(ROOT.kRed)
canvas = ROOT.TCanvas("Positions", "Monte Carlo Vertex Positions")
graph.Draw("P")
canvas.Update()
canvas.SaveAs("../plots/positions_mc.png")

# for layer_index in range(5):
#     graph = ROOT.TGraph2D()
#     r = layer_radii[layer_index]
#     for i in range(len(hits[r])):
#         graph.SetPoint(i, hits[r][i].getPosition().z, hits[r][i].getPosition().x, hits[r][i].getPosition().y)
#     graph.SetTitle(f"Layer {layer_index + 1} Hits;z;x;y")
#     graph.SetMarkerColor(ROOT.kRed)
#     canvas = ROOT.TCanvas("Hits", f"Layer {layer_index + 1}")
#     graph.Draw("P")
#     canvas.Update()
#     canvas.SaveAs(f"../plots/positions_layer{layer_index+1}.png")

# for disk_index in range(6):
#     graph = ROOT.TGraph2D()
#     z = disk_z[disk_index]
#     for i in range(len(hits[z])):
#         graph.SetPoint(i, hits[z][i].getPosition().z, hits[z][i].getPosition().x, hits[z][i].getPosition().y)
#     graph.SetTitle(f"Disk {disk_index + 1} Hits;z;x;y")
#     graph.SetMarkerColor(ROOT.kRed)
#     canvas = ROOT.TCanvas("Hits", f"Disk {disk_index + 1}")
#     graph.Draw("P")
#     canvas.Update()
#     canvas.SaveAs(f"../plots/positions_disk{disk_index+1}.png")
