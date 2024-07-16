from podio import root_io
import ROOT
import numpy as np
import csv

# input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/electronGun/egun_10GeV_10k_Geant4TrackerAction_edep0.root"
input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
layer_radii = [14, 23, 34.5, 141, 316]
true_layer_radii = ["(13.7 mm)", "(23.7 mm)", "(34 mm)", "(130 mm)", "(315 mm)"]
disk_z = [-945, -635, -303, 303, 635, 945]
hit_counter = {coord: [0]*100 for coord in layer_radii + disk_z}
hit_counter["mc_charge"] = [0]*100
hit_counter["mc_neutral"] = [0]*100

def radius(hit):
    """
    Calculates polar radius of particle.
    Inputs: hit, SimTrackerHit object.
    Output: r, float representing polar radius in mm.
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
        if abs(true_z-z) < 30:
            return z
    raise ValueError(f"Not close enough to any of the disks {true_z}")

for i in range(100):
    event = events[i]
    for particle in event.get("MCParticles"):
        if particle.getCharge() != 0:
            hit_counter["mc_charge"][i] += 1
        else:
            hit_counter["mc_neutral"][i] += 1

    for collection in ["VTXIBCollection", "VTXOBCollection", "VTXDCollection"]:
        for hit in event.get(collection):

            if collection != "VTXDCollection":
                hit_counter[radius(hit)][i] += 1
            else:
                hit_counter[z_coord(hit)][i] += 1

print(f"Average layer counts: {[(coord, np.mean(counts)) for coord, counts in hit_counter.items()]}")
print(f"Standard Deviation of layer counts: {[(coord, np.std(counts)) for coord, counts in hit_counter.items()]}")

hist = ROOT.TH1F("counts", "Average Number of Hits per Bunch Crossing;Detector;Number of Hits", 27, 0, 27)
hist.SetBinContent(2, np.mean(hit_counter["mc_charge"]))
hist.GetXaxis().SetBinLabel(2, "MC Charged")
hist.SetBinContent(4, np.mean(hit_counter["mc_neutral"]))
hist.GetXaxis().SetBinLabel(4, "MC Neutral")
# plots layer hits
for i in range(5):
    hist.SetBinContent(6 + 2*i, np.mean(hit_counter[layer_radii[i]]))
    hist.GetXaxis().SetBinLabel(6 + 2*i, f"Layer {i + 1} " + true_layer_radii[i])
# plots disk hits
for j in range(6):
    hist.SetBinContent(16 + 2*j, np.mean(hit_counter[disk_z[j]]))
    hist.GetXaxis().SetBinLabel(16 + 2*j, f"Disk {j + 1}")

hist.SetFillColor(ROOT.kBlue)
canvas = ROOT.TCanvas("counts", "Hit Counts")
canvas.SetLogy(True)
hist.Draw("b")
canvas.Update()
canvas.SaveAs("../plots/monte_carlo/hit_counts_test.png")

# with open("../data_tables/layer_hits.csv", "w") as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerow(["Monte Carlo"] + hit_counter["mc"])
#     for i in range(5):
#         radii = layer_radii[i]
#         writer.writerow([f"Layer {i+1}"] + hit_counter[radii])
#     for j in range(6):
#         z = disk_z[j]
#         writer.writerow([f"Disk {j+1}"] + hit_counter[z])
