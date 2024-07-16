from podio import root_io
import ROOT
import numpy as np
import random

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
e0 = events[0]
inner_barrel = e0.get("VTXIBCollection")
monte_carlo = e0.get("MCParticles")

inner_radii = [np.sqrt((particle.getPosition().x)**2 + (particle.getPosition().y)**2) for particle in inner_barrel]

outer_radii = []
for event in events:
    outer_barrel = event.get("VTXOBCollection")
    for particle in outer_barrel:
        outer_radii.append(np.sqrt((particle.getPosition().x)**2 + (particle.getPosition().y)**2))

def kmeans(radii, k):
    """
    K means clustering function to group determine radii of layers.
    Inputs:
        radii, python list of detector hit radii.
        k, integer representing expected number of layers.
    Outputs:
        centroids, floats representing means of each cluster.
        clusters, python list of lists, one for each cluster.
    """
    centroids = random.sample(radii, k)
    converged = False

    while not converged:
        clusters = [[] for _ in range(k)]

        for r in radii:
            shortest_dist = abs(r-centroids[0])
            best_index = 0

            for i in range(1,k):
                dist = abs(r-centroids[i])

                if dist < shortest_dist:
                    shortest_dist = dist
                    best_index = i

            clusters[best_index].append(r)

        converged = True
        new_centroids = [np.mean(cluster) for cluster in clusters]
        for j in range(k):

            if new_centroids[j] - centroids[j] != 0.0:
                converged = False
                break

        centroids = new_centroids
    return centroids, clusters

for radii, k in [(inner_radii, 3), (outer_radii, 2)]:
    centroids, clusters = kmeans(radii, k)
    print(f"Clusters found with centroids: {centroids}")
    print(f"The size of each cluster is {[len(cluster) for cluster in clusters]}")
    print(f"The standard deviations of each cluster are {[np.std(cluster) for cluster in clusters]}")

mc_radii = [np.sqrt(particle.getVertex().x**2 + particle.getVertex().y**2) for particle in monte_carlo]
hist = ROOT.TH1F("radii", "Monte Carlo Radii", 600, 0, 150)
for particle in monte_carlo:
    hist.Fill(np.sqrt(particle.getVertex().x**2 + particle.getVertex().y**2))
canvas = ROOT.TCanvas("radii", "Monte Carlo Radii")
hist.Draw()
canvas.Update()
canvas.SaveAs("../plots/mc_radii.png")

disk_z = []
events = [events[i] for i in range(100)]
for event in events:
    for hit in event.get("VTXDCollection"):
        disk_z.append(hit.getPosition().z)

centroids, clusters = kmeans(disk_z, 6)
print(f"Clusters found with centroids: {centroids}")
print(f"The size of each cluster is {[len(cluster) for cluster in clusters]}")
print(f"The standard deviations of each cluster are {[np.std(cluster) for cluster in clusters]}")
