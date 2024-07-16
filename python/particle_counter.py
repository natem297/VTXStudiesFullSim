from podio import root_io
import csv

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")

pdg_ids = {11: "Electron", -11: "Positron", 12: "Electron Neutrino", 22: "Photon", 2112: "Neutron", 2212: "Proton", \
           111: "Pion 0", 211: "Pion +", -13: "Muon +", 14: "Muon Neutrino", -14: "Anti Muon Neutrino"}

particle_counter_list = []

for event in events:
    particle_counter = {particle: 0 for particle in pdg_ids.values()}
    other_particles_count = 0

    for particle in event.get("MCParticles"):
        pdg = particle.getPDG()

        if pdg in pdg_ids:
            particle_counter[pdg_ids[pdg]] += 1
        else:
            other_particles_count += 1


    particle_counter_list.append(particle_counter)

with open("data_tables/particle_counts.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    for particle in particle_counter:
        writer.writerow([particle] + [counter[particle] for counter in particle_counter_list])
