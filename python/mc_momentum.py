from podio import root_io
import ROOT
import numpy as np

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)
events = podio_reader.get("events")

def momentum(particle):
    """
    Calculates magnitude of a partricle's momentum.
    Inputs: particle: MCParticle object.
    Outputs: momentum: float representing momentum magnitude.
    """
    return np.sqrt(particle.getMomentum().x**2 + particle.getMomentum().y**2 + particle.getMomentum().z**2)

max_p = 0
min_p = float("Inf")
for event in events:
    for particle in event.get("MCParticles"):
        p = momentum(particle)
        max_p = max(max_p, p)
        min_p = min(min_p, p)

print(f"Maximum momentum: {max_p} GeV")
print(f"Minimum momentum: {min_p} GeV")
