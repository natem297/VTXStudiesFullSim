from podio import root_io
import ROOT
import numpy as np
import math

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/IDEA_01_v03_pairs_all.root"
podio_reader = root_io.Reader(input_file_path)

events = podio_reader.get("events")
layer_radii = [14, 23, 34.5, 141, 316]

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

def match_finder(particle):

    xi = particle.getVertex().x
    yi = particle.getVertex().y
    zi = particle.getVertex().z

    px = particle.getMomentum().x
    py = particle.getMomentum().y
    pz = particle.getMomentum().z

    theta_i = theta(px, py, pz)
    phi_i = phi(px, py)

    for hit in hits[14]:
        if delta_squared(theta_i, phi_i, xi, yi, zi, hit) < 0.01:
            return True
    return False

hist = ROOT.TH1F("momentum", "Guinea Pig MC Particles with Layer 1 Match", 150, 0, 30)
for e in range(100):
    print(f"starting event {e}")
    event = events[e]
    hits = {coord: [] for coord in layer_radii}

    # categorizes all hits by layer
    for collection in ["VTXIBCollection", "VTXOBCollection"]:
        for hit in event.get(collection):
            hits[radius(hit)].append(hit)

    for particle in event.get("MCParticles"):
        if particle.getVertex().x**2 + particle.getVertex().y**2 > 169:
            continue
        if match_finder(particle):
            px = particle.getMomentum().x
            py = particle.getMomentum().y
            pz = particle.getMomentum().z
            hist.Fill(np.sqrt(px**2 + py**2 + pz**2)*1000)

hist.GetXaxis().SetTitle("Momentum (MeV)")
hist.GetYaxis().SetTitle("Number of Particles")
canvas = ROOT.TCanvas("momentum", "MC Particle Momentum")
hist.Draw()
canvas.Update()
canvas.SaveAs("../plots/monte_carlo/gp_mc_momentum_layer1_match.png")
