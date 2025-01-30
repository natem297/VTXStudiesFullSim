from podio import root_io
import numpy as np
import os

##########################################################################################
#  this file is for printing the radius of hits to observe approximate layer radii
##########################################################################################

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/IDEA_wz3p6_ee_qq_ecm91p2_noXing/"
files = os.listdir(folder)

z = []

for filename in files[0:1]:

    print(f"starting {filename}")
    input_file_path = os.path.join(folder, filename)
    podio_reader = root_io.Reader(input_file_path)
    events = podio_reader.get("events")
    for event in events:
        for hit in event.get("VertexEndcapCollection"):
            z.append(hit.getPosition().z)

print(f"z: {sorted(z)}")
