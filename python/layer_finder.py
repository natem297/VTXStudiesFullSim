from podio import root_io
import ROOT
import math
import numpy as np
import os

folder = "/eos/experiment/fcc/users/j/jaeyserm/VTXStudiesFullSim/CLD_wz3p6_ee_qq_ecm91p2"
files = os.listdir(folder)
event_count = len(files)

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
