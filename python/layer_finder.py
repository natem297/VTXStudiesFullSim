from podio import root_io
import numpy as np
import os

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_wz3p6_ee_qq_ecm91p2"
files = os.listdir(folder)
event_count = len(files)

r = []

for filename in files[0:1]:

    print(f"starting {filename}")
    input_file_path = os.path.join(folder, filename)
    podio_reader = root_io.Reader(input_file_path)
    events = podio_reader.get("events")

    for event in events[0:5]:

        for hit in event.get("VertexBarrelCollection"):
            r.append(np.sqrt(hit.getPosition().x**2 + hit.getPosition().y**2))

print(f"r: {sorted(r)}")
