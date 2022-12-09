import ModifyReads as mr
import json
import time
import numpy as np
import os


def TransformReads(reads, numReadGroups, local_path):
    if not os.path.isdir(local_path):
        os.makedirs(local_path)

    numReads = len(reads)
    assert all([len(reads[i]) == len(reads[i+1]) for i in range(numReads-1)])

    st = time.process_time()
    reads = np.array(reads)
    readGroups, readMapAcrossGroups = mr.createGroups(numReads, numReadGroups)
    time_rg = time.process_time() - st

    st = time.process_time()
    readOrder, superReads = mr.createSuperReads(reads, readGroups)
    time_ro_sr = time.process_time() - st

    st = time.process_time()
    np.save(local_path + 'readOrder', readOrder)
    np.save(local_path + 'superReads', superReads)
    np.save(local_path + 'readGroups', readGroups)
    with open(local_path + 'readMapAcrossGroups.json', 'w') as f:
        json.dump(readMapAcrossGroups, f)
    time_load = time.process_time() - st

    data = dict()
    data["time_rg"] = time_rg
    data["time_ro_sr"] = time_ro_sr
    data["time_load"] = time_load

    print("Read report time: ", time_rg + time_ro_sr)

    with open(local_path + "read_times.json", 'w') as outF:
        json.dump(data, outF)

    return
