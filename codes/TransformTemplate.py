import BWTStructure as bwt
import json
import os
import time
import numpy as np


def TransformReferenceTemplate(templateChunk, local_path, chunkingTime):
    if not os.path.isdir(local_path):
        os.makedirs(local_path)

    st = time.process_time()
    Wf = np.append(templateChunk, '$')
    Wl = np.insert(templateChunk, 0, '$')
    q = len(templateChunk) + 1
    time_wf_wl = time.process_time() - st

    st = time.process_time()
    pif = np.random.permutation(q)
    pil = np.random.permutation(q)
    time_pif_pil = time.process_time() - st

    st = time.process_time()
    ipif = bwt.ComputePermInv(pif)
    ipil = bwt.ComputePermInv(pil)
    time_ipif_ipil = time.process_time() - st

    st = time.process_time()
    L = bwt.CreateL(Wl, pil)
    time_l = time.process_time() - st
    F = bwt.CreateF(Wf, pif)
    time_f_l = time.process_time() - st

    st = time.process_time()
    Fp = bwt.CreateFprime(ipil, pif)
    Lp = bwt.CreateLprime(pil, ipif)
    time_fp_lp = time.process_time() - st

    st = time.process_time()
    np.save(local_path + 'SA', pif)
    np.save(local_path + 'templateChunk', templateChunk)
    np.save(local_path + 'L', L)
    np.save(local_path + 'Fp', Fp)
    np.save(local_path + 'Lp', Lp)
    time_load = time.process_time() - st

    np.save(local_path + 'F', F)

    data = dict()  # write to JSON file
    data["time_wf_wl"] = time_wf_wl
    data["time_pif_pil"] = time_pif_pil
    data["time_ipif_ipil"] = time_ipif_ipil
    data["time_f_l"] = time_f_l
    data["time_fp_lp"] = time_fp_lp
    data["time_load"] = time_load
    data["time_chunking"] = chunkingTime

    print("Template report time: ", time_wf_wl + time_ipif_ipil + time_pif_pil + time_l + time_fp_lp)

    with open(local_path + "template_times.json", 'w') as outF:
        json.dump(data, outF)

    return
