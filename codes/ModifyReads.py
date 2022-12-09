import numpy as np


def createGroups(numReads, numGroups):
    # Create one-to-one mappings from one groups of reads to other groups of reads
    # numReads - number of reads
    # numGroups - number of groups to divide the reads
    assert numReads % numGroups == 0

    permReads = np.arange(numReads)
    np.random.shuffle(permReads)
    groupSize = numReads // numGroups

    readGroups = np.reshape(permReads, (numGroups, groupSize))
    readMapAcrossGroups = dict()
    for i in range(groupSize):
        readMapAcrossGroups[str(readGroups[0][i])] = (readGroups[1:numGroups, i]).tolist()

    return readGroups, readMapAcrossGroups


def restrictedPerm(groupSize, readLength):
    res = np.empty(groupSize * readLength, dtype=int)
    perm = np.random.permutation(groupSize * readLength)
    for i in range(groupSize):
        res[np.sort(perm[i * readLength:(i + 1) * readLength])] = np.arange(i * readLength, (i + 1) * readLength)

    return res


def take2s_np(arr, idx_a, idx_b):
    return arr[idx_a].view(("U", 1))[idx_b].view(("U", len(idx_b)))[0]


def createSuperReads(reads, readGroups):
    numReads = len(reads)
    readLength = len(reads[0])
    numGroups = len(readGroups)
    groupSize = numReads // numGroups

    p = restrictedPerm(groupSize, readLength)

    ro = np.repeat(readGroups[0], readLength)
    readOrder = ro[p]

    superReads = []
    for i in range(numGroups):
        superReads.append(take2s_np(reads, readGroups[i], p))

    return readOrder, superReads
