import numpy as np

def ScrambleTemplate(X, numPartitions):
    # Divide and scramble the template
    # X - Reference template
    # numPartitions - number of partitions to divide the template into
    assert len(X) % (2 * numPartitions) == 0

    size_parts = len(X) // (2 * numPartitions)
    template_parts = [X[i: i + size_parts] for i in range(0, len(X), size_parts)]

    # do not scramble the template at the moment
    piTemplate = np.arange(2 * numPartitions)
    templateParts = np.array([list(template_parts[piTemplate[2 * i]]) + list(template_parts[piTemplate[2 * i + 1]]) for i in range(numPartitions)])

    return templateParts, piTemplate


def ComputePermInv(pi):
    # Create inverse of random permutation pi
    inv = np.empty_like(pi)
    inv[pi] = np.arange(len(inv), dtype=inv.dtype)
    return inv


def CreateLprime(pil, ipif):
    # Create links between rows of F and L
    # pil - random permutation 2
    # ipif - inverse of random permutation 1
    pil_minus_one = np.where(pil > 0, pil-1, len(pil)-1)
    return ipif[pil_minus_one]


def CreateFprime(ipil, pif):
    # Create links between rows of L and F
    # pif - random permutation 1
    # ipil - inverse of random permutation 2
    return ipil[pif]


def CreateF(Wf, pif):
    # Create shuffled first column (F)
    # pif - random permutation 1
    # Wf - First column of rotation matrix
    return Wf[pif]


def CreateL(Wl, pil):
    # Create shuffled first column (L)
    # pil - random permutation 2
    # Wl - Last column of rotation matrix
    return Wl[pil]
