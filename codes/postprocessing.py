import json
import os
import time
import numpy as np
import argparse
import Parameters
import editdistance

# Specify the path of the main directory
basePath = os.getcwd() + "/../"

# Specify input sizes for used for preprocessing and alignment
numReadsPerBatch = 100
templateLengthPerChunk = 1000000

parser = argparse.ArgumentParser()
parser.add_argument('--chrnum', type=int, required=True)
parser.add_argument('--tcn', type=int, required=True)
parser.add_argument('--rl', type=int, required=True)
parser.add_argument('--nrg', type=int, required=True)
parser.add_argument('--nrb', type=int)
parser.add_argument('--rbn', type=int)
args = parser.parse_args()

chrLengthFile = "inputs/chrLengths.txt"
with open(basePath + chrLengthFile, 'r') as f:
    chrLengthsString = f.read().splitlines()

argsChrNum = int(args.chrnum)
numTemplateChunks = int(chrLengthsString[int(args.chrnum - 1)][:-1])
templateChunkNum = int(args.tcn)
readLength = int(args.rl)
numReadGroups = int(args.nrg)
numReadBatches = int(args.nrb)
readBatchNum = int(args.rbn)

startTime = time.time()

print("\nStarting postprocessing for template chunk \#" + str(templateChunkNum) + " and read batch \#" + str(
    readBatchNum))

resultsPath = basePath + "results/chr" + str(argsChrNum) + '_' + str(numTemplateChunks) + '_reads_' + str(
    readLength) + "_" + str(numReadGroups) + '_' + str(numReadBatches) + "/"

outputPath = resultsPath + "template_chunk_" + str(templateChunkNum) + "_read_batch_" + str(readBatchNum) + "/"

# Load information
st = time.time()
print("Loading template and read information + alignments...")

with open(outputPath + 'alignments.json') as f:
    alignments = json.load(f)

localPath = resultsPath + "template_chunk_" + str(templateChunkNum) + "/"
SA = np.load(localPath + 'SA.npy')
templateChunk = np.load(localPath + 'templateChunk.npy')

localPath = resultsPath + "read_batch_" + str(readBatchNum) + "/"
with open(localPath + 'readMapAcrossGroups.json') as f:
    readMapAcrossGroups = json.load(f)

primes = np.load(outputPath + 'primes.npy')

localPath = basePath + "inputs/reads.txt"
with open(localPath, "r") as f:
    reads = f.read().splitlines()

currReads = reads[readBatchNum * numReadsPerBatch: (readBatchNum + 1) * numReadsPerBatch]

print("Finished loading information in ", round(time.time() - st, 2))

# Create primes array for group encoding
numPrimes = Parameters.numPrimes

# Create easy way to identify groups
groupDict = dict()
LStart = Parameters.LStart
numLs = Parameters.numLs
l_vals = np.arange(LStart, LStart + numLs)
r_vals = np.empty([numReadGroups, numPrimes * numLs], dtype=np.int64)
for i in range(numLs):
    r_vals[:, i * numPrimes:(i + 1) * numPrimes] = primes * l_vals[i]

AStart = Parameters.AStart
numAs = Parameters.numAs
a_vals = np.arange(AStart, AStart + numAs)
m_vals = np.empty([numReadGroups, numPrimes * numLs, numAs], dtype=np.int64)
for i in range(numAs):
    m_vals[:, :, i] = r_vals * a_vals[i]
for i in range(numReadGroups):
    for j in range(numPrimes * len(l_vals)):
        for k in range(len(a_vals)):
            groupDict[m_vals[i, j, k]] = i

##############################################
# Modules required for post-processing
##############################################

def getGroups(Fbar):
    result = np.empty_like(Fbar)
    result[:, 0] = SA[Fbar[:, 0]]
    tmp = Fbar[:, 1]
    u, inv = np.unique(tmp, return_inverse=True)
    result[:, 1] = np.array([groupDict[x] for x in u])[inv].reshape(tmp.shape)
    return result


def oneSNPAlignment(string1, string2):
    assert len(string1) == len(string2)
    if editdistance.eval(string1, string2) > 2:
        return False
    return True


def originalAlignments():
    origAlignInChunk = dict()
    valid_range_low = templateChunkNum * templateLengthPerChunk
    valid_range_up = (templateChunkNum + 1) * templateLengthPerChunk
    for i in range(readBatchNum * numReadsPerBatch, (readBatchNum + 1) * numReadsPerBatch):
        tmpAlign = int(origAlignInChunk[i])
        if valid_range_low <= tmpAlign < valid_range_up:
            origAlignInChunk[i - readBatchNum * numReadsPerBatch] = [tmpAlign - valid_range_low]
    print(origAlignInChunk)
    return origAlignInChunk


def postprocessing():
    finalAlignments = dict()
    for readNum in alignments:
        Fbar = np.array(alignments[readNum])
        if not np.any(Fbar):
            continue
        tmpAlign = getGroups(Fbar)
        candidateAlign = tmpAlign[np.where(tmpAlign[:, 1] == 0)][:, 0].tolist()
        for align in candidateAlign:
            if align + readLength <= templateLengthPerChunk:
                if oneSNPAlignment(templateChunk[align: align + readLength], currReads[int(readNum)]):
                    if int(readNum) in finalAlignments:
                        finalAlignments[int(readNum)].append(align)
                    else:
                        finalAlignments[int(readNum)] = [align]
        for g in range(1, numReadGroups):
            grp = np.where(tmpAlign[:, 1] == g)
            candidateAlign = tmpAlign[grp][:, 0].tolist()
            for align in candidateAlign:
                if align + readLength <= templateLengthPerChunk:
                    readNum2 = int(readMapAcrossGroups[readNum][g - 1])
                    if oneSNPAlignment(templateChunk[align: align + readLength], currReads[readNum2]):
                        if readNum2 in finalAlignments:
                            finalAlignments[readNum2].append(align)
                        else:
                            finalAlignments[readNum2] = [align]

    sorted_final_alignments = dict(sorted(finalAlignments.items()))
    return sorted_final_alignments


############################################
# BEGIN POSTPROCESSING
############################################

final_alignments = postprocessing()
with open(outputPath + "final_alignments.json", 'w') as outF:
    json.dump(final_alignments, outF)

print("Completed postprocessing in ", round(time.time() - startTime, 2))
