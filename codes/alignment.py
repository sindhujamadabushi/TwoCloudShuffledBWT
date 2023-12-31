import time
import json
import numpy as np
import argparse
import os
import Parameters
import errno
from mpi4py import MPI

basePath = os.getcwd() + "/../"

# extract total number of template chunks
chrLengthFile = "inputs/chrLengths.txt"
with open(basePath + chrLengthFile, 'r') as f:
    chrLengthsString = f.read().splitlines()

parser = argparse.ArgumentParser()
parser.add_argument('--chrnum', type=int, required=True)
parser.add_argument('--tcn', type=int, required=True)
parser.add_argument('--rl', type=int, required=True)
parser.add_argument('--nrg', type=int, required=True)
parser.add_argument('--nrb', type=int, required=True)
parser.add_argument('--rbn', type=int, required=True)
args = parser.parse_args()

# get script arguments
chrNum = int(args.chrnum)
numTemplateChunks = int(chrLengthsString[chrNum-1][:-1])
templateChunkNum = int(args.tcn)
readLength = int(args.rl)
numReadGroups = int(args.nrg)
numReadBatches = int(args.nrb)
readBatchNum = int(args.rbn)

# specify input path and create output path
inputPath = basePath + "results/chr" + str(chrNum) + '_' + str(numTemplateChunks) + '_reads_' + str(readLength) + "_" + str(numReadGroups) + '_' + str(numReadBatches) + "/"
outputRelPath = "template_chunk_" + str(templateChunkNum) + "_read_batch_" + str(readBatchNum) + "/"
outputPath = inputPath + outputRelPath
try:
    os.mkdir(outputPath)
except OSError as e:
    if e.errno == errno.EEXIST:
        print('output path not created.')
    else:
        raise

# Load information for template chunk and read batch
startLoadTime = time.process_time()
print("\nLoading read and template information...")

# Read information
localPath = inputPath + "read_batch_" + str(readBatchNum) + "/"
readOrder = np.load(localPath + "readOrder.npy")
superReads = np.load(localPath + "superReads.npy")
readGroups = np.load(localPath + "readGroups.npy")

# Template information
localPath = inputPath + "template_chunk_" + str(templateChunkNum) + "/"
Lp = np.load(localPath + "Lp.npy")
Fp = np.load(localPath + "Fp.npy")
L = np.load(localPath + "L.npy")

timeLoad = time.process_time() - startLoadTime
print("Finished loading read and template information in ", round(timeLoad, 2))

# Start preprocessing operations in CloudL
stPreprocessingTime = time.process_time()
print("\nPerforming preprocessing operations for CloudL...")

# Create Ldict
Ldict = {'A': np.where(L == 'A')[0], 'C': np.where(L == 'C')[0], 'G': np.where(L == 'G')[0], 'T': np.where(L == 'T')[0]}

# Load primes
primes = np.load(localPath + 'primes.npy')

# Create easy way to identify groups
groupDict = dict()
LStart = Parameters.LStart
numLs = Parameters.numLs
numPrimes = Parameters.numPrimes
l_vals = np.arange(LStart, LStart + numLs)
r_vals = np.empty([numReadGroups, numPrimes*numLs], dtype=np.int64)
for i in range(numLs):
     r_vals[:, i*numPrimes:(i+1)*numPrimes] = primes*l_vals[i]

AStart = Parameters.AStart
numAs = Parameters.numAs
a_vals = np.arange(AStart, AStart + numAs)
m_vals = np.empty([numReadGroups, numPrimes*numLs, numAs], dtype=np.int64)
for i in range(numAs):
     m_vals[:, :, i] = r_vals*a_vals[i]

for i in range(numReadGroups):
     for j in range(numPrimes*len(l_vals)):
         for k in range(len(a_vals)):
             groupDict[m_vals[i, j, k]] = i

timePreprocessing = time.process_time() - stPreprocessingTime
print("Finished performing preprocessing operations for CloudL in ", round(timePreprocessing, 2))

CloudLNumRandomIndices = Parameters.CloudLNumRandomIndices

############################################
# MODULES REQUIRED FOR ALIGNMENT
############################################


def getGroups(Fbar):
    result = np.empty_like(Fbar)
    result[:, 0] = Fp[Fbar[:, 0]]
    tmp = Fbar[:, 1]
    u, inv = np.unique(tmp, return_inverse=True)
    result[:, 1] = np.array([groupDict[x] for x in u])[inv].reshape(tmp.shape)
    return result


def convertLbar2Fbar(Lbar):
    a = np.random.randint(AStart, AStart + numAs, len(Lbar))
    p = Lbar[:, 1]
    q = Lbar[:, 0]
    m = (a * (p - 1))
    t = Lp[q]
    result = np.column_stack((t, m))
    return result


def generateGroupEncoding(groupNum, numIndices):
    ls = np.random.randint(low=LStart, high=LStart + numLs, size=numIndices)
    qs = np.random.choice(primes[groupNum], numIndices, replace=True)
    result = (qs*ls)+1
    return result


def firstIterationInCloudL(iterNum):
    bases = [superReads[g][iterNum] for g in range(numReadGroups)]
    groupSizes = [len(Ldict[bases[g]]) for g in range(numReadGroups)]
    startPos = 0
    Lbar = np.empty([sum(groupSizes), 2], dtype=np.int64)
    for g in range(numReadGroups):
        if g > 0:
            startPos += groupSizes[g-1]
        Lbar[startPos:startPos+groupSizes[g], 0] = Ldict[bases[g]]
        Lbar[startPos:startPos+groupSizes[g], 1] = generateGroupEncoding(g, len(Ldict[bases[g]]))
    return Lbar


def convertFbar2Lbar_and_obfuscate(Fbar, iterNum, mismatchStatus, cutoffStatus):
    intResult = np.empty((0, 2), np.int64)
    currReadMismatchStatus = np.full(numReadGroups, 0)
    readGroupPerm = np.random.permutation(numReadGroups)

    for g in readGroupPerm:
        if not cutoffStatus[g]:
            if not mismatchStatus[g]:
                ind1 = np.where(Fbar[:, 1] == g)[0]
                tmp1 = np.where(L[Fbar[ind1, 0]] == superReads[g][iterNum])[0]

            else:  # mismatch in last iteration, so ignore all indices
                tmp1 = np.array([])

            if np.size(tmp1):  # only happens if no mismatch in previous iteration and current iteration
                tmpResult = np.empty((len(tmp1), 2), np.int64)
                tmpResult[:, 0] = Fbar[ind1[tmp1], 0]
                tmpResult[:, 1] = generateGroupEncoding(g, len(tmp1))
                intResult = np.append(intResult, tmpResult, axis=0)

            else:
                if not mismatchStatus[g]:
                    currReadMismatchStatus[g] = 1

                else:     # only happens if mismatch in previous iter, so didn't check previous iter indices
                    tmp1 = Ldict[superReads[g][iterNum]]
                    tmpResult = np.empty((len(tmp1), 2), np.int64)
                    tmpResult[:, 0] = tmp1
                    tmpResult[:, 1] = generateGroupEncoding(g, len(tmp1))

                    intResult = np.append(intResult, tmpResult, axis=0)

        else:
            currReadMismatchStatus[g] = 1
    return intResult, currReadMismatchStatus


def checkCutoff(countMismatch, cutoff):
    result = np.full(countMismatch.shape, 0)
    ind = np.where(countMismatch >= cutoff)
    result[ind] = 1
    return result


def saveResults(resultsDict):
    timesDict = dict()
    timesDict["Time in CloudF"] = resultsDict["Time in CloudF"]
    timesDict["Time in CloudL"] = resultsDict["Time in CloudL"]
    timesDict["Total Time"] = resultsDict["Total Time"]

    with open(outputPath + "alignment_times.json", 'w') as outF:
        json.dump(timesDict, outF)

    with open(outputPath + 'Communication_overhead_F.txt', 'w') as f:
        f.write(str(resultsDict["commOverheadF"]))

    with open(outputPath + 'Communication_overhead_L.txt', 'w') as f:
        f.write(str(resultsDict["commOverheadL"]))

    finalAlignments = resultsDict["finalAlignments"]
    finalAlignmentsList = {key: finalAlignments[key].tolist() for key in list(finalAlignments.keys())}
    with open(outputPath + "alignments.json", "w") as outF:
        json.dump(finalAlignmentsList, outF)
    return


############################################
# CORE ALIGNMENT ALGORITHM
############################################

alignments = {str(readGroups[0][i]): np.array([]) for i in range(len(readGroups[0]))}
countMismatches = {str(readGroups[0][i]): np.zeros(numReadGroups, dtype=np.int64) for i in range(len(readGroups[0]))}
readMismatchStatus = {str(readGroups[0][i]): np.full(numReadGroups, 0) for i in range(len(readGroups[0]))}

CloudFReadMismatchStatus = np.full(numReadGroups, 0)
CloudFCutoffStatus = np.full(numReadGroups, 0)
CloudL_readMismatchStatus = np.full(numReadGroups, 0)
inexactAlignmentCutoff = Parameters.inexactAlignmentCutoff

readOrderLength = len(readOrder)

print("\nStarting the alignment algorithm...")

commOverheadF = []
commOverheadL = []

timeInCloudF = 0.0
timeInCloudL = 0.0
timeCloudL = 0.0

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

Fbar = np.empty((0, numReadGroups), np.int64)
Lbar = np.empty((0, numReadGroups), np.int64)

# TEMPORARY
CloudLResult = []

# stTotalCommunicationTime = 0.0
stTimeInCloudF = 0.0
stTotalTime = 0.0
endTotaltime = 0.0

stTotalCommunicationTime = time.process_time()
for iterNum in range(readOrderLength - 1, -2, -1):

    # CloudF
    if rank == 0:
        if iterNum == readOrderLength - 1:
            pass
        else:
            CloudLResult = comm.recv(source=1)
            stTimeInCloudF = time.time()
            CloudL_readMismatchStatus = CloudLResult[0].astype(np.int64)
            timeInCloudL = float(CloudLResult[1][0])
            Lbar = CloudLResult[2:].astype(np.int64)

            commOverheadL.append(len(Lbar))

            alignments[str(readOrder[iterNum + 1])] = convertLbar2Fbar(Lbar)
            readMismatchStatus[str(readOrder[iterNum + 1])] = CloudL_readMismatchStatus
            for g in range(numReadGroups):
                if CloudL_readMismatchStatus[g] and (not CloudFCutoffStatus[g]):
                    countMismatches[str(readOrder[iterNum + 1])][g] += 1
            timeInCloudF += time.time() - stTimeInCloudF

        if iterNum == -1:
            resultsDict = dict()
            resultsDict["Time in CloudF"] = timeInCloudF
            resultsDict["Time in CloudL"] = timeInCloudL
            resultsDict["Total Time"] = time.process_time() - stTotalCommunicationTime
            resultsDict["commOverheadF"] = commOverheadF
            resultsDict["commOverheadL"] = commOverheadL
            resultsDict["finalAlignments"] = alignments
            saveResults(resultsDict)
            comm.Abort()
            break

        # Prepare for transferring information to CloudL
        stTimeInCloudF = time.time()
        Fbar = alignments[str(readOrder[iterNum])]
        CloudFReadMismatchStatus = readMismatchStatus[str(readOrder[iterNum])]
        CloudFCutoffStatus = checkCutoff(countMismatches[str(readOrder[iterNum])], inexactAlignmentCutoff)
        CloudLTime = [timeInCloudL, 0.0]
        if Fbar.size:
            CloudFResult = np.vstack((CloudFReadMismatchStatus, CloudFCutoffStatus, CloudLTime, Fbar))
        else:
            CloudFResult = np.vstack((CloudFReadMismatchStatus, CloudFCutoffStatus, CloudLTime))
        commOverheadF.append(len(Fbar))
        comm.send(CloudFResult, dest=1)

    # CloudL
    elif rank == 1:
        CloudFResult = comm.recv(source=0)
        stTimeInCloudL = time.time()
        CloudL_readMismatchStatus = np.full(numReadGroups, 0)
        CloudFReadMismatchStatus = CloudFResult[0].astype(np.int64)
        CloudFCutoffStatus = CloudFResult[1].astype(np.int64)
        timeInCloudL = float(CloudFResult[2][0])
        if len(CloudFResult) > 3:
            Fbar = CloudFResult[3:].astype(np.int64)
        else:
            Fbar = np.empty((0, numReadGroups), np.int64)

        if len(Fbar) == 0 and (not np.all(CloudFCutoffStatus)):
            Lbar = firstIterationInCloudL(iterNum)
        else:
            Fbar = getGroups(Fbar)
            Lbar, CloudL_readMismatchStatus = convertFbar2Lbar_and_obfuscate(Fbar, iterNum,
                                                                             CloudFReadMismatchStatus,
                                                                             CloudFCutoffStatus)
        timeInCloudL += time.time() - stTimeInCloudL
        CloudLTime = [timeInCloudL, 0.0]
        CloudLResult = np.vstack((CloudL_readMismatchStatus, CloudLTime, Lbar))
        comm.send(CloudLResult, dest=0)
