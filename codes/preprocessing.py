import TransformReads as tr
import TransformTemplate as tt
import BWTStructure as bwt
import Parameters
import numpy as np
import time
import argparse
import os

basePath = os.getcwd() + "/../inputs/"

parser = argparse.ArgumentParser()
parser.add_argument('--chrnum', type=int, required=True)
parser.add_argument('--nr', type=int, required=True)
parser.add_argument('--rl', type=int, required=True)
parser.add_argument('--nrg', type=int, required=True)
parser.add_argument('--nrb', type=int, required=True)
args = parser.parse_args()

# get script arguments
chrNum = int(args.chrnum)
numReads = int(args.nr)
readLength = int(args.rl)
numReadGroups = int(args.nrg)
numReadBatches = int(args.nrb)
numReadsPerBatch = numReads//numReadBatches

startTime = time.process_time()

# extract total number of template chunks
chrLengthFile = "chrLengths.txt"
with open(basePath + chrLengthFile, 'r') as f:
    chrLengthsString = f.read().splitlines()

numTemplateChunks = int(chrLengthsString[int(args.chrnum-1)][:-1])

# create a path for storing results
resultsPath = basePath + "../results/chr" + str(chrNum) + '_' + str(numTemplateChunks) + '_reads_' + str(readLength) + "_" + str(numReadGroups) + '_' + str(numReadBatches) + "/"
if not os.path.isdir(resultsPath):
    os.makedirs(resultsPath)

# extract the template
templateFile = basePath + "chr" + str(chrNum) + ".fa"
with open(templateFile, "r") as f:
    template = f.read().splitlines()
    template = template[1]

# extract reads
readFile = basePath + "reads.txt"
with open(readFile, "r") as f:
    reads = f.read().splitlines()

# Create template chunks
print("\nStarting to process chromosome ", chrNum, "...")
st = time.process_time()
scrambledTemplateChunks, scramblePermutation = bwt.ScrambleTemplate(template, numTemplateChunks)
templateChunkTime = time.process_time() - st
print("\nTime for creating all template chunks for chromosome ", chrNum, ": ", round(time.process_time() - st, 2), "\n")

# Load primes
with open(basePath + 'primes.txt', 'r') as f:
    primesList = list(map(int, f.readlines()))

# Transform template
for templateChunk in range(numTemplateChunks):
    localPath = resultsPath + "template_chunk_" + str(templateChunk) + "/"
    st = time.process_time()
    tt.TransformReferenceTemplate(scrambledTemplateChunks[templateChunk], localPath, templateChunkTime/numTemplateChunks)
    templateTime = time.process_time() - st
    print("Template Chunk #", templateChunk, " time: ", round(templateTime, 2))

    # Create primes
    numPrimes = Parameters.numPrimes
    primes = np.reshape(np.random.choice(primesList, numPrimes*numReadGroups, replace=False), (numReadGroups, numPrimes))
    np.save(localPath + 'primes', primes)

# Transform reads
print("\nStarting to process reads...")
for readBatch in range(numReadBatches):
    localPath = resultsPath + "read_batch_" + str(readBatch) + "/"
    st = time.process_time()
    currReads = reads[readBatch*numReadsPerBatch: (readBatch+1)*numReadsPerBatch]
    tr.TransformReads(currReads, numReadGroups, localPath)
    readsTime = time.process_time() - st
    print('Read Batch #', readBatch, ' time: ', round(readsTime, 2))

print('\nTotal preprocessing time: ', round(time.process_time() - startTime, 2))
