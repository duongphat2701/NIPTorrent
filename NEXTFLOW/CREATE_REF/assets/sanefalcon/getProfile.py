# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (r.straver@vumc.nl)
#
# This file is part of SANEFALCON
# SANEFALCON is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.



# 1 Peak file as reference for nucleosome positions
# 2 Read start pos file for a sample
# 3 Reverse: 0 or 1
# 4 Output basename

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

minSide = 25
minCenter = 100

areaSize = 73.0
sideSize = 20.0

filtOut = 0


def nuclFilt(prop):
    global filtOut
    if min(prop[0], prop[2]) / sideSize > prop[1] / 147.0:
        return True
    filtOut += 1
    return False


def loadNucl(nuclLine):
    return [int(nuclLine[0]), float(nuclLine[1])]


peaks = [loadNucl(line.split()) for line in open(sys.argv[1]) if nuclFilt([int(x) for x in line.split()[2:]])]
reads = [int(line.strip()) for line in open(sys.argv[2]) if line.strip()]
shift = int(sys.argv[3])
print(len(peaks), filtOut)

rev = False
if shift != 0:
    rev = True
    peaks.reverse()
    reads.reverse()

peaks.append([-1, 0])

maxDist = 147
sumPeak = [0.0 for _ in range(maxDist)]
read = reads[0] if reads else 0
j = 0
nuclHit = []

if rev:
    for peakPair in peaks:
        peak = peakPair[0]
        thisPeak = [0.0 for _ in range(maxDist)]
        while j < len(reads) and read >= peak:
            if read < peak + maxDist:
                thisPeak[read - peak] += 1
            j += 1
            if j >= len(reads):
                break
            read = reads[j]
        thisSum = float(sum(thisPeak))
        if thisSum == 0:
            continue
        nuclHit.append(peak)
        sumPeak = [sumPeak[x] + thisPeak[x] for x in range(len(thisPeak))]
    nuclHit.reverse()
else:
    for peakPair in peaks:
        peak = peakPair[0]
        thisPeak = [0.0 for _ in range(maxDist)]
        while j < len(reads) and read <= peak:
            if read > peak - maxDist:
                thisPeak[peak - read] += 1
            j += 1
            if j >= len(reads):
                break
            read = reads[j]
        thisSum = float(sum(thisPeak))
        if thisSum == 0:
            continue
        nuclHit.append(peak)
        sumPeak = [sumPeak[x] + thisPeak[x] for x in range(len(thisPeak))]

with open(sys.argv[4] + ".np", "w") as outPeaks:
    outPeaks.write(",".join([str(x) for x in sumPeak]))
