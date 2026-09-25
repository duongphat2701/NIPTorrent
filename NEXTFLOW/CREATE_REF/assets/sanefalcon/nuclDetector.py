# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (r.straver@vumc.nl)
#
# This file is part of SANEFALCON
# SANEFALCON is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.



#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import operator

areaSize = 73  # 60
sideSize = 20  # 25
padding = areaSize + sideSize
outerLen = 2 * sideSize
innerLen = areaSize * 2 + 1

if len(sys.argv) > 3:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt


def plotRegion(start, end, area, bpScores, centers):
    if len(centers) < 6 or len(area) < 2000:
        return

    plt.clf()

    left = len(area) // 2 - 1000
    right = len(area) // 2 + 1000

    smooth = []
    for i in range(len(area)):
        smooth.append(sum(area[max(0, i - 10): i + 10]) / 20.0)

    for center in centers:
        plt.axvspan(center[0] - 73, center[0] + 73, color="gray", alpha=0.3)
        plt.axvline(center[0], color="black")
    plt.plot(smooth)
    plt.plot(bpScores, "red")
    plt.title(f"Smoothed Read Start Count (chr21:{start + left}-{start + right})")
    plt.xlabel("Relative BP position")
    plt.ylabel("Count or Score")
    plt.xlim([left, right])
    plt.ylim([0, 10])
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(16, 2)
    plt.savefig(sys.argv[3] + "_" + str(start) + "-" + str(end) + ".plot.pdf", dpi=100)


def flush(area, endPoint):
    bins = [0] * padding + area + [0] * padding
    areaLen = len(area)
    startPoint = endPoint - areaLen + 1

    bpScores = []
    extra = []

    for i in range(padding, areaLen + padding):
        leftVal = sum(bins[i - areaSize - sideSize: i - areaSize])
        rightVal = sum(bins[i + areaSize + 1: i + areaSize + sideSize + 1])
        innerVal = sum(bins[i - areaSize: i + areaSize + 1])
        outerVal = leftVal + rightVal

        score = ((outerVal + 1) / float(outerLen)) / ((innerVal + 1) / float(innerLen))
        bpScores.append(score)
        extra.append([leftVal, innerVal, rightVal])

    def findCenters(start, end):
        if end - start < 1:
            return []
        maxIndex, bpMax = max(enumerate(bpScores[start:end]), key=operator.itemgetter(1))
        if bpMax > 1:
            tmpIndex = maxIndex
            while tmpIndex < (end - start - 1) and bpScores[tmpIndex + 1] == bpMax:
                tmpIndex += 1
            maxIndex = (maxIndex + tmpIndex) // 2

            left = start + maxIndex - innerLen
            right = start + maxIndex + innerLen + 1
            leftList = findCenters(start, left)
            rightList = findCenters(right, end)
            thisList = [start + maxIndex, bpMax]
            return leftList + [thisList] + rightList
        return []

    newCenters = findCenters(0, areaLen)
    allNucl.extend([[x[0] + startPoint] + x[1:] + extra[x[0]] for x in newCenters])

    if len(sys.argv) > 3:
        plotRegion(startPoint, endPoint, area, bpScores, newCenters)


curArea = [0]
lastPos = 0
maxDist = 190  # slightly over sliding window size
allNucl = []

with open(sys.argv[1], "r") as inFile:
    for line in inFile:
        line = line.strip()
        if not line:
            continue
        splitLine = line.split()
        position = int(splitLine[0])
        distance = position - lastPos

        if distance > maxDist:
            flush(curArea, lastPos)
            curArea = [1]
        else:
            curArea += [0 for _ in range(distance)]
            curArea[-1] += 1

        lastPos = position

flush(curArea, lastPos)

with open(sys.argv[2], "w") as output_file:
    for nucl in allNucl:
        output_file.write("\t".join([str(x) for x in nucl]) + "\n")
