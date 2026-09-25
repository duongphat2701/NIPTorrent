# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (r.straver@vumc.nl)
#
# This file is part of SANEFALCON
# SANEFALCON is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.



# Arguments:
# 1 Training nucleosome profiles
# 2 Training reference data
# 3 Output basename
# 4 Test nucleosome profiles
# 5 Test reference data

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import glob
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import numpy
import random
import pickle

# ---------------------------------------------------------------------------- #
# Data loading functions
# ---------------------------------------------------------------------------- #
ignoredSeries = []  # e.g. ['serie23_140707','serie12_140428']


def loadNuclFile(nuclFile):
    samples = dict()
    coverages = dict()
    with open(nuclFile) as sampleFile:
        for line in sampleFile:
            splitLine = line.strip().split(",")
            if len(splitLine) < 3:
                continue
            sampleName = splitLine[0].split("/")[-1].split(".")[0].split("-")[-1]
            # Keep original behavior: uses argv[6], argv[7] as slice bounds
            values = [float(x) for x in splitLine[int(sys.argv[6]): int(sys.argv[7])]]
            valSum = sum(values)
            if valSum == 0:
                continue
            samples[sampleName] = [x / valSum for x in values]
            coverages[sampleName] = valSum
    return samples, coverages


def loadRefFile(refFile):
    series = dict()
    reference = dict()
    girls = dict()
    bads = dict()
    with open(refFile) as referenceFile:
        for line in referenceFile:
            splitLine = line.strip().split(" ")
            splitLine[0] = splitLine[0].split("-")[-1]

            if len(splitLine) < 3:
                continue

            if float(splitLine[-2]) > 30 or float(splitLine[-2]) < 3:
                continue
            if splitLine[0] in ignoredSeries:
                continue

            if splitLine[2] == "Male":
                reference[splitLine[0]] = float(splitLine[1])
                series[splitLine[0]] = "Training"
            elif splitLine[2] == "Female":
                girls[splitLine[0]] = float(splitLine[1])
            elif splitLine[2] == "BAD":
                bads[splitLine[0]] = float(splitLine[1])
            else:
                print(splitLine[2])
    return reference, series, girls, bads


def splitByReference(samples, reference):
    overlap = [ref for ref in reference if ref in samples]
    overlap.sort()
    print(len(overlap), "samples overlap")

    noOverlap = [sample for sample in samples if sample not in overlap]
    noOverlap.sort()
    print(len(noOverlap), "samples noOverlap")

    return overlap, noOverlap


# ---------------------------------------------------------------------------- #
# Data preparation functions
# ---------------------------------------------------------------------------- #
def getCorrelationProfile(samples, reference, trainingSet, yVals):
    correlations = []
    for i in range(len(samples[trainingSet[0]])):
        bpVals = [samples[x][i] for x in trainingSet]
        correlations.append(pearsonr(bpVals, yVals)[0])
    return correlations


def getNuclRatio(sample, correlations):
    sampleVal = 0.0
    for i, val in enumerate(sample):
        sampleVal += val * correlations[i]
    return sampleVal, sample


def getNuclRatios(names, sampleSet, correlations):
    selectedRegions = []
    sampleScores = []
    for sample in names:
        a, b = getNuclRatio(sampleSet[sample], correlations)
        sampleScores.append(a)
        selectedRegions.append(b)
    return sampleScores, selectedRegions


def getAreaScores(samples):
    scores = []
    for thisSample in samples:
        neg = sum(thisSample[30:60])
        pos = sum(thisSample[80:125])
        scores.append([pos, neg])
    return scores


# ---------------------------------------------------------------------------- #
# Data processing functions
# ---------------------------------------------------------------------------- #
def getErrorRate(prediction, reference):
    errors = []
    for i, val in enumerate(prediction):
        error = (val - reference[i])
        errors.append(abs(error))
    return numpy.median(errors)


def testPolyFit(samples, reference, p, prefix):
    fitSamples = []
    for val in samples:
        fitSamples.append(p[0] * val + p[1])
    print(prefix + " polyFit: Pearson:", pearsonr(fitSamples, reference))
    print(prefix + " polyFit: errorRate:", getErrorRate(fitSamples, reference))
    return fitSamples


def trainPolyFit(samples, reference):
    p = numpy.polyfit(samples, reference, 1)
    return p


def testLinearModel(samples, reference, clf, prefix):
    predicted = clf.predict(samples)
    print(prefix + " linearModel: Pearson:", pearsonr(predicted, reference))
    print(prefix + (" linearModel: Residual sum of squares: %.2f" % numpy.mean((predicted - reference) ** 2)))
    print(prefix + (" linearModel: Variance score: %.2f" % clf.score(samples, reference)))
    print(prefix + " linearModel: errorRate:", getErrorRate(predicted, reference))
    return predicted


def trainLinearModel(samples, reference):
    from sklearn import linear_model
    clf = linear_model.LinearRegression()
    clf.fit(samples, reference)
    return clf


# ---------------------------------------------------------------------------- #
# Plotting functions (unchanged except Py3-safe integer divisions)
# ---------------------------------------------------------------------------- #
def plotCorrelation(correlations, outFile):
    plt.figure(figsize=(16, 2))
    plt.title("Correlation per BP in Artificial Nucleosome on Training Set")
    plt.ylabel("Pearson-Correlation Score")
    plt.xlabel("Aligned Nucleosome BP Position")
    plt.plot(correlations)
    plt.xlim([0, 292])
    center = 147 - 1
    plt.xticks(
        [0, center - 93, center - 73, center, center + 73, center + 93, len(correlations) - 1],
        ["\nUpstream", "93", "73\nStart", "0\nCenter", "73\nEnd", "93", "\nDownstream"],
    )
    for x in [center - 93, center - 73, center, center + 73, center + 93]:
        plt.axvline(x=x, linewidth=1, ls="--", color="k")
    plt.savefig(outFile, dpi=100)


def plotScatter(trainX, trainY, testX, testY, outFile, plotName, recolor):
    recoloredX = [testX[x] for x in recolor]
    recoloredY = [testY[x] for x in recolor]
    normalX = [testX[x] for x in range(len(testX)) if x not in recolor]
    normalY = [testY[x] for x in range(len(testY)) if x not in recolor]

    plt.figure()
    plt.scatter(trainX, trainY, color="blue")
    plt.scatter(normalX, normalY, color="red")
    plt.scatter(recoloredX, recoloredY, color="green")

    plt.xlim([0, 25])
    plt.ylim([0, 25])

    plt.title(plotName)
    plt.xlabel("FF using Nucleosomes")
    plt.ylabel("FF using Y-Chrom")

    plt.savefig(outFile + ".pdf", dpi=100)


def plotProfiles(training, testing, outFile, correlations=None):
    correlations = correlations or []
    train = [sum(x) for x in map(list, zip(*training))]
    test = [sum(x) for x in map(list, zip(*testing))]

    train = [x / sum(train) for x in train]
    test = [x / sum(test) for x in test]

    plt.figure(figsize=(16, 2))
    plt.plot(train)
    plt.plot(test, color="red")

    plt.xlim([0, 292])
    center = 147 - 1
    plt.xticks(
        [0, center - 93, center - 73, center, center + 73, center + 93, len(train) - 1],
        ["\nUpstream", "93", "73\nStart", "0\nCenter", "73\nEnd", "93", "\nDownstream"],
    )
    for x in [center - 93, center - 73, center, center + 73, center + 93]:
        plt.axvline(x=x, linewidth=1, ls="--", color="k")

    plt.title("Aligned Nucleosome Profile")
    plt.xlabel("Aligned Nucleosome BP Position")
    plt.ylabel("Frequency")

    plt.savefig(outFile + ".profiles.pdf", dpi=100)


# ---------------------------------------------------------------------------- #
# Main
# ---------------------------------------------------------------------------- #
def main(argv=None):
    if argv is None:
        argv = sys.argv

    print("- Training stage:")

    samples, coverages = loadNuclFile(argv[1])
    reference, series, girls, bads = loadRefFile(argv[2])

    trainingSet, trainingSetNo = splitByReference(samples, reference)
    covs = [coverages[x] for x in trainingSet]
    yVals = [reference[x] for x in trainingSet]

    correlations = getCorrelationProfile(samples, reference, trainingSet, yVals)
    smoothedCorrelations = [numpy.mean(correlations[max(i - 4, 0): i + 5]) for i in range(len(correlations))]

    sampleScores, selectedRegions = getNuclRatios(trainingSet, samples, correlations)
    profiles = selectedRegions

    selectedRegions = getAreaScores(selectedRegions)

    polyFit = trainPolyFit(sampleScores, yVals)
    linearModel = trainLinearModel(selectedRegions, yVals)

    trPolyFit = testPolyFit(sampleScores, yVals, polyFit, "Train")
    trLinearModel = testLinearModel(selectedRegions, yVals, linearModel, "Train")

    fittedVals = [polyFit[0] * val + polyFit[1] for val in sampleScores]
    print(" ".join([str(x) for x in polyFit]))

    plt.scatter(fittedVals, yVals)
    plt.xlim([0, 25])
    plt.savefig(sys.argv[3] + ".direct2.pdf", dpi=100)

    with open(sys.argv[3] + ".model", "w") as modelFile:
        modelFile.write(" ".join([str(x) for x in correlations]))
        modelFile.write("\n")
        modelFile.write(" ".join([str(x) for x in polyFit]))

    if len(argv) >= 6:
        print("\n- Testing Stage:")
        newSamples, newCoverages = loadNuclFile(argv[4])
        newReference, newSeries, newGirls, newBads = loadRefFile(argv[5])

        testSet, testSetNo = splitByReference(newSamples, newReference)
        newCovs = [newCoverages[x] for x in testSet]
        newYVals = [newReference[x] for x in testSet]
        print(testSet)

        singleRun = ["15P0008B", "15P0135A", "15P0136A", "15P0137A", "15P0138A", "15P0139A",
                     "15P0140A", "15P0148A", "15P0149A", "15P0150A", "15P0151A", "15P0152A"]
        recolor = [i for i, val in enumerate(testSet) if val in singleRun]

        newSampleScores, newSelectedRegions = getNuclRatios(testSet, newSamples, correlations)
        newProfiles = newSelectedRegions
        newSelectedRegions = getAreaScores(newSelectedRegions)

        tePolyFit = testPolyFit(newSampleScores, newYVals, polyFit, "Test")
        teLinearModel = testLinearModel(newSelectedRegions, newYVals, linearModel, "Test")

        plotScatter(trPolyFit, yVals, tePolyFit, newYVals, argv[3] + ".polyfit", "Nucleosome based prediction", recolor)
        plotScatter(trLinearModel, yVals, teLinearModel, newYVals, argv[3] + ".linearmodel", "linearmodel", recolor)
        plotScatter(sampleScores, yVals, newSampleScores, newYVals, argv[3] + ".direct", "direct", recolor)
        plotScatter(sampleScores, covs, newSampleScores, newCovs, argv[3] + ".cov", "coverages", recolor)

        plotProfiles(profiles, newProfiles, argv[3], correlations)

    plotCorrelation(correlations, argv[3] + ".correlations.pdf")
    plotCorrelation(smoothedCorrelations, argv[3] + ".smoothedcorrelations.pdf")


if __name__ == "__main__":
    main()
