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
import glob


def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if length == 0:
        return float("nan")
    if length % 2 == 0:
        return (sorts[length // 2] + sorts[length // 2 - 1]) / 2.0
    return float(sorts[length // 2])


def main():
    if len(sys.argv) < 2:
        raise SystemExit(f"Usage: {sys.argv[0]} <sample_prefix_path>  (expects <prefix>.*.start.* files)")

    inFiles = glob.glob(sys.argv[1] + ".*.start.*")

    # Prepare bins per chromosome by their sizes
    binSize = 50000
    chrSizes = [
        16571, 249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
        159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878,
        107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520,
        48129895, 51304566, 155270560, 59373566
    ]

    chrArrays = []
    for c in chrSizes:
        chrArrays.append([0] * (c // binSize + 1))

    for inFile in inFiles:
        chrom = inFile.split(".")[-3]
        # Turn X into the right chromosome number in our arrays
        if chrom == "X":
            numChrom = 23
        else:
            numChrom = int(chrom)
        curChrArr = chrArrays[numChrom]
        with open(inFile, "r") as reader:
            for line in reader:
                line = line.strip()
                if not line:
                    continue
                pos = int(line) // binSize
                if 0 <= pos < len(curChrArr):
                    curChrArr[pos] += 1

    autosomals = []
    for chrom in chrArrays[1:23]:
        autosomals.extend(chrom)
    medA = median(autosomals)
    medX = median(chrArrays[23])

    sample_prefix = sys.argv[1]
    sample_name = sample_prefix.split("/")[-1]
    ff_metric = ((medA - medX) / medA * 2) * 100 if medA and medA == medA else float("nan")

    print(sample_prefix, sample_name, ff_metric, medA, medX)


if __name__ == "__main__":
    main()
