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


def main():
    if len(sys.argv) < 3:
        raise SystemExit(f"Usage: {sys.argv[0]} <model_file> <profile_file>")

    with open(sys.argv[1], "r") as referenceFile:
        correlations = [float(x) for x in referenceFile.readline().split()]
        scalars = [float(x) for x in referenceFile.readline().split()]

    with open(sys.argv[2], "r") as profileFile:
        profile = [float(x) for x in profileFile.readline().split(",")[:-1]]
        profile.reverse()
        profile.extend([float(x) for x in profileFile.readline().split(",")[1:-1]])

    totalReads = sum(profile)
    normProfile = [x / totalReads for x in profile] if totalReads else [0.0 for _ in profile]

    summed = 0.0
    for i, val in enumerate(normProfile):
        summed += val * correlations[i]
    fetalFraction = scalars[0] * summed + scalars[1]

    if len(correlations) == len(normProfile):
        print("Fetal Fraction:\t", fetalFraction)
    else:
        print("ERROR: correlation and sample profiles are not aligned")
    print("Nucleosome Profile:\t" + "\t".join([str(x) for x in normProfile]))


if __name__ == "__main__":
    main()
