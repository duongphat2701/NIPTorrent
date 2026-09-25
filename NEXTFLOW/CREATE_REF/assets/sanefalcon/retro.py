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
import argparse

parser = argparse.ArgumentParser(
    description="Convert any stream of reads to a filtered stream (RETRO-like). Defaults assume SAM format.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument("-binsize", type=int, default=1000000, help="binsize used for samples")
parser.add_argument("-retdist", type=int, default=4, help="max bp difference between sequential reads to consider same tower")
parser.add_argument("-retthres", type=int, default=4, help="threshold when a group of reads is a tower and will be removed")
parser.add_argument("-colchr", type=int, default=2, help="column containing chromosome, default is for sam format")
parser.add_argument("-colpos", type=int, default=3, help="column containing read start position, default is for sam format")

args = parser.parse_args()

binsize = args.binsize
minShift = args.retdist
threshold = args.retthres
chrColumn = args.colchr
startColumn = args.colpos


def flush(fullBuff):
    """
    Print the buffered lines if buffer size is not considered a tower.
    Keeps original behavior: if stairSize <= threshold OR threshold < 0, print.
    """
    stairSize = len(fullBuff)
    if stairSize <= threshold or threshold < 0:
        for line in fullBuff:
            sys.stdout.write(line.rstrip("\n") + "\n")


prevWords = ["0"] * 10
fullBuff = []
readBuff = []  # kept for conceptual parity (not strictly needed)

for line in sys.stdin:
    curWords = line.split()
    if not curWords:
        continue

    # Not ndup, flush and start new stair
    if not (
        (curWords[chrColumn] == prevWords[chrColumn])
        and (minShift >= (int(curWords[startColumn]) - int(prevWords[startColumn])))
    ):
        flush(fullBuff)
        readBuff = []
        fullBuff = []

    # Ignore reads starting at the exact same position: likely PCR dups
    if len(readBuff) == 0 or int(curWords[startColumn]) != readBuff[-1][1]:
        readBuff.append([curWords[chrColumn], int(curWords[startColumn])])
        fullBuff.append(line)

    prevWords = curWords

# Flush after we're done
flush(fullBuff)
