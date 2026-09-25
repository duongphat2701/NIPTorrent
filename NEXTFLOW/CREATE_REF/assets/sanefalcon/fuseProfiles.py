#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse


def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--version", action="version", version="1.0.0")
    parser.add_argument("-u", dest="up", type=argparse.FileType("r"), required=True, help="upstream csv file")
    parser.add_argument("-d", dest="down", type=argparse.FileType("r"), required=True, help="downstream csv file")
    return parser.parse_args()


def splitcsv(csv_fh):
    csvdict = {}
    for line in csv_fh:
        line = line.rstrip("\n")
        if not line:
            continue
        elem = line.split(",")
        sid = elem[0]
        pos = elem[1:-1]  # keep original behavior (drop last empty element if present)
        csvdict[sid] = pos
    return csvdict


def main(args):
    up = splitcsv(args.up)
    down = splitcsv(args.down)

    for sample in sorted(up.keys()):
        upPos = list(up[sample])
        upPos.reverse()
        downPos = down[sample]
        profiles = [sample] + upPos[:-1] + downPos  # drop duplicate center (same as Py2)
        print(",".join(profiles))


if __name__ == "__main__":
    main(getArgs())
