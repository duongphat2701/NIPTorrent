#!/usr/bin/env python

import pysam
import os
import pandas as pd

def read_file(_file):
    with open(_file) as f1:
        data = [elm.strip() for elm in f1.readlines()]

    return data


def write_df(_dict,columns,output):

    """
    :param _dict: A dictionary {"sample_name" : [value], ...} or {"sample_name" : value, ...}
    :param columns: A list of colname
    :param output: Output file name
    :return: A .csv file

    Writes a dict into a csv file
    """

    _df = pd.DataFrame.from_dict(_dict, orient='index', columns=columns[1:])
    _df.reset_index(level=0, inplace=True)
    _df.rename(columns={'index': 'sample'}, inplace=True)
    _df.to_csv(output, sep='\t', index=False, na_rep='NA')

def read_number(bamfile):
    """
    Calcultate the total number of reads in the bam file.
    Discard read aligned on chrM and unaligned reads (annotated *).

    :param bamfile: path to bamfile. e.g. /path/to/bamfile.bam
    :return: a dict with the total readnumber associated to the sample name.
    """
    _dict = dict()
    _sample = os.path.basename(bamfile).split(".")[0]
    f = pysam.AlignmentFile(bamfile, 'rb')
    d = f.get_index_statistics()
    tot_read = 0

    for stat in d:
        if stat.contig != 'ChrM':
            tot_read += int(stat.mapped)

    _dict[_sample] = tot_read

    return _dict

def main():
    import configparser
    import argparse
    parser = argparse.ArgumentParser(
        description='Compute readnumber.')

    parser.add_argument('file', type=str, help='List of bam files to use. Should contain absolute path of bam files.')

    args = parser.parse_args()

    _dict = dict()

    _list = read_file(args.file)

    for sample in _list:
        _dict.update(read_number(sample))

    _outdir = os.getcwd()
    my_output = os.path.join(_outdir, "read_count.tsv")
    columns = ['sample', 'read number']
    write_df(_dict, columns, my_output)


if __name__ == '__main__':
    main()
