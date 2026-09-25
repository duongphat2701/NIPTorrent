#!/usr/bin/env python

import subprocess
import os
import argparse

wx_exec = 'WisecondorX'

class Sao:

    def __init__(self, reference, resultdir):
        self.reference = reference
        self.resdir = resultdir

    def execute_wisecondorx(self, fname, blacklist=None, gender=None):
        """
        Wrapper to execute WisecondorX predict.

        :param fname: path to .npz file
        :param blacklist: optional path to blacklist bed file
        :param gender: optional gender override ('M' or 'F')
        :return: path to result prefix
        """
        sample_name = os.path.basename(fname)
        result = os.path.join(self.resdir, sample_name.split(".npz")[0])

        cmd = [wx_exec, 'predict', fname, self.reference, result, '--bed', '--plot']

        if blacklist:
            cmd.extend(['--blacklist', blacklist])

        if gender:
            cmd.extend(['--gender', gender])

        subprocess.Popen(cmd).wait()
        return result

    def summarize_wisecondorx(self, file, output_file, _sep='\t'):
        result = {}
        sample_name = os.path.basename(file).split("_statistics.txt")[0]

        with open(file) as f:
            for line in f:
                exclude_line = ('chr', 'Standard', 'Median', 'Gender', 'Number', 'Copy')
                if not line.startswith(exclude_line):
                    parts = line.strip().split('\t')
                    chr, zscore = [parts[i] for i in (0, 3)]
                    result[chr] = float(zscore)
            if 'Y' not in result:
                result['Y'] = 'NA'

        chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
        zscore_values = [str(result.get(chr, 'NA')) for chr in chromosomes]

        with open(output_file, 'w') as out:
            header = ['sample'] + [f'chr{chr}' for chr in chromosomes]
            out.write(_sep.join(header) + '\n')
            out.write(f"{sample_name}{_sep}{_sep.join(zscore_values)}\n")

def main():
    resultdir = os.getcwd()
    table_dir = os.getcwd()

    parser = argparse.ArgumentParser(description='Execute WisecondorX prediction and summarize.')
    parser.add_argument('file', type=str, help='Path to input .npz file.')
    parser.add_argument('--ref', type=str, required=True, help='Path to reference .npz file.')
    parser.add_argument('--blacklist', type=str, help='Optional blacklist .bed file.')
    parser.add_argument('--gender', type=str, choices=['M', 'F'], help='Force gender (M or F).')
    args = parser.parse_args()

    sample = args.file
    sample_name = os.path.basename(sample).split(".npz")[0]
    my_output_file = os.path.join(table_dir, f"{sample_name}_abn.tsv")

    s = Sao(args.ref, resultdir)
    print(f"Running WisecondorX for: {sample_name}")
    s.execute_wisecondorx(sample, args.blacklist, args.gender)

    file_zscore = os.path.join(resultdir, f"{sample_name}_statistics.txt")
    print(f"Summarizing results for: {sample_name}")
    s.summarize_wisecondorx(file_zscore, my_output_file)
    print(f"Done. Output saved to: {my_output_file}")

if __name__ == '__main__':
    main()