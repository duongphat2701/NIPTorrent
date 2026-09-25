#!/usr/bin/env python

import subprocess
import os

samtools = 'samtools'
python = 'python3'

current_dir = os.path.dirname(os.path.abspath(__file__))
assets_dir = os.path.join(current_dir, "../assets")
consam_script = os.path.join(assets_dir, "binning.py")
gcc_script = os.path.join(assets_dir, "gc_correct.py")
ref_gcc = os.path.join(assets_dir, "gccount_hg19.txt")
wx_exec = 'WisecondorX'

def execute_binning(_input,binSizePickle, binSizeNpz):
    """

    :param _input: bam file
    :param _outdir: directory path where three files will be created with file name (*.pickle, *.gcc, *.npz)
    :return: 3 files : 1 file name with .gcc extension, 1 file with .pickle extension and 1 file with .npz extension

    This function apply a lowess function in order to correct GC content compared to reference GC content

    e.g.

    input = list of bam file, should contain the absolute path:
        /home/etc/my_file.bam

    output = my_file.gcc, my_file.pickle and my_file.npz

    """

    _outdir = os.getcwd()
    _pickle = os.path.join(_outdir,os.path.basename(_input).replace('.bam','.pickle'))
    _gcc = os.path.join(_outdir,os.path.basename(_input).replace('.bam','.gcc'))
    _npz = os.path.join(_outdir,os.path.basename(_input).replace('.bam','.npz'))

    list_cmd = list()


    _cmd_pickle = 'samtools view -q 1 {} | {python} {script} -binsize {binsize} -outfile {}'.format(_input,_pickle,python = python, script = consam_script, binsize = binSizePickle)
    _cmd_gcc = "{python} {script} {} -binsize {binsize} {ref_gcc} {}".format(_pickle, _gcc,python = python, script = gcc_script, ref_gcc = ref_gcc, binsize = binSizePickle)
    _cmd_npz = '{wisex} convert --binsize {binsize} {input} {output}'.format(wisex = wx_exec, input = _input, output = _npz, binsize = binSizeNpz)


    list_cmd.append(_cmd_pickle)
    list_cmd.append(_cmd_gcc)
    list_cmd.append(_cmd_npz)

    print('Running...')
    for cmd in list_cmd:
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        out,err = p.communicate()

        if out != 0:
            print(out.decode('utf-8'))

        if err != 0:
            print(err.decode('utf-8'))

def main():
    import configparser
    import argparse


    parser = argparse.ArgumentParser(description='Create bin files .pickle, .gcc and .npz from bam files.')
    parser.add_argument('file', type=str, help='List of bam files to use. Should contain absolute path of bam files.')
    parser.add_argument('-p', '--binSizePickle', type = int, default = 1e6, help = 'Bin size (bp) used for for .pickle and .gcc file. Default is 1Mb.')
    parser.add_argument('-n', '--binSizeNpz', type = int, default = 1e6, help = 'Bin size (bp) used for for .pickle and .gcc file. Default is 1Mb.')

    args = parser.parse_args()

    execute_binning(args.file, int(args.binSizePickle), args.binSizeNpz)


if __name__ == '__main__':
    main()