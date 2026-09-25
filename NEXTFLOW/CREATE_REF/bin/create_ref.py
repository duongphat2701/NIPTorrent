#!/usr/bin/env python

import os
import random
import subprocess
import pandas as pd

def read_csv(_file,_sep):
    _df = pd.read_csv(_file, sep=_sep, keep_default_na=False, na_values=[""])
    return _df

def read_file(_file):
    with open(_file) as f1:
        data = [elm.strip() for elm in f1.readlines()]
    return data

def gather_files(_dir, _extension, _list, isdir=False):
    check_list = []
    if not isdir:
        for prefix in _list:
            if os.path.isfile(os.path.join(_dir,prefix) + _extension):
                check_list.append(os.path.join(_dir,prefix) + _extension)
    else:
        for prefix in _list:
            if os.path.isdir(os.path.join(_dir, prefix) + _extension):
                check_list.append(os.path.join(_dir, prefix) + _extension)
    return check_list

def write_list(_list,output):
    with open(output,'w') as f1:
        for elm in _list:
            f1.write(f"{elm}\n")

def list_sampling(_list,n):
    if len(_list) < n:
        print(f"Warning: Requested sample size {n} exceeds available samples {len(_list)}. Sampling all available samples.")
        return _list.copy()
    sampling = random.sample(_list, n)
    return sampling

def list_gender(dict_gender):
    list_girl = []
    list_boy = []
    for sample, gender in dict_gender.items():
        if gender == "F":
            list_girl.append(sample)
        elif gender == "M":
            list_boy.append(sample)
        else:
            print(f"Warning: Sample '{sample}' has undefined gender '{gender}'. Skipping.")
    print(f"Total girls: {len(list_girl)}, Total boys: {len(list_boy)}")
    return list_girl, list_boy

def csv_to_dict(_file, _sep):
    _df = read_csv(_file, _sep)
    _dict = dict(zip(list(_df.iloc[:,0]), list(_df.iloc[:,1])))
    return _dict

def create_link(input_dir, output_dir, my_ext, my_list):
    my_list2 = gather_files(input_dir, my_ext, my_list)
    for elm in my_list2:
        my_hardlink = os.path.join(output_dir, os.path.basename(elm))
        if not os.path.isfile(my_hardlink):
            os.link(elm, my_hardlink)
        else:
            print(f"File exists: {elm} -> {my_hardlink}")

def create_wisex_ref(_cpu, binSizeNpz, _input, _npz, refsize=300):
    # Prepare command line for WisecondorX newref with given parameters
    cmd_wisex = (
        f'WisecondorX newref --nipt --cpus {_cpu} --binsize {binSizeNpz} --refsize {refsize} '
        f'{_input} {_npz}'
    )
    p = subprocess.Popen(cmd_wisex, shell=True, stdout=subprocess.PIPE)
    out, err = p.communicate()

def main():
    import argparse

    gender_table = "gender_prediction.csv"

    my_outdir = os.getcwd()
    my_outdir_npz = os.path.join(my_outdir, "npz")
    my_girldir = os.path.join(my_outdir, "girldir")
    my_boydir = os.path.join(my_outdir, "boydir")
    input_dir = os.getcwd()
    wisex_ref = os.path.join(my_outdir, "wisecondorx_reference.npz")

    list_dir = [my_outdir, my_outdir_npz, my_girldir, my_boydir]
    for my_dir in list_dir:
        if not os.path.isdir(my_dir):
            print(f"creates {my_dir}")
            os.mkdir(my_dir)

    parser = argparse.ArgumentParser(description='Create a list of files to use as a reference')

    # Positional argument: input file containing list of sample prefixes to process
    parser.add_argument('file', type=str,
                        help='Path to the text file containing list of sample prefixes to process.')

    # Optional argument: number of samples to include in the reference list
    parser.add_argument('-r', '--ref-size', type=int, default=100,
                        help='Total number of reference samples to select (default: 100). '
                             'The program will split this number approximately evenly between female and male samples.')

    # Optional argument: refsize parameter passed to WisecondorX's newref command
    parser.add_argument('-R', '--wise-refsize', type=int, default=300,
                        help='The --refsize parameter used in WisecondorX newref command (default: 300). '
                             'Controls internal reference size setting for WisecondorX.')

    # Optional argument: binsize parameter for scaling samples
    parser.add_argument('-b', '--binSizeNpz', type=int, default=1000000,
                        help='Binsize used for scaling samples (must be multiples of existing binsize). Default: 1000000.')

    # Optional argument: number of CPUs to use for parallel processing
    parser.add_argument('-c', '--cpus', type=int, default=1,
                        help='Number of CPUs to utilize for processing (default: 1).')

    args = parser.parse_args()

    my_list = read_file(args.file)

    if not os.path.isfile(gender_table):
        print(f"{gender_table} is missing in directory. Please Use Preci_gender module first to compute gender prediction.")
        return

    dict_gender = csv_to_dict(gender_table, "\t")
    list_girl, list_boy = list_gender(dict_gender)

    # Get intersection of input list and gender lists to filter valid samples
    list_girl = list(set(my_list).intersection(list_girl))
    list_boy = list(set(my_list).intersection(list_boy))

    # Split reference size approximately evenly between girls and boys
    if args.ref_size % 2 == 0:
        k = args.ref_size // 2
    else:
        k = (args.ref_size - 1) // 2

    sampling_girl = list_sampling(list_girl, k)
    sampling_boy = list_sampling(list_boy, k)
    sampling = sampling_girl + sampling_boy

    # Write the selected samples to files for later use
    write_list(sampling, os.path.join(my_outdir_npz, "reference_files.txt"))
    write_list(sampling_girl, os.path.join(my_girldir, "reference_files.txt"))
    write_list(sampling_boy, os.path.join(my_boydir, "reference_files.txt"))

    # Create hardlinks for reference sample files with appropriate extensions
    create_link(input_dir, my_outdir_npz, ".npz", sampling)
    create_link(input_dir, my_girldir, ".gcc", sampling_girl)
    create_link(input_dir, my_girldir, ".pickle", sampling_girl)
    create_link(input_dir, my_boydir, ".gcc", sampling_boy)
    create_link(input_dir, my_boydir, ".pickle", sampling_boy)

    # Gather .npz files for WisecondorX reference creation and run the command
    list_npz = gather_files(my_outdir_npz, ".npz", sampling)
    create_wisex_ref(args.cpus, int(args.binSizeNpz), " ".join(list_npz), wisex_ref, args.wise_refsize)


if __name__ == '__main__':
    main()