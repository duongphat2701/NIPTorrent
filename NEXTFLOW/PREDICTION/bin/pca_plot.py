#!/usr/bin/env python

import argparse
import os
import pandas as pd
import subprocess
import pickle

current_dir = os.path.dirname(os.path.abspath(__file__))
assets_dir = os.path.join(current_dir, "../assets")
r_script = os.path.join(assets_dir, "pca.R")
rmd_script = os.path.join(assets_dir,"generate_table_html.Rmd")

def read_csv(_file,_sep):
    _df = pd.read_csv(_file, sep = _sep,keep_default_na = False, na_values = [""])
    return _df


def read_pickle(_file):
    with open(_file, 'rb') as f1:
        _pickle = pickle.load(f1)

    return _pickle

def read_file(_file):
    with open(_file) as f1:
        data = [elm.strip() for elm in f1.readlines()]

    return data

def gather_files(_dir, _extension, _list, isdir = False):
    """

    :param _dir: Directory where file of interest are stored.
    :param _extension: Extension of file of interest.
    :param _list: List of files of interest. Look for those files in _dir.
    :param pattern: In some case we need to look for a pattern and not for an extension.
    :return: A list of files.

    Regroup files of interest from a directory and all the subdirectories.
    """
    check_list = list()
    """list_files = [os.path.join(root, file) for root, dirs, files in os.walk(_dir) for file in files if
                   file.endswith(_extension)]
    for file in list_files:
         if pattern is False:
             prefix = os.path.basename(file).split('.')[0]
         else :
             prefix = os.path.basename(file).split(_extension)[0]
         if prefix in _list:
             check_list.append(file)"""
    # if pattern is False:
    #     files_found = [os.path.basename(file).split('.')[0] for file in list_files]
    # else:
    #     files_found = [os.path.basename(file).split(_extension)[0] for file in list_files]
    #
    # # check_list = []

    if isdir is False:
        for prefix in _list:
            if os.path.isfile(os.path.join(_dir,prefix) + _extension):
                check_list.append(os.path.join(_dir,prefix) + _extension)
    else:
        for prefix in _list:
            if os.path.isdir(os.path.join(_dir, prefix) + _extension):
                check_list.append(os.path.join(_dir, prefix) + _extension)
    # return list_files
    return check_list

def load_pickle(my_file):
    """
    Load a gcc file in pickle format.
    The file input file come from the binning step.
    The output consist in a dict of list of bin.

    :param my_file: a file in pickle format
    :return: a dict of list
    """
    my_bin = []
    dict_bin = dict()
    filename = os.path.basename(my_file).split('.gcc')[0]
    data = read_pickle(my_file)
    for k, v in data.items():
        for elm in v:
            my_bin.append(str(elm))

    if filename not in dict_bin:
        dict_bin[filename] = my_bin

    else:
        # logger.error("Redondant file found during the PCA preparation")
        exit("Redondant file found during the PCA preparation")

    return dict_bin


def prepare_dict_pca(bin_dir, my_list):

    """
    Use load_pickle() to create a dict containing a list of bin for each sample.

    :param bin_dir: input dir where are the file we want to use i.e. : ./tests/data/
    :param my_list: list of prefix files to use in the analyze i.e. : ./tests/data/list_prefix.txt
    :return: a dict of all the bin list
    """

    dict_pca = dict()

    liste_gcc = gather_files(bin_dir, '.gcc', my_list)

    for file_gcc in liste_gcc:
        tempmy_dict = load_pickle(file_gcc)

        for k, v in tempmy_dict.items():
            if k not in dict_pca:
                dict_pca[k] = v

    return dict_pca


def res_to_df(dict_pca, df_path):

    """
    Write a df and a PCA plot in outdir.

    :param dict_pca: dict from preparemy_dict_pca()
    :param outdir: output dir where results will be written
    :return: a pandas df
    """
    for key, value in dict_pca.items():
        print(f"Key: {key}, Length: {len(value)}")
    lengths = {key: len(value) for key, value in dict_pca.items()}
    min_length = min(lengths.values())

    print(f"Min length found: {min_length}")

    filtered_dict_pca = {key: value[:min_length] for key, value in dict_pca.items()}

    my_df = pd.DataFrame(filtered_dict_pca).transpose()
    my_df.index.name = 'sample'

    my_df.to_csv(df_path, sep=';')

    return my_df

def plot_PCA(df_path, outdir, sample = None):

    """

    Call an R script to generate a PCA plot

    :param df_path: path of the df to use as input
    :param outdir: output dire where are stored the pca fig
    :param sample: name of sample to color in red
    :return:  plot in .pdf format
    """

    # fig_output = os.path.join(outdir,"".join(["PCA_",sample,".pdf"]))
    fig_output = os.path.abspath(os.path.join(outdir,"".join(["PCA_",sample,".svg"])))


    # TODO improve the fig output (axis, title, size, etc.)
    cmd_r = 'Rscript {} {} {} {}'.format(r_script, df_path, fig_output, sample)
    p = subprocess.Popen(cmd_r, shell=True, stdout=subprocess.PIPE)
    out, err = p.communicate()


def generate_html(input, output_html, tmp_dir):

    """
    Call a rmarkdown to generate an html table with 2 column:
        - col 1 = sample name
        - col 2 = PCA plot

    :param input: a table in csv format containing the data
    :param output_html: a table in html format
    :return:
    """
    cmd = """ Rscript -e \"rmarkdown::render(\'{rmd_script}\',knit_root_dir = {tmp_dir},intermediates_dir = {tmp_dir}, params=list(data = {input_csv}, method = {method}), output_file = {output_html})\" """.format(tmp_dir = tmp_dir, rmd_script = rmd_script, input_csv = input,method = "\'proteus\'",  output_html = output_html)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    print(err.decode('utf-8'))

def main():
    #argparse
    parser = argparse.ArgumentParser(
        description='Prepare the .gcc or .pickle files for the PCA by printing all the bins, the output correspond to the input file name, the extention ".pickle" or ".gcc" is replaced by ".bin"')
    parser.add_argument('file', type=str, help='List of prefix files to use.')
    parser.add_argument('-r', '--routine', default=False, action='store_true', help='Default = False. Set this option to compute PCA in routine.')

    args = parser.parse_args()



    ##define folders where data are stored and output
    table_name = "PCA_bins.csv"
    table_path = "pca.tsv"
    html_path = "pca.html"
    reference_bins = "bin_plot"

    bin_dir = os.getcwd()

    my_outdir = os.getcwd()
    outdir_pca_plots = "pca_plots"
    outdir_pca_bins = "pca_bins"
    my_home= os.getcwd()
    
    os.makedirs(outdir_pca_plots, exist_ok=True)
    os.makedirs(outdir_pca_bins, exist_ok=True)

    my_list = read_file(args.file)

    my_dict = prepare_dict_pca(bin_dir,my_list)

    df_path = os.path.join(outdir_pca_bins, table_name)

    if args.routine is True:
        tmp_dict = dict()
        for k, v in my_dict.items():
            tmp_dict = dict()
            tmp_dict[k] = v

            my_df = res_to_df(tmp_dict, df_path)
            tmp_df = df_path.replace(".csv", ".tmp")
            df_ref = read_csv(os.path.join(reference_bins,"PCA_bins.csv"), ";")

            ## read df from a file otherwise there is data type issues
            my_df = read_csv(df_path, ";")
            my_df = pd.concat([df_ref, my_df])

            my_df.to_csv(tmp_df, index=False, sep=';')

            plot_PCA(tmp_df, outdir_pca_plots, k)

            os.remove(tmp_df)
            print("Delelete temporary file {}".format(tmp_df))

    my_df = res_to_df(my_dict, df_path)
    ## plot PCA and write abs path to a df
    dict_path = dict()

    for sample in my_list:
        # fig_output = os.path.join(outdir_pca_plots, "".join(["PCA_", sample, ".pdf"]))
        fig_output = os.path.join(outdir_pca_plots, "".join(["PCA_", sample, ".svg"]))
        abs_fig_output = os.path.abspath(fig_output)
        if args.routine is False:
            plot_PCA(df_path, outdir_pca_plots, sample)


        if sample not in dict_path:
            dict_path[sample] = abs_fig_output
        else:
            print("Sample duplicated. Already present in the table.")

    output_html = os.path.abspath(os.path.join(my_outdir, html_path))

    df_abs_fig_path = pd.DataFrame.from_dict(dict_path, orient='index',columns=["path"])
    df_abs_fig_path.reset_index(level=0, inplace=True)
    df_abs_fig_path.rename(columns={'index': 'sample'}, inplace=True)

    output = os.path.abspath(os.path.join(my_outdir, table_path))
    df_abs_fig_path.to_csv(output, sep='\t', index=False, na_rep='NA')

    generate_html("'" + output + "'", "'" + output_html + "'", "'" + my_home + "'")




if __name__ == '__main__':

    main()


