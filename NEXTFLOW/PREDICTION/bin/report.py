#!/usr/bin/env python

import pandas as pd
import configparser
import os
from functools import reduce
import subprocess

def read_csv(_file,_sep):
    _df = pd.read_csv(_file, sep = _sep,keep_default_na = False, na_values = [""])
    return _df

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

current_dir = os.path.dirname(os.path.abspath(__file__))
assets_dir = os.path.join(current_dir, "../assets")
rmd_script = os.path.join(assets_dir,"generate_table_html.Rmd")

def summary_table(list_df, output):

    """
    Take a list of df and merge them into one table.

    :param list_df: a list of df. e.g. [df1, df2]
    :param output: output file name
    :return:
    """


    list_df2 = list()
    for table in list_df:
        table_name = os.path.basename(table).split('.')[0]
        if table_name == "abnormal_results":
            df_sao = read_csv(table, "\t")
            df_sao = df_sao[["sample", "chr13", "chr18", "chr21"]]
            df_sao.columns = ['sample', 'Z-score chr13', 'Z-score chr18', 'Z-score chr21']
            list_df2.append(df_sao)

        else:
            list_df2.append(read_csv(table, "\t"))


    df_final = reduce(lambda left, right: pd.merge(left, right, on='sample'), list_df2)
    df_final.to_csv(output, sep='\t', index=False, na_rep='NA')
    print("Final table in {}".format(output))

def generate_html(input, output_html, tmp_dir):

    """
     Call a rmarkdown to generate an html table with x column:

     :param input: a table in csv format containing the data
     :param output_html: a table in html format
     :return:
     """

    # cmd = """ Rscript -e \"rmarkdown::render(\'{rmd_script}\',params=list(data = {input_csv}, method = {method}), output_file = {output_html})\" """.format(rmd_script = rmd_script, input_csv = input, method = "\'thalassa\'", output_html = output_html)
    cmd = """ Rscript -e \"rmarkdown::render(\'{rmd_script}\',knit_root_dir = {tmp_dir},intermediates_dir = {tmp_dir}, params=list(data = {input_csv}, method = {method}), output_file = {output_html})\" """.format(tmp_dir = tmp_dir, rmd_script = rmd_script, input_csv = input, method = "\'thalassa\'", output_html = output_html)

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    print(err.decode('utf-8'))

def main():
    table_dir = os.getcwd()
    my_home = os.getcwd()
    my_output_file = "report.tsv"
    my_output_file = os.path.join(table_dir,my_output_file)

    output_html = "report.html"
    output_html = os.path.join(table_dir, output_html)



    my_list = ['read_count','preci_pca','preci_gender','ff_male', 'ff_all', 'abnormal_results', 'path_plots_wisecondorx']
    my_list = gather_files(table_dir, ".tsv", my_list)

    if my_list == []:
        print("No table found in {}".format(table_dir))
    else:
        print("Report found the following tables: ")
        for elm in my_list:
            print(elm)

        summary_table(my_list, my_output_file)
        generate_html("'" + my_output_file + "'", "'" + output_html + "'", "'" + my_home + "'")


if __name__ == '__main__':
    main()






