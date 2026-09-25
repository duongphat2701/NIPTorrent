#!/usr/bin/env python

import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import pysam
import os
import configparser
import pandas as pd
import sys

def read_file(_file):
    with open(_file) as f1:
        data = [elm.strip() for elm in f1.readlines()]

    return data

def y_fraction_specific(bam_file):

    '''
    Compute y fraction on specific gene present ONLY on the chrY

    :param bam_file: path to a file in bam format with corresponding index file. e.g. /path/to/bamfile.bam
    :return: a single y fraction value in dict. e.g. {"sample_name" : my_value}
    '''
    Y = {"HSFY1": "chrY:20708557-20750849",
         "BPY2": "chrY:25119966-25151612",
         "BPY2B": "chrY:26753707-26785354",
         "BPY2C": "chrY:27177048-27208695",
         "XKRY ": "chrY:19880860-19889280",
         "PRY": "chrY:24636544-24660784",
         "PRY2": "chrY:24217903-24242154"
         }

    dict_yfraction = dict()
    count_spec = 0
    count_Y = 0

    samfile = pysam.AlignmentFile(bam_file, "rb")
    for gene,coord in Y.items():

        for read in samfile.fetch(region = coord):
            count_spec += 1

    for read in samfile.fetch('chrY'):
        count_Y += 1

    samfile.close()
    if count_Y == 0:
        count_Y = 1

    fraction_specific_to_Y = (count_spec * 100) / count_Y
    dict_yfraction[get_prefix(bam_file)] = fraction_specific_to_Y
    # list_yfraction.append(fraction_specific_to_Y)

    # return fraction_specific_to_Y
    return dict_yfraction

def get_prefix(bamfile):

    return os.path.basename(bamfile).split(".")[0]



def calculate_threshold(y_fractions):
    """
    Compute a cut off with a Gaussian Mixture Model representing a bimodal gaussian distribution.

    :param y_fractions: list of y_fractions
    :return: a threshohld for a specific method
    """
    y_fractions = np.array(y_fractions)


    gmm = GaussianMixture(n_components=2, covariance_type='full', reg_covar=1e-99, max_iter=10000, tol=1e-99)
    gmm.fit(X=y_fractions.reshape(-1, 1))
    # gmm_x = np.linspace(0, 0.02, 5000)
    gmm_x = np.linspace(min(y_fractions), max(y_fractions), 5000)
    gmm_y = np.exp(gmm.score_samples(gmm_x.reshape(-1, 1)))

    sort_idd = np.argsort(gmm_x)
    sorted_gmm_y = gmm_y[sort_idd]

    local_min_i = argrelextrema(sorted_gmm_y, np.less)

    cut_off = gmm_x[local_min_i][0]

    #plot gmm fig
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(y_fractions, bins=50, density=True)
    ax.plot(gmm_x, gmm_y, 'r-', label='Gaussian mixture fit')
    # ax.set_xlim([0.001, 0.01])
    ax.set_xlim([min(y_fractions), max(y_fractions)])
    plt.axvline(x=cut_off,linestyle='dashed',color = 'green', label = 'threshold')
    ax.legend(loc='best')
    # plt.show()

    return cut_off


def gender_prediction(my_threshold,dict_y_fractions):

    """
    Use the cut off computed with function calculate_threshold() to determine the sample gender.
    :param my_threshold:
    :param dict_y_fractions:
    :return: dict with sample name associated with gender.
    """
    dict_gender_prediction = dict()

    for sample,y_fraction in dict_y_fractions.items():
        if y_fraction >= my_threshold:
            dict_gender_prediction[sample] = "M"

        elif y_fraction < my_threshold:
            dict_gender_prediction[sample] = "F"


    return dict_gender_prediction


def res_to_df(_pred,outfile,colname = None):
    _df = pd.DataFrame.from_dict(_pred, orient='index', columns=colname)
    _df.index.name = 'sample'
    _df.to_csv(outfile, sep='\t')

def main():
    import argparse


    parser = argparse.ArgumentParser(description='Compute gender prediction.')

    parser.add_argument('file', type=str, help='List of bam files to use. Should contain absolute path of bam files.')
    parser.add_argument('-s','--setThreshold', default=False, action='store_true', help='Default = False. Compute a threshold.')
    parser.add_argument('-t','--thresholdValue', default = False, type = float, help='You can use your own threshold. When absent use the threshold from niptune_results/gender_prediction/threshold.txt.')
    args = parser.parse_args()

    ####

    gender_pred_dir = os.getcwd()
    yfrac_dir = os.getcwd()
    table_dir = os.getcwd()

    table_output = "gender_prediction.csv"
    table_output = os.path.join(table_dir,table_output)

    table_thr = "threshold.txt"
    table_thr = os.path.join(gender_pred_dir,table_thr)

    table_yfrac = "gender_prediction_Yfrac.tsv"
    table_yfrac = os.path.join(yfrac_dir,table_yfrac)




    my_list = read_file(args.file)

    dict_magicY = dict()
    for sample in my_list:
        dict_magicY.update(y_fraction_specific(sample))

    if args.setThreshold == True and args.thresholdValue == True:
        print("You can only use 1 option at the time, either -s or -t.")
        sys.exit(1)

    my_thr = None


    if args.setThreshold == True:
        my_thr = float(calculate_threshold(list(dict_magicY.values())))
        with open (table_thr,"w") as f1:
            f1.write(str(my_thr) + "\n")

    elif args.setThreshold == False:
        if args.thresholdValue is not None:
            my_thr = args.thresholdValue
        else:
            try:
                with open(table_thr,"r") as f1:
                    my_thr = float([elm.strip() for elm in f1.readlines()][0])
            except FileNotFoundError:
                print("File {} is Missing".format(table_thr))
                print("Please use -s option first")
                sys.exit(1)

    
    dict_gender_pred = gender_prediction(my_thr, dict_magicY)
    res_to_df(dict_gender_pred, table_output, ["Gender Prediction (MagicY)"])
    res_to_df(dict_magicY, table_yfrac, ["MagicY"])





if __name__ == '__main__':

    main()