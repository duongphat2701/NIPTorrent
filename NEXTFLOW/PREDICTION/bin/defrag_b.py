#!/usr/bin/env python

###############################################################################
#                                                                             #
#    DEFRAG	(DEtection of fetal FRaction And Gender)        				  #
#    Copyright(C) 2014  VU University Medical Center Amsterdam    			  #
#    Authors: 																  #
#	 Daphne van Beek, d.vanbeek@vumc.nl  									  #
#	 Roy Straver, r.straver@vumc.nl                                 		  #
#                                                                             #
#    This script is supplementary to WISECONDOR.                       		  #
#                                                                             #
#    WISECONDOR is free software: you can redistribute it and/or 	  		  #
#	 modify it under the terms of the GNU General Public License as 		  #
#	 published by the Free Software Foundation, either version 3 of the 	  #
# 	 License, or (at your option) any later version.                          #
#                                                                             #
#    WISECONDOR is distributed in the hope that it will be useful,     		  #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with WISECONDOR. If not, see <http://www.gnu.org/licenses/>.       #
#                                                                             #
###############################################################################
#
# Modified by :
# David Pratella, david.pratella@univ-cotedazur.fr

import glob
import matplotlib
import pickle
matplotlib.use('Agg')
from pylab import *
import argparse
import os
import pandas as pd

def write_df(_dict, columns, output):
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

def gather_files(_dir, _extension, _list, isdir=False):
    """
    :param _dir: Directory where file of interest are stored.
    :param _extension: Extension of file of interest.
    :param _list: List of files of interest. Look for those files in _dir.
    :param pattern: In some case we need to look for a pattern and not for an extension.
    :return: A list of files.

    Regroup files of interest from a directory and all the subdirectories.
    """
    check_list = list()
    
    if isdir is False:
        for prefix in _list:
            if os.path.isfile(os.path.join(_dir, prefix) + _extension):
                check_list.append(os.path.join(_dir, prefix) + _extension)
    else:
        for prefix in _list:
            if os.path.isdir(os.path.join(_dir, prefix) + _extension):
                check_list.append(os.path.join(_dir, prefix) + _extension)
    return check_list

def read_file(_file):
    with open(_file) as f1:
        data = [elm.strip() for elm in f1.readlines()]
    return data

def read_pickle(_file):
    with open(_file, 'rb') as f1:
        _pickle = pickle.load(f1)
    return _pickle

def getCoverage(sample):
    return sum([sum(sample[chrom]) for chrom in sample])

def getYPerc(sample):
    return sum(sample[testChrom]) / sum([sum(sample[chrom]) for chrom in sample])

def getYPercMean(sampleList):
    values = []
    for sample in sampleList:
        values.append(getYPerc(sampleList[sample]))
    return mean(values)

def getYPercGrand(sampleList):
    values = []
    for sample in sampleList:
        values.append(getYPerc(sampleList[sample]))
    return values

def solveFetalFraction(percYMales, percYFemales, percYSample):
    # Based on: %chrY sample = meanY% males * FF + meanY% women with female fetusses * (1 - FF)
    # Taken from: Chiu et al, Non-invasive prenatal assessment of trisomy 21 by multiplexed maternal plasma DNA sequencing: large scale validity study, 2011
    return (percYSample - percYFemales) / (percYMales - percYFemales)

testChrom = 'Y'

def list_ref(my_dir):
    dict_samples = dict()
    dict_samples_pickle = dict()
    list_files = glob.glob(my_dir + '/*.gcc')

    for my_file in list_files:
        curFile = read_pickle(my_file)
        dict_samples[my_file] = curFile
        if len(curFile[testChrom]) > 1:
            my_pickle = my_file.replace('.gcc', '.pickle')
            curFile = read_pickle(my_pickle)
            dict_samples_pickle[my_pickle] = curFile
    return dict_samples, dict_samples_pickle

def bin_to_keep(boy_dir, girl_dir, gc_corr=True):
    girlData = []
    boyData = []
    removables = []
    keepers = []

    boySamples = list_ref(boy_dir)[0]
    girlSamples = list_ref(girl_dir)[0]

    if gc_corr is False:
        boySamples = list_ref(boy_dir)[1]
        girlSamples = list_ref(girl_dir)[1]

    minLen = min([len(girlSamples[girlSample][testChrom]) for girlSample in girlSamples])
    print(minLen)
    
    for i in range(minLen):
        girlData.append([girlSamples[girlSample][testChrom][i] for girlSample in girlSamples])    
    print(girlData)

    for i in range(minLen):
        boyData.append([boySamples[boySample][testChrom][i] for boySample in boySamples])

    for pos, values in enumerate(girlData):
        if median(values) != 0 or sum(boyData[pos]) == 0:
            removables.append(pos)
        else:
            keepers.append(pos)

    for i in reversed(removables):
        boyData.pop(i)
        girlData.pop(i)

    print(keepers)
    return keepers

def compute_defrag(fname, percYMales, girl_dir, gender_dict):
    defrag_dict = dict()

    my_samplename = os.path.basename(fname)
    if my_samplename.endswith('.gcc'):
        my_samplename = os.path.basename(fname).split('.gcc')[0]
    elif my_samplename.endswith('.pickle'):
        my_samplename = os.path.basename(fname).split('.pickle')[0]

    # Check gender from the provided dictionary
    sample_gender = gender_dict.get(my_samplename, None)
    
    # If female, set ff to 0
    if sample_gender == 'F':
        defrag_dict[my_samplename] = [0.0]
        return defrag_dict
    
    # Otherwise proceed with normal calculation
    girlSamplesPickle = list_ref(girl_dir)[1]
    percYGirls = getYPercMean(girlSamplesPickle)

    pickle_file = read_pickle(fname.replace('.gcc', '.pickle'))
    
    # Only calculate Defrag B (percentage-based method)
    daphGender = solveFetalFraction(percYMales, percYGirls, getYPerc(pickle_file))
    res_defrag_b = daphGender * 100

    defrag_dict[my_samplename] = [res_defrag_b]
    return defrag_dict

def compute_corrMalesMedian(male_dir, keepers):
    maleSamples, maleSamplesPickle = list_ref(male_dir)
    percYMales = getYPercMean(maleSamplesPickle)
    maleCorMedian = []
    for male in maleSamples:
        corrMales = [maleSamples[male][testChrom][pos] for pos in keepers]
        maleCorMedian.append(median(corrMales))
    corrMalesMedian = mean(maleCorMedian) * 10
    return percYMales, corrMalesMedian

def read_gender_file(gender_file):
    """Read gender prediction TSV file and return dictionary of {sample: gender}"""
    gender_dict = {}
    if gender_file and os.path.isfile(gender_file):
        with open(gender_file) as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        sample = parts[0]
                        gender = parts[1][0].upper()  # Take first character and uppercase
                        gender_dict[sample] = gender
    return gender_dict

def main():
    import configparser
    parser = argparse.ArgumentParser(description='Execute fetal fraction calcultation with Defrag b.')
    parser.add_argument('file', type=str, help='pickle file')
    parser.add_argument('girldir', type=str, help="Path to female reference directory")
    parser.add_argument('output_tsv', type=str, help="Output TSV file for the sample")
    parser.add_argument('gender_file', type=str, help='Path to TSV file with gender predictions (sample\tgender)')
    args = parser.parse_args()

    # Read gender predictions if file provided
    gender_dict = read_gender_file(args.gender_file)

    # We only need percYMales for Defrag B calculation
    percYMales = 0.00278246251169

    ff_results = compute_defrag(args.file, percYMales, args.girldir, gender_dict)

    columns = ['sample', 'ff (Defrag b)']
    write_df(ff_results, columns, args.output_tsv)

if __name__ == '__main__':
    main()