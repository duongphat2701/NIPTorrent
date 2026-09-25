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

# logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
#                     filename='NiPTUNE.log', filemode='w', level=logging.DEBUG)
#logger = logging.getLogger('NiPTUNE')


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
    dict_samples=dict()
    dict_samples_pickle = dict()
    list_files = glob.glob(my_dir + '/*.gcc')

    for my_file in list_files:
        # curFile = pickle.load(open(my_file,'rb'))
        curFile = read_pickle(my_file)
        dict_samples[my_file]=curFile
        if len(curFile[testChrom]) > 1:
            my_pickle = my_file.replace('.gcc','.pickle')
            # curFile = pickle.load(open(m:y_pickle, 'rb'))
            curFile = read_pickle(my_pickle)
            dict_samples_pickle[my_pickle] = curFile
    return dict_samples,dict_samples_pickle



def bin_to_keep(boy_dir,girl_dir, gc_corr=True):

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
    print (minLen)
      ##repose sur les echantillons femelle
    # minLen correspond à l'échantillon fille avec le chrY le plus court ?

    for i in range(minLen):
        girlData.append([girlSamples[girlSample][testChrom][i] for girlSample in girlSamples])    
    print(girlData)
    #logger.debug('girlData : {}'.format(girlData[4]))
    #logger.debug('len girlData = {}'.format([len(elm) for elm in girlData]))

    for i in range(minLen):
        boyData.append([boySamples[boySample][testChrom][i] for boySample in boySamples])

    #logger.debug('boyData : {}'.format(boyData[4]))
    #logger.debug('len boyData = {}'.format([len(elm) for elm in boyData]))

    for pos, values in enumerate(girlData):
        # keepers
        #logger.debug('median = {}, sum = {}'.format(median(values),sum(boyData[pos])))
        if median(values) != 0 or sum(boyData[pos]) == 0:
            # sum or median, delete sample for which chrY is either different from 0 for girl or equal to 0 for boys
            removables.append(pos)
            #logger.debug('removed append : {}'.format(pos))

        else:
            #logger.debug('keepers append : {}'.format(pos))
            keepers.append(pos)

    for i in reversed(removables):
        boyData.pop(i)
        girlData.pop(i)

    #logger.debug('list of keepers : {}'.format(keepers))
    print (keepers)
    return keepers
    

def compute_defrag(fname,keepers,corrMalesMedian, percYMales, girl_dir,boy_dir,gc_corr = True, output=False):
    defrag_dict=dict()

    girlSamplesPickle = list_ref(girl_dir)[1]
    percYGirls = getYPercMean(girlSamplesPickle)


    gcc_file = read_pickle(fname)
    pickle_file = read_pickle(fname.replace('.gcc','.pickle'))

    my_samplename = os.path.basename(fname)

    if my_samplename.endswith('.gcc') :
        my_samplename = os.path.basename(fname).split('.gcc')[0]
    elif my_samplename.endswith('.pickle') :
        my_samplename = os.path.basename(fname).split('.pickle')[0]


    result = 0

    if gc_corr is True:
        result = [gcc_file[testChrom][pos] for pos in keepers]

    elif gc_corr is False:
        result = [pickle_file[testChrom][pos] for pos in keepers]


    daphGender = solveFetalFraction(percYMales, percYGirls, getYPerc(pickle_file))

    #logger.debug('results : {}'.format(result))

    res_defrag_a = (median(result)/corrMalesMedian)*100
    res_defrag_b = daphGender*100



    ## Build trainingset for gender determination
    girlSamplesPickle = list_ref(girl_dir)[1]
    boySamplesPickle = list_ref(boy_dir)[1]

    from sklearn.neighbors import KNeighborsClassifier

    training = getYPercGrand(girlSamplesPickle)[:]
    training.extend(getYPercGrand(boySamplesPickle)[:])
    training = np.array([[x] for x in training])

    targets = [0] * len(getYPercGrand(girlSamplesPickle))
    targets.extend([1] * len(getYPercGrand(boySamplesPickle)))
    targets = np.array(targets)

    gnb = KNeighborsClassifier(5)
    gnb.fit(training, targets)
    # y_pred = gnb.predict(training)
    # print ("Testing classifier on trainingset.\nNumber of mislabeled points : %d" % (targets != y_pred).sum())

    # updated this line is not working like this...
    # prediction = gnb.predict(getYPerc(pickle_file))
    prediction = gnb.predict(getYPerc(pickle_file).reshape(1, -1))


    if median(result) / corrMalesMedian == 0.0 and getGender(prediction) == 'Male':
        cluster = "BAD"
    elif median(result) / corrMalesMedian == 0.0:
        cluster = "Girls"
    else:
        cluster = "Boys"

    # defrag_dict[my_samplename]=[res_defrag_a,res_defrag_b,cluster]
    defrag_dict[my_samplename]=[res_defrag_a]

    return defrag_dict



def compute_corrMalesMedian(male_dir,keepers):

    maleSamples,maleSamplesPickle = list_ref(male_dir)
    percYMales = getYPercMean(maleSamplesPickle)
    maleCorMedian = []
    for male in maleSamples:
        corrMales = [maleSamples[male][testChrom][pos] for pos in keepers]
        maleCorMedian.append(median(corrMales))
    corrMalesMedian = mean(maleCorMedian) * 10

    return percYMales,corrMalesMedian

def getYPercGrand(sampleList):
	values=[]
	for sample in sampleList:
		values.append(getYPerc(sampleList[sample]))

	return values

def getGender(prediction):

	if prediction == [1]:
		return "Male"
	elif prediction == [0]:
		return "Female"
	else:
		return None

def main():
    import configparser
    parser = argparse.ArgumentParser(description='Execute fetal fraction calcultation with Defrag a.')
    parser.add_argument('file', type=str, help='List of prefix files to use.')
    parser.add_argument('-g','--GCcorrection', type=bool, default=True, help='Apply GC correction. True/False. Default is True')



    args = parser.parse_args()

    table_dir = os.getcwd()
    # ref_dir = os.path.expanduser(os.path.join(my_home,ref_dir))
    ref_dir = os.getcwd()
    converted_dir = os.getcwd()

    girldir = os.path.expanduser(os.path.join(ref_dir, "girldir"))
    boydir = os.path.expanduser(os.path.join(ref_dir, "boydir"))
    maledir=""
    """
    If maledir isn't set in the conf file then we don't look for it in the reference folder and use the corrMaleMedian et percYMales parameters from the conf file.
    """
    # if maledir != "":
    #     maledir = os.path.expanduser(os.path.join(ref_dir, maledir))

    my_output = os.path.join(table_dir, "defrag_a_summary.tsv")


    corrMalesMedian = 0.412516803449
    percYMales = 0.00278246251169

    keepers = bin_to_keep(boydir, girldir)


    if maledir != "":
        percYMales, corrMalesMedian = compute_corrMalesMedian(maledir, keepers)

    _dict = dict()

    _list = read_file(args.file)
    _list = gather_files(converted_dir, ".gcc", _list)

    if not _list:
        pass

    else:


        for sample in _list:
            _dict.update(compute_defrag(sample, keepers, corrMalesMedian, percYMales, girldir, boydir,args.GCcorrection, output=True))

        # columns = ['sample', 'defrag_a', 'defrag_b', 'defrag_gender']
        columns = ['sample', 'ff (Defrag a)']
        write_df(_dict,columns,my_output)


if __name__ == '__main__':

    main()







