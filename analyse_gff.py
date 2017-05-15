#!/usr/bin/python

# analyse_gff.py
# Analyses features in gff files
# 1- histogram and summary of features in only_chr gff files and cdna-ncrna gff files
# Sophie Lemoine
# Genomicpariscentre



# Usage : python analyse_gff.py -o Ensembl organism name -e Ensembl version -f feature_list -v
# Arguments :
# -o or --organism    	 Ensembl organism name Ex: Mus_musculus
# -e or --ensemblversion Ensembl version Ex: 84
#-f  or --feature        List of features, separated by commas Ex: ncrna,cdna
# -v or --verbose

# Exemple :
# python analyse_gff.py -o Mus_musculus -e 84 -f cdna,ncrna -v




import argparse
import os.path
import fileinput
from pylab import *
import subprocess
import logging
import pandas as pd
from validannot_env import gff3_path
from validannot_env import log_path

#############
# Functions #
#############


# Checks if the path and file exist
# Argv = path and file
# Returns a var that can be a boolean or an integer
# If true, path and gz file exist, if false, path and gunzipped file exist, if 2, path or file do not exist
def make_sure_file_exists(p,f):
    if os.path.exists(p):
        if os.path.isfile(p+"/"+f[:-3]):
            logging.info("Found gunzipped file: "+ f[:-3] + "\n No need to gunzip...")
            gz = False
        elif os.path.isfile(p+"/"+f):
            logging.info("Found gzipped file: "+ f +"...")
            gz = True
        else:
            gz=2
    else:
            gz=2
    return gz

# Gets the length of specific features in a gff file
# Argv = gff file and a feature list
# Returns dictionary of ID:Length
def get_length(gff_in,feature_list):
    length ={}
    for line in fileinput.input(gff_in):
        if line.startswith('#'):
            continue
        gff_fields = line.strip().split('\t')
        attribute_field = gff_fields[8].split(';')
        for i,a in enumerate(attribute_field):
            if a.startswith("ID"):
                id_feature=a
        if gff_fields[2] in feature_list:
            i = id_feature
            l= int(gff_fields[4]) -int(gff_fields[3])
            length.update({i:l})

    fileinput.close()

    return length


# Gets the total length of feature from gff file
# Argv = gff file and a feature list
# Returns dictionary of total length
def get_total_length(gff_in,feature_list):
    length = {}
    for line in fileinput.input(gff_in):
        if line.startswith('#'):
            continue
        gff_fields = line.strip().split('\t')
        attribute_field = gff_fields[8].split(';')
        for i,a in enumerate(attribute_field):
            if a.startswith("ID"):
                id_feature=a
        if gff_fields[2] in feature_list:
            i = id_feature
            l = int(gff_fields[4]) - int(gff_fields[3])
            length.update({i: l})
    fileinput.close()
    return sum(length.values())



# Prints length data.frame summary
# Argv = length data.frame and filehandle
def length_summary(df,report):
    min = df['Length'].min()
    index_min = df['Length'].idxmin()
    ID_min = df['Ensembl ID'][index_min]
    max = df['Length'].max()
    index_max = df['Length'].idxmax()
    ID_max= df['Ensembl ID'][index_max]
    std = df['Length'].std()
    mean = df['Length'].mean()
    median = df['Length'].median()
    q1 = df['Length'].quantile(0.25)
    q2 = df['Length'].quantile(0.5)
    q3 = df['Length'].quantile(0.75)
    report.write("Min length: " + str(min) + " (" + str(ID_min)+")\n")
    report.write("Max length: " + str(max) + " (" + str(ID_max)+")\n")
    report.write("Mean length: " + str(std) + "\n")
    report.write("Length std: " + str(mean) + "\n")
    report.write("Median length: " + str(median) + "\n")
    report.write("Median length, First to third quantile: " + str(q1) + ", " + str(q2) + ", " + str(q3) + "\n")

# Generates an configured histogram of a feature length
# Argv = length data.frame, bin number, list of features, Ensembl organism name and Ensembl version
def configured_histogram(df,b,features,o,v,path):
    # Counts of feature inside each bin
    bins = np.linspace(df['Length'].min(), df['Length'].max(), b)
    grouped = df.groupby(np.digitize(df['Length'], bins))
    size = grouped.size()
    mean = grouped.mean()
    okbins = size[size > 20]
    lastbin = okbins.keys()[-1]
    lastbin_mean = mean.loc[lastbin].values[0]
    f="_".join(features)
    plt.hist(df['Length'], b, facecolor='green')
    plt.xlabel('Length')
    plt.xlim(0, lastbin_mean)
    plt.ylabel('Number')
    plt.title("Length of "+f+" "+o + "_ens" + v)
    savefig(path+f +"_"+ o + "_ens" + v + "_size_histogram.png")
    plt.close()

#############
# Arguments and usage

parser = argparse.ArgumentParser(description='A script to build an annotation file describing rna (cdna and ncrna from Ensembl) from Ensembl fasta files.')
parser.add_argument('-o', '--organism', help='Ensembl organism name (Ex: Mus_musculus)', required=True)
parser.add_argument('-e', '--ensemblversion', help='Ensembl version Ex: 84', required=True)
parser.add_argument('-f', '--feature', help='List of features, separated by commas Ex: ncrna,cdna', required=True)
parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
args = vars(parser.parse_args())


organism = args['organism']
version = args['ensemblversion']
feature = args['feature'].split(',')
if args['verbose']:
    logging.basicConfig(filename=log_path+'analyse_gff.log',level=logging.INFO,format='%(asctime)s %(message)s')

## only_chr_gff
only_chr_feature = feature[:]

if 'cdna' in only_chr_feature:
    only_chr_feature.remove('cdna')
if 'ncrna' in only_chr_feature:
    only_chr_feature.remove('ncrna')

only_chr_gff = "only_chr_" + organism + "_ens" + version+"_sgdb.gff.gz"
if len(only_chr_feature) != 0:
    isgz_only_chr_gff = make_sure_file_exists(gff3_path,only_chr_gff)
    # gunzip or not only_chr_gff
    if type(isgz_only_chr_gff) is int:
     logging.info(os.strerror(isgz_only_chr_gff))
     sys.exit(1)
    else:
        if not isgz_only_chr_gff:
           only_chr_gff = only_chr_gff[:-3]
           only_chr_gff_report = only_chr_gff[:-6] + "_report.txt"
        else:
            inF = gff3_path + only_chr_gff
            proc = subprocess.Popen(["gunzip", gff3_path + only_chr_gff], stderr=subprocess.PIPE)
            for line in proc.stderr:
                logging.info(line)
                only_chr_gff = only_chr_gff[:-3]
                only_chr_gff_report = only_chr_gff[:-6] + "_report.txt"

    report = open(gff3_path+only_chr_gff_report,"w")
    length = get_length(gff3_path+only_chr_gff,only_chr_feature)
    df_only_chr = pd.DataFrame(length.items(), columns = ['Ensembl ID','Length'])
    length_summary(df_only_chr,report)
    configured_histogram(df_only_chr, 500,only_chr_feature,organism,version,gff3_path)
    report.close()
else:
    logging.info("Some features don't exist in file " +only_chr_gff+"\n: Can't test...\n")



## cdna_ncrna_gff
cdna_ncrna_feature=[]
if 'cdna' in feature:
    cdna_ncrna_feature.append('cdna')
if 'ncrna' in feature:
    cdna_ncrna_feature.append('ncrna')

cdna_ncrna_gff = organism + "_ens" + version +"_cdna_ncrna.gff.gz"
if len(cdna_ncrna_feature) != 0:
    isgz_cdna_ncrna_gff = make_sure_file_exists(gff3_path,cdna_ncrna_gff)
    # gunzip or not cdna_ncrna_gff
    if type(isgz_cdna_ncrna_gff) is int:
        logging.info(os.strerror(isgz_cdna_ncrna_gff))
        sys.exit(1)
    else:
        if not isgz_cdna_ncrna_gff:
            cdna_ncrna_gff = cdna_ncrna_gff[:-3]
            cdna_ncrna_gff_report = cdna_ncrna_gff[:-3] + "_report.txt"
        else:
            inF = gff3_path + cdna_ncrna_gff
            proc = subprocess.Popen(["gunzip", gff3_path + cdna_ncrna_gff], stderr=subprocess.PIPE)
            for line in proc.stderr:
                logging.info(line)
                cdna_ncrna_gff = cdna_ncrna_gff[:-3]
                cdna_ncrna_gff_report = cdna_ncrna_gff[:-3] + "_report.txt"

    report = open(gff3_path+cdna_ncrna_gff_report,"w")
    length = get_length(gff3_path+cdna_ncrna_gff,feature)
    df_cdna_ncrna = pd.DataFrame(length.items(), columns = ['Ensembl ID','Length'])
    length_summary(df_cdna_ncrna,report)
    total_length=get_total_length(gff3_path+cdna_ncrna_gff,feature)
    report.write("Total length of features:"+str(total_length)+"\n" )
    configured_histogram(df_cdna_ncrna, 500,feature,organism,version,gff3_path)
    report.close()
else:
    logging.info("cdna and ncrna are the only features available in file " + cdna_ncrna_gff + "\n: Can't test any other feature....\n")
