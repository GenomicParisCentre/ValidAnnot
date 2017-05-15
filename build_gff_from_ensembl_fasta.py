#!/usr/bin/python


# build_gff_from_ensembl_fasta.py
# Builds a a fasta and a gff file from cdna and ncrna fasta files from Ensembl versionned files
# Sophie Lemoine
# Genomicpariscentre



# Usage : python build_gff_from_ensembl_fasta.py -o Ensembl organism name -e Ensembl version -v
# Arguments :
# -o or --organism    	 Ensembl organism name Ex: Mus_musculus
# -e or --ensemblversion Ensembl version Ex: 84
# -v or --verbose

# Exemple :
# python build_gff_from_ensembl_fasta.py -o Mus_musculus -e 84 -v





from validannot_env import gff3_path
from validannot_env import log_path
from validannot_env import cdna_fasta_path
from validannot_env import ncrna_fasta_path
import argparse
import subprocess
import fileinput
import os
import time
import sys
import logging
from Bio import SeqIO


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

def merge_fasta(fasta1,fasta2,fasta_out):
    fout = open(fasta_out, 'w')
    for line1 in fileinput.input(fasta1):
        fout.write(line1)
    fileinput.close()
    for line2 in fileinput.input(fasta2):
        fout.write(line2)
    fileinput.close()
    fout.close()


def build_gff(fasta,gff):
    fin = open(fasta, "rU")
    gffout = open(gff,"w")
    gffout.write("##gff-version\t3\n")
    gffout.write('##created by SGDB on the ' + time.strftime("%Y-%m-%d") + '\n')
    gffout.write("##source file is " + fasta +"\n")
    for record in SeqIO.parse(fin, "fasta"):
        ID = record.id
        Desc = record.description
        Seq=str(record.seq)
        fields = Desc.split(' ')
        f = fields[1].split(':')
        feature = f[0]
        ref = fields[2].split(':')
        source = "Ensembl_"+ref[1]
        chromosome = "chromosome="+ref[2]
        chr_start = "start="+ref[3]
        chr_end = "end="+ref[4]
        chr_str = "strand="+ref[5]
        p = fields[3].split(':')
        id="ID="+ID
        parent = "Parent="+p[1]
        start = '1'
        end = str(len(Seq)+1)
        strand = '1'
        score = '.'
        phase = '.'
        attributes=";".join([id,parent,chromosome,chr_start,chr_end,chr_str])
        line = "\t".join([ID,source,feature,start,end,strand,score,phase,attributes])
        gffout.write(line+"\n")
    gffout.close()
    fin.close()

#############
# Arguments and usage

parser = argparse.ArgumentParser(description='A script to build an annotation file describing rna (cdna and ncrna from Ensembl) from Ensembl fasta files.')
parser.add_argument('-o', '--organism', help='Ensembl organism name (Ex: Mus_musculus)', required=True)
parser.add_argument('-e', '--ensemblversion', help='Ensembl version Ex: 84', required=True)
parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
args = vars(parser.parse_args())


organism = args['organism']
version = args['ensemblversion']
if args['verbose']:
    logging.basicConfig(filename=log_path+'build_ensembl_cdna_ncrna_gff_from_fasta.log',level=logging.INFO,format='%(asctime)s %(message)s')


cdna_fasta_in = organism + "_ens" + version + "_cdna.fa.gz"
# Check if cdna_fasta coming in exists
isgz_cdna_fasta = make_sure_file_exists(cdna_fasta_path,cdna_fasta_in)
# gunzip or not cdna_fasta
if type(isgz_cdna_fasta) is int:
    logging.info(os.strerror(isgz_cdna_fasta))
    sys.exit(1)
else:
    if not isgz_cdna_fasta:
        cdna_fasta_in = cdna_fasta_in[:-3]
        cdna_gff_out = cdna_fasta_in[:-6]+"_sgdb.gff"
    else:
        inF = cdna_fasta_path + cdna_fasta_in
        proc = subprocess.Popen(["gunzip", cdna_fasta_path + cdna_fasta_in], stderr=subprocess.PIPE)
        for line in proc.stderr:
            logging.info(line)
            cdna_fasta_in = cdna_fasta_in[:-3]
            cdna_gff_out = cdna_fasta_in[:-6]+"_sgdb.gff"


ncrna_fasta_in = organism + "_ens" + version + "_ncrna.fa.gz"
# Check if ncrna_fasta coming in exists
isgz_ncrna_fasta = make_sure_file_exists(ncrna_fasta_path,ncrna_fasta_in)
# gunzip or not ncrna_fasta
if type(isgz_ncrna_fasta) is int:
    logging.info(os.strerror(isgz_ncrna_fasta))
    sys.exit(1)
else:
    if not isgz_ncrna_fasta:
        ncrna_fasta_in = ncrna_fasta_in[:-3]
        ncrna_gff_out = ncrna_fasta_in[:-6] + "_sgdb.gff"
    else:
        inF = ncrna_fasta_path + ncrna_fasta_in
        proc = subprocess.Popen(["gunzip", ncrna_fasta_path + ncrna_fasta_in], stderr=subprocess.PIPE)
        for line in proc.stderr:
            logging.info(line)
            ncrna_fasta_in = ncrna_fasta_in[:-3]
            ncrna_gff_out = ncrna_fasta_in[:-6] + "_sgdb.gff"

# merge cdna and ncrna fasta files in a single file saved in the same directory than cdna fasta files
cdna_ncrna_fasta = cdna_fasta_path + organism + "_ens" + version + "_cdna_ncrna.fa"
cdna_ncrna_gff = gff3_path + organism + "_ens" + version + "_cdna_ncrna.gff"
merge_fasta(cdna_fasta_path+cdna_fasta_in, ncrna_fasta_path+ncrna_fasta_in,cdna_ncrna_fasta)
logging.info("New fasta file "+cdna_ncrna_fasta+" written.")
# build gff file from merged fasta file and save it in the gff3 directory
build_gff(cdna_ncrna_fasta,cdna_ncrna_gff)
logging.info("New gff file "+cdna_ncrna_gff+" written.")

