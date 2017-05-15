#!/usr/bin/python


# select_ensembl_fasta_from_gffid_ensembl.py
# Selects only the official chromosome sequences in a fasta file from an only_chr_Xxxx_xxxx_ensNN_sgdb.gff file
# Sophie Lemoine
# Genomicpariscentre



# Usage : python select_ensembl_fasta_from_gffid_ensembl.py -o Ensembl organism name -e Ensembl version -v
# Arguments :
# -o or --organism    	 Ensembl organism name Ex: Mus_musculus
# -e or --ensemblversion Ensembl version Ex: 84
# -v or --verbose

# Exemple :
# python select_ensembl_fasta_from_gffid_ensembl.py -o Mus_musculus -e 84 -v

from validannot_env import gff3_path
from validannot_env import log_path
from validannot_env import dna_fasta_path
import argparse
import subprocess
import fileinput
import os
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


# Gets the chromosomes names in the gff3 file to format the fasta file
# Argv = the gff file
# Returns a chromosome list
def get_chromosome_from_ensemblgff(gff_in):
    wanted = []
    for line in fileinput.input(gff_in):
            if line.startswith('#'):
                continue
            gff_fields = line.strip().split('\t')
            if gff_fields[2] == 'chromosome' and not '_' in gff_fields[0]:
                wanted.append(gff_fields[0])
    logging.info("chromosome list for "+gff_in+" : "+','.join(wanted))
    fileinput.close()
    return wanted




#############
# Arguments and usage

parser = argparse.ArgumentParser(description='A script to build a fasta file containing only the chromosomes described in a gff3 file.')
parser.add_argument('-o', '--organism', help='Ensembl organism name (Ex: Mus_musculus)', required=True)
parser.add_argument('-e', '--ensemblversion', help='Ensembl version Ex: 84', required=True)
parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
args = vars(parser.parse_args())


organism = args['organism']
version = args['ensemblversion']
if args['verbose']:
    logging.basicConfig(filename=log_path+'select_ensembl_fasta_from_gffid.log',level=logging.INFO,format='%(asctime)s %(message)s')

gff_in = "only_chr_" + organism + "_ens" + version + "_sgdb.gff.gz"
# Check if gff coming in exists
isgz_gff = make_sure_file_exists(gff3_path,gff_in)
# gunzip or not gff
if type(isgz_gff) is int:
    logging.info(os.strerror(isgz_gff))
    sys.exit(1)
else:
    if not isgz_gff:
        gff_in = gff_in[:-3]
    else:
        inF = gff3_path + gff_in
        proc = subprocess.Popen(["gunzip", gff3_path + gff_in], stderr=subprocess.PIPE)
        for line in proc.stderr:
            logging.info(line)
        gfffile_in = gff_in[:-3]


fasta_in = organism + "_ens" + version + ".fa.gz"
# Check if fasta coming in exists
isgz_fasta = make_sure_file_exists(dna_fasta_path,fasta_in)
# gunzip or not fasta
if type(isgz_fasta) is int:
    logging.info(os.strerror(isgz_fasta))
    sys.exit(1)
else:
    if not isgz_fasta:
        fasta_in = fasta_in[:-3]
    else:
        inF = dna_fasta_path + fasta_in
        proc = subprocess.Popen(["gunzip", dna_fasta_path + fasta_in], stderr=subprocess.PIPE)
        for line in proc.stderr:
            logging.info(line)
        fasta_in = fasta_in[:-3]


# Create chromosome list from gff file
wanted = get_chromosome_from_ensemblgff(gff3_path+gff_in)

logging.info("Writing new fasta file....")
# Open fasta coming in
fin = open(dna_fasta_path+fasta_in, "rU")
# Create output fasta file and handle it
fasta_out = dna_fasta_path + "only_chr_" + fasta_in
fout = open(fasta_out, "w")

for record in SeqIO.parse(fin, "fasta") :
    if record.id in wanted:
        logging.info(record.id)
        SeqIO.write(record, fout, "fasta")
        wanted.remove(record.id)
if wanted:
    logging.info("Missing chromosomes in new fasta file: "+','.join(wanted))
else:
    logging.info("All chromosomes found and written in new fasta file")

fin.close()
fout.close()
logging.info("New fasta file "+fasta_out+" written.")