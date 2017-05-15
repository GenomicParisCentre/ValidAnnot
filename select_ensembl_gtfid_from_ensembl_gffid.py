#!/usr/bin/python

# select_ensembl_gtfid_from_ensembl_gffid.py
# Format Ensembl gtf files to cope with the formatted gff3 files
# Sophie Lemoine
# Genomicpariscentre

# Usage : python select_ensembl_gtfid_from_ensembl_gffid.py -o Ensembl organism name -e Ensembl version -v
# Arguments :
# -o or --organism    	 Ensembl organism name Ex: Mus_musculus
# -e or --ensemblversion Ensembl version Ex: 84
# -v or --verbose

# Exemple :
# python select_ensembl_gtfid_from_ensembl_gffid.py -o Mus_musculus -e 84 -v



from validannot_env import gff3_path
from validannot_env import gtf_path
from validannot_env import log_path
import argparse
import subprocess
import time
import fileinput
import os
import sys
import logging

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


# Add lines in the gtf file header
# Argv = the gtf file handle
# Returns a new header
def format_new_header(gtfin_h):
    header_fields = []
    for line in gtfin_h:
        if line.startswith('#!genome-') or line.startswith('#!genebuild-'):
            header_fields.append(line)
    header_fields.append('#!modified by SGDB on the '+time.strftime("%Y-%m-%d")+'\n')
    header_fields.append('#!regular chromosome only gtf file\n')
    return header_fields

#############
# Arguments and usage


parser = argparse.ArgumentParser(description='A script to build a gtf file containing only the chromosomes described in a gff3 file.')
parser.add_argument('-o', '--organism', help='Ensembl organism name (Ex: Mus_musculus)', required=True)
parser.add_argument('-e', '--ensemblversion', help='Ensembl version Ex: 84', required=True)
parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
args = vars(parser.parse_args())


organism = args['organism']
version = args['ensemblversion']
if args['verbose']:
    logging.basicConfig(filename=log_path +'select_ensembl_gtfid_from_gffid.log',level=logging.INFO,format='%(asctime)s %(message)s')


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

gtf_in = organism + "_ens" + version + ".gtf.gz"
# Check if gtf coming in exists
isgz_gtf = make_sure_file_exists(gtf_path, gtf_in)
# gunzip or not gtf
if type(isgz_gtf) is int:
    logging.info(os.strerror(isgz_gtf))
    sys.exit(1)
else:
    if not isgz_gtf:
        gtf_in = gtf_in[:-3]
    else:
        inF = gtf_path + gtf_in
        proc = subprocess.Popen(["gunzip", gtf_path + gtf_in], stderr=subprocess.PIPE)
        for line in proc.stderr:
            logging.info(line)
        gtf_in = gtf_in[:-3]


# GFF analysis
wanted = get_chromosome_from_ensemblgff(gff3_path+gff_in)

logging.info("writing new gtf file....")
# Open gtf coming in
gtfin = open(gtf_path+gtf_in, "rU")
# Create output gtf file and handle it
gtf_out = gtf_path + "only_chr_" + gtf_in
gtfout = open(gtf_out, "w")



# Manage header
header = format_new_header(gtfin)
for h in header:
    gtfout.write(h)
gtfin.seek(0)
# Read and select data
for line in gtfin:
    if not line.startswith('#!'):
        gtf_fields = line.strip().split('\t')
        if gtf_fields[0] in wanted:
            gtfout.write(line)

gtfout.close()
gtfin.close()
logging.info("New gtf file written.")
