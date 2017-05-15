#!/usr/bin/python

# retrieve_files_from_ensembl.py
# Retrieve fasta and annotation files from Ensembl databases using ftp
# Sophie Lemoine
# Genomicpariscentre

# Usage : python retrieve_files_from_ensembl.py -o organism_ensembl -e ensemblversion -f file type -t generic -v

# Arguments :
# -o or --organism         = Bos_taurus, Mus_musculus, Homo_sapiens
# -t or type               = plants, fungi, metazoa, bacteria, protists, generic
# -e or --ensemblversion   = 83
# -f or --files            = all, gff3, gtf, dna_fasta, cdna_fasta, ncrna_fasta (default=all)
# -v or --verbose

# Exemple :
# python retrieve_files_from_ensembl.py -o Bos_taurus -e 83 -f all -t generic -v

# Results :

from validannot_env import gff3_path
from validannot_env import log_path
from validannot_env import dna_fasta_path
from validannot_env import gtf_path
from validannot_env import cdna_fasta_path
from validannot_env import ncrna_fasta_path
import argparse
import logging
from ftplib import FTP


#############
# Functions #
#############
# Ftp connexion
def ftp_connect(url):
    ftp = FTP(url)
    ftp.login()
    logging.info(ftp.getwelcome())
    return ftp

def ftp_defdir(t):
    if t == 'generic':
        d ='pub'
    elif t == 'plants':
        d ='pub/plants/'
    elif t == 'metazoa':
        d = 'pub/metazoa/'
    elif t == 'fungi':
        d = 'pub/fungi/'
    elif t == 'protists':
        d = 'pub/protists/'
    elif t == 'bacteria':
        d= 'pub/bacteria/'
    return d

# Retrieve gff3 file from ftp://ftp.ensembl.org/pub/release-xx/gff3
def retrieve_gff3(ftp ,s, v, t):
    # Paths and filename preparation
    releasenumber = "release" + "-" + str(v)
    gfffile_out = gff3_path + s + "_ens" + str(v) + ".gff.gz"
    root=ftp.pwd()
    ftpdir = ftp_defdir(t)
    ftp.cwd(ftpdir)
    ftp.cwd(releasenumber)
    ftp.cwd('gff3/' + s.lower())
    logging.info(ftp.pwd())
    files = []
    ftp.retrlines('NLST', files.append)
    for f in files:
        if f.endswith(v + ".gff3.gz"):
            logging.info(f)
            ftp.retrbinary('RETR ' + f, open(gfffile_out, 'wb').write)
    ftp.cwd(root)


# Retrieve gtf file from ftp://ftp.ensembl.org/pub/release-xx/gtf
def retrieve_gtf(ftp, s, v, t):
    # Paths and filename preparation
    releasenumber = "release" + "-" + str(v)
    gtffile_out = gtf_path + s + "_ens" + str(v) + ".gtf.gz"
    ftpdir = ftp_defdir(t)
    ftp.cwd(ftpdir)
    ftp.cwd(releasenumber)
    ftp.cwd('gtf/' + s.lower())
    logging.info(ftp.pwd())
    files = []
    ftp.retrlines('NLST', files.append)
    for f in files:
        if f.endswith(v + ".gtf.gz"):
            logging.info(f)
            ftp.retrbinary('RETR ' + f, open(gtffile_out, 'wb').write)


# ftp://ftp.ensembl.org/pub/release-xx/fasta/
def retrieve_dna_fasta(ftp, s, v, t):
    # Paths and filename preparation
    releasenumber = "release" + "-" + str(v)
    fastafile_out = dna_fasta_path + s + "_ens" + str(v) + ".fa.gz"
    root = ftp.pwd()
    ftpdir = ftp_defdir(t)
    ftp.cwd(ftpdir)
    ftp.cwd(releasenumber)
    ftp.cwd('fasta/' + s.lower() + '/dna')
    logging.info(ftp.pwd())
    files = []
    ftp.retrlines('NLST', files.append)
    for f in files:
        if f.endswith(".dna.toplevel.fa.gz"):
            logging.info (f)
            ftp.retrbinary('RETR ' + f, open(fastafile_out, 'wb').write)
    ftp.cwd(root)


# ftp://ftp.ensembl.org/pub/release-xx/fasta/
def retrieve_cdna_fasta(ftp, s, v, t):
    # Paths and filename preparation
    releasenumber = "release" + "-" + str(v)
    fastafile_out = cdna_fasta_path + s + "_ens" + str(v) + "_cdna.fa.gz"
    root = ftp.pwd()
    ftpdir = ftp_defdir(t)
    ftp.cwd(ftpdir)
    ftp.cwd(releasenumber)
    ftp.cwd('fasta/' + s.lower() + '/cdna')
    logging.info(ftp.pwd())
    files = []
    ftp.retrlines('NLST', files.append)
    for f in files:
        if f.endswith(".cdna.all.fa.gz"):
            logging.info (f)
            ftp.retrbinary('RETR ' + f, open(fastafile_out, 'wb').write)
    ftp.cwd(root)


# ftp://ftp.ensembl.org/pub/release-xx/fasta/
def retrieve_ncrna_fasta(ftp, s, v, t):
    # Paths and filename preparation
    releasenumber = "release" + "-" + str(v)
    fastafile_out = ncrna_fasta_path + s + "_ens" + str(v) + "_ncrna.fa.gz"
    root = ftp.pwd()
    ftpdir = ftp_defdir(t)
    ftp.cwd(ftpdir)
    ftp.cwd(releasenumber)
    ftp.cwd('fasta/' + s.lower() + '/ncrna')
    logging.info(ftp.pwd())
    files = []
    ftp.retrlines('NLST', files.append)
    for f in files:
        if f.endswith(".ncrna.fa.gz"):
            logging.info (f)
            ftp.retrbinary('RETR ' + f, open(fastafile_out, 'wb').write)
    ftp.cwd(root)



#############

# Arguments and usage

parser = argparse.ArgumentParser(description='A script to retrieve gff3, gtf and fasta files from the ensembl genome database')
parser.add_argument('-o', '--organism', help='Ensembl organism name (Ex: Mus musculus)', required=True)
parser.add_argument('-e', '--ensemblversion', help='Ensembl version (Ex:84)', required=True)
parser.add_argument('-t', '--type',help='Type of organism you want (plants, fungi, metazoa, bacteria, protists, generic)',default='generic')
parser.add_argument('-f', '--files', help='gff3, gtf, dna_fasta, cdna_fasta, ncrna_fasta or all of them', default='all')
parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
args = vars(parser.parse_args())

organism = args['organism']
version = args['ensemblversion']
files = args['files']
type = args['type']
if type == 'generic':
    url = 'ftp.ensembl.org'
else:
    url = 'ftp.ensemblgenomes.org'
if args['verbose']:
    logging.basicConfig(filename=log_path +'retrieve_ensembl.log',level=logging.INFO,format='%(asctime)s %(message)s')

ftp = ftp_connect(url)

if files == 'gff3':
    retrieve_gff3(ftp, organism, version, type)
elif files == 'gtf':
    retrieve_gtf(ftp, organism, version, type)
elif files == 'dna_fasta':
    retrieve_dna_fasta(ftp, organism, version, type)
elif files == 'cdna_fasta':
    retrieve_cdna_fasta(ftp, organism, version, type)
elif files == 'ncrna_fasta':
    retrieve_ncrna_fasta(ftp, organism, version, type)
elif files == 'all':
    retrieve_dna_fasta(ftp, organism, version, type)
    retrieve_cdna_fasta(ftp, organism, version, type)
    retrieve_ncrna_fasta(ftp, organism, version, type)
    retrieve_gff3(ftp, organism, version, type)
    retrieve_gtf(ftp, organism, version, type)
else:
    logging.info(files +" unknown. Check your parameters.")

ftp.quit()
ftp.close()