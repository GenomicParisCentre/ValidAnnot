#!/usr/bin/python


# modify_ensembl_gff.py
# Format Ensembl gff3 files to make real and clean gff3 files with chromosome only if desired

# Sophie Lemoine
# Genomicpariscentre

# Usage : python modify_ensembl_gff.py -o Ensembl organism name -e Ensembl version -c y -v

# Arguments :
# -o or --organism    	 Ensembl organism name Ex: Mus_musculus
# -e or --ensemblversion Ensembl version Ex: 84
# -c or --chronly        Chromosome only : y
# -v or --verbose

# Exemple :
# python modify_ensembl_gff.py -o Mus_musculus -e 84 -c y -v

from validannot_env import gff3_path
from validannot_env import log_path
import sys
import argparse
import subprocess
import os
import os.path
import fnmatch
import fileinput
import time
import logging


#############
# Functions #
#############

# Finds a file in a directory
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
            return result

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


# Gets the chromosomes names in the gff3 file
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
    return wanted

def convert_ensemblgff_attributes(gff_fields):

    attributes = gff_fields[8].split(';')
    for a in attributes:
        if a.startswith('ID=gene:') and gff_fields[2].endswith('_gene'):
            gff_fields[2] = 'gene'
        if a.startswith('ID=gene:') and gff_fields[2] == 'RNA':
            gff_fields[2] = 'gene'
        if a.startswith('ID=chromosome:'):
            i = attributes.index(a)
            attributes[i] = a.replace('ID=chromosome:', 'ID=')
        if a.startswith('ID=gene:'):
            i = attributes.index(a)
            attributes[i] = a.replace('ID=gene:', 'ID=')
        if a.startswith('ID=transcript:'):
            i = attributes.index(a)
            attributes[i] = a.replace('ID=transcript:', 'ID=')
        if a.startswith('ID=CDS:'):
            i = attributes.index(a)
            attributes[i] = a.replace('ID=CDS:', 'ID=')
        if a.startswith('Parent=gene:'):
            i = attributes.index(a)
            attributes[i] = a.replace('Parent=gene:', 'Parent=')
        if a.startswith('Parent=transcript'):
            i = attributes.index(a)
            attributes[i] = a.replace('Parent=transcript:', 'Parent=')

    gff_fields[8] = ";".join(attributes)
    n = "\t".join(gff_fields)
    return n


#############

# Arguments and usage
parser = argparse.ArgumentParser(description='A script to format Ensembl a raw gff3 file into a usable gff3 file')
parser.add_argument('-o', '--organism', help='Ensembl organism name (Ex: Mus_musculus)', required=True)
parser.add_argument('-e', '--ensemblversion', help='Ensembl version Ex: 84', required=True)
parser.add_argument('-c', '--chronly', help='Chromosome only', default='y')
parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
args = vars(parser.parse_args())

organism = args['organism']
version = args['ensemblversion']
chronly= args['chronly']
if chronly == 'y':
    wanted = []
if args['verbose']:
    logging.basicConfig(filename=log_path +'modify_ensembl_gff.log',level=logging.INFO,format='%(asctime)s %(message)s')

gfffile_in = organism + "_ens" + version + ".gff.gz"
logging.info(gfffile_in )


isgz = make_sure_file_exists(gff3_path,gfffile_in)
# gunzip or not
if type(isgz) is int:
    logging.info(os.strerror(isgz))
    sys.exit(1)
else:
    if not isgz:
        gfffile_in = gfffile_in[:-3]
    else :
        inF = gff3_path + gfffile_in
        proc = subprocess.Popen(["gunzip", gff3_path + gfffile_in],stderr = subprocess.PIPE)
        for line in proc.stderr:
            logging.info(line)
        gfffile_in = gfffile_in[:-3]


# format gff to build a new gff and retrieve chr if necessary
gfffile_out = gfffile_in[:-4]
if chronly == 'y':
    gfffile_out = "only_chr_"+gfffile_out + '_sgdb.gff'
    # GFF analysis
    wanted = get_chromosome_from_ensemblgff(gff3_path + gfffile_in)
    logging.info("Analysing raw Ensembl gff3 file and writing new only chromosome Ensembl sgdb gff3 file....")
else:
    gfffile_out = gfffile_out + '_sgdb.gff'
    logging.info("Analysing raw Ensembl gff3 file and writing new Ensembl sgdb gff3 file....")

gffout = open(gff3_path + gfffile_out, 'w')


for line in fileinput.input(gff3_path + gfffile_in):
    if not line.startswith('#'):
        gff_fields = line.strip().split('\t')
        if chronly == 'y':
            if gff_fields[0] in wanted:
                newline = convert_ensemblgff_attributes(gff_fields)
                gffout.write(newline + "\n")
            elif gff_fields[0] not in wanted:
                pass
        else:
            newline = convert_ensemblgff_attributes(gff_fields)
            gffout.write(newline + "\n")
    elif line.startswith('##gff-version'):
        gffout.write(line)
        gffout.write('##modified by SGDB on the ' + time.strftime("%Y-%m-%d") + '\n')
        if chronly == 'y':
            gffout.write('##regular chromosome only gff file\n')
    elif line.startswith('##sequence-region'):
        if chronly == 'y':
                header_fields = line.strip().split(' ')
                header_fields = filter(None, header_fields)
                if header_fields[1] in wanted:
                  gffout.write(line)
        else:
            gffout.write(line)


fileinput.close()
gffout.close()
if chronly == 'y':
    logging.info("Only chromosome Ensembl_sgdb gff3 file written.")
else:
    logging.info("Ensembl_sgdb gff3 file written.")