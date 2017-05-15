#!/usr/bin/python

# format_ncbi_fasta_from_gffid.py
# Format fasta and annotation files retrieved from the NCBI database.
# Sophie Lemoine
# Genomicpariscentre

# Usage : python format_ncbi_fasta_from_gffid.py -i ncbi_gff3_file -f ncbi_merged_fasta_file -v

# Arguments :
# -i or --gffin        = 20160203_Citrus_sinensis_ncbi.gff.gz or .gff
# -f or --fastain      = 20160203_Citrus_sinensis_ncbi_genome.fa.gz or .fa
# -v or --verbose

# Exemple :
# python format_ncbi_fasta_from_gffid.py -i 20160203_Citrus_sinensis_ncbi.gff.gz -f 20160203_Citrus_sinensis_ncbi_genome.fa.gz -v


from validannot_env import gff3_path
from validannot_env import log_path
from validannot_env import dna_fasta_path
import argparse
import os
import os.path
from os.path import basename
import gzip
import time
import logging
import fileinput
import sys
from Bio import SeqIO


#############
# Functions #
#############

# Checks if fasta and gff files exists
# Argv = fasta and gff complete paths
# Returns nothing
def make_sure_files_exist(fasta, gff):
    if not os.path.isfile(gff):
        logging.info("No such file: " + gff)
        sys.exit(1)
    elif not os.path.isfile(fasta):
        logging.info("No such file: " + fasta)
        sys.exit(1)
    else:
        logging.info("Found " + gff + " and " + fasta + " files.\nAnalyzing files...")

# Gets the chromosomes names to format from the gff file in the fasta file
# Argv = the gff file
# Returns a chromosome list
def get_chromosome_from_ncbigff(gff):
    wanted = []
    if gff.endswith(".gz"):
        gffin = gzip.open(gff, "rU")
    else:
        gffin = open(gff, "rU")
    for line in gffin:
            if line.startswith('#'):
                continue
            gff_fields = line.strip().split('\t')
            attributes = gff_fields[8].split(';')
            if gff_fields[2] == 'region':
                attributes = gff_fields[8].split(';')
                for a in attributes:
                    if a.startswith("chromosome=") and a != "chromosome=Unknown" and gff_fields[0] not in wanted:
                        wanted.append(gff_fields[0])
                    elif a.startswith("chromosome=") and a != "chromosome=Unknown" and gff_fields[0] in wanted:
                        continue
                    else:
                        continue
    logging.info("chromosome list from " + basename(gff) + " : " + ','.join(wanted))
    gffin.close()
    return wanted

# Builds a new gff file restricted to features referenced by a real chromosome ncbi id
# Argv = the gff input file and the wanted chromosome id list
# Returns the new output gff file
def build_newgff(gff,wanted):
    # Format output name
    if gff.endswith(".gz"):
        gff_out = gff3_path + "only_chr_" + gff[:-3]  # Output gff file
        gff_in = gzip.open(gff3_path + gff, "rU")
    else:
        gff_out = gff3_path + "only_chr_" + gff  # Output gff file
        gff_in = open(gff3_path + gff, "rU")

    # Sort and write new GFF
    logging.info("writing new gff3 file....")
    gffout = open(gff_out, "w")


    for line in gff_in:
        if not line.startswith('#'):
            gff_fields = line.strip().split('\t')
            if gff_fields[0] in wanted:
                gffout.write(line)
        elif line.startswith('##gff-version'):
            gffout.write(line)
            gffout.write('##modified by SGDB on the '+time.strftime("%Y-%m-%d")+'\n')
            gffout.write('##regular chromosome only gff file\n')
        elif line.startswith('#!'):
            gffout.write(line)
        elif line.startswith('##species'):
            gffout.write(line)
        elif line.startswith('##sequence-region'):
            header_fields = line.strip().split(' ')
            header_fields = filter(None, header_fields)
            if header_fields[1] in wanted:
                gffout.write(line)
        else:
            gffout.write(line)
    gffout.close()
    fileinput.close()
    return gff_out


#############

# Arguments and usage
parser = argparse.ArgumentParser(description='A script to format the fasta file so that the sequence ids (chromosomes) cope to the gff3 files references.')
parser.add_argument('-i', '--gff', help='input gff file name', required=True)
parser.add_argument('-f', '--fasta', help='input fasta file name', required=True)
parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
args = vars(parser.parse_args())

gff_in = args['gff']
fasta_in = args['fasta']
if args['verbose']:
    logging.basicConfig(filename=log_path +'format_ncbi_fasta_from_gffid.log',level=logging.INFO,format='%(asctime)s %(message)s')
# Check if input files exist
make_sure_files_exist(gff3_path+gff_in,dna_fasta_path+fasta_in)

# GFF analysis
# Parse the gff file to get the chromosome list.
# In ncbi files, chromosomes are defined like region features.
# To ensure the features are really chromosomes, ckeck the attributes field to find the chromosome id
wanted = get_chromosome_from_ncbigff(gff3_path+gff_in)

# Sort and write new GFF
gff_out = build_newgff(gff_in,wanted)
logging.info(gff_out+" : New gff3 file written.")

# Sort and write new fasta
# Open the input fasta file differently if it's a gz or not
if fasta_in.endswith(".gz"):
    fasta_out = dna_fasta_path + "only_chr_sgdb_" + fasta_in[:-3] # Output fasta file
    fastain = gzip.open(dna_fasta_path+fasta_in,"rU")
else:
    fasta_out = dna_fasta_path + "only_chr_sgdb_" + fasta_in  # Output fasta file
    fastain = open(dna_fasta_path+fasta_in, "rU")

fout = open(fasta_out, "w")
# Parse the fasta file to get the ref id and check if it's in the wanted list, following the gff
for record in SeqIO.parse(fastain, "fasta") :
    fasta_id = record.id
    fasta_id = fasta_id.split('|')
    ref = fasta_id[3]
    if ref in wanted:
        logging.info(record.id+"=>"+ref)
        record.id = ref
        SeqIO.write(record, fout, "fasta")
        wanted.remove(record.id)
if wanted:
    logging.info("Missing chromosomes in new fasta file: "+','.join(wanted))
else:
    logging.info("All chromosomes found and written in new fasta file")

fastain.close()
fout.close()
logging.info("New fasta file "+fasta_out+" written.")
