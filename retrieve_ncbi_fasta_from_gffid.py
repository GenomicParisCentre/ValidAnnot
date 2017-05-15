#!/usr/bin/python


# retrieve_ncbi_fasta_from_gffid.py
# Retrieve fasta and annotation files from NCBI database using ftp
# Sophie Lemoine
# Genomicpariscentre

# Usage : python retrieve_files_from_ensembl.py -o organism_ncbi -v

# Arguments :
# -o or --organism         = Capra_hircus
# -v or --verbose

# Exemple :
# python retrieve_ncbi_fasta_from_gffid.py -o Capra_hircus -v


from validannot_env import gff3_path
from validannot_env import dna_fasta_path
from validannot_env import log_path
import argparse
import os
import datetime
import os.path
import hashlib
import gzip
import urllib2
import logging
import glob
from ftplib import FTP


#############
# Functions #
#############

# Retrieves server date link to file name in a ftp connexion
# Argv = the ftp connexion
# Returns a dictionary containing filename:date
def get_ftp_file_datetime(ftp):
    data = []
    filetodate = {}
    corresponding_month = {'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04', 'May': '05', 'Jun': '06', 'Jul': '07',
                           'Aug': '08', 'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'}
    # Get files and file parameters
    ftp.dir(data.append)
    for line in data:
        print line
        col = line.split()
        # Format date
        if ":" in col[7]:
            now = datetime.datetime.now()
            year = str(now.year)
        else:
            year = col[7]
        month =  corresponding_month[col[5]]
        if len(col[6]) == 1:
            day = "0"+col[6]
        else:
            day = col[6]
        datestr = year + month + day
        # Link file name and its date
        filetodate.update({col[8]: datestr})
    return filetodate

# Retrieves server file sizes in a ftp connexion
# Argv = the ftp connexion
# Returns a dictionary containing filename:size
def get_ftp_file_size(ftp):
    data = []
    filetosize = {}
    # Get files and file parameters
    ftp.dir(data.append)
    for line in data:
        col = line.split()
        # Format date
        size = col[4]
        # Link file name and its date
        filetosize.update({col[8]: size})
    return filetosize

# Computes the md5 key for a remote file without copying it locally
# Would not work for big files
# Argv = remote file path
# Returns the md5sum key
def remote_md5Checksum(remotefile):
    remote = urllib2.urlopen(remotefile)
    m = hashlib.md5()
    while True:
        data = remote.read(8192)
        if not data:
            break
        m.update(data)
    return m.hexdigest()

# Computes the md5 key for a local file
# Would not work for big files
# Argv = file path
# Returns the md5sum key
def md5Checksum(filePath):
    with open(filePath, 'rb') as fh:
        m = hashlib.md5()
        while True:
            data = fh.read(8192)
            if not data:
                break
            m.update(data)
        return m.hexdigest()


# Retrieves gff3 *top_level.gff3.gz file from a ftp connexion on the ncbi genome site
# Reformats the name following the server date of the file
# Checks if already exists locally with its date and its size
# Does not retrieve it if it's the same
# Argv = the root path where to save the files and the organism name following the ncbi convention ex: Mus_musculus
# Gets the distant file and returns the reformated gff3 file path and name
# ftp://ftp.ncbi.nih.gov/genomes/xxx_xxx/GFF/
def retrieve_gff3(p, s):
    ncbi = 'ftp.ncbi.nih.gov'
    # Paths and filename preparation
    ftp = FTP(ncbi)
    ftp.login()
    ftp.getwelcome()
    ftp.cwd('genomes')
    ftp.cwd(s)
    ftp.cwd('GFF')
    filetodate = get_ftp_file_datetime(ftp)
    filetosize = get_ftp_file_size(ftp)
    for f in ftp.nlst("."):
        if f.startswith("ref") and f.endswith("top_level.gff3.gz"):
            gfffile_out = p + filetodate[f] + "_" + s + "_ncbi.gff.gz"
            # Looks for local gff3.gz file
            if os.path.exists(gfffile_out): # Finds a local file
                oldfile_size=str(os.path.getsize(gfffile_out))
                newfile_size = filetosize[f]
                if newfile_size != oldfile_size:
                    logging.info("retrieve gff3 file: " + f)
                    logging.info("write to file: " + gfffile_out)
                    outf = open(gfffile_out, 'wb')
                    ftp.retrbinary('RETR ' + f, outf.write)
                    outf.close()
                else: # Server and local files are the same
                    logging.info(gfffile_out+" already exist and is the same as server file\nThere is no need to get it again")
            else: # No local file
                logging.info("retrieve gff3 file: " + f)
                logging.info("write to file: " + gfffile_out)
                outf = open(gfffile_out, 'wb')
                ftp.retrbinary('RETR ' + f, outf.write)
                outf.close()
    return gfffile_out
    ftp.quit()

# Retrieves chromosome fasta files from a ftp connexion on the ncbi genome site following a chromosome list
# Reformats the name following the server date of the files
# Checks if already exists locally with its date and its size
# Does not retrieve it if it's the same
# Gets the distant file
# Argv = the root path where to save the files, the organism name following the ncbi convention ex: Mus_musculus and the chromosome list
# Returns a prefix following this pattern = ncbiftpfiledate_latinspecies_name__ncbi
# ftp://ftp.ncbi.nih.gov/genomes/xxx_xxx/Assembled_chromosomes/seq/
def retrieve_fasta(p, s, chromo):
    # Paths and filename preparation
    ftp = FTP('ftp.ncbi.nih.gov')
    ftp.login()
    ftp.getwelcome()
    ftp.cwd('genomes')
    ftp.cwd(s)
    ftp.cwd('Assembled_chromosomes/seq')
    filetodate = get_ftp_file_datetime(ftp)
    filetosize = get_ftp_file_size(ftp)
    prefix = ""
    for f in ftp.nlst("."):
        for c in chromo:
            if 'ref' in f and f.endswith("chr" + c + ".fa.gz"):
                newfile_size = filetosize[f]
                fastafile_out = p + filetodate[f] + "_" + s + "_ncbi_" + f
                prefix = p + filetodate[f] + "_" + s + "_ncbi_"
                if os.path.exists(fastafile_out):
                    oldfile_size = str(os.path.getsize(fastafile_out))
                    if newfile_size != oldfile_size:
                        logging.info("retrieve fasta file lines: " + f)
                        ouf = open(fastafile_out, 'w')
                        ftp.retrbinary('RETR ' + f, ouf.write)
                        ouf.close()
                    else: # Server and local files are the same
                        logging.info(fastafile_out + " already exist locally and is the same as the server file: there is no need to get it again")
                else: # No local file
                    logging.info("retrieve fasta file lines: " + f)
                    ouf = open(fastafile_out, 'w')
                    ftp.retrbinary('RETR ' + f, ouf.write)
                    ouf.close()

    ftp.quit()
    return prefix

# Analyses a gff3.gz file to build a defined chromosome list found in the gff3
# Gets the chromosome number or name in the attribute column, the reference column being the ncbi id for the chromosome
# Argv = the path (not the root path) where to find the gff3 file and the organism name following the ncbi convention ex: Mus_musculus
# Returns the chromosome list
def get_reference_and_chromosome_from_ncbigff(p,s):
    wanted = []
    with gzip.open(p,'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if fields[2] == 'region':
                attributes = fields[8].split(';')
                for a in attributes:
                    if a.startswith("chromosome=") and a != "chromosome=Unknown" and a[11:] not in wanted:
                        chromosome = a[11:]
                        wanted.append(chromosome)
                    elif a.startswith("chromosome=") and a != "chromosome=Unknown" and a[11:] in wanted:
                        continue
                    elif a.startswith("chromosome=") and a == "chromosome=Unknown" and "unplaced" not in wanted:
                        wanted.append("unplaced")
                    else:
                        continue
        logging.info ("chromosome list for " + s + " : " + ','.join(wanted))
        return wanted


# Gets the fa.gz files defined by a specific pattern from a directory and concatenate them into a single fa.gz file
# Argv = pattern to find the fa.gz files, the pattern is the concatenation between the path and the prefix of the file names
# Returns the output file name
def concatenate_fasta(pattern):
    filenames = glob.glob(pattern+"*.fa.gz")
    outfileName = pattern+"genome.fa.gz"
    with open(outfileName, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                outfile.write(infile.read())
    outfile.close()
    return outfileName


#############

# Arguments and usage

parser = argparse.ArgumentParser(description='A script to retrieve and format gff3 and fasta files from the ncbi genome database')
parser.add_argument('-o', '--organism', help='NCBI organism name (Ex: Mus_musculus)', required=True)
parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
args = vars(parser.parse_args())

organism = args['organism']
if args['verbose']:
    logging.basicConfig(filename=log_path +'retrieve_ncbi.log',level=logging.INFO,format='%(asctime)s %(message)s')

# Retrieves gff3 file
logging.info ("Get " + organism + " gff3 file from the NCBI")
gfffile_in = retrieve_gff3(gff3_path, organism)
# Gets chromosome numbers from the gff3 file
logging.info ("Get chromosome fasta file needed by gff3 file "+gfffile_in)
ref_chr = get_reference_and_chromosome_from_ncbigff(gfffile_in,organism)


# Gets the fasta files following the chromosome list
logging.info ("Get chromosome fasta files")
prefix = retrieve_fasta(dna_fasta_path, organism, ref_chr)
logging.info ("Merge all the chromosome fasta files into one")
genomeFile = concatenate_fasta(prefix)
logging.info (genomeFile+" created")