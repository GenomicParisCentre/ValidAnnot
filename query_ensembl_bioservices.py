#!/usr/bin/python

# query_ensembl_bioservices.py
# Retrieve biomart annotations from Ensembl gene database using bioservices web services (https://pythonhosted.org/bioservices/quickstart.html) 
# Sophie Lemoine 
# Genomicpariscentre

# docker run -t -i -v /.../Scripts/ValidAnnot/Scripts:/scripts -v /..../validation_genomeannot/biomart/:/biomart -v /..../validation_genomeannot/log/:/log --rm genomicpariscentre/bioservices bash
# Usage : python /scripts/query_ensembl_bioservices.py -o organism_ensembl_name -e ensemblversion -f filenb

# Arguments :
# -o or --organism         = Bos_taurus, Mus_musculus, Homo_sapiens
# -e or --ensemblversion   = 83
# -f or --filenb           = 2 for genes and transcripts separated, 1 means together
# -v or --verbose

# Exemple :
# python query_ensembl_bioservices.py -o Bos_taurus -e 83 -f 1 -v

# Results :
# 2 files if file_number_for_gene_and_transcript = 2   btaurus_ens83_transcriptid.tsv and btaurus_ens83_geneid.tsv
# 1 file  if file_number_for_gene_and_transcript = 1   btaurus_ens83.tsv

# To be done :
# Issue about scerevisiae gene ids are mostly the same as transcrit ids (nearly no intron), building 2 files is probably better 


import os
import argparse
import httplib2
from bioservices import *
import logging


#############
# Functions #
#############

# Get the redirected url of an Ensembl version (The real url is referenced by date and not by version number)
def verified_url(version):
    initial_url = "e" + version + ".ensembl.org"
    h = httplib2.Http()
    h.follow_redirects = False
    (response, body) = h.request("http://" + initial_url)
    if response.status == 200:
        url = initial_url
    elif response.status == 301:
        url = response['location'][7:]
        u=url.split("/")
        url=u[0]
        #if url.endswith('/'):
        #   url=url[:-1]
    logging.info(initial_url + "  >>>>  " + url + " \n")
    return url


# Replace the header line of the biomart query in case transcript and gene annotations are mixed in the same file, to suit both
def replace_header(src_filename, target_filename, replacement_line):
    f = open(src_filename)
    first_line, remainder = f.readline(), f.read()
    t = open(target_filename, "w")
    t.write(replacement_line)
    t.write(remainder)
    t.close()


def format_organism_name(organism):
    gender,species=organism.split('_')
    biomart_orga=str(gender[0].lower()+species)
    return biomart_orga

#############

# Arguments and usage
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--organism', help='Ensembl organism name (Ex: Mus_musculus)', required=True)
parser.add_argument('-e', '--ensemblversion', help='Ensembl version (Ex:84)', required=True)
parser.add_argument('-f', '--filenb',
                    help='1 if gene and transcript annotations in the same file, 2 if annotations in separated files',
                    default='1')
parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
args = vars(parser.parse_args())

if args['verbose']:
    logging.basicConfig(filename='/log/query_ensembl_bioservices.log',level=logging.INFO,format='%(asctime)s %(message)s')

ensembl_organism = args['organism']
version = args['ensemblversion']
filenb = args['filenb']





# Format organism name
organism=format_organism_name(ensembl_organism)

if organism == "scerevisiae" and filenb == "1":
    logging.info("scerevisiae gene ids are mostly the same as transcrit ids (nearly no intron).\nBuilding distinct files for gene annotations and transcript annotations is probably better.\n")

url = verified_url(version)
logging.info("url="+url+"\n")
print url
s = BioMart(verbose=False, host=url)
s.datasets("ENSEMBL_MART_ENSEMBL")

# Retrieve transcript linked annotations
s.new_query()
s.add_dataset_to_xml(organism + "_gene_ensembl")
s.add_attribute_to_xml("ensembl_transcript_id")
s.add_attribute_to_xml("external_gene_name")
s.add_attribute_to_xml("description")
s.add_attribute_to_xml("chromosome_name")
s.add_attribute_to_xml("transcript_start")
s.add_attribute_to_xml("transcript_end")
s.add_attribute_to_xml("strand")
if "refseq_ncrna" in s.attributes(organism + "_gene_ensembl"):
    s.add_attribute_to_xml("refseq_ncrna")
if organism == "hsapiens":
    s.add_attribute_to_xml("ucsc")
    s.add_attribute_to_xml("entrezgene")
elif organism == "mmusculus":
    s.add_attribute_to_xml("mgi_id")
    s.add_attribute_to_xml("mgi_symbol")
elif organism == "drerio":
    s.add_attribute_to_xml("zfin_id_id")
    #s.add_attribute_to_xml("zfin_id")
    #s.add_attribute_to_xml("zfin_symbol")
    s.add_attribute_to_xml("zfin_id_symbol")
elif organism == "btaurus":
    s.add_attribute_to_xml("hgnc_id")
    s.add_attribute_to_xml("hgnc_symbol")
elif organism == "scerevisiae":
    s.add_attribute_to_xml("sgd_gene")
elif organism == "rnorvegicus":
    s.add_attribute_to_xml("rgd_id")
    s.add_attribute_to_xml("rgd_symbol")
elif organism == "celegans":
    s.add_attribute_to_xml("external_transcript_name")
    s.add_attribute_to_xml("external_transcript_source_name")
    s.add_attribute_to_xml("transcript_count")
else:
    logging.info("Undefined organism, write only standard transcript annotations to tsv files\n")

s.add_attribute_to_xml("transcription_start_site")
s.add_attribute_to_xml("transcript_length")
s.add_attribute_to_xml("ensembl_gene_id")

xmlq = s.get_xml()
xmlq = xmlq.replace("header = \"0\"", "header = \"1\"")
xmlq = xmlq.replace("uniqueRows = \"0\"", "uniqueRows = \"1\"")
transcript = s.query(xmlq)

# Write annotations to file : 2 cases
if filenb == "2":
    logging.info("writing transcript related annotations in a transcript file.\n")
    tsvfile_out = "/biomart/" + organism + "_ens" + version + "_transcriptid.tsv"
    #tsvfile_out = "./" + organism + "_ens" + version + "_transcriptid.tsv"
    tsvout = open(tsvfile_out, 'w')
    tsvout.write(transcript)
    tsvout.close()
else:
    logging.info("writing transcript related annotations in an intermediate file.\n")
    tsvfile_out = "/biomart/" + organism + "_ens" + version + "_raw.tsv"
    #tsvfile_out = "./" + organism + "_ens" + version + "_raw.tsv"
    tsvout = open(tsvfile_out, 'a')
    tsvout.write(transcript)

# Retrieve gene linked annotations
s.new_query()
s.add_dataset_to_xml(organism + "_gene_ensembl")
s.add_attribute_to_xml("ensembl_gene_id")
s.add_attribute_to_xml("external_gene_name")
s.add_attribute_to_xml("description")
s.add_attribute_to_xml("chromosome_name")
s.add_attribute_to_xml("start_position")
s.add_attribute_to_xml("end_position")
s.add_attribute_to_xml("strand")
if "refseq_ncrna" in s.attributes(organism + "_gene_ensembl"):
    s.add_attribute_to_xml("refseq_ncrna")
if organism == "hsapiens":
    s.add_attribute_to_xml("ucsc")
    s.add_attribute_to_xml("entrezgene")
elif organism == "mmusculus":
    s.add_attribute_to_xml("mgi_id")
    s.add_attribute_to_xml("mgi_symbol")
elif organism == "drerio":
    s.add_attribute_to_xml("zfin_id_id")
    #s.add_attribute_to_xml("zfin_id")
    s.add_attribute_to_xml("zfin_symbol")
    s.add_attribute_to_xml("zfin_id_symbol")
elif organism == "btaurus":
    s.add_attribute_to_xml("hgnc_id")
    s.add_attribute_to_xml("hgnc_symbol")
elif organism == "scerevisiae":
    s.add_attribute_to_xml("sgd_gene")
elif organism == "rnorvegicus":
    s.add_attribute_to_xml("rgd_id")
    s.add_attribute_to_xml("rgd_symbol")
elif organism == "celegans":
    s.add_attribute_to_xml("external_gene_name")
    s.add_attribute_to_xml("external_gene_source")
    s.add_attribute_to_xml("transcript_count")
else:
    logging.info("Undefined organism, write only standard gene annotations to tsv files\n")

xmlq = s.get_xml()
xmlq = xmlq.replace("uniqueRows = \"0\"", "uniqueRows = \"1\"")

# Write annotations to file : 2 cases
if filenb == "2":
    logging.info("writing gene related annotations in a gene file.\n")
    xmlq = xmlq.replace("header = \"0\"", "header = \"1\"")
    tsvfile_out = "/biomart/" + organism + "_ens" + version + "_geneid.tsv"
    #tsvfile_out = "./" + organism + "_ens" + version + "_geneid.tsv"
    tsvout = open(tsvfile_out, 'w')

gene = s.query(xmlq)
tsvout.write(gene)
tsvout.close()

# Change file header when transcript and gene annotations are written in the same file
# The file header needs to suit both cases, default Ensembl header does not.
if filenb == "1":
    logging.info("writing gene related annotations in an intermediate file.\n")
    newtsvfile_out = "/biomart/" + organism + "_ens" + version + ".tsv"
    #newtsvfile_out = "./" + organism + "_ens" + version + ".tsv"
    #prepare new header and replace
    f = open(tsvfile_out, 'r')
    header = f.readline()
    header = header.replace("Ensembl Transcript ID", "Ensembl ID")
    header = header.replace("Transcript ID", "ID")
    header = header.replace("Transcript stable ID", "ID")
    header = header.replace("Transcript Start (bp)", "Start (bp)")
    header = header.replace("Transcript start (bp)", "Start (bp)")
    header = header.replace("Transcript End (bp)", "End (bp)")
    header = header.replace("Transcript end (bp)", "End (bp)")
    f.close()
    logging.info("merging transcript and gene annotations in a final single file.\n")
    replace_header(tsvfile_out, newtsvfile_out, header)
    #tsvout.close()
    f.close()
    os.remove(tsvfile_out)
