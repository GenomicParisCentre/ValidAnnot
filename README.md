# ValidAnnot

These scripts aim at facilitating all the annotation formatting steps before going into an RNASeq analysis.
It is an on going suite.

## Getting Started

To be defined before launching the scripts:
=> validannot_env.py

# validannot_env.py
Define paths to all the directories needed
Needs to be included in all the validannot scripts

These paths are required :
gff3_path = ".../gff3/"
dna_fasta_path = ".../dna_fasta/"
cdna_fasta_path = ".../cdna_fasta/"
ncrna_fasta_path = ".../ncrna_fasta/"
gtf_path = ".../gtf/"
log_path = ".../log/"

### Prerequisites

python 2.7
biopython
docker


### retrieve_files_from_ensembl.py


# retrieve_files_from_ensembl.py
Retrieve fasta and annotation files from Ensembl ftp  (ensembl.org or ensemblgenomes.org)

``
python retrieve_files_from_ensembl.py -o organism_ensembl -e ensemblversion -f file type -t generic -v`

```
# Arguments :
 -o or --organism         = Bos_taurus, Mus_musculus, Homo_sapiens
 -t or type               = plants, fungi, metazoa, bacteria, protists, generic
 -e or --ensemblversion   = 83
 -f or --files            = all, gff3, gtf, dna_fasta, cdna_fasta, ncrna_fasta (default=all)
 -v or --verbose

``
python retrieve_files_from_ensembl.py -o Bos_taurus -e 83 -f all -t generic -v`

```



### modify_ensembl_gff.py


# modify_ensembl_gff.py

Format Ensembl gff3 files to make real and clean gff3 files with chromosome only if desired

```
python modify_ensembl_gff.py -o Ensembl organism name -e Ensembl version -c y -v

```
# Arguments :
-o or --organism    	 Ensembl organism name Ex: Mus_musculus
-e or --ensemblversion Ensembl version Ex: 84
-c or --chronly        Chromosome only : y
-v or --verbose

```
python modify_ensembl_gff.py -o Mus_musculus -e 84 -c y -v

```

### select_ensembl_fasta_from_gffid_ensembl.py


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

#####################
select_ensembl_gtfid_from_ensembl_gffid.py
#####################

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

#####################
build_gff_from_ensembl_fasta.py
#####################

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

#####################
analyse_gff.py
#####################

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


#####################
query_ensembl_bioservices.py
#####################

# query_ensembl_bioservices.py
# Retrieve biomart annotations from Ensembl gene database using bioservices web services (https://pythonhosted.org/bioservices/quickstart.html)
# Sophie Lemoine
# Genomicpariscentre

# docker run -t -i -v /.../Scripts/ValidAnnot/:/test --rm genomicpariscentre/bioservices bash
# Usage : python query_ensembl_bioservices.py -o organism_ensembl_name -e ensemblversion -f filenb

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

#####################
retrieve_ncbi_fasta_from_gffid.py
#####################

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

#####################
format_ncbi_fasta_from_gffid.py
#####################

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
