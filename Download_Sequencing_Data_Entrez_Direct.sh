# Turn on strict error checking.
set -uex

# Bioproject number
PNUM=PRJNA272617

# Accession number
ACC=AF086833

# Accession number AF086833 in Genbank format.
efetch -db nuccore -format gb -id $ACC > $ACC.gb

# Accession number AF086833 in Fasta format.
efetch -db nuccore -format fasta -id $ACC > $ACC.fa

# Select a subsequence from the genome.
efetch -db nuccore -format fasta -id $ACC -seq_start=1000 -seq_stop=1030

# Get the sequencing run information deposited for a project.
esearch  -db sra -query $PNUM | efetch -format runinfo > runinfo.csv

# What is inside of the runinfo file?
cat runinfo.csv | cut -d , -f 1,2,16 | head -3

# Get all genomic sequences deposited for a project.
esearch -db nucleotide -query $PNUM | efetch -format fasta > genomes.fa

# Get all protein sequences deposited for a project.
esearch -db protein -query $PNUM | efetch -format fasta > proteins.fa