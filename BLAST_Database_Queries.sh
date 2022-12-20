# Stop script on any error.
set -ux

# Make a directory to store the blast database
mkdir -p db

# Download the proteins deposited at NCBI for project PRJNA257197
esearch -db protein -query PRJNA257197 | efetch -format fasta > db/proteins.fa

# Download the nucleotides deposited at NCBI for project PRJNA257197
esearch -db nuccore -query PRJNA257197 | efetch -format fasta > db/genomes.fa

# Build the blast index for proteins.
makeblastdb -in db/proteins.fa -dbtype prot -out db/ebov_prots -parse_seqids

# Build the blast index for proteins.
makeblastdb -in db/genomes.fa -dbtype nucl -out db/ebov_genomes -parse_seqids

# Get information on the database
blastdbcmd -db db/ebov_prots -info

# Make a mystery query
echo ">mystery" > foo.fa
echo "ATGGACTCTCGTCCTCAGAAAGTCTGGATGACGCCGAGTCTCACTGAATCTGACATGGAT" >> foo.fa
echo "TACCACAAGATCTTGACAGCAGGTCTGTCCGTTCAACAGGGGGTTGTTCGGCAAAGAGTC" >> foo.fa
echo "ATCCCAGTGTATCAAGTAAACAATCTTGAG" >> foo.fa

# The format string for each query.
F='6 qacc sacc pident qlen length frames'

# Perform a search in nucleotide space.
blastn -db db/ebov_genomes  -query foo.fa -outfmt "$F" | sort -k 3rn > results_blastn.txt

# Perform a search in protein space.
blastx -db db/ebov_prots  -query foo.fa -outfmt "$F" | sort -k 3rn > results_blastx.txt

# Perform a search in translated nucleotide spaces.
tblastx -db db/ebov_genomes -query foo.fa -outfmt "$F" | sort -k 3rn > results_tblastn.txt