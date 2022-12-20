# Stop on errors.
set -uex

# Reference genome accession number.
ACC=AF086833

# The SRR number for the sequencing data.
SRR=SRR1972739

# How many reads to unpack
N=10000

# The reference genome stored locally.
REF=refs/$ACC.fa

# The directory that store the reference.
mkdir -p refs

# Get the reference genome in FASTA format.
efetch -db nuccore -format fasta -id $ACC > $REF

# Build the bwa index for the reference genome.
bwa index $REF  2>> bwa.log.txt

# Build IGV index for the reference genome.
samtools faidx $REF

# Obtain the FASTQ sequences for the SRR number.
fastq-dump -X $N --split-files $SRR  > fastqdump.log.txt

# The name for the read pairs.
R1=${SRR}_1.fastq
R2=${SRR}_2.fastq

# Run the bwa aligner. Creates a SAM file.
bwa mem $REF $R1 $R2 > $SRR.sam 2>> bwa.log.txt

# Convert the SAM file to BAM format.
cat $SRR.sam  | samtools sort > $SRR.bam

# Index the BAM file.
samtools index $SRR.bam

#Generate an alignment report.
samtools flagstat $SRR.bam > report.txt