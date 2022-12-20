# Stop on errors.
set -uex

# Reference genome accession number.
ACC=AF086833

# The reference genome stored locally.
REF=refs/$ACC.fa

# The "real" genome that we will simulate reads from.
GENOME=genome.fa

# How many reads to simulate.
N=1000

# Output bam file.
BAM=results.bam

# The directory that store the reference.
mkdir -p refs

# Delete the log file if it exists.
rm -f log.txt

# If the reference file does not exist.
if [ ! -f $REF ]; then
	# Get the reference genome in FASTA format.
    efetch -db nuccore -format fasta -id $ACC > $REF
	
	# Build the bwa index for the reference genome.
	bwa index $REF  2>> log.txt
	
	# Build IGV index for the reference genome.
	samtools faidx $REF
fi

# Copy the reference to genome only once first time around.

# If the genome does not exist, create it.
if [ ! -f $GENOME ]; then
	cp $REF $GENOME
fi

# Edit the genome between runs.
# Introduce changes into it then rerun this recipe.
# Visualize the BAM file.

# The read pair names.
R1=read1.fq
R2=read2.fq

# Simulate reads from the genome.
# No sequencing errors. Don't mutate the genome.
wgsim -N $N -e 0 -r 0 -R 0 $GENOME $R1 $R2

# Run the bwa aligner to create a BAM file.
bwa mem $REF $R1 $R2| samtools sort > $BAM

# Index the BAM file.
samtools index $BAM