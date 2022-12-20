# Stop script on error.
set -uex

# The SRR BioProject number for the sequencing data.
PROJECT=PRJNA257197

# The number of datasets to subselect from the project.
N=5

# Get the project run information.
esearch -db sra -query $PROJECT  | efetch -format runinfo > runinfo.txt

# Select the first N elements. Keep only valid SRR numbers.
cat runinfo.txt | cut -f 1 -d , | grep SRR | head -$N > selected.txt

# Store the data in the reads folder.
mkdir -p reads

# Download the SRR data for each
cat selected.txt | parallel fastq-dump -O reads -X 1000 --split-files {}

# Create a directory for bam files
mkdir -p bam

# Generate a separate BAM file for each SAMPLE.
cat selected.txt | parallel "picard FastqToSam F1=reads/{}_1.fastq F2=reads/{}_1.fastq O=bam/{}.bam  RG=GROUP-{} LB=LIB-{} SM=SAMPLE_{} QUIET=true 2>> log.txt"

# Merge all the BAM files into one.
samtools merge -f all.bam bam/*.bam

# Investigate the readgroups in the header.
echo ""
echo "SAM file header:"
samtools view -H all.bam

echo ""
echo "Number of alignments with read group: GROUP-SRR1972919"
samtools view -c -r GROUP-SRR1972919 all.bam

# Reverting the process is to extract reads, tagged with readgroups to paired files.
samtools fastq -t -1 all1.fq -2 all2.fq all.bam

# To convert just one specific read group.
samtools view -r GROUP-SRR1972919 all.bam | samtools fastq -t -1 all_SRR1972919_1.fq -2 all_SRR1972919_2fq -