# Set strict error checking.
set -eux

# The accession of the data to be downloaded from SRA. Must be a valid SRA id.

SRR=SRR519926

# The number of reads to unpack.
N=1000

# ------------ No changes needed below this line --------

# Directory containing the reads.
mkdir -p reads

# Directory containing the fastqc reports.
mkdir -p reports

# First in pair.
READ1=reads/${SRR}_1.fastq

# The second in pair.
READ2=reads/${SRR}_2.fastq

# First in pair, passed QC.
GOOD1=reads/good1.fq

# Second in pair, passed QC.
GOOD2=reads/good2.fq

# First in pair, dropped by the QC.
BAD1=reads/bad1.fq

# Second in pair, dropped by the QC.
BAD2=reads/bad2.fq

# Get fastq files from SRA. Limit to 1000 reads.
fastq-dump -X $N -O reads --split-files $SRR

# Run the initial fastqc report.
fastqc $READ1 $READ2

# Trim reads by quality.
trimmomatic PE $READ1 $READ2 $GOOD1 $BAD1 $GOOD2 $BAD2 SLIDINGWINDOW:4:30

# Run fastqc again to evaluate improvement.
fastqc -o reports $GOOD1 $BAD1 $GOOD2 $BAD2