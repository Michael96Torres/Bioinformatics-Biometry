# Stop on errors
set -uexo pipefail

# Install the bio package
# pip install bio --upgrade

# The reference genome.
ACC=NZ_CP008918

# Select sequencing run.
SRR=SRR4124989

# The directory that stores the reference sequence.
mkdir -p ref

# The file that stores the reference genome as FASTA.
REF=ref/${ACC}.fa

# The file that stores the reference genome as GFF.
GFF=ref/${ACC}.gff

# The resulting alignment file.
BAM=$ACC.bam

# Download the reference accession number.
bio fetch -format fasta  $ACC > $REF

# Download the accessions as intervals.
bio fetch -format gff  $ACC > $GFF

# The directory that stores the sequencing reads.
mkdir -p reads

# Download the sequencing data.
fastq-dump --origfmt  --split-files -O reads $SRR

# Generate a report on the data size
seqkit stat reads/${SRR}_?.fastq

# The one-liner code comes from:
#
# https://thegenomefactory.blogspot.com/2018/10/a-unix-one-liner-to-call-bacterial.html

# The number of CPUs to use.
CPUS=4

# The names for the read pairs.
R1=reads/${SRR}_1.fastq
R2=reads/${SRR}_2.fastq

# Unix trick: A single command line that is too long may be broken
# into several lines as long as each ends with a \ symbol.
#
# Generate the alignments.
minimap2 -a -x sr -t $CPUS $REF $R1 $R2 \
 | samtools sort -l 0 --threads $CPUS > $BAM
 
# Call the variants.
cat $BAM | bcftools mpileup -Ou -B --min-MQ 60 -f $REF - \
 | bcftools call -Ou -v -m - \
 | bcftools norm -Ou -f $REF -d all - \
 | bcftools filter -Ov -e 'QUAL<40 || DP<10 || GT!="1/1"' \
 > variants.vcf

# Generate statistics on variants
bcftools stats variants.vcf > variants-stats.txt

# Plot the variant statistics.
plot-vcfstats -P -p plots variants-stats.txt
