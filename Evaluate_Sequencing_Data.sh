# Reference genome accession number.
ACC=MN996532

# The SRR number for the sequencing data.
SRR=SRR11085797

# Alignment file.
BAM=${SRR}-${ACC}.bam

# Readcount
N=10000

# The name for the read pairs.
R1=reads/${SRR}_1.fastq
R2=reads/${SRR}_2.fastq

# The reference genome stored locally.
REF=refs/${ACC}.fa

genome:
	mkdir -p refs
	bio fetch -format fasta ${ACC} > ${REF}
	samtools faidx ${REF}
	bwa index ${REF}

data:
	# Obtain the FASTQ sequences for the SRR number.
	mkdir -p reads
	fastq-dump -F -X ${N} --split-files ${SRR} -O reads

align:
	# Run the bwa aligner. 
	bwa mem -t 4 ${REF} ${R1} ${R2} | samtools sort > ${BAM}
	samtools index ${BAM}
	samtools flagstat ${BAM} > report.txt
	cat report.txt