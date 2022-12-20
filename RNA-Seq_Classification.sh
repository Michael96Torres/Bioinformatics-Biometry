# Stop on errors. Print the commands.
set -uex

# Download the reference genome.
wget -q -nc http://data.biostarhandbook.com/books/rnaseq/data/golden.genome.tar.gz

# Unpack the reference genome.
tar xzvf golden.genome.tar.gz

# Download the data
wget -q -nc http://data.biostarhandbook.com/books/rnaseq/data/golden.reads.tar.gz

# Unpack the data
tar zxvf golden.reads.tar.gz

# Set the name for the reference.
REF=refs/transcripts.fa

# Set the name of the index.
IDX=refs/transcript.idx

# Index the genome with kallisto.
kallisto index -i $IDX $REF

# Index the reference genome with samtools
samtools faidx $REF

# Let's make gnu parallel nagging go away.
echo 'will cite' | parallel --citation 2> /dev/null

# Create the root ids of the data layout.
parallel -j 1 echo {1}_{2} ::: BORED EXCITED ::: 1 2 3 > ids

# Make a directory for the results
mkdir -p output

# Run Kallisto to classify the reads.
cat ids | parallel kallisto quant -i $IDX -o output/{} reads/{}_R1.fq reads/{}_R2.fq

# Download the custom script to combine kallisto outputs.
curl http://data.biostarhandbook.com/books/rnaseq/code/combine.py > combine.py

# Combine the outputs created by kallisto.
cat ids | python combine.py output > counts.txt

# Download the edger R script.
curl http://data.biostarhandbook.com/books/rnaseq/code/edger.r > edger.r

# Perform the differential expression detection with edger.
cat counts.txt | Rscript edger.r 3x3 > results.csv

# Draw the heatmap from the results.
cat results.csv | Rscript heatmap.r > results.pdf
