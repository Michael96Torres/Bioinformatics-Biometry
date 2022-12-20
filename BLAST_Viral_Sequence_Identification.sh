# Stop script on any error.
set -uex

# The database name.
DBNAME=ref_viruses_rep_genomes

# If you set the BLASTDB variable blast will look in that directory to find databases.
# Setting the BLASTDB variable is handy as it makes for shorter commands.
# It also allows for taxonomy information to be loaded if you have the taxdb installed there.

# Set the directory for storing BLAST databases in.
# We may repeatedly use databases, hence we store them
# outside the current project directory.
export BLASTDB=~/blastdb

# Make the BLAST database directory if it does not exist..
mkdir -p $BLASTDB

# The update_blast.pl code is not following proper UNIX standards.
# It exits with a nonzero status code even when it completes correctly (Thanks NCBI!).
# Hence we have to turn off the error watch for this line of code.
set +e

# Temporarily switch to the BLASTDB directory and run the blast database updater
# to download a viral reference database.
# The command will only update the blast database if it is incomplete or has changed.
(cd $BLASTDB && update_blastdb.pl --quiet --decompress $DBNAME)

# When one or more commands are enclosed with parenthesis as above then a new bash shell
# is opened for the commands inside the parenthesis. This is why we don't have to change
# the directory back once we are done with the commands. Only the newly opened "child" bash
# shell will see the change directory command. The "parent" bash stays where it started.
# The && symbol is treated as a logical AND. Use it to list multiple commands on a single line.

# Turn the error watch back on.
set -e

# Get the "mystery sequence" from the book site.
wget -q -nc http://data.biostarhandbook.com/fasta/mystery1.fa

# We can list the short database name because BLASTDB is set properly
blastn -db ref_viruses_rep_genomes -query mystery1.fa > blast.results.1.txt

# If you did not set the BLASTB the same command would be:
# blastn -db ~/blastdb/ref_viruses_rep_genomes -query mystery1.fa > blast.results.1.txt

# We can format the output differently.
# Look at the blast -help for information on what each field means.
blastn -db ref_viruses_rep_genomes -query mystery1.fa -outfmt 7 > blast.results.2.txt

# We can also selectively output only certain columns:
# subject accession, query accession and percent identity.
# Look at blast --help for information on how to encode the format string.
blastn -db ref_viruses_rep_genomes -query mystery1.fa -outfmt "7 sacc qacc pident nident mismatch btop" > blast.results.3.txt