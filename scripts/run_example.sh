#!/bin/sh -x

# Assume user has done this.
#
###  Build Celera Assembler in directory called wgs-assembler
#    $ generic_build_ca.sh
###  Clone Alpaca in directory called alpaca
#    $ git clone https://github.com/VicugnaPacos/ALPACA/ alpaca
###  Run the example
#    $ mkdir test
#    $ cd test
#    $ ../alpaca/scripts/run_example.sh

echo "ALPACA NEEDS CELERA ASSEMBLER"
which gatekeeper
echo "ALPACA NEEDS BOWTIE"
which bowtie2
echo "ALPACA NEEDS SAMTOOLS"
which samtools

# absolute or relative path to the git directory containing this script
ALPACA=../alpaca
# absolute or relative path to Celera Assembler binaries
export CABIN=../wgs-assembler/Linux-amd64/bin/
# make sure commands like 'gatekeeper' use the Celera Assembler version
export PATH=${CABIN}:${PATH}
# provided with alpaca is this guide (i.e. ALLPATHS scaffolds of target genome)
GUIDE=yeast.guide.fasta
# output filename
EXAMPLE=YEAST
# generate this file of simulated long reads from the guide
SIMULATED_LONG_READS=${EXAMPLE}.s.fastq

echo "REMOVE ANY PREVIOUS RUNS"
rm -v *.frg
rm -v *.spec
rm -rfv [0-9]-*
rm -v filtered_scaffolds.*.fasta
rm -v pair_coverage*.txt
${ALPACA}/scripts/cleanup.sh

cp ${ALPACA}/example_data/*.fasta .

echo
echo "BEFORE CREATE EXAMPLE"
ls -l

echo
echo "RUN CREATE EXAMPLE..."
${ALPACA}/scripts/create_example.sh ${GUIDE} ${EXAMPLE}

echo
echo "AFTER CREATE EXAMPLE"
ls -l

echo "RUN ALPACA"
${ALPACA}/scripts/run_alpaca.sh ${GUIDE} ${SIMULATED_LONG_READS}
# The longest running step by far is Celera Assembler utgcns: 30 cpu min.
# To speed this up, use the Celera Assembler spec file to parallelize consensus.
