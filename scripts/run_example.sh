#!/bin/sh -x

HERE=.
ALPACA=../gitrepo
export CABIN=../wgs-assembler/Linux-amd64/bin/
export PATH=${CABIN}:${PATH}
GUIDE=yeast.guide.fasta
EXAMPLE=YEAST
SIMULATED_LONG_READS=${EXAMPLE}.s.fastq

echo "REMOVE ANY PREVIOUS RUNS"
rm -v *.frg
rm -v *.spec
rm -rfv [0-9]-*
rm -v filtered_scaffolds.*.fasta
rm -v pair_coverage*.txt
${ALPACA}/scripts/cleanup.sh

cp ${ALPACA}/example_data/*.fasta .

echo "BEFORE CREATE EXAMPLE"
ls -l

echo "RUN CREATE EXAMPLE..."
${ALPACA}/scripts/create_example.sh ${GUIDE} ${EXAMPLE}

echo "AFTER CREATE EXAMPLE"
ls -l

echo "RUN ALPACA"
${ALPACA}/scripts/run_alpaca.sh ${GUIDE} ${SIMULATED_LONG_READS}
# The longest running step by far is utgcns: 30 min?
# To speed this up, use the Celera Assembler spec file to parallelize consensus.
