#!/bin/sh

#  Create an example for testing the ALPACA assembler.
#  Create simulated corrected Long Reads given a genome fasta.

#  Find the input

if [ $# -ne 2 ] ; then
    echo "Usage: $0 <guide_scaffold_fasta_input> <simulated_long_read_fastq_output>"
    exit 1
fi
INPUT_SCAFFOLDS=$1
OUTPUT_COREADS=$2
echo INPUT_SCAFFOLDS ${INPUT_SCAFFOLDS}
echo OUTPUT_COREADS ${OUTPUT_COREADS}
if [ ! -f ${INPUT_SCAFFOLDS} ] ; then
    echo "ERROR: Cannot read INPUT_SCAFFOLDS file"
    exit 1
fi
#  Find the assembler

if   [ ! -z "${CABIN}" ]; then
    echo "Using assembler in CABIN environment variable = $CABIN"

elif [ ! -z `which gatekeeper` ] ; then
    CABIN=`which gatekeeper`
    CABIN=`dirname $CABIN`
    echo "Using assembler in \$PATH = $CABIN"

else
    echo "Didn't find Celera Assembler in your \$PATH.  Set \$CABIN to the path."
    exit 1
fi

#  Make some long reads.

echo BEGIN at `date`

echo ""
echo "Creating simulated long reads in"
echo "    ${OUTPUT_COREADS}.s.fastq"
echo ""

${CABIN}/fastqSimulate \
	-f ${INPUT_SCAFFOLDS} \
	-o ${OUTPUT_COREADS} \
	-l 8000 \
	-x 10 \
	-em 0 \
	-ei 0 \
	-ed 0 \
	-se

echo ""
echo "Run ALPACA with:"
echo "  run_alpaca.sh \\"
echo "    ${INPUT_SCAFFOLDS} \\"
echo "    ${OUTPUT_COREADS}.s.fastq"
echo ""

echo "Line count, word count, char count:"
wc ${OUTPUT_COREADS}*

echo DONE at `date`
