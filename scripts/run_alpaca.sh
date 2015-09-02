#!/bin/bash

#  Run the ALPACA assembler.
#
#  It assumes that:
#   - botwie2 is in your path.
#   - samtools is in your path.
#
#  The example can be run as:
#    scripts/run_alpaca.sh yeast.guide.scf.fasta Yeast.correctedLongRead.s.fastq

#  Paths to executables.  If scripts not set, it should be discovered based on the loction of this script.

bowtie2=`which bowtie2`
bowtie2build=`which bowtie2-build`
samtools=`which samtools`
gatekeeper=`which gatekeeper`    #  Used only for dirname below
cabin=`dirname $gatekeeper`      #  Probably fails if gatekeeper isn't in your path
scripts=""

#  Try a bit harder to find the assembler locally

if [ -z $cabin ] ; then
    gatekeeper=`echo wgs-assembler/*/bin/gatekeeper`
    cabin=`dirname $gatekeeper`
fi

#  Some more contortions to get the scripts directory.

if [ -z $scripts ] ; then
    scripts=`dirname $0`
    if [ -z $scripts ] ; then
        scripts=`which alpaca`
        scripts=`dirname $scripts`
    fi

    before=`pwd`
    cd $scripts
    scripts=`pwd`
    cd $before
fi

#  Check that stuff is found

if [ ! -x $bowtie2 ] ; then
    echo "bowtie2 in '$bowtie2' is not found or not executable"
    exit 1
fi

if [ ! -x $samtools ] ; then
    echo "samtools in '$samtools' is not found or not executable"
    exit 1
fi

if [ ! -x ${cabin}/gatekeeper ] ; then
    echo "Celera Assembler not found in $cabin"
    exit 1
fi

if [ ! -x ${scripts}/run_alpaca.sh ] ; then
    echo "alpaca.sh in '$scripts/' is not found"
    exit 1
fi

if [ $# -ne 2 ] ; then
    echo "Usage: $0 <guide_scaffold_fasta> <corrected_reads_fastq>"
    exit 1
fi


wrk=${PWD}
asm="alpaca"

INPUT_SCAFFOLDS=$1
INPUT_COREADS=$2

echo ""
echo "Run the ALPACA assembler starting at `date`"
echo ""
echo " - bowtie2  version `$bowtie2 --version | 2>&1 grep version | grep -v Comp` in $bowtie2"
echo " - samtools version `$samtools --version | 2>&1 grep samtools` in $samtools"
echo " - wgs-assembler in $cabin/"
echo " - alpaca in $scripts/"
echo ""
echo " - working in directory   $wrk"
echo " - file name prefix       $asm"
echo ""
echo " - input scaffolds        $INPUT_SCAFFOLDS"
echo " - input corrected reads  $INPUT_COREADS"

if [ -z "${INPUT_SCAFFOLDS}" ]; then
    echo "ERROR - input scaffolds not found"
fi

if [ -z "${INPUT_COREADS}" ]; then
    echo "ERROR - input corrected reads not found"
fi

if [ -z "${INPUT_SCAFFOLDS}" -o -z ${INPUT_COREADS} ]; then
    exit 1
fi


SIMU_MATE_READ_COVER=20
SIMU_MATE_CLONE_COVER=200
SIMU_MATE_READ_LENGTH=2000
SIMU_MATE_INSERT_SIZE=0

if [ ! -e 1-simumates ] ; then
  echo ""
  echo "ALPACA STEP: Simulator(GuideScaffolds) -> SimuMates"
  echo ""

    echo  > ${asm}.utg.spec "# Assemble Unitigs then stop"
    echo >> ${asm}.utg.spec "useGrid=0"
    echo >> ${asm}.utg.spec "scriptOnGrid=0"
    echo >> ${asm}.utg.spec "ovlThreads=2"
    echo >> ${asm}.utg.spec "batThreads=4"
    echo >> ${asm}.utg.spec "ovlConcurrency=8"
    echo >> ${asm}.utg.spec "cnsConcurrency=16"
    echo >> ${asm}.utg.spec "merSize=22 "
    echo >> ${asm}.utg.spec "merylMemory=48000  "
    echo >> ${asm}.utg.spec "#merylMemory=8192  "
    echo >> ${asm}.utg.spec "merylThreads=4"
    echo >> ${asm}.utg.spec "merDistinct=0.99"
    echo >> ${asm}.utg.spec "doOBT=1                 # use 1 for unitigs, 0 for scaffolds"
    echo >> ${asm}.utg.spec "obtErrorRate=0.04"
    echo >> ${asm}.utg.spec "ovlErrorRate=0.03"
    echo >> ${asm}.utg.spec "ovlStoreMemory=48000    # consumes 130GB"
    echo >> ${asm}.utg.spec "ovlMinLen=500"
    echo >> ${asm}.utg.spec "doFragmentCorrection=0"
    echo >> ${asm}.utg.spec "unitigger=bogart"
    echo >> ${asm}.utg.spec "utgBubblePopping=1"
    echo >> ${asm}.utg.spec "utgGraphErrorRate=0.01"
    echo >> ${asm}.utg.spec "utgGraphErrorLimit=1.25"
    echo >> ${asm}.utg.spec "utgMergeErrorRate=0.02"
    echo >> ${asm}.utg.spec "utgMergeErrorLimit=2.25"
    echo >> ${asm}.utg.spec "computeInsertSize=0     # Leave insert sizes as is"
    echo >> ${asm}.utg.spec "cnsReuseUnitigs=1       # Speed scaffolds by re-using consensus of unitigs"
    echo >> ${asm}.utg.spec "doExtendClearRanges=0   # Use 2 for scaffolds, 0 for unitigs"
    echo >> ${asm}.utg.spec "cgwDemoteRBP=0          # This is important for large unitigs"
    echo >> ${asm}.utg.spec "cnsMinFrags=10000"
    echo >> ${asm}.utg.spec "cnsMaxCoverage=100      # Save time on consens of really deep unitigs"
    echo >> ${asm}.utg.spec ""

    ${cabin}/fastqToCA \
      -libraryname corPacBio \
      -technology pacbio-corrected \
      -type sanger \
    	-reads ${INPUT_COREADS} \
    > ${PWD}/${asm}.corLongReads.frg

    echo >> ${asm}.utg.spec "${PWD}/${asm}.corLongReads.frg"

    mkdir 1-simumates

    #  If using other than three values here, update the allfragsdeleted lines below.

    for SIMU_MATE_INSERT_SIZE in 010000 040000 160000 ; do
	      file_name=${asm}.${SIMU_MATE_INSERT_SIZE}bpPE.${SIMU_MATE_READ_LENGTH}bp 

	      echo "  simulate mates ${file_name}"

	      # Question of whether to use -alowns or -alowgaps

	      ${cabin}/fastqSimulate \
                -em 0 -ei 0 -ed 0 \
                -allowns \
                -l ${SIMU_MATE_READ_LENGTH} \
                -pe ${SIMU_MATE_INSERT_SIZE} 0 \
		            -x ${SIMU_MATE_READ_COVER} \
                -X ${SIMU_MATE_CLONE_COVER} \
                -f ${INPUT_SCAFFOLDS} \
                -o 1-simumates/${file_name} \
        > 1-simumates/${file_name}.fastqSimulate.err 2>&1

	      rm -f 1-simumates/${file_name}.i.fastq
	      rm -f 1-simumates/${file_name}.c.fastq

        echo "  find reads with too many N's"
        cat 1-simumates/${file_name}.?.fastq \
        | \
        ${scripts}/list_pairs_with_n.pl 5 \
        >  1-simumates/${file_name}.mates_with_5N \
        2> 1-simumates/${file_name}.mates_with_5N.err

	      # indicate an assembly standard deviation of 20% insert size
	      assemble_sd=`expr 2 \* ${SIMU_MATE_INSERT_SIZE} / 10`

	      ${cabin}/fastqToCA \
		            -libraryname ${file_name} \
		            -insertsize ${SIMU_MATE_INSERT_SIZE} ${assemble_sd} constant \
		            -technology none \
                -type sanger \
		            -mates ${PWD}/1-simumates/${file_name}.1.fastq,${PWD}/1-simumates/${file_name}.2.fastq \
		            -nonrandom > ${PWD}/1-simumates/${file_name}.frg

        echo >> ${asm}.utg.spec "${PWD}/1-simumates/${file_name}.frg"
    done
fi






if [ ! -e 2-unitigs ] ; then
    echo 
    echo "ALPACA STEP: Assemble(correctedLongReads) to unitigs"
    echo ""

    mkdir  2-unitigs

    echo "  launch Celera Assembler, but stop after gatekeeper (logging to 2-unitigs/${asm}.utg.1.err)"
    ${cabin}/runCA -d assembly -p ${asm} stopAfter=initialStoreBuilding -s ${asm}.utg.spec > 2-unitigs/${asm}.utg.1.err 2>&1

    echo "  save mate info for restoration later"
    ${cabin}/gatekeeper \
      -dumpfragments \
      -tabular \
      assembly/${asm}.gkpStore \
	  | awk '{if (c++>0 && $4>0) print "frg iid",$2,"mateiid",$4;}' \
    > 2-unitigs/gatekeeper.restore_mates.edit
    
    echo "  delete the SimuMates to keep them out of unitig assembly"
    echo "lib iid 2 allfragsdeleted t"  > 2-unitigs/gatekeeper.delete_mates.edit
    echo "lib iid 3 allfragsdeleted t" >> 2-unitigs/gatekeeper.delete_mates.edit
    echo "lib iid 4 allfragsdeleted t" >> 2-unitigs/gatekeeper.delete_mates.edit
    ${cabin}/gatekeeper \
      --edit \
        gatekeeper.delete_mates.edit \
        assembly/${asm}.gkpStore \
	  2> 2-unitigs/gatekeeper.delete_mates.log

    echo "  restart Celera Assembler, but stop after unitig consensus (logging to 2-unitigs/${asm}.utg.2.err)"
    ${cabin}/runCA -d assembly -p ${asm} stopAfter=utgcns -s ${asm}.utg.spec > 2-unitigs/${asm}.utg.2.err 2>&1

    echo "  extract unitig consensus"
    ${cabin}/tigStore \
      -g assembly/${asm}.gkpStore \
      -t assembly/${asm}.tigStore 5 \
	    -U \
      -d consensus \
    > 2-unitigs/${asm}.corLongReads.unitigs.fasta
fi



#  It is ok if STDERR contains "Invalid CIGAR remnant" errors.  These are SimuMates with Ns and we discard them.

if [ ! -e 3-bowtie ] ; then
    echo ""
    echo "ALPACA STEP: Map SimuMates to unitigs"
    echo ""

    mkdir 3-bowtie

    echo "  bowtie2 build"
    ${bowtie2build} 2-unitigs/${asm}.corLongReads.unitigs.fasta 2-unitigs/${asm}.corLongReads.unitigs > 2-unitigs/${asm}.corLongReads.unitigs.bowtie2build.log 2>&1

    echo ""
    echo "  bowtie2 align for 10Kbp mates"
    ${bowtie2} \
      -p 16 \
      --end-to-end --sensitive --fr \
      -q --phred33 \
      -x 2-unitigs/${asm}.corLongReads.unitigs \
      --minins 7000 --maxins 13000 --un-conc 3-bowtie/bowtie.PairNotConcordant.10K.fasta \
      -1 1-simumates/${asm}.010000bpPE.2000bp.1.fastq \
      -2 1-simumates/${asm}.010000bpPE.2000bp.2.fastq \
      -S 3-bowtie/unitigs_10K.sam \
    > 3-bowtie/unitigs_10K.bowtie2.log 2>&1

    echo "  samtools convert for 10Kbp mates"
    ${samtools} view -S -b 3-bowtie/unitigs_10K.sam > 3-bowtie/unitigs_10K.bam

    echo "  coords convert for 10Kbp mates"
    ${scripts}/bam-to-nucmercoords.pl 2-unitigs/${asm}.corLongReads.unitigs.fasta 3-bowtie/unitigs_10K.bam > 3-bowtie/${asm}.010000bpPE.2000bp.coords 2> 3-bowtie/${asm}.010000bpPE.2000bp.err

    echo "  reformat for 10Kbp mates"
    cat 3-bowtie/${asm}.010000bpPE.2000bp.coords \
    | \
    ${scripts}/filter_coords_file.pl 1-simumates/${asm}.010000bpPE.2000bp.mates_with_5N 2> 3-bowtie/${asm}.010000bpPE.2000bp.mates_with_5N.err \
    | \
    ${scripts}/convert_to_iid.pl 2-unitigs/${asm}.gkpStore.fastqUIDmap 2> 3-bowtie/${asm}.010000bpPE.2000bp.reads-to-unitigs.err \
    > 3-bowtie/${asm}.010000bpPE.2000bp.reads-to-unitigs


    echo ""
    echo "  bowtie2 align for 40Kbp mates"
    ${bowtie2} \
      -p 16 \
      --end-to-end --sensitive --fr \
      -q --phred33 \
      --minins 28000 --maxins 58000 --un-conc 3-bowtie/bowtie.PairNotConcordant.40K.fasta \
      -x 2-unitigs/${asm}.corLongReads.unitigs \
      -1 1-simumates/${asm}.040000bpPE.2000bp.1.fastq \
      -2 1-simumates/${asm}.040000bpPE.2000bp.2.fastq \
      -S 3-bowtie/unitigs_40K.sam \
    > 3-bowtie/bowtie2.40K.log 2>&1

    echo "  samtools convert for 40Kbp mates"
    ${samtools} view -S -b 3-bowtie/unitigs_40K.sam > 3-bowtie/unitigs_40K.bam

    echo "  coords convert for 40Kbp mates"
    ${scripts}/bam-to-nucmercoords.pl 2-uinitigs/${asm}.corLongReads.unitigs.fasta 3-bowtie/unitigs_40K.bam > 3-bowtie/${asm}.040000bpPE.2000bp.coords 2> 3-bowtie/${asm}.040000bpPE.2000bp.err

    echo "  reformat for 40Kbp mates"
    cat 3-bowtie/${asm}.040000bpPE.2000bp.coords \
    | \
    ${scripts}/filter_coords_file.pl 1-simumates/${asm}.040000bpPE.2000bp.mates_with_5N 2> 3-bowtie/${asm}.040000bpPE.2000bp.mates_with_5N.err \
    | \
    ${scripts}/convert_to_iid.pl 2-unitigs/${asm}.gkpStore.fastqUIDmap 2> 3-bowtie/${asm}.040000bpPE.2000bp.reads-to-unitigs.err \
    > 3-bowtie/${asm}.040000bpPE.2000bp.reads-to-unitigs

    echo ""
    echo "  bowtie2 align for 160Kbp mates"
    ${bowtie2} \
      -p 16 \
      --end-to-end --sensitive --fr \
      -q --phred33 \
      --minins 112000 --maxins 208000 --un-conc 3-bowtie/bowtie.PairNotConcordant.160K.fasta \
      -x ${asm}.corLongReads.unitigs \
      -1 ${asm}.160000bpPE.2000bp.1.fastq \
      -2 ${asm}.160000bpPE.2000bp.2.fastq \
      -S 3-bowtie/unitigs_160K.sam \
    > 3-bowtie/bowtie2.160K.log 2>&1

    echo "  samtools convert for 160Kbp mates"
    ${samtools} view -S -b 3-bowtie/unitigs_160K.sam > 3-bowtie/unitigs_160K.bam

    echo "  coords convert for 160Kbp mates"
    ${scripts}/bam-to-nucmercoords.pl 2-unitigs/${asm}.corLongReads.unitigs.fasta 3-bowtie/unitigs_160K.bam > 3-bowtie/${asm}.160000bpPE.2000bp.coords 2> 3-bowtie/${asm}.160000bpPE.2000bp.err

    echo "  reformat for 160Kbp mates"
    cat 3-bowtie/${asm}.160000bpPE.2000bp.coords \
    | \
    ${scripts}/filter_coords_file.pl 1-simumates/${asm}.160000bpPE.2000bp.mates_with_5N 2> 3-bowtie/${asm}.160000bpPE.2000bp.mates_with_5N.err \
    | \
    ${scripts}/convert_to_iid.pl 2-unitigs/${asm}.gkpStore.fastqUIDmap 2> 3-bowtie/${asm}.160000bpPE.2000bp.reads-to-unitigs.err \
    > 3-bowtie/${asm}.160000bpPE.2000bp.reads-to-unitigs
fi



if [ ! -e 4-inject ] ; then
    echo ""
    echo "ALPACA STEP: Inject mappings into unitigs"
    echo ""

    mkdir 4-inject

    ${cabin}/gatekeeper -tabular -dumpinfo assembly/${asm}.gkpStore > 4-inject/gkp-1-before-edit.stats

    echo "  restore SimuMates in gkpStore"
    ${cabin}/gatekeeper --edit 2-unitigs/gatekeeper.restore_mates.edit assembly/${asm}.gkpStore > 2-unitigs/gatekeeper.restore_mates.log 2>&1

    ${cabin}/gatekeeper -tabular -dumpinfo assembly/${asm}.gkpStore > 4-inject/gkp-2-after-edit-before-inject

    echo "  inject 10 Kbp SimuMate-to-Unitig mappings into tigStore"
    ${cabin}/addReadsToUnitigs  \
      -g assembly/${asm}.gkpStore \
      -t assembly/${asm}.tigStore 5 \
      -lookup assembly/${asm}.gkpStore.fastqUIDmap \
      -m 3-bowtie/${asm}.010000bpPE.2000bp.reads-to-unitigs  \
    > 4-inject/addReadsToUnitigs.${asm}.010000bpPE.2000bp.err 2>&1

    echo "  inject 40 Kbp SimuMate-to-Unitig mappings into tigStore"
    ${cabin}/addReadsToUnitigs  \
      -g assembly/${asm}.gkpStore \
      -t assembly/${asm}.tigStore 5 \
      -lookup assembly/${asm}.gkpStore.fastqUIDmap \
      -m 3-bowtie/${asm}.040000bpPE.2000bp.reads-to-unitigs  \
    > 4-inject/addReadsToUnitigs.${asm}.040000bpPE.2000bp.err 2>&1

    echo "  inject 160 Kbp SimuMate-to-Unitig mappings into tigStore"
    ${cabin}/addReadsToUnitigs  \
      -g assembly/${asm}.gkpStore \
      -t assembly/${asm}.tigStore 5 \
      -lookup assembly/${asm}.gkpStore.fastqUIDmap \
      -m 3-bowtie/${asm}.160000bpPE.2000bp.reads-to-unitigs  \
    > 4-inject/addReadsToUnitigs.${asm}.160000bpPE.2000bp.err 2>&1

    ${cabin}/gatekeeper -tabular -dumpinfo assembly/${asm}.gkpStore > 4-inject/gkp-3-after-inject
fi


if [ ! -e 5-scaffold ] ; then
    echo ""
    echo "ALPACA STEP: finish assembly (logging to 5-scaffold/${asm}.scf.err)"
    echo ""

    mkdir 5-scaffold

    echo  > ${asm}.scf.spec "# Assemble scaffolds from unitigs"
    echo >> ${asm}.scf.spec "useGrid=0"
    echo >> ${asm}.scf.spec "scriptOnGrid=0"
    echo >> ${asm}.scf.spec "ovlThreads=2"
    echo >> ${asm}.scf.spec "batThreads=4"
    echo >> ${asm}.scf.spec "ovlConcurrency=8"
    echo >> ${asm}.scf.spec "cnsConcurrency=16"
    echo >> ${asm}.scf.spec "merSize=22"
    echo >> ${asm}.scf.spec "merylMemory=48000"
    echo >> ${asm}.scf.spec "#merylMemory=8192"
    echo >> ${asm}.scf.spec "merylThreads=4"
    echo >> ${asm}.scf.spec "merDistinct=0.99"
    echo >> ${asm}.scf.spec "doOBT=0                 # use 1 for unitigs, 0 for scaffolds"
    echo >> ${asm}.scf.spec "obtErrorRate=0.04"
    echo >> ${asm}.scf.spec "ovlErrorRate=0.03"
    echo >> ${asm}.scf.spec "ovlStoreMemory=48000    # consumes 130GB"
    echo >> ${asm}.scf.spec "ovlMinLen=500"
    echo >> ${asm}.scf.spec "doFragmentCorrection=0"
    echo >> ${asm}.scf.spec "unitigger=bogart"
    echo >> ${asm}.scf.spec "utgBubblePopping=1"
    echo >> ${asm}.scf.spec "utgGraphErrorRate=0.01"
    echo >> ${asm}.scf.spec "utgGraphErrorLimit=1.25"
    echo >> ${asm}.scf.spec "utgMergeErrorRate=0.02"
    echo >> ${asm}.scf.spec "utgMergeErrorLimit=2.25"
    echo >> ${asm}.scf.spec "computeInsertSize=0     # Leave insert sizes as is"
    echo >> ${asm}.scf.spec "cnsReuseUnitigs=1       # Speed scaffolds by re-using consensus of unitigs"
    echo >> ${asm}.scf.spec "doExtendClearRanges=2   # Use 2 for scaffolds, 0 for unitigs"
    echo >> ${asm}.scf.spec "cgwDemoteRBP=0          # This is important for large unitigs"
    echo >> ${asm}.scf.spec "cnsMinFrags=10000"
    echo >> ${asm}.scf.spec "cnsMaxCoverage=100      # Save time on consens of really deep unitigs"

    ${cabin}/runCA -d assembly -p ${asm} -s ${asm}.scf.spec > 5-scaffold/${asm}.scf.err 2>& 1
fi


echo ""
echo "ALPACA STEP: Analysis"
echo ""

CA_POSMAP_FRGSCF=${wrk}/assembly/9-terminator/${asm}.posmap.frgscf
CA_POSMAP_SCFLEN=${wrk}/assembly/9-terminator/${asm}.posmap.scflen
CA_SCAFFOLD_FASTA=${wrk}/assembly/9-terminator/${asm}.scf.fasta

WINDOW_SIZE=1
LONG_READ_HOLE_COVERAGE=1  # Define a hole as at least this many uncovered bases 
LONG_READ_TMP1=long_read.covered.intervals
LONG_READ_TMP2=long_read.coverage.txt
LONG_READ_TMP3=long_read.coverage.histogram
LONG_READ_HOLES_FILE=long_read.coverage_holes.txt

if [ ! -e LongReadCoverage ] ; then
    echo ""
    echo "Long Read Coverage"
    echo ""

    mkdir LongReadCoverage
    cd LongReadCoverage

    # Search for 10 as read UID prefix for long reads.
    # CA was given long reads first so they make up library #1 in gkpStore. 

    echo "  extract Long Reads"
    grep '^10' $CA_POSMAP_FRGSCF \
	      | awk '{print $2,$3,$4,$1;}' \
	      | sort -k1,1 -k2,3n > $LONG_READ_TMP1

    # Lengthy run
    # This warns "inputs with bad coords" for mappings where start==end. That is ok.

    echo "  compute coverage"
    ${scripts}/scaffold_intervals_to_coverage.pl $WINDOW_SIZE < $LONG_READ_TMP1 > $LONG_READ_TMP2

    # This output is informational not required

    echo "  report histogram"
    cat $LONG_READ_TMP2 | \
	      awk '{c=$3; A[c]++; if (c>M)M=c;} END {for (i=0; i<=M; i++) print i,A[i];}' \
	          > $LONG_READ_TMP3

    # Lengthy run

    echo "  list internal coverage holes"
    cat $LONG_READ_TMP2 | \
	      ${scripts}/list_coverage_holes.pl $LONG_READ_HOLE_COVERAGE \
		              > $LONG_READ_HOLES_FILE    

    cd ..
fi

SIMU_MATE_TMP1=simu_mate.tmp.1
SIMU_MATE_TMP2=simu_mate.tmp.2
SIMU_MATE_TMP3=simu_mate.tmp.3
SIMU_MATE_TMP4=simu_mate.tmp.4
SIMU_MATE_HOLE_COVERAGE=0  # Define a hole as at least this many uncovered bases
SIMU_MATE_COVERAGE_FILE=pair_coverage_per_window.txt
SIMU_MATE_HOLES_FILE=pair_coverage_holes.txt
SIMU_MATE_HOLES_PER_SCAF=holes_per_scaffold.txt

if [ ! -e SimuMateCoverage ] ; then
    echo ""
    echo "SimuMate Coverage"
    echo ""

    mkdir SimuMateCoverage
    cd SimuMateCoverage

    # Characters 2, 3, 4 at start of line correspond to 3 libraries of SimuMate
    # in the Celera Assembler gkpStore after the Long Reads which are library #1.

    echo "  generate the file of read 1 mapped intervals"
    grep '^[234]1' $CA_POSMAP_FRGSCF \
	      | awk '{delim="_"; print $2 delim substr($1,3),$3,$4,$5;}' \
	      | sort -k1,1 -k2,3n > $SIMU_MATE_TMP1

    echo "  generate the file of read 2 mapped intervals"
    grep '^[234]2' $CA_POSMAP_FRGSCF \
	      | awk '{delim="_"; print $2e delim substr($1,3),$3,$4,$5,$2,$1;}' \
	      | sort -k1,1 -k2,3n > $SIMU_MATE_TMP2

    # In Celera Assembler gkpStore, reads of a pair share same suffix
    # but differ in their 2-char prefix e.g. 2100003 and 2200003 are third pair in second library.
    # Exploit this for the sort.

    echo "  join on pair ID"
    join $SIMU_MATE_TMP1 $SIMU_MATE_TMP2 > $SIMU_MATE_TMP3

    #Exploit the fact that posmap always lists lo before hi in an interval.
    #Still, we must check which interval is lower and which is higher.
    #Check relative mate orientations are not equal (should be f/r or r/f).

    echo "  report union of intervals of two reads of each pair"
    cat $SIMU_MATE_TMP3 \
	      | awk '{lo=$2; if ($5<lo)lo=$5; hi=$3; if ($6>hi) hi=$6; if ($4!=$7) print $8,lo,hi,$9;}' \
	            > $SIMU_MATE_TMP4

    # Note we are missing data for zero coverage at scaffold start and end.
    # We infer zero coverage at start by position offset.
    # We could infer zero coverage at end using scaffold length but we do not.

    echo "  print average pair coverage per window"
    ${scripts}/scaffold_intervals_to_coverage.pl $WINDOW_SIZE \
	            < $SIMU_MATE_TMP4 \
	            > $SIMU_MATE_COVERAGE_FILE

    # List all internal coverage holes.
    # These are holes in SimuMate coverage excluding leading and trailing regions.

    echo "  list coverage holes (hole defined as $HOLE_COVERAGE or less)"
    cat $SIMU_MATE_COVERAGE_FILE | \
	      ${scripts}/list_coverage_holes.pl \
		              $SIMU_MATE_HOLE_COVERAGE \
		              > $SIMU_MATE_HOLES_FILE

    # This is a non-required informational report

    echo "  report holes per scaffold with length"
    cut -d ' ' -f 1 $SIMU_MATE_HOLES_FILE | sort | uniq -c | sort -k 2 \
	      | join -1 2 -2 1 - $CA_POSMAP_SCFLEN > $SIMU_MATE_HOLES_PER_SCAF    

    cd ..

    ln -s SimuMateCoverage/${SIMU_MATE_COVERAGE_FILE}  .
    ln -s SimuMateCoverage/${SIMU_MATE_HOLES_FILE}     .
fi

CHIMER_READS_IN_HOLES=long_reads_in_holes.txt

if [ ! -e IdentifyChimer ] ; then
    echo ""
    echo "Identify Chimer Suspects"
    echo ""

    mkdir IdentifyChimer
    cd IdentifyChimer

    echo "  find reads in holes"
    ${scripts}/find_reads_in_holes.pl \
		          ../LongReadCoverage/${LONG_READ_HOLES_FILE} \
		          ../SimuMateCoverage/${SIMU_MATE_HOLES_FILE} \
		          ${CA_POSMAP_FRGSCF} \
		          > ${CHIMER_READS_IN_HOLES}

    cd ..
fi

# This is minimum read UID of any SimuMate.
# This is true because long reads are loaded into gkpStore as first library.
# Assume all long reads have UIDs less than this.
FIRST_MATE=200000000000

FILTERED_SCAFFOLDS_PREFIX=filtered_scaffolds
TRUSTED_SCAFFOLDS=${FILTERED_SCAFFOLDS_PREFIX}.trusted.fasta
SUSPECT_SCAFFOLDS=${FILTERED_SCAFFOLDS_PREFIX}.suspect.fasta

if [ ! -e FilterTrusted ] ; then
    echo ""
    echo "Filter trusted scaffolds"
    echo ""

    mkdir FilterTrusted
    cd FilterTrusted

    ${scripts}/filter_scaffolds_by_read_counts.pl \
		          ${FIRST_MATE} \
		          ${CA_POSMAP_FRGSCF} \
		          ${CA_SCAFFOLD_FASTA} \
		          ${FILTERED_SCAFFOLDS_PREFIX}

    cd ..

    ln -s FilterTrusted/${TRUSTED_SCAFFOLDS} .
    ln -s FilterTrusted/${SUSPECT_SCAFFOLDS} .
fi

BREAKER_HOLES=scaffold_holes.txt
BREAKER_SCAFFOLDS_TMP1=trusted.scf.fasta
BREAKER_SCAFFOLDS_TMP2=oneline.scf.fasta
BREAKER_SCAFFOLDS_OUTPUT=alpaca.scf.fasta

if [ ! -e BreakScaffolds ] ; then
    echo ""
    echo "Break scaffolds"
    echo ""

    mkdir BreakScaffolds
    cd BreakScaffolds

    # For scaffolds generated by breaking, start UID higher than any existing scaffold.

    LAST_SCAFFOLD_ID=`cut -f 1 ${CA_POSMAP_SCFLEN} | sort -n | tail -n 1`
    NEXT_SCAFFOLD_ID=$((${LAST_SCAFFOLD_ID}+1000))

    echo LAST_SCAFFOLD_ID ${LAST_SCAFFOLD_ID}
    echo NEXT_SCAFFOLD_ID ${NEXT_SCAFFOLD_ID}

    echo "  prepare inputs"
    cut -d ' ' -f 3,5 ../IdentifyChimer/${CHIMER_READS_IN_HOLES} | tr '-' ' ' | sort -k2n | uniq > ${BREAKER_HOLES}

    echo "  convert FASTA to one line per scaffold"
    ${scripts}/fasta-to-one-line.pl < ../FilterTrusted/${TRUSTED_SCAFFOLDS} > ${BREAKER_SCAFFOLDS_TMP2}

    echo "  start breaking"
    ${scripts}/break_scaffolds.pl \
		          ${NEXT_SCAFFOLD_ID} \
		          ${BREAKER_SCAFFOLDS_TMP2} \
		          ${BREAKER_HOLES} \
		          > ${BREAKER_SCAFFOLDS_OUTPUT}

    cd ..

    ln -s BreakScaffolds/${BREAKER_SCAFFOLDS_OUTPUT} .
fi

echo ""
echo ""
echo "SIMU_MATE_COVERAGE        ${SIMU_MATE_COVERAGE_FILE}"
echo "SIMU_MATE_HOLES           ${SIMU_MATE_HOLES_FILE}"
echo "TRUSTED_SCAFFOLDS         ${TRUSTED_SCAFFOLDS}"
echo "SUSPECT_SCAFFOLDS         ${SUSPECT_SCAFFOLDS}"
echo "BREAKER_SCAFFOLDS_OUTPUT  ${BREAKER_SCAFFOLDS_OUTPUT}"
