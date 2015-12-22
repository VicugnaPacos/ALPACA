# ALPACA
A hybrid strategy for assembly of genomic DNA shotgun sequencing reads.

## Install ##
### Install samtools ###
Install samtools with a package mananger like brew or apt-get.
Alternately, build samtools with source from [htslib].
Put samtools on your PATH and run 'samtools' as a test.

[htslib]: http://www.htslib.org/download/

### Install bowtie2 ###
Install bowtie2 with a package manager like brew or apt-get.
Alternately, build bowtie2 from source or install binaries from Source Forge.

[bowtie2]: http://sourceforge.net/projects/bowtie-bio/

### Install Celera Assembler ###
It is recommended to build Celera Assembler from source. 
(ALPACA will work with the pre-build binaries for Celera Assembler version [CA_8.3rc2].
However, one Celera Assembler component will fail on the yeast example data provided with ALPACA.
The CA main branch includes a patched [buildPosMap].
The fixed source is also provided with ALPACA merely for documentation.) 
To build CA from source, follow instructions for [CA_check_out_and_compile].
Add the Celera Assembler binaries directory to your PATH.
Try running 'gatekeeper' to test whether the program is on your path.
The Celera Assembler binaries get installed to <HOME>/<ENV>/bin
where HOME is the top level directory that contains src and kmer 
and ENV describes your environment (e.g. Linux_x86_64).

[CA_check_out_and_compile]: http://wgs-assembler.sourceforge.net/wiki/index.php/Check_out_and_Compile
[CA_8.3rc2]: http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.3/
[buildPosMap]: https://github.com/VicugnaPacos/ALPACA/blob/master/patch/wgs-assembler/src/AS_TER/buildPosMap.C

### Install ALPACA ###

Create a local copy of the ALPACA source. This command creates and populates the ALPACA direcdtory.
```
git clone https://github.com/VicugnaPacos/ALPACA/
```

## Test ##
The ALPACA/example_data contains a small test set. 
The file yeast.reference.fasta contains a single sequence representing the finished chromosome.
The file yeast.guide.fasta contains two sequences representing a guide assembly with two scaffolds.
The guide assembly could have been generated, for example, with ALLPATHS-LG from Illumina reads.
These commands run the test in the directory that contains ALPACA.
The first step creates a FASTQ file of corrected long reads sampled from the reference chromosome.
The second step runs ALPACA with the corrected long reads and the guide assembly.
```
mkdir TEST
cd TEST
../ALPACA/scripts/create_example.sh ../ALPACA/example_data/yeast.reference.fasta yeast
../ALPACA/scripts/run_alpaca.sh ../ALPACA/example_data/yeast.guide.fasta yeast.s.fastq
```


