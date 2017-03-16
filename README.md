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
ALPACA requires Celera Assembler 8.3 or later.
It is recommended to build Celera Assembler from source. 
(Why?
The pre-built binaries CA_8.3rc1 and CA8.3rc2 will work for any large data set.
They will fail on tiny data sets such as the ALPACA examples.
The CA source code main branch includes a patched [buildPosMap].
The patched source is also provided with ALPACA merely for documentation.
End of why.) 
To build CA from source, use scripts/generic_build_ca.sh
(after modifying the path to your C++ compiler and library).
Or, to build CA from source, follow instructions for [CA_check_out_and_compile].
Add the Celera Assembler binaries directory to your PATH.
Try running 'gatekeeper' to test whether the program is on your path.
The Celera Assembler binaries get installed to <HOME>/<ENV>/bin
where HOME is the top level directory that contains src and kmer 
and ENV describes your environment (e.g. Linux_x86_64).
Not all Celera Assembler dependencies get tested during installation.
For example, caqc.pl requires the Statistics::Descriptive perl module.
If run without this, CA fails silently after logging the following error.
```
Can't locate Statistics/Descriptive.pm in @INC (you may need to install the Statistics::Descriptive module)
```

[CA_check_out_and_compile]: http://wgs-assembler.sourceforge.net/wiki/index.php/Check_out_and_Compile
[CA_8.3rc2]: http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.3/
[buildPosMap]: https://github.com/VicugnaPacos/ALPACA/blob/master/patch/wgs-assembler/src/AS_TER/buildPosMap.C

### Install ALPACA ###

Create a local copy of the ALPACA source. This command creates and populates the ALPACA direcdtory.
```
git clone https://github.com/VicugnaPacos/ALPACA/ alpaca
```

## Test ##
The ALPACA/example_data contains a small test set. 
The file yeast.reference.fasta contains a single sequence representing the finished chromosome.
The file yeast.guide.fasta contains two sequences representing a guide assembly with two scaffolds.
The guide assembly could have been generated, for example, with ALLPATHS-LG from Illumina reads.
The commands below run the test in the directory that contains ALPACA.
The first step creates a FASTQ file of corrected long reads sampled from the given FASTA.
The second step runs ALPACA with the corrected long reads and mates sampled from the given FASTA.

### Test 0 ###
Run the scripted test called scripts/run_example.sh
(after editing the script to adjust any paths).

### Test 1 ###
Generate long reads from the one-chromosome reference.
Expect one contig.
Generate mates from the two-scaffold reference.
Expect one scaffold, identical to the contig.
```
mkdir TEST
cd TEST
../ALPACA/scripts/create_example.sh ../ALPACA/example_data/yeast.reference.fasta yeast
../ALPACA/scripts/run_alpaca.sh ../ALPACA/example_data/yeast.guide.fasta yeast.s.fastq
```

### Test 2 ###
Generate long reads from the two-scaffold assembly.
Expect two contigs.
Generate mates from the one-chromosome reference.
Expect one scaffold that joins the two contigs.
```
mkdir TEST
cd TEST
../ALPACA/scripts/create_example.sh ../ALPACA/example_data/yeast.reference.fasta yeast
../ALPACA/scripts/run_alpaca.sh ../ALPACA/example_data/yeast.guide.fasta yeast.s.fastq
```


