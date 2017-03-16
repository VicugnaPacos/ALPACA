#!/bin/sh

echo "CELERA ASSEMBLER CHECKOUT"

export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64:${LD_LIBRARY_PATH}
export LD_RUN_PATH=/usr/local/packages/gcc/lib64

CVS_RSH=ssh
export CVS_RSH

svn checkout svn://svn.code.sf.net/p/wgs-assembler/svn/trunk wgs-assembler

cd wgs-assembler

echo "KMER CHECKOUT"
svn checkout svn://svn.code.sf.net/p/kmer/code/trunk kmer

echo "KMER BUILD"
cd kmer
sh configure.sh
make
make install
cd ..

echo "CELERA ASSEMBLER BUILD"
cd src
make

echo "DONE"
