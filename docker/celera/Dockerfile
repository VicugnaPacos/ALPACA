FROM centos:centos6
MAINTAINER JasonRafeMiller <Jason.Rafe.Miller@gmail.com>

# this script requires internet access
# confirm the docker-machine can ping internet hosts
# some configuration may be required to work with a VPN or proxy

# the yum package manager is ideal for centOS
# the -y means assume yes to interactive prompts
RUN yum -y update

# compilers
RUN yum -y install gcc gcc-c++

# compression library for C programs
RUN yum -y install zlib-devel
RUN yum -y install tar

# meryl, part of kmer, is needed for celera
# the kmer build requires many other dependencies not actually needed for celera
# python dev library for C programs needed for atac (not needed for celera)
# which needed for configure.sh
RUN yum -y install python
RUN yum -y install python-devel
RUN python --version
RUN yum -y install which

# celera pipeline implemented in perl (runCA)
RUN yum -y install perl
# statistics library used by caqc.pl, part of celera
# will use cpan to get statistics library
RUN yum -y install cpan
# use yes to preclude stops for interactive prompts
# cpan likes to use YAML (though it does work without YAML, eventually)
RUN yes | cpan YAML
# get the statistics module for perl
RUN yes | cpan Statistics::Descriptive

# source control to get Celera Assembler latest source
# side effect: installs perl 5 if not already there
RUN yum -y install subversion

# get the Celera Assembler source
# get the latest (version 8.3 after rc2)
RUN svn checkout svn://svn.code.sf.net/p/wgs-assembler/svn/trunk wgs-assembler

# get the KMER source
# we only need the meryl program but must build it all
RUN svn checkout svn://svn.code.sf.net/p/kmer/code/trunk kmer && \
    mv kmer wgs-assembler

# compile
RUN cd /wgs-assembler/kmer && sh configure.sh 
RUN cd /wgs-assembler/kmer && make && \
    make install
RUN cd /wgs-assembler/src && make
RUN cd /

# test for successful installation: expect help text
# perl wgs-assembler/Linux-amd64/bin/caqc.pl

