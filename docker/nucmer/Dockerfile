FROM centos:centos6
MAINTAINER JasonRafeMiller <Jason.Rafe.Miller@gmail.com>

# this script requires internet access
# confirm the docker-machine can ping internet hosts
# some configuration may be required to work with a VPN or proxy

# the yum package manager is ideal for centOS
# the -y means assume yes to interactive prompts
RUN yum -y update
# required to download the package
RUN yum -y install wget
RUN yum -y install tar gzip
# required to build or run the package
RUN yum -y install gcc gcc-c++
RUN yum -y install perl 
RUN yum -y install csh

# nucmer source code from SourceForge
RUN wget -O /mummer.tar.gz \
    http://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz/download
RUN mkdir -p /usr/local/packages && \
    mv /mummer.tar.gz /usr/local/packages
RUN cd /usr/local/packages && \
    tar -xvzf mummer.tar.gz && \
    rm mummer.tar.gz
RUN cd /usr/local/packages/MUMmer3.23 && \
    make check
RUN cd /usr/local/packages/MUMmer3.23 && \
    make install && make
RUN cd /usr/local/bin && \
    ln -s `find /usr/local/packages/MUMmer3.23 -maxdepth 1 -type f -perm "-g=x"` .
RUN cd / && \
    type -a nucmer

