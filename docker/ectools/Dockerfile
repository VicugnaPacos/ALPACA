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
RUN yum -y install python

# ectools source code from SourceForge, zip of code from arxive paper
RUN wget -O /ectools.tar.gz \
    http://www.github.com/jgurtowski/ectools/archive/v0.1.tar.gz
RUN mkdir -p /usr/local/packages && \
    mv /ectools.tar.gz /usr/local/packages
RUN cd /usr/local/packages && \
    tar -xvzf ectools.tar.gz && \
    rm ectools.tar.gz
RUN cd /usr/local/packages/ectools-0.1 && \
    chmod 775 *.py *.sh
RUN cd /usr/local/bin && \
    ln -s /usr/local/packages/ectools-0.1/*.py . && \
    ln -s /usr/local/packages/ectools-0.1/*.sh .

