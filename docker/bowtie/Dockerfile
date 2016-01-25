FROM centos:centos6
MAINTAINER JasonRafeMiller <Jason.Rafe.Miller@gmail.com>

# this script requires internet access
# confirm the docker-machine can ping internet hosts
# some configuration may be required to work with a VPN or proxy

# the yum package manager is ideal for centOS
# the -y means assume yes to interactive prompts
RUN yum -y update
RUN yum -y install wget
RUN yum -y install zip unzip

# bowtie linux binaries from SourceForge
RUN wget -O /bowtie.zip \
    http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip/download
RUN mkdir -p /usr/local/packages && \
    mv /bowtie.zip /usr/local/packages
RUN cd /usr/local/packages && \
    unzip bowtie.zip && \
    rm bowtie.zip
RUN cd /usr/local/bin && \
    ln -s /usr/local/packages/bowtie2-2.2.3/bowtie* .
RUN type -a bowtie2

