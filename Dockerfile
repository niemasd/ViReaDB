# Minimal Docker image for samtools v1.16.1 using Alpine base
FROM alpine:3.13.5
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>

# install samtools + minimap2 + vireadb
RUN apk update && \

    # install packages
    apk add --no-cache bash bzip2-dev gcc g++ make musl-dev python3 py3-pip xz-dev zlib-dev && \

    # install samtools
    wget -qO- "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2" | tar -xj && \
    cd samtools-* && \
    ./configure --without-curses && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-* && \

    # install minimap2
    wget -qO- "https://github.com/lh3/minimap2/archive/refs/tags/v2.24.tar.gz" | tar -zx && \
    cd minimap2-* && \
    make && \
    chmod a+x minimap2 && \
    mv minimap2 /usr/local/bin/minimap2 && \
    cd .. && \
    rm -rf minimap2-* &&\

    # install vireadb
    pip3 install --no-cache-dir --update vireadb
