# Multi-stage Dockerfile for PHLAWD and phlawd_db_maker
FROM ubuntu:22.04 AS builder

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    gcc \
    make \
    cmake \
    wget \
    git \
    autoconf \
    automake \
    libtool \
    libsqlite3-dev \
    zlib1g-dev \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /build

# Copy phlawd source and dependencies
COPY . /build/phlawd/

# Build MAFFT
WORKDIR /build/phlawd/deps
RUN tar -xzf mafft-6.864-without-extensions-src.tgz && \
    if [ -d "mafft-6.864-without-extensions/core" ]; then \
        cd mafft-6.864-without-extensions/core && \
        sed -i 's/PREFIX = \/usr\/local/PREFIX = \/usr\/local/' Makefile && \
        make && \
        make install || (mkdir -p /usr/local/bin && cp mafft /usr/local/bin/); \
    else \
        cd mafft-6.864-without-extensions && \
        make && \
        mkdir -p /usr/local/bin && \
        find . -name mafft -executable -exec cp {} /usr/local/bin/ \; || \
        find . -name "mafft*" -type f -executable -exec cp {} /usr/local/bin/ \; ; \
    fi

# Build MUSCLE
WORKDIR /build/phlawd/deps
RUN tar -xzf muscle3.8.31_src.tar.gz && \
    cd muscle3.8.31/src && \
    mkdir -p /usr/local/bin && \
    make && \
    cp muscle /usr/local/bin/ && \
    chmod +x /usr/local/bin/muscle

# Build QuickTree
WORKDIR /build/phlawd/deps
RUN tar -xzf quicktree.tar.gz && \
    cd quicktree_1.1 && \
    make && \
    cp bin/quicktree /usr/local/bin/ && \
    chmod +x /usr/local/bin/quicktree

# Build sqlitewrapped for Linux
WORKDIR /build/phlawd/deps
RUN tar -xzf sqlitewrapped-1.3.1.tar.gz && \
    cd sqlitewrapped-1.3.1 && \
    make clean || true && \
    make && \
    mkdir -p ../linux && \
    cp libsqlitewrapped.a ../linux/

# Build PHLAWD
WORKDIR /build/phlawd/src
RUN autoreconf -fi && \
    ./configure && \
    make && \
    mkdir -p /usr/local/bin && \
    cp PHLAWD /usr/local/bin/ && \
    chmod +x /usr/local/bin/PHLAWD

# Clone and build phlawd_db_maker
WORKDIR /build
RUN git clone https://github.com/jonchang/phlawd_db_maker.git && \
    cd phlawd_db_maker && \
    cmake . && \
    make && \
    find . -maxdepth 1 -type f -executable -name "phlawd_db_maker*" -exec cp {} /usr/local/bin/ \;

# Stage all binaries in a single directory for easy copying
RUN mkdir -p /build/binaries && \
    cp /usr/local/bin/mafft* /build/binaries/ 2>/dev/null || true && \
    cp /usr/local/bin/muscle /build/binaries/ && \
    cp /usr/local/bin/quicktree /build/binaries/ && \
    cp /usr/local/bin/PHLAWD /build/binaries/ && \
    cp /usr/local/bin/phlawd_db_maker* /build/binaries/ 2>/dev/null || true

# Final stage - runtime image
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libsqlite3-0 \
    libgomp1 \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy all binaries from builder staging directory
COPY --from=builder /build/binaries/ /usr/local/bin/

# Set working directory
WORKDIR /workspace

# Default command
CMD ["/bin/bash"]

