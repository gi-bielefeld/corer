FROM debian:12-slim AS build-stage
# Download dependencies for compilation
RUN apt update && apt -y install build-essential cmake zlib1g-dev git

# Download Bifrost, compile and install it
ARG BIFROST_VER=v1.2.1
ADD --keep-git-dir=false https://github.com/pmelsted/bifrost.git#$BIFROST_VER /usr/src/bifrost
WORKDIR /usr/src/bifrost/build
# Set maximal kmer length as build parameter
ARG MAX_K=64
RUN echo $MAX_K
RUN cmake -DMAX_KMER_SIZE=$MAX_K .. && make && make install

# Copy source code and compile corer
WORKDIR /usr/src/corer
COPY ./src ./
# Set max kmer length
RUN sed -i 's/-O3/-O3 -DMAX\_KMER\_SIZE='"$MAX_K"'/' makefile
# Compile Corer
RUN make

# Prod stage without compiler, etc.
FROM debian:12-slim AS prod-stage
# Install runtime dependency
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt update && apt -y install zlib1g procps

# Set env variables, see https://github.com/pmelsted/bifrost#troubleshooting
ENV LD_LIBRARY_PATH=/usr/local/lib/
ENV LIBRARY_PATH=/usr/local/lib/
ENV PATH=$PATH:/usr/local/lib/

# Copy binary and .so files
COPY --from=build-stage /usr/local/lib /usr/local/lib
COPY --from=build-stage /usr/src/corer/Corer /usr/local/bin
COPY --from=build-stage /usr/local/bin/Bifrost /usr/local/bin

#ENTRYPOINT [ "Corer" ]
