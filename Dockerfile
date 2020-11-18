FROM centos/devtoolset-7-toolchain-centos7:7
USER root

MAINTAINER KNMI <info@knmi.nl>

# Install system libraries
RUN yum update -y && yum install -y \
    epel-release deltarpm

# Install dependencies from repository
RUN yum install -y \
    hdf5 \
    netcdf \
    make \
    netcdf-fortran-devel \
    lapack-devel.x86_64 \
    m4
    
# Set environment
ENV PKG_CONFIG_PATH /build/lib/pkgconfig/
ENV FORTRAN_FLAGS "-I/usr/include -I/build/include/ -I/build/include/fgsl/ -L/build/lib -L/usr/lib64"
ENV CPPFLAGS "-I/usr/include -I/build/include/ -I/build/include/fgsl/ -L/build/lib -L/usr/lib64"
ENV LD_LIBRARY_PATH "/build/lib:$LD_LIBRARY_PATH"


# Add test data
ADD debiltdata /testdata

# Install gsl
WORKDIR /src
RUN curl -L "ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz" > gsl-2.5.tar.gz && tar -xzvf gsl-2.5.tar.gz 
RUN cd /src/gsl-2.5 && ./configure --prefix /build && make -j4 && make install

# Install fgsl
WORKDIR /src

# Original fgfs URL from Dockerfile is dead:
# RUN curl -L "https://doku.lrz.de/download/attachments/28051060/fgsl-1.2.0.tar.gz" > fgsl.tar.gz && tar -xzvf fgsl.tar.gz 

# Try fgsl1.2 release from github - Error:
# /bin/sh: ./configure: No such file or directory
# The command '/bin/sh -c cd /src/fgsl-1.2.0 && ./configure --prefix /build && make && make install' returned a non-zero code: 127
#RUN curl -L "https://github.com/reinh-bader/fgsl/archive/v1.2.0.tar.gz" > fgsl.tar.gz && tar -xzvf fgsl.tar.gz 

# Try new ugly URL from LRZ
RUN curl -L "https://doku.lrz.de/download/attachments/43321199/fgsl-1.2.0.tar.gz?version=1&modificationDate=1594902950587&api=v2&download=true" > fgsl.tar.gz && tar -xzvf fgsl.tar.gz 
RUN cd /src/fgsl-1.2.0 && ./configure --prefix /build && make && make install

# Install climate explorer
WORKDIR /src
COPY . climexp
ENV PVM_ARCH build
WORKDIR /src/climexp/${PVM_ARCH}
COPY ./Docker/Makefile.docker /src/climexp/${PVM_ARCH}/Makefile

RUN make

CMD bash
