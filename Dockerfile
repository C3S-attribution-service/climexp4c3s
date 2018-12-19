FROM centos/devtoolset-7-toolchain-centos7:7
USER root

MAINTAINER KNMI <info@knmi.nl>

# Install system libraries
RUN yum update -y && yum install -y \
    epel-release deltarpm

# Install dependencies from repository
RUN yum install -y \
    make \
    netcdf-devel \
    lapack-devel.x86_64 \
    m4

# Set environment
ENV PKG_CONFIG_PATH /build/lib/pkgconfig/
ENV FORTRAN_FLAGS "-I/usr/include -I/build/include/ -I/build/include/fgsl/ -L/build/lib -L/usr/lib64"
ENV CPPFLAGS "-I/usr/include -I/build/include/ -I/build/include/fgsl/ -L/build/lib -L/usr/lib64"
ENV LD_LIBRARY_PATH "/build/lib:$LD_LIBRARY_PATH"

# Install HDF5
WORKDIR /src
RUN curl -L https://s3.amazonaws.com/hdf-wordpress-1/wp-content/uploads/manual/HDF5/HDF5_1_10_4/hdf5-1.10.4.tar.gz > /src/hdf5-1.10.4.tar.gz && tar -xzvf hdf5-1.10.4.tar.gz
RUN cd /src/hdf5-1.10.4 && ./configure --prefix /build && make -j4 && make install

# Install NetCDF
WORKDIR /src
RUN curl -L https://github.com/Unidata/netcdf-c/archive/v4.6.1.tar.gz > /src/netcdf-4.6.1.tar.gz && tar -xzvf netcdf-4.6.1.tar.gz
RUN cd /src/netcdf-c-4.6.1 && ./configure --prefix /build && make -j4 && make install
# 
# # Install NetCDF fortran
# WORKDIR /src
# RUN curl -L ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz > netcdf-fortran-4.4.4.tar.gz && tar -xzvf netcdf-fortran-4.4.4.tar.gz
# RUN cd /src/netcdf-fortran-4.4.4 && ./configure --prefix /build && make -j4 && make install
#     
# # Install gsl
# WORKDIR /src
# RUN curl -L "ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz" > gsl-2.5.tar.gz && tar -xzvf gsl-2.5.tar.gz 
# RUN cd /src/gsl-2.5 && ./configure --prefix /build && make -j4 && make install
# 
# # Install fgsl
# WORKDIR /src
# RUN curl -L "https://doku.lrz.de/download/attachments/28051060/fgsl-1.2.0.tar.gz" > fgsl.tar.gz && tar -xzvf fgsl.tar.gz 
# RUN cd /src/fgsl-1.2.0 && ./configure --prefix /build && make && make install
#  
# # Install climate explorer
# WORKDIR /src
# COPY . climexp
# ENV PVM_ARCH build
# WORKDIR /src/climexp/${PVM_ARCH}
# COPY ./Docker/Makefile.docker /src/climexp/${PVM_ARCH}/Makefile

#RUN make

CMD bash



# docker build -t climexp_numerical .
# mkdir ./data/
# wget "http://opendap.knmi.nl/knmi/thredds/fileServer/climate_explorer/nino3.nc" -O ./data/nino3.nc
# wget "http://opendap.knmi.nl/knmi/thredds/fileServer/climate_explorer/cru_ts3.22.1901.2013.pre.dat.nc" -O ./data/cru_ts3.22.1901.2013.pre.dat.nc
# docker run -m ./data:/data -it climexp_numerical /src/climexp/build/correlatefield
# docker run -v `pwd`/data:/data -it climexp_numerical
# ./correlatefield /data/cru_ts3.22.1901.2013.pre.dat.nc /data/nino3.nc mon 1:12 ave 3 out.nc
